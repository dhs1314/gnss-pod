#!/usr/bin/env python3
"""
Sequential EKF Reduced-Dynamic POD for GRACE-FO.

Processes epochs one at a time with STM-based state/covariance propagation,
avoiding the S-matrix cross-segment instability of batch LSQ.

State vector (ECI): [r(3), v(3), aR, aT, aN, zwd, clk, amb_1..amb_Nsv]

Usage:
  py -3.12 run_sequential_pod.py --date 2024-04-29 --hours 0.17 --interval 30
  py -3.12 run_sequential_pod.py --date 2024-04-29 --hours 0.5 --interval 30
  py -3.12 run_sequential_pod.py --date 2024-04-29 --hours 2 --interval 30
"""
import sys, os, pickle, math, argparse
from pathlib import Path
from datetime import datetime, timedelta
from collections import defaultdict
import numpy as np

sys.path.insert(0, str(Path(__file__).parent / 'src'))
from gps1b_rnx_loader import load_gps1b_rnx
from sp3_loader import get_gps_pos_from_sp3 as _sp3_get
from coordinates import ecef_to_eci, eci_to_ecef
from gravity_model import read_icgem_gfc
from sequential_filter import SequentialEKF, EKFState, N_BASE
from troposphere import ecef_to_geodetic

C_LIGHT = 299792458.0
SEC_PER_DAY = 86400.0
OMEGA_E = 7.2921151467e-5
J2000 = datetime(2000, 1, 1, 12, 0, 0)
MJD_J2000 = 51544.5
GPS_UTC_OFFSET = 18.0  # seconds, valid since 2017-01-01


def load_gnv1b(filepath):
    """Load GNV1B reference orbit (position + velocity).

    Returns (pos_dict, vel_dict) where each maps gps_sod -> np.array([x,y,z]) or [vx,vy,vz].
    """
    pos_orbit, vel_orbit = {}, {}
    with open(filepath, encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 12:
                continue
            try:
                t = float(parts[0])
                flag = parts[2]
                if flag not in ('C', 'E'):
                    continue
                X, Y, Z = float(parts[3]), float(parts[4]), float(parts[5])
                VX, VY, VZ = float(parts[9]), float(parts[10]), float(parts[11])
                if abs(X) < 1e3:
                    continue
                pos_orbit[t] = np.array([X, Y, Z])
                vel_orbit[t] = np.array([VX, VY, VZ])
            except (ValueError, IndexError):
                continue
    return pos_orbit, vel_orbit


def interpolate_ref(ref_dict, gps_sod):
    """Linear interpolation in reference dictionary (gps_sod -> vector)."""
    ts = sorted(ref_dict.keys())
    t0 = t1 = None
    for i, ti in enumerate(ts):
        if ti >= gps_sod:
            t1 = ti
            t0 = ts[i - 1] if i > 0 else None
            break
        t0 = ti
    if t1 is None:
        t0 = t1 = ts[-1]
    if t0 is None:
        t0 = ts[0]
    if t0 == t1:
        return ref_dict[t0]
    a = (gps_sod - t0) / (t1 - t0)
    return ref_dict[t0] * (1 - a) + ref_dict[t1] * a


def get_sat_geometry(sp3, sv, utc_dt, rcv_pos, clk_data=None):
    """Compute satellite geometry with light-time iteration.

    Returns (sat_pos_ecef, sat_clk, rho_corr) or (None, None, None).

    If clk_data is provided (from read_rinex_clk), satellite clock is taken
    from the RINEX CLK file instead of the SP3-interpolated value.
    """
    pos, clk, vel = _sp3_get(sp3, sv, utc_dt)
    if pos is None or abs(clk) > 0.1 * C_LIGHT:
        return None, None, None
    rho = float(np.linalg.norm(pos - rcv_pos))
    for _ in range(5):
        travel_time = rho / C_LIGHT
        tx_dt = utc_dt - timedelta(seconds=travel_time)
        pos_tx, clk_tx, vel_tx = _sp3_get(sp3, sv, tx_dt)
        if pos_tx is None:
            break
        rho_new = float(np.linalg.norm(pos_tx - rcv_pos))
        if abs(rho_new - rho) < 1e-8:
            pos, clk = pos_tx, clk_tx
            rho = rho_new
            break
        pos, clk = pos_tx, clk_tx
        rho = rho_new
    if not (1.8e7 < rho < 2.8e7):
        return None, None, None
    # Override clock with precise CLK file if available
    if clk_data is not None:
        from src.precision_products import get_clock_from_rinex_clk
        clk_precise = get_clock_from_rinex_clk(clk_data, sv, tx_dt if pos_tx is not None else utc_dt)
        if abs(clk_precise) < 0.1 * C_LIGHT:
            clk = clk_precise
    sag = (OMEGA_E / C_LIGHT) * (pos[0] * rcv_pos[1] - pos[1] * rcv_pos[0])
    return pos, clk, rho + sag


def compute_epoch_geometry(gps_sod, gps1b_raw, sp3, rcv_pos_ecef, clk_data=None):
    """Compute satellite geometry for all SVs at one epoch.

    Uses the given receiver position (filter estimate, not GNV1B reference).
    """
    utc_dt = J2000 + timedelta(seconds=gps_sod)
    ep_data = []
    recs = gps1b_raw.get(gps_sod, {})
    for sv_id, rec in recs.items():
        if 'L_if' not in rec:
            continue
        sat_pos, sat_clk, rho_corr = get_sat_geometry(sp3, sv_id, utc_dt, rcv_pos_ecef, clk_data)
        if sat_pos is None:
            continue
        # Elevation angle (approximate, using vertical component)
        el = math.asin(abs(sat_pos[2] - rcv_pos_ecef[2]) / rho_corr)
        if el < 0.087:  # ~5 deg
            continue
        # Raw L1/L2/P1/P2 for cycle slip detection (TurboEdit needs raw, not IF)
        F1, F2, C_LIGHT = 1575.42e6, 1227.60e6, 299792458.0
        L1_m = float(rec.get('L1_cyc', 0)) * C_LIGHT / F1
        L2_m = float(rec.get('L2_cyc', 0)) * C_LIGHT / F2
        p1 = float(rec.get('P1', 0))
        p2 = float(rec.get('P2', 0))
        has_slip_lli = rec.get('slip_L1', False) or rec.get('slip_L2', False)

        ep_data.append({
            'sv': sv_id,
            'sat_pos': sat_pos,     # ECEF [m]
            'sat_clk': sat_clk,     # [m]
            'rho_corr': rho_corr,   # [m] (with Sagnac)
            'el': el,               # [rad]
            'L_if_raw': float(rec['L_if']),
            'P_if_raw': float(rec['P_if']),
            'L1_raw': L1_m,         # L1 phase [m] for MW/GF
            'L2_raw': L2_m,         # L2 phase [m] for MW/GF
            'L1_cyc': float(rec.get('L1_cyc', 0)),  # L1 phase [cycles] for MW
            'L2_cyc': float(rec.get('L2_cyc', 0)),  # L2 phase [cycles] for MW
            'P1_raw': p1,           # P1/C1 pseudorange [m] for MW
            'P2_raw': p2,           # P2 pseudorange [m] for MW
            'slip_lli': has_slip_lli,  # RINEX LLI-based slip flag
        })
    return ep_data


def main():
    parser = argparse.ArgumentParser(description='Sequential EKF POD')
    parser.add_argument('--date', required=True, help='YYYY-MM-DD')
    parser.add_argument('--hours', type=float, default=0.5)
    parser.add_argument('--interval', type=float, default=30.0)
    parser.add_argument('--data-dir', default='./data')
    parser.add_argument('--grace-id', default='C')
    parser.add_argument('--gravity-nmax', type=int, default=90)
    parser.add_argument('--sigma-acc', type=float, default=1e-3,
                        help='Process noise: unmodeled acceleration per epoch [m/s^2] '
                             '(1e-3 for GGM05C-only, use 1e-7 with full dynamics)')
    parser.add_argument('--tau-emp', type=float, default=600.0,
                        help='Empirical RTN correlation time [s]')
    parser.add_argument('--sigma-emp-ss', type=float, default=1e-8,
                        help='Empirical RTN steady-state sigma [m/s^2]')
    parser.add_argument('--sigma-zwd-rw', type=float, default=1e-9,
                        help='ZWD random walk sigma [m/鈭歴]')
    parser.add_argument('--sigma-phase', type=float, default=0.20,
                        help='Phase measurement sigma [m]')
    parser.add_argument('--sigma-code', type=float, default=0.30,
                        help='Code measurement sigma [m]')
    parser.add_argument('--chi2-threshold', type=float, default=25.0,
                        help='Chi-square innovation test threshold (5-sigma=25, '
                             'alpha=0.001=10.828)')
    # 鈹€鈹€ Phase 2.3: Measurement corrections 鈹€鈹€
    parser.add_argument('--sp3-file', default=None,
                        help='Override SP3 pickle file (e.g. CODE final SP3)')
    parser.add_argument('--clk-file', default=None,
                        help='RINEX CLK file for precise satellite clocks (e.g. CODE 30s CLK)')
    parser.add_argument('--antex-file', default=None,
                        help='IGS ANTEX file path (e.g. igs14.atx)')
    parser.add_argument('--dcb-file', default=None,
                        help='CODE P1-P2 DCB file path')
    parser.add_argument('--dcb-p1c1-file', default=None,
                        help='CODE P1-C1 DCB file path (optional, auto-detected if --dcb-file given)')
    parser.add_argument('--iers-c04', default=None,
                        help='IERS C04 ERP file path (eopc04_IAU2000.txt)')
    parser.add_argument('--enable-phase-windup', action='store_true', default=False,
                        help='Enable phase wind-up correction (Wu et al. 1993)')
    parser.add_argument('--enable-relativity', action='store_true', default=False,
                        help='Enable relativistic Shapiro delay correction')
    parser.add_argument('--enable-cycle-slip', action='store_true', default=False,
                        help='Enable MW+GF cycle slip detection (TurboEdit)')
    parser.add_argument('--ar-min-epochs', type=int, default=10,
                        help='Min epochs before WL ambiguity fixing (default: 10)')
    parser.add_argument('--enable-all-corrections', action='store_true', default=False,
                        help='Enable all Phase 2.3 measurement corrections')
    # 鈹€鈹€ Dynamics mode 鈹€鈹€
    parser.add_argument('--dynamics-mode', default='simplified',
                        choices=['simplified', 'orekit'],
                        help='Dynamics model: simplified (V2.2.1) or orekit')
    # 鈹€鈹€ Batch LSQ (Phase 6.0/10.0/11.0) 鈹€鈹€
    parser.add_argument('--batch-lsq', action='store_true', default=False,
                        help='Phase 10.0: Batch solver SD NL + Pass 2 EKF')
    parser.add_argument('--batch-lsq-v2', action='store_true', default=False,
                        help='Phase 11.0 (V2.3.0): Full-arc batch LSQ (FD Jacobian)')
    parser.add_argument('--batch-ar', action='store_true', default=False,
                        help='Two-pass batch ambiguity resolution (Phase 6.0)')
    parser.add_argument('--arc-ar', action='store_true', default=False,
                        help='Phase 12.0 (V2.3.0): Arc-based ambiguity resolution + Pass 2 EKF')
    parser.add_argument('--clk1b', default=None,
                        help='CLK1B receiver clock pkl for precise clock constraint')
    parser.add_argument('--osb-file', default=None,
                        help='Phase 9.0: CODE OSB SINEX BIAS file for undifferenced NL fixing')
    args = parser.parse_args()

    # Auto-enable all corrections
    if args.enable_all_corrections:
        args.enable_phase_windup = True
        args.enable_relativity = True
        args.enable_cycle_slip = True
        # Auto-detect IERS C04 and DCB if not specified
        if args.iers_c04 is None:
            default_c04 = dp / 'IERS' / 'eopc04_IAU2000.txt'
            if default_c04.exists():
                args.iers_c04 = str(default_c04)

    date_str = args.date
    y, m, d = [int(x) for x in date_str.split('-')]
    doy = (datetime(y, m, d) - datetime(y, 1, 1)).days + 1
    dp = Path(args.data_dir)
    grace_id = args.grace_id

    # -- Load GNV1B reference orbit --
    gnv_path = dp / 'gracefo' / str(y) / date_str / f'GNV1B_{date_str}_{grace_id}_04.txt'
    ref_orbit, ref_vel = load_gnv1b(str(gnv_path))
    print(f"[GNV1B] {len(ref_orbit)} pos, {len(ref_vel)} vel epochs")

    # -- Load GPS1B observations --
    rnx_path = dp / f'GPS1B_{date_str}_{grace_id}_04.rnx'
    pkl_path = dp / 'gracefo' / str(y) / date_str / f'GPS1B_{date_str}_{grace_id}_04.pkl'
    if rnx_path.exists():
        gps1b_raw = load_gps1b_rnx(str(rnx_path))
        print(f"[RINEX] {len(gps1b_raw)} epochs")
    elif pkl_path.exists():
        gps1b_raw = pickle.load(open(str(pkl_path), 'rb'))
        print(f"[GPS1B PKL] {len(gps1b_raw)} epochs")
    else:
        print(f"[FATAL] No GPS data found")
        return

    # -- Load SP3 --
    if args.sp3_file:
        sp3_path = Path(args.sp3_file)
        if sp3_path.suffix == '.pkl':
            sp3 = pickle.load(open(str(sp3_path), 'rb'))
        else:
            # Load from text SP3 file (auto-convert to pickle cache)
            from src.sp3_loader import parse_sp3_text
            print(f"  [SP3] Parsing {sp3_path.name}...")
            with open(str(sp3_path), 'r') as f:
                epochs_dict, ts_list = parse_sp3_text(f.read())
            sp3 = {
                'ts': ts_list,
                'epochs': epochs_dict,
                'source': 'CODE',
                'product': 'FIN',
            }
            # Cache as pickle for next run
            cache_pkl = sp3_path.with_suffix('.pkl')
            pickle.dump(sp3, open(str(cache_pkl), 'wb'))
            print(f"  [SP3] Cached 鈫?{cache_pkl.name}")
    else:
        sp3_pkl = dp / str(y) / f'{doy:03d}' / 'igs_sp3_FIN.pkl'
        if not sp3_pkl.exists():
            sp3_pkl = dp / str(y) / f'{doy:03d}' / 'igs_sp3.pkl'
        sp3 = pickle.load(open(str(sp3_pkl), 'rb'))
    print(f"[SP3] {len(sp3['ts'])} epochs ({sp3.get('source', 'IGS')} {sp3.get('product', '')})")

    # -- Load RINEX CLK (optional) --
    clk_data = None
    if args.clk_file:
        from src.precision_products import read_rinex_clk
        clk_path = args.clk_file
        if not Path(clk_path).exists():
            # Try gunzipped version
            clk_gz = Path(str(clk_path) + '.gz')
            if clk_gz.exists():
                import gzip
                with gzip.open(str(clk_gz), 'rt') as gf:
                    tmp_path = str(clk_gz).replace('.gz', '')
                    with open(tmp_path, 'w') as tf:
                        tf.write(gf.read())
                clk_path = tmp_path
        if Path(clk_path).exists():
            clk_data = read_rinex_clk(clk_path)
            n_sv = len(clk_data)
            print(f"[CLK] {n_sv} SVs, {len(clk_data[list(clk_data.keys())[0]]['epochs']) if clk_data else 0} epochs (CODE 30s)")
        else:
            print(f"[CLK] WARNING: CLK file not found: {args.clk_file}")

    # -- Select epochs --
    gps_sod_start = min(gps1b_raw.keys())
    gps_sod_end = gps_sod_start + args.hours * 3600
    epochs = []
    for gps_sod in sorted(gps1b_raw.keys()):
        if not (gps_sod_start <= gps_sod <= gps_sod_end):
            continue
        dt_ep = gps_sod - gps_sod_start
        nearest = round(dt_ep / args.interval) * args.interval
        if abs(dt_ep - nearest) > max(2.0, args.interval * 0.1):
            continue
        epochs.append(gps_sod)
    print(f"[EPOCH] {len(epochs)} selected (dt={args.interval}s)")

    if len(epochs) < 2:
        print("[FATAL] Need at least 2 epochs")
        return

    # -- Estimate per-SV code biases --
    # Use a short window (60 epochs 鈮?30 min) to keep receiver clock roughly constant.
    # Biases are zero-meaned: the common receiver clock/DCB is absorbed by clk state.
    N_BIAS = min(60, len(epochs))

    # Load PCO correction function if ANTEX is available
    _pco_fn = None
    if args.antex_file:
        from src.precision_products import read_antex, get_satellite_pco
        from src.measurement_corrections import compute_pco_ecef_from_nadir
        _antex = read_antex(args.antex_file)
        _pco_fn = lambda sv, pos: pos + compute_pco_ecef_from_nadir(
            pos, float(get_satellite_pco(_antex, sv, 'L1')[2]))
        print(f"  [ANTEX] PCO corrections enabled for bias estimation")

    # Load CODE DCB data
    _dcb_if = {}
    if args.dcb_file:
        p1c1_path = args.dcb_p1c1_file
        if p1c1_path is None:
            # Auto-detect P1C1 file from P1P2 filename
            p1c1_path = args.dcb_file.replace('P1P2', 'P1C1')
        from src.precision_products import load_code_dcb_pair, compute_dcb_if_correction
        dcb_pair = load_code_dcb_pair(p1c1_path, args.dcb_file)
        for prn in dcb_pair:
            _dcb_if[prn] = compute_dcb_if_correction(dcb_pair, prn)
        print(f"  [DCB] {len(_dcb_if)} GPS satellites loaded")

    print(f"\n-- Estimating per-SV code biases (first {N_BIAS} of {len(epochs)} epochs) --")
    sv_p_residuals = defaultdict(list)
    for gps_sod in epochs[:N_BIAS]:
        ref_pos = interpolate_ref(ref_orbit, gps_sod)
        utc_dt = J2000 + timedelta(seconds=gps_sod)
        recs = gps1b_raw.get(gps_sod, {})
        for sv_id, rec in recs.items():
            if 'L_if' not in rec:
                continue
            sat_pos, sat_clk, rho_corr = get_sat_geometry(sp3, sv_id, utc_dt, ref_pos, clk_data)
            if sat_pos is None:
                continue
            # Apply satellite PCO to ensure consistency with measurement model
            if _pco_fn is not None:
                sat_pos = _pco_fn(sv_id, sat_pos)
                rho_corr = (float(np.linalg.norm(sat_pos - ref_pos))
                            + (OMEGA_E / C_LIGHT) * (sat_pos[0] * ref_pos[1]
                                                     - sat_pos[1] * ref_pos[0]))
            P_if_raw = float(rec['P_if'])
            # Apply CODE DCB IF correction
            dcb_corr = _dcb_if.get(sv_id, 0.0)
            P_r = P_if_raw + dcb_corr + sat_clk - rho_corr
            sv_p_residuals[sv_id].append(P_r)

    sv_bias = {}
    sv_bias_ref = 0.0  # reference mean subtracted from all biases
    for sv, p_vals in sv_p_residuals.items():
        if len(p_vals) >= 10:
            sv_bias[sv] = float(np.median(p_vals))
    if sv_bias:
        sv_bias_ref = float(np.mean(list(sv_bias.values())))
        for sv in sv_bias:
            sv_bias[sv] -= sv_bias_ref
    print(f"  Per-SV code biases: {len(sv_bias)} SVs (zero-meaned, ref={sv_bias_ref:.3f}m)")

    # -- Load gravity model --
    GRAVITY_NMAX = args.gravity_nmax
    gravity_path = dp / 'gravity' / 'GGM05C.gfc'
    if gravity_path.exists():
        Cnm, Snm, _, GM_grav, R_grav = read_icgem_gfc(str(gravity_path))
        print(f"  Gravity: GGM05C Nmax={GRAVITY_NMAX}")
    else:
        print(f"  [FATAL] Gravity model not found: {gravity_path}")
        return

    # -- Initial state from GNV1B --
    r0_ecef = interpolate_ref(ref_orbit, epochs[0])
    v0_ecef = interpolate_ref(ref_vel, epochs[0])
    mjd_start = MJD_J2000 + (epochs[0] - GPS_UTC_OFFSET) / SEC_PER_DAY
    mjd_tt_start = mjd_start + 69.184 / SEC_PER_DAY

    r0_eci, v0_eci = ecef_to_eci(r0_ecef, v0_ecef, mjd_start)
    print(f"\n  Initial ECEF: r0=[{r0_ecef[0]/1000:.1f},{r0_ecef[1]/1000:.1f},{r0_ecef[2]/1000:.1f}]km "
          f"|v0|={np.linalg.norm(v0_ecef)/1000:.3f}km/s")

    # 鈹€鈹€ CLK1B receiver clock (Phase 16.0) 鈹€鈹€
    # Remove constant bias: use 未clk = clk(t) - clk(t0)
    _clk1b_data = None
    _clk1b_ref = 0.0
    if args.clk1b:
        _clk1b_raw = pickle.load(open(args.clk1b, 'rb'))
        # Reference clock at first epoch
        first_sod = min(_clk1b_raw.keys())
        _clk1b_ref = float(_clk1b_raw[first_sod]) * 299792458.0  # s 鈫?m
        # Store bias-removed clock: negative because clock offset
        # in CLK1B: GPS_time = rcv_time + eps_time 鈫?eps_time is negative
        _clk1b_data = {}
        for sod, val in _clk1b_raw.items():
            _clk1b_data[sod] = float(val) * 299792458.0 - _clk1b_ref  # relative to t0
        print(f"  [CLK1B] {len(_clk1b_data)} epochs, ref={_clk1b_ref:.0f}m removed")

    # -- Initialize EKF --
    ekf_config = {
        'clk1b_data': _clk1b_data,
        'clk1b_sigma': 0.50,  # 蟽=0.50m (soft constraint, USO drift)
        'Cd': 2.2, 'CR': 1.3,
        'area_drag': 0.68, 'area_srp': 3.4, 'mass': 580.0,
        'dt_integ': 10.0,
        'bodies': ['Sun', 'Moon'],
        'Cnm': Cnm, 'Snm': Snm, 'GM_grav': GM_grav, 'R_grav': R_grav,
        'gravity_nmax': GRAVITY_NMAX,
        'sigma_acc_process': args.sigma_acc,
        'tau_emp': args.tau_emp,
        'sigma_emp_ss': args.sigma_emp_ss,
        'sigma_zwd_rw': args.sigma_zwd_rw,
        'sigma_phase': args.sigma_phase,
        'sigma_code': args.sigma_code,
        'chi2_threshold': args.chi2_threshold if args.chi2_threshold != 25.0
                          else (100.0 if args.hours >= 0.3 else 25.0),
        'el_min': 0.087,
        # Elevation weighting: adaptive
        'elev_exp_phase': 1.0,
        'elev_exp_code': 0.70 if args.hours >= 0.3 else 1.0,
        # Clock process noise: adaptive
        'clock_rw': 0.0004 if args.hours < 0.3 else 0.001,
        # Pseudo-stochastic pulses (Phase 14.0 experimental, disabled by default)
        'tau_emp': args.tau_emp,  # 600s default
        'pulse_interval': 0,       # 0 = disabled
        'pulse_amplify': 50.0,
        # Phase 2.3 measurement corrections
        'use_phase_windup': args.enable_phase_windup,
        'use_relativity': args.enable_relativity,
        'use_cycle_slip': args.enable_cycle_slip,
        # Phase 3.0 PPP-AR
        'ar_min_epochs': args.ar_min_epochs,
        # Phase 3.1 Dynamics mode
        'dynamics_mode': args.dynamics_mode,
        'orekit_gravity_degree': args.gravity_nmax,  # reuse --gravity-nmax for Orekit
        'gravity_field': str(dp / 'gravity' / 'GGM05C.gfc'),  # GGM05C for Orekit
    }

    # Load ANTEX data if available
    if args.antex_file:
        from src.precision_products import read_antex
        print(f"  [ANTEX] Loading {args.antex_file}...")
        antex_data = read_antex(args.antex_file)
        ekf_config['antex_data'] = antex_data
        print(f"  [ANTEX] {len(antex_data.get('G', {}))} GPS satellites loaded")

    # Load DCB data if available
    if args.dcb_file:
        p1c1_path = args.dcb_p1c1_file
        if p1c1_path is None:
            p1c1_path = args.dcb_file.replace('P1P2', 'P1C1')
        from src.precision_products import load_code_dcb_pair
        print(f"  [DCB] Loading P1C1: {p1c1_path}")
        print(f"  [DCB] Loading P1P2: {args.dcb_file}")
        dcb_pair = load_code_dcb_pair(p1c1_path, args.dcb_file)
        ekf_config['dcb_data'] = dcb_pair
        print(f"  [DCB] {len(dcb_pair)} satellites loaded")

    # 鈹€鈹€ Phase 9.0: CODE OSB integer clock (undifferenced AR) 鈹€鈹€
    if args.osb_file:
        from src.batch_lsq import read_code_osb
        osb_wl, osb_nl = read_code_osb(args.osb_file)
        ekf_config['osb_wl'] = osb_wl
        ekf_config['osb_nl'] = osb_nl
        print(f"  [OSB] {len(osb_wl)} SVs: undifferenced NL fixing enabled")

    # Load IERS C04 Earth Orientation Parameters
    if args.iers_c04:
        from src.precision_products import setup_iers_from_c04
        print(f"  [IERS] Loading {args.iers_c04}...")
        setup_iers_from_c04(args.iers_c04)

    # 鈹€鈹€ Batch LSQ (Phase 6.0) 鈹€鈹€
    ekf = SequentialEKF(ekf_config)
    state = ekf.initialize(r0_eci, v0_eci, mjd_start, epochs[0])

    # Track each SV's B_if at the moment it first gets WL-fixed in pass 1.
    # Used by batch AR to compute clock-consistent absolute ambiguities.
    _sv_first_wl_bif = {}  # sv -> B_if_float
    _wl_fixed_prev = set()
    # Store per-epoch geometry for batch linear solver (Phase 6.0)
    _pass1_geometry = []      # list of ep_data per epoch
    _pass1_r_ecef = []        # list of r_ecef per epoch

    # 鈹€鈹€ Arc-based AR tracker (Phase 12.0) 鈹€鈹€
    if args.arc_ar:
        from src.arc_ambiguity import ArcTracker
        from src.sequential_filter import I_AMB_START
        _arc_tracker = ArcTracker()

    # -- Run sequential filter --
    print(f"\n-- Sequential EKF ({args.hours}h arc, {args.interval}s interval) --")
    results = {
        'epochs': [], 'r_ecef': [], 'v_ecef': [],
        'a_rtn': [], 'zwd': [], 'clk': [],
        'stats': [], 'n_sv': [],
        'r_gnv': [],  # GNV1B position at each epoch
    }

    for i_epoch, gps_sod in enumerate(epochs):
        mjd_utc = MJD_J2000 + (gps_sod - GPS_UTC_OFFSET) / SEC_PER_DAY
        mjd_tt = mjd_utc + 69.184 / SEC_PER_DAY

        if i_epoch > 0:
            # Predict to current epoch
            dt = gps_sod - epochs[i_epoch - 1]
            mjd_utc_prev = MJD_J2000 + (epochs[i_epoch - 1] - GPS_UTC_OFFSET) / SEC_PER_DAY
            mjd_tt_prev = mjd_utc_prev + 69.184 / SEC_PER_DAY
            state = ekf.predict(state, gps_sod, mjd_utc_prev, mjd_tt_prev)

        # Compute receiver position in ECEF for geometry computation
        rcv_pos_ecef, _ = eci_to_ecef(state.r_eci, state.v_eci, mjd_utc)

        # Compute geometry using filter's predicted position
        ep_data = compute_epoch_geometry(gps_sod, gps1b_raw, sp3, rcv_pos_ecef, clk_data)

        if not ep_data:
            continue

        # Process epoch
        state, stats = ekf.process_epoch(state, ep_data, sp3, sv_bias, sv_bias_ref,
                                          mjd_utc, mjd_tt, doy)

        # Record results
        r_ecef, v_ecef = eci_to_ecef(state.r_eci, state.v_eci, mjd_utc)
        r_gnv = interpolate_ref(ref_orbit, gps_sod)

        results['epochs'].append(gps_sod)
        results['r_ecef'].append(r_ecef)
        results['v_ecef'].append(v_ecef)
        results['a_rtn'].append(state.a_rtn.copy())
        results['zwd'].append(state.zwd)
        results['clk'].append(state.clk)
        results['stats'].append(stats)
        results['n_sv'].append(state.n_sv)
        results['r_gnv'].append(r_gnv)

        # Store geometry for batch linear solver (Phase 6.0/11.0)
        if args.batch_ar or args.batch_lsq or args.batch_lsq_v2:
            # Compute geometric range base (without clock/zwd/amb) for each SV.
            # geo_base = |sat_eci - rcv_eci| + sag - sat_clk
            # This allows the batch solver to solve for clock/zwd/amb
            # by subtracting this known part from the observation.
            _OM = 7.2921151467e-5
            _CL = 299792458.0
            # Compute ZHD once per epoch
            lat_rad = math.asin(r_ecef[2] / np.linalg.norm(r_ecef))
            h_m = np.linalg.norm(r_ecef) - 6378137.0
            from src.troposphere import saastamoinen_zhd
            _zhd = saastamoinen_zhd(lat_rad, h_m)
            for d in ep_data:
                sv = d['sv']
                sat_ecef = np.asarray(d['sat_pos'], dtype=float)
                sat_clk = float(d.get('sat_clk', 0))
                sat_eci, _ = ecef_to_eci(sat_ecef, np.zeros(3), mjd_utc)
                rho_eci = np.linalg.norm(sat_eci - state.r_eci)
                sag = (_OM / _CL) * (sat_ecef[0] * r_ecef[1]
                                      - sat_ecef[1] * r_ecef[0])
                el = float(d.get('el', 0.5))
                mf_h = 1.001 / np.sqrt(0.002001 + np.sin(el)**2)
                # Full geometric range without clock/zwd/amb
                # = ECI range + Sagnac - sat_clk + ZHD*mf_h
                geo_full = rho_eci + sag - sat_clk + _zhd * mf_h
                d['_geo_full'] = geo_full
                # Observation corrections (same as EKF applies)
                dcb = ekf._dcb_if.get(sv, 0.0)
                d['_obs_code'] = float(d.get('P_if_raw', 0)) + dcb - sv_bias.get(sv, 0.0)
                d['_obs_phase'] = float(d.get('L_if_raw', 0)) - sv_bias.get(sv, 0.0)
            _pass1_geometry.append(ep_data)
            _pass1_r_ecef.append(r_ecef)

        # Track B_if at first WL-fix for batch AR (Phase 6.0)
        if args.batch_ar or args.batch_lsq or args.batch_lsq_v2:
            new_wl = set(ekf._wl_fixed.keys()) - _wl_fixed_prev
            for sv in new_wl:
                if sv in state.sv_to_idx:
                    idx = state.sv_to_idx[sv]
                    amb_idx = 11 + idx  # I_AMB_START
                    _sv_first_wl_bif[sv] = float(state.x[amb_idx])
            _wl_fixed_prev = set(ekf._wl_fixed.keys())

        # 鈹€鈹€ Arc-based AR tracking (Phase 12.0) 鈹€鈹€
        if args.arc_ar:
            for d in ep_data:
                sv = d['sv']
                if sv not in state.sv_to_idx:
                    continue
                L1c = float(d.get('L1_cyc', 0))
                L2c = float(d.get('L2_cyc', 0))
                P1 = float(d.get('P1_raw', d.get('P_if_raw', 0)))
                P2 = float(d.get('P2_raw', d.get('P_if_raw', 0)))
                idx = state.sv_to_idx[sv]
                B_if = float(state.x[I_AMB_START + idx])
                has_slip = False
                _arc_tracker.add_epoch(sv, L1c, L2c, P1, P2, B_if,
                                        mjd_utc, has_slip)

        if i_epoch % max(1, len(epochs) // 10) == 0 or i_epoch < 3:
            pos_diff = np.linalg.norm(r_ecef - r_gnv)
            print(f"  [{i_epoch:4d}/{len(epochs)}] t={gps_sod:.0f}  "
                  f"|dr|={pos_diff:.3f}m  "
                  f"aR={state.a_rtn[0]:.2e} aT={state.a_rtn[1]:.2e} aN={state.a_rtn[2]:.2e}  "
                  f"zwd={state.zwd:.4f}m  "
                  f"phs={stats['n_phase']}/{stats['n_phase']+stats['n_rej']} "
                  f"rms_phs={stats['rms_phase']:.4f}m  "
                  f"n_sv={state.n_sv}  WL={stats.get('n_wl_fixed', 0)}")

    # -- Summary --
    print(f"\n-- Results --")
    pos_diffs = np.array([np.linalg.norm(results['r_ecef'][i] - results['r_gnv'][i])
                           for i in range(len(results['epochs']))])
    rms_3d = np.sqrt(np.mean(pos_diffs**2))
    print(f"  3D RMS vs GNV1B: {rms_3d:.3f} m  "
          f"(mean={np.mean(pos_diffs):.3f}, max={np.max(pos_diffs):.3f})")

    # Mean empirical accelerations
    a_rtn_all = np.array(results['a_rtn'])
    print(f"  Mean aR={np.mean(a_rtn_all[:,0]):.3e} "
          f"aT={np.mean(a_rtn_all[:,1]):.3e} "
          f"aN={np.mean(a_rtn_all[:,2]):.3e} m/s^2")
    print(f"  Mean ZWD={np.mean(results['zwd']):.4f} m")

    # Innovation stats
    n_phase_total = sum(s['n_phase'] for s in results['stats'])
    n_code_total = sum(s['n_code'] for s in results['stats'])
    n_rej_total = sum(s['n_rej'] for s in results['stats'])
    rms_phase_total = np.sqrt(np.mean([s['rms_phase']**2 for s in results['stats']]))
    rms_code_total = np.sqrt(np.mean([s['rms_code']**2 for s in results['stats']]))
    print(f"  Phase: {n_phase_total} accepted, RMS={rms_phase_total:.4f}m")
    print(f"  Code:  {n_code_total} accepted, RMS={rms_code_total:.4f}m")
    print(f"  Rejected: {n_rej_total}")
    print(f"  Final SV count: {results['n_sv'][-1] if results['n_sv'] else 0}")
    rms_float = rms_3d  # save for comparison

    # 鈹€鈹€ Framework v3: Single-pass Batch Solver (Phase 15.0) 鈹€鈹€
    # Fixed EKF orbit + batch clock/zwd/amb 鈫?one-pass final solution.
    # Compares phase residuals: batch clock/amb vs EKF clock/amb.
    if getattr(args, 'batch_v3', False):
        from src.batch_solver import BatchLinearSolver

        print(f"\n-- Framework v3 Batch Solver ({args.hours}h, "
              f"{len(_pass1_geometry)} epochs) --")

        bls = BatchLinearSolver(_pass1_geometry,
                                 sigma_phase=args.sigma_phase,
                                 sigma_code=args.sigma_code)
        bls_sol = bls.solve()

        clk_b = bls_sol['clk']
        zwd_b = bls_sol['zwd']
        amb_b = bls_sol['amb_dict']

        print(f"  Batch: {len(amb_b)} SVs, "
              f"phase_rms={bls_sol['rms_phase']:.3f}m, "
              f"code_rms={bls_sol['rms_code']:.3f}m")

        # Compute EKF phase residuals (at last epoch, as a comparison)
        ekf_phase_errs = []
        bch_phase_errs = []
        for i_ep, ep_list in enumerate(_pass1_geometry):
            for d in ep_list:
                sv = d['sv']
                if '_obs_phase' not in d: continue
                if '_geo_full' not in d: continue
                geo = float(d['_geo_full'])
                el = float(d.get('el', 0.5))
                mf = 1.001 / np.sqrt(0.002001 + np.sin(el)**2)
                obs = float(d['_obs_phase'])

                # EKF residual (use pass 1 results)
                ekf_clk = float(results['clk'][i_ep]) if i_ep < len(results['clk']) else 0
                ekf_zwd = float(results['zwd'][i_ep]) if i_ep < len(results['zwd']) else 0
                # EKF amb not easily accessible per-epoch... skip

                # Batch residual
                if sv in amb_b:
                    bch_model = geo + float(clk_b[i_ep]) + float(zwd_b[i_ep]) * mf + float(amb_b[sv])
                    bch_phase_errs.append(obs - bch_model)

        if bch_phase_errs:
            bch_rms = np.sqrt(np.mean([e**2 for e in bch_phase_errs]))
            ekf_rms = results['stats'][-1]['rms_phase'] if results['stats'] else 0
            print(f"  Phase residuals: EKF={ekf_rms:.3f}m  Batch={bch_rms:.3f}m")
            if ekf_rms > 0:
                print(f"  Improvement: {(ekf_rms - bch_rms)/ekf_rms*100:+.1f}%")

        print(f"  Note: orbit is EKF trajectory (3D RMS = {rms_float:.3f}m unchanged)")

        # Save
        out_d = Path('results') / 'batch_v3'
        out_d.mkdir(parents=True, exist_ok=True)
        pickle.dump({
            'rms_ekf_orbit': rms_float,
            'clk_batch': clk_b.tolist() if hasattr(clk_b, 'tolist') else clk_b,
            'amb_batch': {k: float(v) for k, v in amb_b.items()},
            'rms_phase_batch': bls_sol['rms_phase'],
            'rms_code_batch': bls_sol['rms_code'],
            'epochs': epochs,
            'config': vars(args),
        }, open(str(out_d / f"batchv3_{date_str}_{grace_id}_{args.hours}h.pkl"), 'wb'))
        print(f"  Saved: {out_d / f'batchv3_{date_str}_{grace_id}_{args.hours}h.pkl'}")

    # 鈹€鈹€ Phase 12.0: Arc-based Ambiguity Resolution (V2.3.0) 鈹€鈹€
    if args.arc_ar:
        from src.arc_ambiguity import ArcAmbiguityResolver, LAM_NL, COEFF_W

        # Finalize arc tracking
        arc_dict = _arc_tracker.finalize()
        n_svs = len(arc_dict)
        n_arcs = sum(len(arcs) for arcs in arc_dict.values())
        print(f"\n-- Phase 12.0 Arc AR ({args.hours}h) --")
        print(f"  Arcs: {n_arcs} across {n_svs} SVs")

        # Load OSB biases if available
        osb_wl = {}; osb_nl = {}
        if args.osb_file:
            from src.batch_lsq import read_code_osb
            osb_wl, osb_nl = read_code_osb(args.osb_file)

        # Extract EKF float B_if for WL-only SVs (not in arc_dict)
        ekf_amb = {}
        for sv in state.sv_list:
            idx = state.sv_to_idx.get(sv)
            if idx is not None:
                ekf_amb[sv] = float(state.x[I_AMB_START + idx])

        # Resolve arc-level WL+NL
        resolver = ArcAmbiguityResolver(
            wl_bias=osb_wl, nl_bias=osb_nl,
            min_epochs=6, max_wl_std=0.30, max_nl_resid=0.30)
        wl_fixed, nl_fixed, amb_fixed = resolver.resolve(arc_dict, ekf_amb)

        n_wl = len(wl_fixed); n_nl = len(nl_fixed)
        print(f"  Arc WL: {n_wl} SVs, Arc NL: {n_nl} SVs "
              f"(OSB: {'yes' if osb_nl else 'no'})")

        if n_wl >= 3:
            print(f"  Pass 2: fresh EKF with {len(amb_fixed)} arc-fixed ambs")
            for sv in sorted(amb_fixed.keys())[:6]:
                w = wl_fixed.get(sv, '-')
                n = nl_fixed.get(sv, '—')
                arc_len = max((a.n_epochs for a in arc_dict.get(sv, [])), default=0)
                print(f"    {sv}: B_if={amb_fixed[sv]:.3f}m "
                      f"(N_w={w}, N1={n}, arc={arc_len}ep)")

            ekf_cfg2 = dict(ekf_config)
            ekf_cfg2['amb_batch_fixed'] = amb_fixed
            ekf_cfg2['amb_batch_var'] = 0.0004
            ekf_cfg2['use_cycle_slip'] = False

            ekf2 = SequentialEKF(ekf_cfg2)
            state2 = ekf2.initialize(r0_eci, v0_eci, mjd_start, epochs[0])

            results2 = {'epochs': [], 'r_ecef': [], 'v_ecef': [],
                        'a_rtn': [], 'zwd': [], 'clk': [],
                        'stats': [], 'n_sv': [], 'r_gnv': []}

            for i_ep, gps_sod in enumerate(epochs):
                mjd_utc = MJD_J2000 + (gps_sod - GPS_UTC_OFFSET) / SEC_PER_DAY
                mjd_tt = mjd_utc + 69.184 / SEC_PER_DAY
                if i_ep > 0:
                    mjd_prev = MJD_J2000 + (epochs[i_ep-1] - GPS_UTC_OFFSET) / SEC_PER_DAY
                    state2 = ekf2.predict(state2, gps_sod, mjd_prev,
                                           mjd_prev + 69.184 / SEC_PER_DAY)
                rcv_e, _ = eci_to_ecef(state2.r_eci, state2.v_eci, mjd_utc)
                ep_d = compute_epoch_geometry(gps_sod, gps1b_raw, sp3, rcv_e, clk_data)
                if ep_d:
                    state2, s2 = ekf2.process_epoch(state2, ep_d, sp3, sv_bias,
                                                     sv_bias_ref, mjd_utc, mjd_tt, doy)
                    r_e, v_e = eci_to_ecef(state2.r_eci, state2.v_eci, mjd_utc)
                    r_g = interpolate_ref(ref_orbit, gps_sod)
                    results2['epochs'].append(gps_sod)
                    results2['r_ecef'].append(r_e); results2['r_gnv'].append(r_g)
                    results2['n_sv'].append(state2.n_sv)
                if i_ep % max(1, len(epochs)//10) == 0 or i_ep < 3:
                    dr2 = (np.linalg.norm(results2['r_ecef'][-1] - results2['r_gnv'][-1])
                           if results2['r_gnv'] else 0)
                    print(f"  [ArcAR {i_ep:4d}/{len(epochs)}] t={gps_sod:.0f}  "
                          f"|dr|={dr2:.3f}m  n_sv={state2.n_sv}")

            rms2 = np.sqrt(np.mean([np.linalg.norm(results2['r_ecef'][i] - results2['r_gnv'][i])**2
                                    for i in range(len(results2['epochs']))]))
            imp = (rms_float - rms2) / rms_float * 100 if rms_float > 0 else 0
            print(f"\n-- Phase 12.0 Results --")
            print(f"  Pass 1 (EKF float):     3D RMS = {rms_float:.3f} m")
            print(f"  Pass 2 (arc-fixed):     3D RMS = {rms2:.3f} m")
            print(f"  Improvement:             {imp:+.1f}%")
            print(f"  Arc WL: {n_wl} SVs, Arc NL: {n_nl} SVs "
                  f"(OSB: {'yes' if osb_nl else 'no'})")

            out_dir = Path('results') / 'arc_ar'
            out_dir.mkdir(parents=True, exist_ok=True)
            out_path = out_dir / f"arcar_{date_str}_{grace_id}_{args.hours}h.pkl"
            pickle.dump({
                'results_pass1': results, 'results_pass2': results2,
                'rms_pass1': rms_float, 'rms_pass2': rms2,
                'wl_fixed': wl_fixed, 'nl_fixed': nl_fixed,
                'arc_info': {sv: len(arcs) for sv, arcs in arc_dict.items()},
                'config': vars(args),
            }, open(str(out_path), 'wb'))
            print(f"\n  Saved: {out_path}")
            return

    # 鈹€鈹€ Framework v3: Single-pass batch solver on EKF orbit 鈹€鈹€
    if args.batch_lsq_v2 and _pass1_geometry:
        from src.batch_solver import BatchLinearSolver

        print(f"\n-- Framework v3 Batch Solver ({args.hours}h, "
              f"{len(_pass1_geometry)} epochs) --")

        bls = BatchLinearSolver(_pass1_geometry,
                                 sigma_phase=args.sigma_phase,
                                 sigma_code=args.sigma_code)
        bls_sol = bls.solve()

        ekf_phase_rms = results['stats'][-1]['rms_phase'] if results['stats'] else 0
        bch_phase_rms = bls_sol['rms_phase']
        bch_code_rms  = bls_sol['rms_code']

        print(f"\n-- Framework v3 Results --")
        print(f"  Geometry: EKF orbit (3D RMS = {rms_float:.3f}m, unchanged)")
        print(f"  Phase RMS:  EKF={ekf_phase_rms:.3f}m  Batch={bch_phase_rms:.3f}m")
        if ekf_phase_rms > 0:
            imp = (ekf_phase_rms - bch_phase_rms)/ekf_phase_rms*100
            print(f"  Phase improvement: {imp:+.1f}%")
        print(f"  Code RMS:   Batch={bch_code_rms:.3f}m")
        print(f"  Ambs: {len(bls_sol['amb_dict'])} SVs")

        out_d = Path('results') / 'batch_v3'
        out_d.mkdir(parents=True, exist_ok=True)
        pickle.dump({
            'rms_ekf_orbit': rms_float,
            'rms_phase_ekf': ekf_phase_rms,
            'rms_phase_batch': bch_phase_rms,
            'rms_code_batch': bch_code_rms,
            'epochs': epochs,
            'config': vars(args),
        }, open(str(out_d / f'batchv3_{date_str}_{grace_id}_{args.hours}h.pkl'), 'wb'))
        print(f"  Saved: {out_d / f'batchv3_{date_str}_{grace_id}_{args.hours}h.pkl'}")
        return
    if args.batch_lsq:
        from src.batch_solver import BatchLinearSolver
        from src.batch_lsq import BatchAmbiguityResolver, COEFF_W, LAM_NL
        from src.sequential_filter import I_AMB_START

        print(f"\n-- Phase 10.0 Batch LSQ ({args.hours}h, {len(_pass1_geometry)} epochs) --")

        # Step 1: Batch linear solver 鈫?jointly-solved amb (self-consistent)
        bls = BatchLinearSolver(_pass1_geometry,
                                 sigma_phase=args.sigma_phase,
                                 sigma_code=args.sigma_code)
        bls_sol = bls.solve()
        amb_bls = bls_sol['amb_dict']
        print(f"  [BS] {len(amb_bls)} SVs, phase_rms={bls_sol['rms_phase']:.3f}m")

        # Step 2: Batch MW 鈫?WL integers (full-arc smoothing)
        all_ep_data = []
        for i_ep, gps_sod in enumerate(epochs):
            r_gnv = interpolate_ref(ref_orbit, gps_sod)
            if r_gnv is None: continue
            ep_data = compute_epoch_geometry(gps_sod, gps1b_raw, sp3, r_gnv, clk_data)
            all_ep_data.extend(ep_data)

        ekf_amb = {}
        for sv in state.sv_list:
            idx = state.sv_to_idx.get(sv)
            if idx is not None: ekf_amb[sv] = float(state.x[I_AMB_START + idx])

        resolver = BatchAmbiguityResolver(all_ep_data, sv_to_idx_map=state.sv_to_idx,
                                          min_epochs_wl=args.ar_min_epochs)
        batch_result = resolver.resolve(ekf_amb)
        wl_batch = batch_result['wl_fixed']
        nl_batch = batch_result['nl_fixed']
        nl_ref   = batch_result.get('nl_ref')
        print(f"  [MW] {len(wl_batch)} WL-fixed, {len(nl_batch)} NL-fixed "
              f"(ref={nl_ref}, b_r_wl={batch_result['b_r_wl']:+.4f}cyc)")

        if len(wl_batch) < 4:
            print("  Too few WL-fixed SVs, skip")
        else:
            # Step 3: Extract SD NL from batch solver (more accurate than MW bootstrapping).
            # Clock-independent: clock cancels in between-satellite difference.
            nl_from_bls = {}
            if nl_ref and nl_ref in amb_bls and nl_ref in wl_batch:
                B_ref = float(amb_bls[nl_ref])
                N_w_ref = wl_batch[nl_ref]
                for sv in wl_batch:
                    if sv == nl_ref or sv not in amb_bls: continue
                    N_w_k = wl_batch[sv]
                    B_k = float(amb_bls[sv])
                    sd_val = B_k - B_ref - COEFF_W * (N_w_k - N_w_ref)
                    dN1 = int(round(sd_val / LAM_NL))
                    if abs(sd_val / LAM_NL - dN1) < 0.40:
                        nl_from_bls[sv] = dN1

            n_sd_bls = len(nl_from_bls)
            print(f"  [BS-SD] {n_sd_bls} SVs from batch solver "
                  f"(vs {len(nl_batch)} from MW bootstrapping)")

            # Step 4: Pass 2 鈥?fresh EKF with WL pre-fill + batch solver SD NL
            ekf_cfg2 = dict(ekf_config)
            ekf_cfg2['use_cycle_slip'] = False

            ekf2 = SequentialEKF(ekf_cfg2)
            ekf2._wl_fixed = dict(wl_batch)
            ekf2._wl_epochs = {sv: 0 for sv in wl_batch}
            ekf2._nl_batch = dict(nl_from_bls)
            ekf2._nl_batch_ref = nl_ref
            state2 = ekf2.initialize(r0_eci, v0_eci, mjd_start, epochs[0])

            results2 = {'epochs': [], 'r_ecef': [], 'v_ecef': [],
                        'a_rtn': [], 'zwd': [], 'clk': [], 'stats': [],
                        'n_sv': [], 'r_gnv': []}

            for i_ep, gps_sod in enumerate(epochs):
                mjd_utc = MJD_J2000 + (gps_sod - GPS_UTC_OFFSET) / SEC_PER_DAY
                mjd_tt = mjd_utc + 69.184 / SEC_PER_DAY
                if i_ep > 0:
                    mjd_up = MJD_J2000 + (epochs[i_ep-1] - GPS_UTC_OFFSET) / SEC_PER_DAY
                    state2 = ekf2.predict(state2, gps_sod, mjd_up, mjd_up + 69.184 / SEC_PER_DAY)
                rcv_e, _ = eci_to_ecef(state2.r_eci, state2.v_eci, mjd_utc)
                ep_d = compute_epoch_geometry(gps_sod, gps1b_raw, sp3, rcv_e, clk_data)
                if ep_d:
                    state2, s2 = ekf2.process_epoch(state2, ep_d, sp3, sv_bias, sv_bias_ref, mjd_utc, mjd_tt, doy)
                    r_e, v_e = eci_to_ecef(state2.r_eci, state2.v_eci, mjd_utc)
                    r_g = interpolate_ref(ref_orbit, gps_sod)
                    if r_g is not None:
                        results2['epochs'].append(gps_sod)
                        results2['r_ecef'].append(r_e); results2['v_ecef'].append(v_e)
                        results2['a_rtn'].append(state2.a_rtn.copy())
                        results2['zwd'].append(state2.zwd); results2['clk'].append(state2.clk)
                        results2['stats'].append(s2); results2['n_sv'].append(state2.n_sv)
                        results2['r_gnv'].append(r_g)
                if i_ep % max(1, len(epochs)//10) == 0 or i_ep < 3:
                    dr = np.linalg.norm(results2['r_ecef'][-1] - results2['r_gnv'][-1]) if results2['r_gnv'] else 0
                    print(f"  [Pass2 {i_ep:4d}/{len(epochs)}] t={gps_sod:.0f}  |dr|={dr:.3f}m  n_sv={state2.n_sv}")

            rms2 = np.sqrt(np.mean([np.linalg.norm(results2['r_ecef'][i] - results2['r_gnv'][i])**2
                                    for i in range(len(results2['epochs']))]))
            imp = (rms_float - rms2) / rms_float * 100 if rms_float > 0 else 0
            print(f"\n-- Phase 10.0 Results --")
            print(f"  Pass 1 (EKF float):    3D RMS = {rms_float:.3f} m")
            print(f"  Pass 2 (batch solver SD): 3D RMS = {rms2:.3f} m")
            print(f"  Improvement:            {imp:+.1f}%")
            print(f"  SD NL from batch solver: {n_sd_bls} SVs")

            out_dir = Path('results') / 'phase10'
            out_dir.mkdir(parents=True, exist_ok=True)
            pickle.dump({
                'results_pass1': results, 'results_pass2': results2,
                'rms_pass1': rms_float, 'rms_pass2': rms2,
                'batch_solver': {k: v for k, v in bls_sol.items() if k in ('rms_phase', 'rms_code')},
                'batch_result': {k: str(v) for k, v in batch_result.items() if k in ('wl_fixed', 'nl_fixed', 'nl_ref', 'b_r_wl')},
                'config': vars(args),
            }, open(str(out_dir / f"phase10_{date_str}_{grace_id}_{args.hours}h.pkl"), 'wb'))
            print(f"  Saved: {out_dir / f'phase10_{date_str}_{grace_id}_{args.hours}h.pkl'}")
            return

    # 鈹€鈹€ Batch AR (Phase 6.0) 鈹€鈹€
    if args.batch_ar:
        from src.batch_lsq import BatchAmbiguityResolver, COEFF_W, LAM_NL
        from src.batch_solver import BatchLinearSolver
        from src.sequential_filter import I_AMB_START
        print(f"\n-- Batch AR (Phase 6.0) --")

        # 鈹€鈹€ Batch linear solver: joint clock+zwd+amb 鈹€鈹€
        print(f"  Building batch linear solver ({len(_pass1_geometry)} epochs)...")
        bls = BatchLinearSolver(_pass1_geometry,
                                 sigma_phase=args.sigma_phase,
                                 sigma_code=args.sigma_code)
        bls_solution = bls.solve()
        amb_dict = bls_solution['amb_dict']  # {sv: B_if}, self-consistent

        print(f"  Batch solver: {len(amb_dict)} SVs, "
              f"rms_phase={bls_solution['rms_phase']:.3f}m, "
              f"rms_code={bls_solution['rms_code']:.3f}m")

        # 鈹€鈹€ Batch WL+NL resolution (full-arc MW smoothing) 鈹€鈹€
        all_ep_data = []
        for i_ep, gps_sod in enumerate(epochs):
            r_gnv = interpolate_ref(ref_orbit, gps_sod)
            if r_gnv is None:
                continue
            ep_data = compute_epoch_geometry(gps_sod, gps1b_raw, sp3, r_gnv, clk_data)
            all_ep_data.extend(ep_data)

        # Extract float ambiguity estimates from EKF final state
        ekf_amb = {}
        for sv in state.sv_list:
            idx = state.sv_to_idx.get(sv)
            if idx is not None:
                ekf_amb[sv] = float(state.x[I_AMB_START + idx])

        resolver = BatchAmbiguityResolver(all_ep_data,
                                          sv_to_idx_map=state.sv_to_idx,
                                          min_epochs_wl=args.ar_min_epochs)
        batch_result = resolver.resolve(ekf_amb)

        n_wl = len(batch_result['wl_fixed'])
        n_nl = len(batch_result['nl_fixed'])
        print(f"  Batch WL: {n_wl} WL-fixed, {n_nl} NL-fixed "
              f"(b_r_wl={batch_result['b_r_wl']:+.4f} cyc)")

        if n_wl >= 4:
            wl_batch = batch_result['wl_fixed']
            nl_batch = batch_result['nl_fixed']
            nl_ref = batch_result.get('nl_ref')

            # 鈹€鈹€ Absolute B_if: use first-WL-epoch B_if (clock-consistent) 鈹€鈹€
            # Per-SV first-WL-fix B_if has the same clock reference as the
            # epoch when the SV first joined the filter. This avoids the
            # ~0.15m clock drift at 0.5h that degrades final-state B_if.
            amb_for_ref = {}
            for sv in wl_batch:
                if sv in _sv_first_wl_bif:
                    amb_for_ref[sv] = _sv_first_wl_bif[sv]
                elif sv in ekf_amb:
                    amb_for_ref[sv] = ekf_amb[sv]
            # Also need reference SV's B_if (may or may not be in _sv_first_wl_bif)
            for sv in ekf_amb:
                if sv not in amb_for_ref:
                    amb_for_ref[sv] = ekf_amb[sv]

            amb_batch = resolver._build_absolute_if(wl_batch, nl_batch,
                                                      nl_ref, amb_for_ref)
            amb_var = max(0.0004, (args.hours * 0.10) ** 2)

            n_first = len(set(wl_batch.keys()) & set(_sv_first_wl_bif.keys()))
            print(f"  Pass 2: {len(amb_batch)} absolute B_if (P_amb={amb_var:.4f}, "
                  f"ref={nl_ref}, {n_first}/{n_wl} from first-WL-epoch)")
            for sv in sorted(amb_batch.keys())[:6]:
                w = wl_batch.get(sv, '?')
                d = nl_batch.get(sv, '?')
                tag = '*' if sv in _sv_first_wl_bif else ' '
                print(f"   {tag} {sv}: B_if={amb_batch[sv]:.3f} m (N_w={w}, dN1={d})")

            ekf_cfg2 = dict(ekf_config)
            ekf_cfg2['amb_batch_fixed'] = amb_batch
            ekf_cfg2['amb_batch_var'] = amb_var
            ekf_cfg2['use_cycle_slip'] = False

            ekf2 = SequentialEKF(ekf_cfg2)
            ekf2._wl_fixed = dict(wl_batch)
            ekf2._wl_epochs = {sv: 0 for sv in wl_batch}
            state2 = ekf2.initialize(r0_eci, v0_eci, mjd_start, epochs[0])

            # Re-process all epochs with locked ambiguities
            results2 = {
                'epochs': [], 'r_ecef': [], 'v_ecef': [],
                'a_rtn': [], 'zwd': [], 'clk': [],
                'stats': [], 'n_sv': [],
                'r_gnv': [],
            }

            for i_ep, gps_sod in enumerate(epochs):
                mjd_utc = MJD_J2000 + (gps_sod - GPS_UTC_OFFSET) / SEC_PER_DAY
                mjd_tt = mjd_utc + 69.184 / SEC_PER_DAY

                if i_ep > 0:
                    dt = gps_sod - epochs[i_ep - 1]
                    mjd_utc_prev = MJD_J2000 + (epochs[i_ep - 1] - GPS_UTC_OFFSET) / SEC_PER_DAY
                    mjd_tt_prev = mjd_utc_prev + 69.184 / SEC_PER_DAY
                    state2 = ekf2.predict(state2, gps_sod, mjd_utc_prev, mjd_tt_prev)

                rcv_pos_ecef2, _ = eci_to_ecef(state2.r_eci, state2.v_eci, mjd_utc)
                ep_data2 = compute_epoch_geometry(gps_sod, gps1b_raw, sp3, rcv_pos_ecef2, clk_data)

                if ep_data2:
                    state2, stats2 = ekf2.process_epoch(
                        state2, ep_data2, sp3, sv_bias, sv_bias_ref, mjd_utc, mjd_tt, doy)

                    r_ecef2, v_ecef2 = eci_to_ecef(state2.r_eci, state2.v_eci, mjd_utc)
                    r_gnv2 = interpolate_ref(ref_orbit, gps_sod)
                    results2['epochs'].append(gps_sod)
                    results2['r_ecef'].append(r_ecef2)
                    results2['v_ecef'].append(v_ecef2)
                    results2['a_rtn'].append(state2.a_rtn.copy())
                    results2['zwd'].append(state2.zwd)
                    results2['clk'].append(state2.clk)
                    results2['stats'].append(stats2)
                    results2['n_sv'].append(state2.n_sv)
                    results2['r_gnv'].append(r_gnv2)

                if i_ep % max(1, len(epochs) // 10) == 0 or i_ep < 3:
                    pos_diff2 = 0
                    if results2['r_gnv']:
                        pos_diff2 = np.linalg.norm(
                            results2['r_ecef'][-1] - results2['r_gnv'][-1])
                    print(f"  [Pass2 {i_ep:4d}/{len(epochs)}] t={gps_sod:.0f}  "
                          f"|dr|={pos_diff2:.3f}m  n_sv={state2.n_sv}")

            # Compute pass 2 RMS
            if results2['r_gnv']:
                pos_diffs2 = np.array([np.linalg.norm(results2['r_ecef'][i] - results2['r_gnv'][i])
                                        for i in range(len(results2['epochs']))])
                rms_3d2 = np.sqrt(np.mean(pos_diffs2**2))
                improvement = (rms_float - rms_3d2) / rms_float * 100 if rms_float > 0 else 0
                print(f"\n-- Batch AR Results --")
                print(f"  Pass 1 (float EKF): 3D RMS = {rms_float:.3f} m")
                print(f"  Pass 2 (warm-start fixed): 3D RMS = {rms_3d2:.3f} m")
                print(f"  Improvement:          {improvement:+.1f}%")

                # Save both passes
                out_dir = Path('results') / 'batch_ar'
                out_dir.mkdir(parents=True, exist_ok=True)
                out_path = out_dir / f"batchar_{date_str}_{grace_id}_{args.hours}h.pkl"
                pickle.dump({
                    'results_pass1': results,
                    'results_pass2': results2,
                    'rms_pass1': rms_float,
                    'rms_pass2': rms_3d2,
                    'config': vars(args),
                    'batch_result': {k: v for k, v in batch_result.items()
                                     if k in ('wl_fixed', 'nl_fixed', 'nl_ref', 'b_r_wl')},
                }, open(str(out_path), 'wb'))
                print(f"\n  Saved (both passes): {out_path}")
                return  # skip normal EKF save

    # -- Save results --
    out_dir = Path('results') / 'sequential_ekf'
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"seq_{date_str}_{grace_id}_{args.hours}h.pkl"
    pickle.dump({
        'results': results,
        'config': vars(args),
        'ekf_config': ekf_config,
    }, open(str(out_path), 'wb'))
    print(f"\n  Saved: {out_path}")


if __name__ == '__main__':
    main()
