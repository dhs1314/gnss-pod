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
                        help='ZWD random walk sigma [m/√s]')
    parser.add_argument('--sigma-phase', type=float, default=0.20,
                        help='Phase measurement sigma [m]')
    parser.add_argument('--sigma-code', type=float, default=0.30,
                        help='Code measurement sigma [m]')
    parser.add_argument('--chi2-threshold', type=float, default=100.0,
                        help='Chi-square innovation test threshold (alpha=0.001 → 10.828; '
                             'tuned: 100 for measurement-driven GGM05C)')
    # ── Phase 2.3: Measurement corrections ──
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
    # ── Dynamics mode ──
    parser.add_argument('--dynamics-mode', default='simplified',
                        choices=['simplified', 'orekit'],
                        help='Dynamics model: simplified (V2.2.1) or orekit')
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
            print(f"  [SP3] Cached → {cache_pkl.name}")
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
    # Use a short window (60 epochs ≈ 30 min) to keep receiver clock roughly constant.
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

    # -- Initialize EKF --
    ekf_config = {
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
        'chi2_threshold': args.chi2_threshold,
        'el_min': 0.087,
        # Phase 2.3 measurement corrections
        'use_phase_windup': args.enable_phase_windup,
        'use_relativity': args.enable_relativity,
        'use_cycle_slip': args.enable_cycle_slip,
        # Phase 3.0 PPP-AR
        'ar_min_epochs': args.ar_min_epochs,
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

    # Load IERS C04 Earth Orientation Parameters
    if args.iers_c04:
        from src.precision_products import setup_iers_from_c04
        print(f"  [IERS] Loading {args.iers_c04}...")
        setup_iers_from_c04(args.iers_c04)

    ekf = SequentialEKF(ekf_config)
    state = ekf.initialize(r0_eci, v0_eci, mjd_start, epochs[0])

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
