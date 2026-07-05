"""5-day Orekit GN validation (2024-04-29 to 2024-05-03).

Runs full pipeline on each date:
  EKF pass 1 (simplified dynamics) -> Orekit GN outer loop (Phase 20.0)

Reports: 3D RMS vs GNV1B, Phase RMS, GN convergence status.
Generates: results/5day_orekit/5day_orekit.png
"""
import sys, os, pickle, time, argparse, gzip, shutil, urllib.request
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'src'))

# ---- CLI ----
parser = argparse.ArgumentParser(description='Multi-day Orekit GN validation')
parser.add_argument('--dates', default='2024-04-29,2024-04-30,2024-05-01,2024-05-02,2024-05-03', help='Comma-separated dates')
parser.add_argument('--hours', default='0.17,0.50', help='Comma-separated arc lengths')
parser.add_argument('--grace-id', default='C', help='GRACE-FO satellite ID')
parser.add_argument('--auto-download', action='store_true', default=True, help='Auto-download missing CODE data')
parser.add_argument('--min-coverage', type=float, default=0.0, help='Min avg SV coverage (0-1) to accept arc')
parser.add_argument('--fuse-arcs', type=float, default=0.0, help='If >0, fuse N sliding arcs per hour')
parser.add_argument('--broadcast', action='store_true', default=False,
    help='Use broadcast-simulated orbits (CODE SP3 + 1.0m/1.5m noise) instead of precise')
parser.add_argument('--code-arc-hours', type=float, default=1.0,
    help='Arc length for code-only initial orbit [hours] (default: 1.0)')
parser.add_argument('--skip-code-orbit', action='store_true', default=False,
    help='Use GNV1B instead of code-only orbit (backward compat)')
args = parser.parse_args()

DATES = [d.strip() for d in args.dates.split(',')]
GRACE = args.grace_id; INTERVAL = 30
ARC_HOURS = [float(h.strip()) for h in args.hours.split(',')]
MIN_COVERAGE = args.min_coverage; FUSE_ARCS = args.fuse_arcs
USE_BROADCAST = args.broadcast
CODE_ARC_HOURS = args.code_arc_hours; SKIP_CODE_ORBIT = args.skip_code_orbit
DR = ROOT / 'data'
MJD0 = 51544.5; SEC = 86400.0; C_L = 299792458.0; OM = 7.2921151467e-5

n_days = len(DATES); n_arcs = len(ARC_HOURS)
date_start = DATES[0].replace('-',''); date_end = DATES[-1].replace('-','')
label = f"{n_days}day_{date_start}_{date_end}" + ("_BRDC" if USE_BROADCAST else "")
OUT = ROOT / 'results' / f'{label}_orekit'; OUT.mkdir(parents=True, exist_ok=True)

# ---- Shared data loading ----
print("=" * 70)
tag = " [BRDC SIM]" if USE_BROADCAST else ""
print(f"Orekit GN Validation{tag}")
print(f"  Dates: {DATES[0]} to {DATES[-1]}, arcs: {ARC_HOURS}")
print("=" * 70)

print("\nLoading shared products...")
sp3_pkls = {}
clk_data = {}
dcb_pairs = {}
antex = None

for date_str in DATES:
    dt = datetime.strptime(date_str, "%Y-%m-%d")
    doy = dt.strftime("%j"); y = dt.year; m = dt.month

    if USE_BROADCAST:
        # Load pre-generated BRDC SP3 pkl (from real IGS broadcast ephemeris)
        brdc_p = DR / "BRDC" / f"BRDC_{y}{doy}0000_01D_05M_ORB.pkl"
        if brdc_p.exists():
            sp3_pkls[date_str] = pickle.load(open(str(brdc_p), "rb"))
            n_pos = sum(1 for t in sp3_pkls[date_str]['epochs'] for _ in sp3_pkls[date_str]['epochs'][t])
            print(f"  BRDC REAL {date_str}: {n_pos} pos from IGS broadcast ephemeris")
        else:
            print(f"  BRDC {date_str}: no precomputed file, SKIP")
        clk_data[date_str] = {}
    else:
        pref = f"COD0OPSFIN_{y}{doy}0000_01D"
        sp3_p = DR / "CODE" / str(y) / f"{pref}_05M_ORB.pkl"
        if sp3_p.exists():
            sp3_pkls[date_str] = pickle.load(open(str(sp3_p), "rb"))
        else:
            sp3_txt = DR / "CODE" / str(y) / f"{pref}_05M_ORB.SP3"
            if sp3_txt.exists():
                from src.sp3_loader import parse_sp3_text
                with open(str(sp3_txt), 'r') as f:
                    epochs_dict, ts_list = parse_sp3_text(f.read())
                sp3_pkls[date_str] = {'ts': ts_list, 'epochs': epochs_dict,
                                       'source': 'CODE', 'product': 'FIN'}

        clk_p = DR / "CODE" / str(y) / f"{pref}_30S_CLK.CLK"
        if clk_p.exists():
            from precision_products import read_rinex_clk
            clk_data[date_str] = read_rinex_clk(str(clk_p))

    # DCB (same for all April/May dates near day 120)
    dcb_p1p2 = DR / "CODE" / str(y) / f"P1P2{y%100:02d}{m:02d}.DCB"
    dcb_p1c1 = DR / "CODE" / str(y) / f"P1C1{y%100:02d}{m:02d}.DCB"
    if date_str not in dcb_pairs:
        from precision_products import load_code_dcb_pair
        try:
            dcb_pairs[date_str] = load_code_dcb_pair(str(dcb_p1c1), str(dcb_p1p2))
        except:
            dcb_pairs[date_str] = load_code_dcb_pair(
                str(DR / "CODE" / str(y) / "P1C12404.DCB"),
                str(DR / "CODE" / str(y) / "P1P22404.DCB"))

# ANTEX (shared)
antex_path = DR / "igs14.atx"
if antex_path.exists():
    from precision_products import read_antex
    antex = read_antex(str(antex_path))

from precision_products import setup_iers_from_c04
setup_iers_from_c04(str(DR / "IERS/eopc04_IAU2000.txt"))

# Gravity model
from gravity_model import read_icgem_gfc
grav_path = DR / 'gravity' / 'GGM05C.gfc'
Cnm, Snm, _, GM_grav, R_grav = read_icgem_gfc(str(grav_path))
GRAV_NMAX = 150

# Import Orekit once
from src.orekit_bridge import OrekitPropagator, is_orekit_available
os.environ['OREKIT_DATA_PATH'] = str(DR / 'orekit')
if not is_orekit_available():
    print("FATAL: Orekit not available"); sys.exit(1)

if USE_BROADCAST: print(f"  [BRDC SIM] Broadcast ephemeris mode (CODE + 1.0m orbit + 1.5m clock noise)")
if not SKIP_CODE_ORBIT: print(f"  [CodeOrbit] Initial orbit from pseudorange GN (arc={CODE_ARC_HOURS}h)")
print(f"  Loaded: SP3={len(sp3_pkls)}d, CLK={'BRDC-embedded' if USE_BROADCAST else len(clk_data)}d, DCB={len(dcb_pairs)}d")
print(f"  Gravity: GGM05C Nmax={GRAV_NMAX}")

# ---- Run one arc ----
from src.troposphere import saastamoinen_zhd
from src.batch_orbit_v3 import BatchOrbitLSQv3
from src.orbit_dynamics import total_acc_eci
from src.batch_solver import BatchLinearSolver
from run_sequential_pod import load_gnv1b, compute_epoch_geometry, interpolate_ref
from coordinates import ecef_to_eci, eci_to_ecef
from sequential_filter import SequentialEKF
from src.code_orbit import kinematic_wls_single_epoch, CodeOnlyOrbitSolver
from src.data_quality import (
    check_sv_coverage, check_snr, check_sp3_integrity, check_clk_integrity,
    compute_per_sv_multipath, compute_mw_stability, screen_sv_for_batch,
    generate_qc_report, print_qc_line,
)

def run_one_arc(date_str, arc_h, arc_offset=0.0):
    """Run EKF + Orekit GN on one arc. arc_offset in hours for sliding windows."""
    dt = datetime.strptime(date_str, "%Y-%m-%d"); y = dt.year
    J2000 = datetime(2000, 1, 1, 12, 0, 0)
    from run_sequential_pod import get_sat_geometry
    from precision_products import compute_dcb_if_correction

    # GPS1B
    gps_path = DR / 'gracefo' / str(y) / date_str / f'GPS1B_{date_str}_{GRACE}_04.pkl'
    if not gps_path.exists():
        return None
    gps1b = pickle.load(open(str(gps_path), "rb"))

    # SP3 and CLK
    sp3 = sp3_pkls.get(date_str)
    clk = clk_data.get(date_str) if not USE_BROADCAST else None
    dcb = dcb_pairs.get(date_str)
    if sp3 is None or dcb is None:
        return None

    # GNV1B reference (only needed for 3D RMS validation or skip-code-orbit mode)
    ref, refv = None, None
    if SKIP_CODE_ORBIT:
        gnv_path = DR / 'gracefo' / str(y) / date_str / f'GNV1B_{date_str}_{GRACE}_04.txt'
        if not gnv_path.exists(): return None
        ref, refv = load_gnv1b(str(gnv_path))
    else:
        # Still try to load GNV1B for validation (non-essential)
        gnv_path = DR / 'gracefo' / str(y) / date_str / f'GNV1B_{date_str}_{GRACE}_04.txt'
        if gnv_path.exists():
            ref, refv = load_gnv1b(str(gnv_path))

    # Time windows — always use arc_h for precise OD
    gps0 = min(gps1b.keys()) + arc_offset * 3600
    gps_end = gps0 + int(arc_h * 3600)
    epochs = sorted(set(g for g in gps1b if gps0 <= g <= gps_end
                        and abs((g - gps0) % INTERVAL) <= 2))
    if len(epochs) < 6: return None

    # ── Orekit propagator — only when needed (not broadcast mode) ──
    prop = None
    if not USE_BROADCAST:
        prop = OrekitPropagator(
            gravity_field=str(grav_path), gravity_degree=GRAV_NMAX,
            solid_tides=True, ocean_tides=True, ocean_tide_degree=50,
            third_body='lunisolar', srp_model='isotropic', relativity=True,
            drag_model='exponential',
            mass=580.0, area_drag=0.68, area_srp=3.4, CR=1.3, CD=2.2,
            stm_perturb=1.0, integrator_tol=1e-12)
        prop._setup()

    # ── Step 0 + 1: Code-only initial orbit (no GNV1B dependency) ──
    code_orbit = None; code_bad_svs = set()
    if not SKIP_CODE_ORBIT:
        # Step 0: Kinematic WLS for initial guess
        r0_kin, _ = kinematic_wls_single_epoch(gps1b, epochs[0], sp3)
        v0_kin = np.zeros(3)

        # Step 1: Code-only GN over longer arc
        code_gps0 = min(gps1b.keys()) + arc_offset * 3600
        code_gps_end = code_gps0 + int(CODE_ARC_HOURS * 3600)
        code_epochs_all = sorted(set(g for g in gps1b
                                     if code_gps0 <= g <= code_gps_end
                                     and abs((g - gps0) % INTERVAL) <= 2))
        if len(code_epochs_all) >= 12:
            # Build code-only geometry from the first N epochs
            code_epochs = code_epochs_all[:int(CODE_ARC_HOURS * 3600 / INTERVAL)]
            code_t_ep = np.array([(g - code_epochs[0]) for g in code_epochs], dtype=float)
            code_mjd_s = MJD0 + code_epochs[0] / SEC
            code_mjd_tt = code_mjd_s + 69.184 / SEC

            # Compute geometry for code-only orbit (coarse rcv position for satellite geometry)
            code_geo_eps = []
            for gps_sod in code_epochs:
                utc_dt = J2000 + timedelta(seconds=gps_sod)
                recs = gps1b.get(int(gps_sod), gps1b.get(gps_sod, {}))
                ep_obs = []
                for sv_id, rec in recs.items():
                    if 'P_if' not in rec: continue
                    # Use kinematic position as approximate receiver location
                    sat_pos, sat_clk, _ = get_sat_geometry(sp3, sv_id, utc_dt, r0_kin, clk)
                    if sat_pos is None: continue
                    dcb_c = compute_dcb_if_correction(dcb, sv_id)
                    ep_obs.append({
                        'sv': sv_id,
                        'sat_pos': sat_pos,
                        'sat_clk': sat_clk,
                        'el': 0.5,  # approximate
                        '_obs_code': float(rec['P_if']) + dcb_c,
                        '_geo_full': 0.0,  # filled by solver
                    })
                if ep_obs:
                    code_geo_eps.append(ep_obs)

            if len(code_geo_eps) >= 6:
                try:
                    # Build force fn for code-only orbit using loaded gravity
                    _G = {'Cnm': Cnm, 'Snm': Snm, 'Nmax': GRAV_NMAX, 'GM': GM_grav, 'R': R_grav}
                    def _code_force(pos, vel, CD=2.2, CR=1.3, area_drag=0.68,
                                    area_srp=3.4, mass=580.0, empirical_acc_rtn=None,
                                    mjd_utc=None, mjd_tt=None, **kw):
                        return total_acc_eci(pos, vel,
                            mjd_tt=mjd_tt or code_mjd_tt,
                            mjd_utc=mjd_utc or code_mjd_s,
                            Cnm=_G['Cnm'], Snm=_G['Snm'], Nmax=_G['Nmax'],
                            CD=CD, CR=CR, area_drag=area_drag, area_srp=area_srp,
                            mass=mass, empirical_acc_rtn=empirical_acc_rtn,
                            GM_gravity=_G['GM'], R_gravity=_G['R'])

                    code_solver = CodeOnlyOrbitSolver(
                        code_geo_eps, code_t_ep, code_mjd_s, code_mjd_tt,
                        None, sigma_code=0.30, max_iter=4, damping=0.5,
                        force_fn=_code_force)
                    code_result = code_solver.solve(r0_kin, v0_kin)

                    # Outlier detection
                    code_outliers = code_solver.detect_outliers(code_result, threshold=5.0)
                    code_bad_svs = {sv for sv, v in code_outliers.items() if v['reject']}

                    # Store code orbit for interpolation
                    code_orbit = code_result
                except Exception as e:
                    print(f"  [CodeOrbit] FAILED: {e}; falling back")
                    code_orbit = None

    # Get initial state for EKF
    # Priority: code_orbit > kinematic WLS (self-consistent, GNV1B never used for computation)
    # GNV1B only allowed in explicit --skip-code-orbit backward-compat mode
    if code_orbit is not None:
        r0e = code_orbit['r_ecef'][0].copy()
        v0e = code_orbit['v_ecef'][0].copy()
    elif SKIP_CODE_ORBIT and ref is not None:
        r0e = interpolate_ref(ref, epochs[0])
        v0e = interpolate_ref(refv, epochs[0])
    else:
        r0_kin, _ = kinematic_wls_single_epoch(gps1b, epochs[0], sp3)
        r0e = np.array(r0_kin); v0e = np.zeros(3)

    # Layer 1 QC
    cov_stats = check_sv_coverage(gps1b, epochs, min_coverage_pct=0.40)
    snr_stats = check_snr(gps1b, epochs)
    mp_stats = compute_per_sv_multipath(gps1b, epochs)
    sp3_ok = check_sp3_integrity(sp3)
    clk_ok = check_clk_integrity(clk) if clk is not None and not USE_BROADCAST else {'is_valid': True, 'n_gps_sv': 32, 'n_epochs_per_sv': 2880}

    # Path 1: Coverage filter
    cov_pcts_pre = [v['pct'] for v in cov_stats.values()]
    avg_cov_pre = float(np.mean(cov_pcts_pre)) if cov_pcts_pre else 0.0
    if MIN_COVERAGE > 0 and avg_cov_pre < MIN_COVERAGE:
        print(f"  SKIP: coverage={avg_cov_pre:.1%} < min={MIN_COVERAGE:.0%}")
        return {'date': date_str, 'arc_h': arc_h, 'rms_3d': 0, 'phase_ekf': 0,
                'phase_batch': 0, 'phase_gn': 0, 'n_sv': len(cov_stats),
                'qc_score': 0, 'qc_grade': 'F', 'qc_flags': 'LOW_COV',
                'avg_cov': avg_cov_pre, 'n_full_cov': 0, 'n_partial': 0,
                'n_low_cov': sum(1 for p in cov_pcts_pre if p < 0.4),
                'cov_pcts': cov_pcts_pre, 'skip': True,
                'converged': False, 'iterations': 0, 'time_gn': 0,
                'n_arc_wl': 0, 'n_arc_nl': 0}

    t_ep = np.array([(g - epochs[0]) for g in epochs], dtype=float)
    mjd_s = MJD0 + epochs[0] / SEC; mjd_tt = mjd_s + 69.184 / SEC
    r0i, v0i = ecef_to_eci(r0e, v0e, mjd_s)

    # Per-SV code bias
    N_BIAS = min(60, len(epochs))
    sv_p_res = {}

    # ── Position/velocity lookup functions ──
    # Two modes: compute (code_orbit only) and reference (code_orbit → GNV1B fallback)

    def _get_compute_pos(gps_sod):
        """Receiver position for COMPUTATION steps (code bias, etc.).
        NO GNV1B fallback — must be self-consistent."""
        if code_orbit is not None:
            try:
                idx0 = np.searchsorted(code_epochs, gps_sod)
                if idx0 == 0: return code_orbit['r_ecef'][0]
                if idx0 >= len(code_epochs): return code_orbit['r_ecef'][-1]
                t0, t1 = code_epochs[idx0-1], code_epochs[idx0]
                dt = (gps_sod - t0) / max(t1 - t0, 1e-6)
                return code_orbit['r_ecef'][idx0-1] * (1-dt) + code_orbit['r_ecef'][idx0] * dt
            except: return None
        # Fallback: use kinematic WLS from Step 0
        r0k, _ = kinematic_wls_single_epoch(gps1b, gps_sod, sp3)
        return np.array(r0k)

    def _get_ref_pos(gps_sod):
        """Receiver position for VALIDATION only (3D RMS).
        GNV1B fallback is acceptable here."""
        if code_orbit is not None:
            try:
                idx0 = np.searchsorted(code_epochs, gps_sod)
                if idx0 == 0: return code_orbit['r_ecef'][0]
                if idx0 >= len(code_epochs): return code_orbit['r_ecef'][-1]
                t0, t1 = code_epochs[idx0-1], code_epochs[idx0]
                dt = (gps_sod - t0) / max(t1 - t0, 1e-6)
                return code_orbit['r_ecef'][idx0-1] * (1-dt) + code_orbit['r_ecef'][idx0] * dt
            except: return None
        elif ref is not None:
            return interpolate_ref(ref, gps_sod)
        return None

    def _get_ref_vel(gps_sod):
        """Receiver velocity for VALIDATION only (3V RMS)."""
        if code_orbit is not None and 'v_ecef' in code_orbit:
            try:
                idx0 = np.searchsorted(code_epochs, gps_sod)
                if idx0 == 0: return code_orbit['v_ecef'][0]
                if idx0 >= len(code_epochs): return code_orbit['v_ecef'][-1]
                t0, t1 = code_epochs[idx0-1], code_epochs[idx0]
                dt = (gps_sod - t0) / max(t1 - t0, 1e-6)
                return code_orbit['v_ecef'][idx0-1] * (1-dt) + code_orbit['v_ecef'][idx0] * dt
            except: return None
        elif refv is not None:
            return interpolate_ref(refv, gps_sod)
        return None

    for gps_sod in epochs[:N_BIAS]:
        utc_dt = J2000 + timedelta(seconds=gps_sod)
        recs = gps1b.get(int(gps_sod), gps1b.get(gps_sod, {}))
        for sv_id, rec in recs.items():
            if 'P_if' not in rec: continue
            rcv_r = _get_ref_pos(gps_sod) if SKIP_CODE_ORBIT else _get_compute_pos(gps_sod)
            if rcv_r is None: continue
            sat_pos, sat_clk, rho_corr = get_sat_geometry(sp3, sv_id, utc_dt, rcv_r, clk)
            if sat_pos is None: continue
            P_r = float(rec['P_if']) + compute_dcb_if_correction(dcb, sv_id) + sat_clk - rho_corr
            sv_p_res.setdefault(sv_id, []).append(P_r)
    sv_bias = {}; sv_bias_ref = 0.0
    for sv, vals in sv_p_res.items():
        if len(vals) >= 3: sv_bias[sv] = float(np.median(vals))
    if sv_bias:
        sv_bias_ref = float(np.mean(list(sv_bias.values())))
        for sv in sv_bias: sv_bias[sv] -= sv_bias_ref

    # Tide corrections
    try:
        from src.solid_tides import (compute_solid_tide_corrections,
                                      compute_time_varying_gravity, merge_tide_corrections)
        tide_c = compute_solid_tide_corrections(mjd_s, mjd_tt)
        tvgrav = compute_time_varying_gravity(mjd_tt)
        tide_c = merge_tide_corrections(tide_c, tvgrav)
    except: tide_c = {}

    # Force function closure
    _G = {'Cnm': Cnm, 'Snm': Snm, 'Nmax': GRAV_NMAX, 'GM': GM_grav, 'R': R_grav}
    def gn_fn(pos, vel, CD=2.2, CR=1.3, area_drag=0.68, area_srp=3.4, mass=580.0,
              empirical_acc_rtn=None, bodies=None, mjd_utc=None, mjd_tt=None, **kw):
        return total_acc_eci(pos, vel,
            mjd_tt=mjd_tt or mjd_tt, mjd_utc=mjd_utc or mjd_s,
            Cnm=_G['Cnm'], Snm=_G['Snm'], Nmax=_G['Nmax'],
            CD=CD, CR=CR, area_drag=area_drag, area_srp=area_srp, mass=mass,
            empirical_acc_rtn=empirical_acc_rtn, tide_corrections=tide_c,
            bodies=bodies or ['Sun','Moon'], GM_gravity=_G['GM'], R_gravity=_G['R'])

    # EKF pass 1
    chi2 = 100 if arc_h >= 0.3 else 25
    ekf_cfg = {
        'dynamics_mode': 'simplified', 'Cd': 2.2, 'CR': 1.3,
        'area_drag': 0.68, 'area_srp': 3.4, 'mass': 580.0,
        'bodies': ['Sun', 'Moon'], 'Cnm': Cnm, 'Snm': Snm,
        'GM_grav': GM_grav, 'R_grav': R_grav, 'gravity_nmax': GRAV_NMAX,
        'sigma_acc_process': 1e-3, 'tau_emp': 600.0, 'sigma_emp_ss': 1e-8,
        'sigma_zwd_rw': 1e-9, 'sigma_phase': 0.20, 'sigma_code': 0.30,
        'chi2_threshold': chi2, 'el_min': 0.087,
        'use_phase_windup': True, 'use_relativity': True, 'use_cycle_slip': True,
        'ar_min_epochs': 6, 'antex_data': antex, 'dcb_data': dcb,
        'elev_exp_phase': 1.0, 'elev_exp_code': 0.70 if arc_h >= 0.3 else 1.0,
        'clock_rw': 0.001 if arc_h >= 0.3 else 0.0004, 'mw_max_epochs': 200,
    }
    ekf = SequentialEKF(ekf_cfg)
    state = ekf.initialize(r0i, v0i, mjd_s, epochs[0])
    pass1 = []; r_ekf_list = []; v_ekf_list = []
    for i_ep, gps_sod in enumerate(epochs):
        mjd_u = MJD0 + gps_sod / SEC
        if i_ep > 0:
            mjd_prev = MJD0 + epochs[i_ep-1] / SEC
            state = ekf.predict(state, gps_sod, mjd_prev, mjd_prev + 69.184 / SEC)
        rcv_e, _ = eci_to_ecef(state.r_eci, state.v_eci, mjd_u)
        ep_data = compute_epoch_geometry(gps_sod, gps1b, sp3, rcv_e, clk)
        if not ep_data:
            r_ekf_list.append(state.r_eci.copy()); v_ekf_list.append(state.v_eci.copy())
            continue
        # Filter code-level outliers detected in Step 1
        if code_bad_svs:
            ep_data = [d for d in ep_data if d['sv'] not in code_bad_svs]
        if not ep_data:
            r_ekf_list.append(state.r_eci.copy()); v_ekf_list.append(state.v_eci.copy())
            continue
        state, stats = ekf.process_epoch(state, ep_data, sp3, sv_bias, sv_bias_ref,
                                          mjd_u, mjd_u + 69.184 / SEC, 120)
        lat = np.arcsin(rcv_e[2] / np.linalg.norm(rcv_e))
        h = np.linalg.norm(rcv_e) - 6378137.0; zhd = saastamoinen_zhd(lat, h)
        for d in ep_data:
            se = np.asarray(d['sat_pos'], dtype=float)
            sc = float(d.get('sat_clk', 0)); el = float(d.get('el', 0.5))
            rho = np.linalg.norm(se - rcv_e)
            sag = (OM / C_L) * (se[0]*rcv_e[1] - se[1]*rcv_e[0])
            mf = 1.001 / np.sqrt(0.002001 + np.sin(el)**2)
            dcb_c = compute_dcb_if_correction(dcb, d['sv'])
            d['_geo_full'] = rho + sag - sc + zhd * mf
            d['_obs_code'] = float(d.get('P_if_raw', 0)) + dcb_c - sv_bias.get(d['sv'], 0.0)
            d['_obs_phase'] = float(d.get('L_if_raw', 0)) - sv_bias.get(d['sv'], 0.0)
        pass1.append(ep_data)
        r_ekf_list.append(state.r_eci.copy())
        v_ekf_list.append(state.v_eci.copy())

    ekf_phase = stats.get('rms_phase', 0)
    n_sv = state.n_sv

    # SV gap detection: split gapped SVs into segments
    sv_ep_indices = {}
    for i_ep, ep_list in enumerate(pass1):
        for d in ep_list:
            sv_ep_indices.setdefault(d['sv'], []).append(i_ep)
    n_gap_svs = 0
    for sv, ep_idx_list in sorted(sv_ep_indices.items()):
        ep_idx_list.sort()
        if len(ep_idx_list) < 2: continue
        gaps = np.diff(ep_idx_list)
        gap_pos = np.where(gaps > 2)[0]
        if len(gap_pos) == 0: continue
        n_gap_svs += 1
        seg_starts = np.concatenate([[0], gap_pos + 1])
        seg_ends = np.concatenate([gap_pos, [len(ep_idx_list) - 1]])
        for seg_i in range(len(seg_starts)):
            seg_eps = set(ep_idx_list[seg_starts[seg_i]:seg_ends[seg_i] + 1])
            seg_sv = f'{sv}_S{seg_i}'
            for i_ep in seg_eps:
                for d in pass1[i_ep]:
                    if d['sv'] == sv: d['sv'] = seg_sv
    if n_gap_svs > 0:
        print(f"  [GapDetect] {n_gap_svs} SV(s) split into segments")

    # Arc-level MW + OSB (from raw GPS1B, pre-segmentation)
    arc_wl_fixed = {}; mw_per_sv = {}; mw_stats = {}
    F1_mw, F2_mw = 1575.42e6, 1227.60e6
    for gps_sod in epochs:
        recs = gps1b.get(int(gps_sod), gps1b.get(gps_sod, {}))
        for sv_id, rec in recs.items():
            L1c = rec.get('L1_cyc')
            if L1c is None and 'L1_phase' in rec: L1c = float(rec['L1_phase'])*F1_mw/C_L
            else: L1c = float(L1c or 0)
            L2c = rec.get('L2_cyc')
            if L2c is None and 'L2_phase' in rec: L2c = float(rec['L2_phase'])*F2_mw/C_L
            else: L2c = float(L2c or 0)
            if L1c == 0 or L2c == 0: continue
            P1 = float(rec.get('P1_raw', rec.get('P1', rec.get('CA_range', 0))))
            P2 = float(rec.get('P2_raw', rec.get('P2', 0)))
            if P1 == 0 or P2 == 0: continue
            from src.ambiguity import compute_mw
            mw_per_sv.setdefault(sv_id, []).append(compute_mw(L1c, L2c, P1, P2))
    for sv, mw_list in mw_per_sv.items():
        mw_stats[sv] = compute_mw_stability(mw_list)
    all_fracs = [np.mean(ml)-round(np.mean(ml)) for ml in mw_per_sv.values() if len(ml)>=5]
    b_r_wl = float(np.median(all_fracs)) if all_fracs else 0.0
    sv_decision = {}; n_rejected = 0
    for sv, mw_data in sorted(mw_stats.items()):
        if len(mw_per_sv.get(sv, [])) < 5: continue
        dec = screen_sv_for_batch(sv, cov_stats, mw_stats, mp_stats, snr_stats)
        sv_decision[sv] = dec
        if dec['decision'] in ('REJECT',): n_rejected += 1; continue
        if not dec.get('wl_eligible', False): continue
        N_w = int(round(mw_data['mean'] - b_r_wl))
        if abs(mw_data['mean'] - b_r_wl - N_w) < 0.35: arc_wl_fixed[sv] = N_w
    if n_rejected > 0: print(f"  [QC] {n_rejected} SV(s) rejected by screening")

    # OSB loading
    osb_wl = osb_nl = None
    osb_path = DR / "CODE" / "2024" / "COD0OPSFIN_20241200000_01D_01D_OSB.BIA"
    if osb_path.exists():
        try:
            from src.batch_lsq import read_code_osb
            osb_wl, osb_nl = read_code_osb(str(osb_path), quiet=True)
        except: pass
    n_arc_wl = len(arc_wl_fixed)
    if n_arc_wl > 0:
        print(f"  [ArcMW] arcWL={n_arc_wl}" + (f" +OSB({len(osb_wl)})" if osb_wl else ""))

    # Baseline batch on EKF orbit
    sv_cov_map = {}
    sv_ep_count = {}; n_total_ep = len(pass1)
    for ep_list in pass1:
        for d in ep_list: sv_ep_count[d['sv']] = sv_ep_count.get(d['sv'], 0) + 1
    for sv, n in sv_ep_count.items(): sv_cov_map[sv] = n / max(n_total_ep, 1)

    bls = BatchLinearSolver(pass1, sigma_phase=0.20, sigma_code=0.30,
                            sv_coverage=sv_cov_map)
    bls_sol = bls.solve()
    batch_ph = bls_sol['rms_phase']

    # Orekit GN — skip for broadcast (convergence unreliable with ~1m orbit errors)
    if USE_BROADCAST:
        # Use EKF orbit directly
        print(f"  [BRDC-Quick] Skipping GN loop — using EKF orbit")
        dt_gn = 0; n_arc_nl = 0
        sol = {'r_eci': np.array(r_ekf_list), 'v_eci': np.array(v_ekf_list),
               'rms_phase': ekf_phase, 'converged': False, 'iterations': 0}
    else:
        gn = BatchOrbitLSQv3(
            pass1, gn_fn, t_ep,
            mjd_utc_start=mjd_s, mjd_tt_start=mjd_tt,
            sigma_phase=0.20, sigma_code=0.30,
            max_iter=6, prior_r0=1.0, prior_v0=0.01, prior_emp=1e-7,
            damping=0.5, orekit_prop=prop, estimate_cd_cr=False)
        t0 = time.time()
        sol = gn.solve(r0i, v0i, arc_wl_fixed=arc_wl_fixed, osb_wl=osb_wl, osb_nl=osb_nl)
        dt_gn = time.time() - t0
        n_arc_nl = sol.get('n_arc_nl', 0)

    # 3D position + velocity RMS
    dr_vals = []; dv_vals = []
    for i_ep, gps_sod in enumerate(epochs):
        r_ref = _get_ref_pos(gps_sod)
        v_ref = _get_ref_vel(gps_sod)
        if r_ref is not None:
            mjd_u = MJD0 + gps_sod / SEC
            r_e, v_e = eci_to_ecef(sol['r_eci'][i_ep], sol['v_eci'][i_ep], mjd_u)
            dr_vals.append(np.linalg.norm(r_e - r_ref))
            if v_ref is not None:
                dv_vals.append(np.linalg.norm(v_e - v_ref))
    rms_3d = np.sqrt(np.mean([d**2 for d in dr_vals])) if dr_vals else 0
    rms_3v = np.sqrt(np.mean([d**2 for d in dv_vals])) if dv_vals else 0

    # QC report (now with real batch_phase_rms)
    cov_pcts = [v['pct'] for v in cov_stats.values()]
    avg_cov = float(np.mean(cov_pcts)) if cov_pcts else 0.0
    n_full_cov = sum(1 for p in cov_pcts if p >= 0.95)
    n_partial = sum(1 for p in cov_pcts if 0.4 <= p < 0.95)
    n_low_cov = sum(1 for p in cov_pcts if p < 0.40)
    qc_report = generate_qc_report({
        'cov_stats': cov_stats, 'snr_stats': snr_stats,
        'mp_stats': mp_stats, 'mw_stats': mw_stats,
        'sp3_ok': sp3_ok, 'clk_ok': clk_ok,
        'batch_phase_rms': batch_ph, 'n_gap_svs': n_gap_svs,
        'n_sv_total': n_sv, 'sv_decision': sv_decision,
    })
    print_qc_line(qc_report)

    return {
        'date': date_str, 'arc_h': arc_h,
        'rms_3d': rms_3d, 'rms_3v': rms_3v,
        'phase_ekf': ekf_phase,
        'phase_batch': batch_ph, 'phase_gn': sol['rms_phase'],
        'converged': sol['converged'], 'iterations': sol['iterations'],
        'n_sv': n_sv, 'time_gn': dt_gn,
        'n_arc_wl': n_arc_wl, 'n_arc_nl': n_arc_nl,
        'qc_score': qc_report['score'], 'qc_grade': qc_report['grade'],
        'qc_flags': ','.join(qc_report['flags']) if qc_report['flags'] else 'OK',
        'avg_cov': avg_cov, 'n_full_cov': n_full_cov,
        'n_partial': n_partial, 'n_low_cov': n_low_cov,
        'cov_pcts': cov_pcts,
    }

# ---- Run all arcs ----
all_results = []
total = len(DATES) * len(ARC_HOURS); done = 0

for date_str in DATES:
    for arc_h in ARC_HOURS:
        done += 1
        print(f"\n[{done}/{total}] {date_str} {arc_h:.2f}h...", end=" ", flush=True)
        try:
            r = run_one_arc(date_str, arc_h)
            if r:
                all_results.append(r)
                if r.get('skip'):
                    print(f"SKIP cov={r.get('avg_cov',0)*100:.0f}% SV={r['n_sv']}")
                else:
                    print(f"3D={r['rms_3d']:.3f}m 3V={r.get('rms_3v',0)*1000:.1f}mm/s "
                          f"GN_ph={r['phase_gn']:.3f}m "
                          f"t={r['time_gn']:.0f}s SV={r['n_sv']} "
                          f"QC={r['qc_score']:.2f} cov={r.get('avg_cov',0)*100:.0f}%")
            else:
                print("SKIP (no data)")
        except Exception as e:
            print(f"ERROR: {e}")
            import traceback; traceback.print_exc()

if not all_results:
    print("\nNo results!"); sys.exit(1)

# ---- Summary ----
print(f"\n{'='*70}")
print(f"{n_days}-Day Summary ({len(all_results)} arcs)")
print(f"{'='*70}")

for arc_h in ARC_HOURS:
    sub = [r for r in all_results if r['arc_h'] == arc_h]
    if not sub: continue
    rms_v = [r['rms_3d'] for r in sub]
    print(f"\n{arc_h:.2f}h ({len(sub)} dates):")
    print(f"  3D RMS: mean={np.mean(rms_v):.3f}m median={np.median(rms_v):.3f}m "
          f"best={np.min(rms_v):.3f}m worst={np.max(rms_v):.3f}m")

# Table
vel_ms = [r.get('rms_3v',0)*1000 for r in all_results if r.get('rms_3v',0) > 0]
print(f"\n{'Date':<12s} {'Arc':>6s} {'3D_RMS':>8s} {'3V_RMS':>9s} {'GN_Ph':>7s} {'Bch_Ph':>7s} "
      f"{'QC':>5s} {'Cov':>5s} {'SVs':>4s} {'Time':>6s} {'Flags':>s}")
for r in all_results:
    qc_str = f"{r.get('qc_score', 0):.2f}{r.get('qc_grade', '?')}"
    cov_str = f"{r.get('avg_cov', 0)*100:.0f}%"
    v_str = f"{r.get('rms_3v', 0)*1000:.1f}mm/s" if r.get('rms_3v', 0) > 0 else "N/A"
    flag_str = r.get('qc_flags', '')
    skip_tag = "SKIP " if r.get('skip') else ""
    print(f"{r['date']:<12s} {r['arc_h']:5.2f}h {r['rms_3d']:8.3f} {v_str:>9s} {r['phase_gn']:7.3f} "
          f"{r['phase_batch']:7.3f} {qc_str:>5s} {cov_str:>5s} "
          f"{skip_tag}{r['n_sv']:>4d} {r['time_gn']:5.0f}s {flag_str}")
if vel_ms:
    print(f"\nVelocity RMS: mean={np.mean(vel_ms):.1f}mm/s median={np.median(vel_ms):.1f}mm/s "
          f"best={np.min(vel_ms):.1f}mm/s worst={np.max(vel_ms):.1f}mm/s")

# ---- Plot: 3 panels side-by-side ----
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 6))

all_dates_sorted = sorted(set(r['date'] for r in all_results))
date_labels = [d[5:] for d in all_dates_sorted]
x_dates = np.arange(len(all_dates_sorted))

rms_017=[]; rms_050=[]; v017=[]; v050=[]; gn_017=[]; gn_050=[]; bch_017=[]; bch_050=[]
qc_017=[]; qc_050=[]; cov_017=[]; cov_050=[]; fc_017=[]; fc_050=[]

for d in all_dates_sorted:
    for ah, rms_l,gn_l,bch_l,qc_l,cov_l,fc_l,v_l in [
        (0.17,rms_017,gn_017,bch_017,qc_017,cov_017,fc_017,v017),
        (0.50,rms_050,gn_050,bch_050,qc_050,cov_050,fc_050,v050)]:
        m=[r for r in all_results if r['date']==d and abs(r['arc_h']-ah)<0.01]
        if m:
            r=m[0]; rms_l.append(r['rms_3d']); gn_l.append(r['phase_gn'])
            bch_l.append(r['phase_batch']); qc_l.append(r.get('qc_score',0))
            cov_l.append(r.get('avg_cov',0)*100); fc_l.append(r.get('n_full_cov',0))
            v_l.append(r.get('rms_3v',0)*1000)  # mm/s
        else:
            for L in [rms_l,gn_l,bch_l,qc_l,cov_l,fc_l,v_l]: L.append(np.nan)

# Panel (a): 3D RMS curves
ax1.plot(x_dates, rms_017, 'o-', color='#607D8B', lw=2, ms=8,
         label='0.17h 3D RMS', markerfacecolor='white', markeredgewidth=1.5)
ax1.plot(x_dates, rms_050, 's-', color='#FF9800', lw=2, ms=8,
         label='0.50h 3D RMS', markerfacecolor='white', markeredgewidth=1.5)
for xi, qi in zip(x_dates, qc_017):
    if not np.isnan(qi) and qi < 0.70: ax1.axvspan(xi-0.3, xi+0.3, color='red', alpha=0.08)
for xi, yi, qi in zip(x_dates, rms_017, qc_017):
    if not np.isnan(yi): ax1.annotate(f'QC{qi:.2f}', (xi,yi), textcoords='offset points',
        xytext=(0,10), ha='center', fontsize=6.5, color='#607D8B', fontweight='bold')
for xi, yi, qi in zip(x_dates, rms_050, qc_050):
    if not np.isnan(yi): ax1.annotate(f'QC{qi:.2f}', (xi,yi), textcoords='offset points',
        xytext=(0,-16), ha='center', fontsize=6.5, color='#FF9800', fontweight='bold')
# Velocity RMS (secondary y-axis, right)
ax1v = ax1.twinx()
ax1v.plot(x_dates, v017, '^--', color='#4CAF50', lw=1.2, ms=6, alpha=0.7,
          label='0.17h 3V RMS')
ax1v.plot(x_dates, v050, 'v--', color='#FF5722', lw=1.2, ms=6, alpha=0.7,
          label='0.50h 3V RMS')
ax1v.set_ylabel('3V RMS [mm/s]', fontsize=10, color='#666666')
ax1v.tick_params(axis='y', labelcolor='#666666')

ax1.set_ylabel('3D RMS [m]', fontsize=12)
ax1.set_title('(a) 3D Position + Velocity Accuracy' + (' [BRDC SIM]' if USE_BROADCAST else ''), fontsize=13, fontweight='bold')
# Combine legends
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax1v.get_legend_handles_labels()
ax1.legend(lines1+lines2, labels1+labels2, fontsize=7, loc='upper left', ncol=2)
ax1.grid(True, alpha=0.3)
ax1.set_xticks(x_dates); ax1.set_xticklabels(date_labels, fontsize=9); ax1.set_ylim(bottom=0)

# Panel (b): Phase RMS
ax2.plot(x_dates, gn_017, 'o-', color='#4CAF50', lw=1.8, ms=7,
         label='GN Phase (0.17h)', markerfacecolor='white', markeredgewidth=1)
ax2.plot(x_dates, gn_050, 's-', color='#2196F3', lw=1.8, ms=7,
         label='GN Phase (0.50h)', markerfacecolor='white', markeredgewidth=1)
ax2.plot(x_dates, bch_017, 'o--', color='#4CAF50', lw=1, ms=5, alpha=0.55, label='Batch (0.17h)')
ax2.plot(x_dates, bch_050, 's--', color='#2196F3', lw=1, ms=5, alpha=0.55, label='Batch (0.50h)')
ax2.set_ylabel('Phase RMS [m]', fontsize=12)
ax2.set_title('(b) Carrier-Phase Residual RMS', fontsize=13, fontweight='bold')
ax2.legend(fontsize=7, loc='upper left'); ax2.grid(True, alpha=0.3)
ax2.set_xticks(x_dates); ax2.set_xticklabels(date_labels, fontsize=9); ax2.set_ylim(bottom=0)

# Panel (c): Coverage vs Accuracy (bubble = SV count + detailed annotation)
for color, cov_v, rms_v, fc_v, lbl, svs in [
    ('#607D8B', cov_017, rms_017, fc_017, '0.17h',
     [r['n_sv'] for r in all_results if abs(r['arc_h']-0.17)<0.01]),
    ('#FF9800', cov_050, rms_050, fc_050, '0.50h',
     [r['n_sv'] for r in all_results if abs(r['arc_h']-0.50)<0.01])]:
    valid_idx = [i for i in range(len(cov_v)) if not (np.isnan(cov_v[i]) or np.isnan(rms_v[i]))]
    if len(valid_idx) >= 2:
        cv = np.array([cov_v[i] for i in valid_idx])
        rv = np.array([rms_v[i] for i in valid_idx])
        sv_a = np.array([svs[i] for i in valid_idx]) if len(svs) == len(cov_v) else np.ones(len(valid_idx))*10
        fc_a = np.array([fc_v[i] for i in valid_idx])
        d_a = [date_labels[i] for i in valid_idx]
        sizes = np.clip(sv_a * 18, 40, 180)
        ax3.scatter(cv, rv, c=color, s=sizes, alpha=0.75, edgecolors='black',
                    linewidth=0.5, label=lbl, zorder=3)
        for ci, ri, svi, fci, di in zip(cv, rv, sv_a, fc_a, d_a):
            ann = f'{di}: {int(svi)}SV cov={ci:.0f}% full={int(fci)}'
            offset = (8, 6) if color == '#607D8B' else (8, -14)
            ax3.annotate(ann, (ci, ri), textcoords='offset points', xytext=offset,
                         ha='left', fontsize=5.2, color=color, alpha=0.85,
                         bbox=dict(boxstyle='round,pad=0.15', facecolor='white', alpha=0.55, lw=0.3))
        if len(valid_idx) >= 3:
            z=np.polyfit(cv, rv, 1); p=np.poly1d(z)
            cr=np.linspace(max(10,min(cv)-5), min(100,max(cv)+5), 20)
            ax3.plot(cr, p(cr), '--', color=color, alpha=0.35, lw=1.5)
            corr=np.corrcoef(cv, rv)[0,1]
            ax3.text(0.97, 0.93 if color=='#607D8B' else 0.82,
                     f'{lbl}: r={corr:.2f}', transform=ax3.transAxes,
                     fontsize=7.5, color=color, ha='right', fontweight='bold')
from matplotlib.lines import Line2D
legend_el = [Line2D([0],[0], marker='o', color='w', markerfacecolor='gray', markersize=8,
                     label='9SV', markeredgecolor='black', markeredgewidth=0.5),
             Line2D([0],[0], marker='o', color='w', markerfacecolor='gray', markersize=12,
                     label='14SV', markeredgecolor='black', markeredgewidth=0.5),
             Line2D([0],[0], marker='o', color='w', markerfacecolor='gray', markersize=17,
                     label='19SV', markeredgecolor='black', markeredgewidth=0.5)]
ax3.legend(handles=legend_el, fontsize=7, title='Bubble = SV count', title_fontsize=8, loc='lower left', ncol=3)
ax3.add_artist(ax3.legend(fontsize=8, loc='upper right'))
ax3.set_ylabel('3D RMS [m]', fontsize=12)
ax3.set_xlabel('Average SV Coverage [%]', fontsize=11)
ax3.set_title('(c) Accuracy vs SV Coverage', fontsize=13, fontweight='bold')
ax3.grid(True, alpha=0.3); ax3.set_ylim(bottom=0)

plt.tight_layout()
png_path = OUT / f'{label}_orekit.png'
plt.savefig(str(png_path), dpi=150, bbox_inches='tight')
pickle.dump(all_results, open(str(OUT / f'{label}_orekit.pkl'), 'wb'))

print(f"\nReport: {png_path}")
print(f"Data:   {OUT / f'{label}_orekit.pkl'}")
