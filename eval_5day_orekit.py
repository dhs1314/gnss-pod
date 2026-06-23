"""5-day Orekit GN validation (2024-04-29 to 2024-05-03).

Runs full pipeline on each date:
  EKF pass 1 (simplified dynamics) -> Orekit GN outer loop (Phase 20.0)

Reports: 3D RMS vs GNV1B, Phase RMS, GN convergence status.
Generates: results/5day_orekit/5day_orekit.png
"""
import sys, os, pickle, time, numpy as np
from pathlib import Path
from datetime import datetime, timedelta
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'src'))

DATES = ['2024-04-29', '2024-04-30', '2024-05-01', '2024-05-02', '2024-05-03']
GRACE = 'C'; INTERVAL = 30; ARC_HOURS = [0.17, 0.50]
DR = ROOT / 'data'
MJD0 = 51544.5; SEC = 86400.0; C_L = 299792458.0; OM = 7.2921151467e-5

OUT = ROOT / 'results' / '5day_orekit'
OUT.mkdir(parents=True, exist_ok=True)

# ---- Shared data loading ----
print("=" * 70)
print("Orekit GN 5-Day Validation (Phase 20.0)")
print("=" * 70)

print("\nLoading shared products...")
sp3_pkls = {}
clk_data = {}
dcb_pairs = {}
antex = None

for date_str in DATES:
    dt = datetime.strptime(date_str, "%Y-%m-%d")
    doy = dt.strftime("%j"); y = dt.year; m = dt.month
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

print(f"  Loaded: SP3={len(sp3_pkls)}d, CLK={len(clk_data)}d, DCB={len(dcb_pairs)}d")
print(f"  Gravity: GGM05C Nmax={GRAV_NMAX}")

# ---- Run one arc ----
from src.troposphere import saastamoinen_zhd
from src.batch_orbit_v3 import BatchOrbitLSQv3
from src.orbit_dynamics import total_acc_eci
from src.batch_solver import BatchLinearSolver
from run_sequential_pod import load_gnv1b, compute_epoch_geometry, interpolate_ref
from coordinates import ecef_to_eci, eci_to_ecef
from sequential_filter import SequentialEKF

def run_one_arc(date_str, arc_h):
    """Run EKF + Orekit GN on one arc. Returns dict or None."""
    dt = datetime.strptime(date_str, "%Y-%m-%d"); y = dt.year

    # GPS1B
    gps_path = DR / 'gracefo' / str(y) / date_str / f'GPS1B_{date_str}_{GRACE}_04.pkl'
    if not gps_path.exists():
        return None
    gps1b = pickle.load(open(str(gps_path), "rb"))

    # SP3 and CLK
    sp3 = sp3_pkls.get(date_str)
    clk = clk_data.get(date_str)
    dcb = dcb_pairs.get(date_str)
    if sp3 is None or clk is None or dcb is None:
        return None

    # GNV1B reference
    gnv_path = DR / 'gracefo' / str(y) / date_str / f'GNV1B_{date_str}_{GRACE}_04.txt'
    if not gnv_path.exists(): return None
    ref, refv = load_gnv1b(str(gnv_path))

    # Time windows
    gps0 = min(gps1b.keys())
    gps_end = gps0 + int(arc_h * 3600)
    epochs = sorted(set(g for g in gps1b if gps0 <= g <= gps_end
                        and abs((g - gps0) % INTERVAL) <= 2))
    if len(epochs) < 6: return None
    t_ep = np.array([(g - epochs[0]) for g in epochs], dtype=float)
    mjd_s = MJD0 + epochs[0] / SEC; mjd_tt = mjd_s + 69.184 / SEC
    r0e, v0e = interpolate_ref(ref, epochs[0]), interpolate_ref(refv, epochs[0])
    if r0e is None or v0e is None: return None
    r0i, v0i = ecef_to_eci(r0e, v0e, mjd_s)

    # Per-SV code bias
    N_BIAS = min(60, len(epochs))
    sv_p_res = {}
    J2000 = datetime(2000, 1, 1, 12, 0, 0)
    for gps_sod in epochs[:N_BIAS]:
        utc_dt = J2000 + timedelta(seconds=gps_sod)
        recs = gps1b.get(int(gps_sod), gps1b.get(gps_sod, {}))
        from run_sequential_pod import get_sat_geometry
        for sv_id, rec in recs.items():
            if 'P_if' not in rec: continue
            rcv_r = interpolate_ref(ref, gps_sod)
            if rcv_r is None: continue
            sat_pos, sat_clk, rho_corr = get_sat_geometry(sp3, sv_id, utc_dt, rcv_r, clk)
            if sat_pos is None: continue
            from precision_products import compute_dcb_if_correction
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
        'use_phase_windup': True, 'use_relativity': True, 'use_cycle_slip': False,
        'ar_min_epochs': 6, 'antex_data': antex, 'dcb_data': dcb,
        'elev_exp_phase': 1.0, 'elev_exp_code': 0.70 if arc_h >= 0.3 else 1.0,
        'clock_rw': 0.001 if arc_h >= 0.3 else 0.0004, 'mw_max_epochs': 200,
    }
    ekf = SequentialEKF(ekf_cfg)
    state = ekf.initialize(r0i, v0i, mjd_s, epochs[0])
    pass1 = []
    for i_ep, gps_sod in enumerate(epochs):
        mjd_u = MJD0 + gps_sod / SEC
        if i_ep > 0:
            mjd_prev = MJD0 + epochs[i_ep-1] / SEC
            state = ekf.predict(state, gps_sod, mjd_prev, mjd_prev + 69.184 / SEC)
        rcv_e, _ = eci_to_ecef(state.r_eci, state.v_eci, mjd_u)
        ep_data = compute_epoch_geometry(gps_sod, gps1b, sp3, rcv_e, clk)
        if not ep_data: continue
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

    ekf_phase = stats.get('rms_phase', 0)
    n_sv = state.n_sv

    # Baseline batch on EKF orbit
    bls = BatchLinearSolver(pass1, sigma_phase=0.20, sigma_code=0.30)
    bls_sol = bls.solve()
    batch_ph = bls_sol['rms_phase']

    # Orekit GN
    prop = OrekitPropagator(
        gravity_field=str(grav_path), gravity_degree=GRAV_NMAX,
        solid_tides=True, ocean_tides=True, ocean_tide_degree=50,
        third_body='lunisolar', srp_model='isotropic', relativity=True,
        drag_model='exponential',
        mass=580.0, area_drag=0.68, area_srp=3.4, CR=1.3, CD=2.2,
        stm_perturb=1.0, integrator_tol=1e-12)
    prop._setup()

    gn = BatchOrbitLSQv3(
        pass1, gn_fn, t_ep,
        mjd_utc_start=mjd_s, mjd_tt_start=mjd_tt,
        sigma_phase=0.20, sigma_code=0.30,
        max_iter=6, prior_r0=1.0, prior_v0=0.01, prior_emp=1e-7,
        damping=0.5, orekit_prop=prop, estimate_cd_cr=False)

    t0 = time.time()
    sol = gn.solve(r0i, v0i)
    dt_gn = time.time() - t0

    # 3D RMS
    dr_vals = []
    for i_ep, gps_sod in enumerate(epochs):
        r_gnv = interpolate_ref(ref, gps_sod)
        if r_gnv is not None:
            mjd_u = MJD0 + gps_sod / SEC
            r_e, _ = eci_to_ecef(sol['r_eci'][i_ep], sol['v_eci'][i_ep], mjd_u)
            dr_vals.append(np.linalg.norm(r_e - r_gnv))
    rms_3d = np.sqrt(np.mean([d**2 for d in dr_vals])) if dr_vals else 0

    return {
        'date': date_str, 'arc_h': arc_h,
        'rms_3d': rms_3d, 'phase_ekf': ekf_phase,
        'phase_batch': batch_ph, 'phase_gn': sol['rms_phase'],
        'converged': sol['converged'], 'iterations': sol['iterations'],
        'n_sv': n_sv, 'time_gn': dt_gn,
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
                print(f"3D={r['rms_3d']:.3f}m GN_ph={r['phase_gn']:.3f}m "
                      f"t={r['time_gn']:.0f}s SV={r['n_sv']}")
            else:
                print("SKIP (no data)")
        except Exception as e:
            print(f"ERROR: {e}")
            import traceback; traceback.print_exc()

if not all_results:
    print("\nNo results!"); sys.exit(1)

# ---- Summary ----
print(f"\n{'='*70}")
print(f"5-Day Summary ({len(all_results)} arcs)")
print(f"{'='*70}")

for arc_h in ARC_HOURS:
    sub = [r for r in all_results if r['arc_h'] == arc_h]
    if not sub: continue
    rms_v = [r['rms_3d'] for r in sub]
    print(f"\n{arc_h:.2f}h ({len(sub)} dates):")
    print(f"  3D RMS: mean={np.mean(rms_v):.3f}m median={np.median(rms_v):.3f}m "
          f"best={np.min(rms_v):.3f}m worst={np.max(rms_v):.3f}m")

# Table
print(f"\n{'Date':<12s} {'Arc':>6s} {'3D_RMS':>8s} {'GN_Ph':>7s} {'Bch_Ph':>7s} "
      f"{'EKF_Ph':>7s} {'Conv':>5s} {'SVs':>4s} {'Time':>6s}")
for r in all_results:
    conv = "Y" if r['converged'] else "N"
    print(f"{r['date']:<12s} {r['arc_h']:5.2f}h {r['rms_3d']:8.3f} {r['phase_gn']:7.3f} "
          f"{r['phase_batch']:7.3f} {r['phase_ekf']:7.3f} {conv:>5s} "
          f"{r['n_sv']:>4d} {r['time_gn']:5.0f}s")

# ---- Plot ----
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8),
    gridspec_kw={'height_ratios': [2, 1]})
w = 0.3

for i, arc_h in enumerate(ARC_HOURS):
    sub = [r for r in all_results if r['arc_h'] == arc_h]
    labels = [r['date'][5:] for r in sub]
    x = np.arange(len(sub)) + i * (len(sub) + 0.5)
    rms_v = [r['rms_3d'] for r in sub]
    gn_v  = [r['phase_gn'] for r in sub]
    bch_v = [r['phase_batch'] for r in sub]
    ekf_v = [r['phase_ekf'] for r in sub]

    ax1.bar(x - w, rms_v, w*0.9, color='#607D8B', alpha=0.85,
            label=f'{arc_h}h 3D RMS')
    ax1.bar(x, gn_v, w*0.9, color='#4CAF50', alpha=0.85,
            label=f'{arc_h}h GN Phase')
    ax1.bar(x + w, bch_v, w*0.9, color='#2196F3', alpha=0.85,
            label=f'{arc_h}h Batch Ph')

    for xi, rv, gv in zip(x, rms_v, gn_v):
        imp = (r['phase_ekf'] - r['phase_batch']) / r['phase_ekf'] * 100 if r['phase_ekf'] > 0 else 0
        ax2.bar(xi + w*0.5, [imp], w*0.8,
                color='#FF9800' if imp > 0 else '#F44336', alpha=0.7)

ax1.set_ylabel('RMS [m]', fontsize=12)
ax1.set_title('Orekit GN 5-Day Validation (Phase 20.0, 0.043m best)', fontsize=14)
ax1.legend(fontsize=8, loc='upper left', ncol=3)
ax1.grid(True, alpha=0.3)
all_dates = sorted(set(r['date'][5:] for r in all_results))
ax1.set_xticks(range(len(all_dates))); ax1.set_xticklabels(all_dates, fontsize=10)

ax2.set_ylabel('Batch Ph Impr [%]', fontsize=11)
ax2.set_xlabel('Date (MM-DD)', fontsize=11)
ax2.axhline(y=0, color='gray', linewidth=0.5)
ax2.grid(True, alpha=0.3)
ax2.set_xticks(range(len(all_dates))); ax2.set_xticklabels(all_dates, fontsize=10)

plt.tight_layout()
png_path = OUT / '5day_orekit.png'
plt.savefig(str(png_path), dpi=120)
pickle.dump(all_results, open(str(OUT / '5day_orekit.pkl'), 'wb'))

print(f"\nReport: {png_path}")
print(f"Data:   {OUT / '5day_orekit.pkl'}")
