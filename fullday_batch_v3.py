"""Full-day batch solver assessment (Framework v3).

Runs EKF + BatchLinearSolver on 0.17h arcs every 2h across 24h.
Compares EKF phase RMS vs batch solver phase RMS.
"""
import sys, os, pickle, math, numpy as np
from pathlib import Path
from datetime import datetime, timedelta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'src'))

DATE_STR = "2024-04-29"; GRACE_ID = "C"; INTERVAL = 30
DATA_ROOT = ROOT / 'data'
GPS_SOD_START = 767620800; OUT_DIR = ROOT / 'results' / 'batch_v3_fullday'
OUT_DIR.mkdir(parents=True, exist_ok=True)

print("Loading shared data...")
gps1b = pickle.load(open(str(DATA_ROOT / f"GPS1B_{DATE_STR}_{GRACE_ID}_04.pkl"), "rb"))
sp3 = pickle.load(open(str(DATA_ROOT / "CODE/2024/COD0OPSFIN_20241200000_01D_05M_ORB.pkl"), "rb"))
from precision_products import *
clk_data = read_rinex_clk(str(DATA_ROOT / "CODE/2024/COD0OPSFIN_20241200000_01D_30S_CLK.CLK"))
antex = read_antex(str(DATA_ROOT / "igs14.atx"))
dcb_pair = load_code_dcb_pair(str(DATA_ROOT / "CODE/2024/P1C12404.DCB"),
                               str(DATA_ROOT / "CODE/2024/P1P22404.DCB"))
setup_iers_from_c04(str(DATA_ROOT / "IERS/eopc04_IAU2000.txt"))
from gravity_model import read_icgem_gfc
Cnm, Snm, Nmax, GM_grav, R_grav = read_icgem_gfc(str(DATA_ROOT / 'gravity' / 'GGM05C.gfc'))
Nmax = min(Nmax, 90)
print(f"GPS1B={len(gps1b)}ep SP3={len(sp3['ts'])}ep grav Nmax={Nmax}")


def run_one_segment(gps_start, arc_hours, chi2_thresh):
    """Run EKF + BatchLinearSolver for one arc segment.
    Returns (rms_3d, ekf_phase_rms, batch_phase_rms, batch_code_rms, n_sv) or None."""
    from run_sequential_pod import load_gnv1b, compute_epoch_geometry, interpolate_ref
    from sequential_filter import SequentialEKF
    from coordinates import ecef_to_eci, eci_to_ecef

    gps_end = gps_start + int(arc_hours * 3600)
    epochs = sorted(set(g for g in gps1b.keys()
                        if gps_start <= g <= gps_end
                        and abs((g - gps_start) % INTERVAL) <= 2.0))
    if len(epochs) < 6: return None

    gnv_path = DATA_ROOT / 'gracefo' / '2024' / DATE_STR / f'GNV1B_{DATE_STR}_{GRACE_ID}_04.txt'
    ref_orbit, ref_vel = load_gnv1b(str(gnv_path))

    # Per-SV code bias
    N_BIAS = min(60, len(epochs))
    sv_p_res = {}
    J2000 = datetime(2000, 1, 1, 12, 0, 0)
    for gps_sod in epochs[:N_BIAS]:
        ref_pos = interpolate_ref(ref_orbit, gps_sod)
        if ref_pos is None: continue
        utc_dt = J2000 + timedelta(seconds=gps_sod)
        recs = gps1b.get(int(gps_sod), gps1b.get(gps_sod, {}))
        from run_sequential_pod import get_sat_geometry
        for sv_id, rec in recs.items():
            if 'P_if' not in rec: continue
            sat_pos, sat_clk, rho_corr = get_sat_geometry(sp3, sv_id, utc_dt, ref_pos, clk_data)
            if sat_pos is None: continue
            dcb_corr = compute_dcb_if_correction(dcb_pair, sv_id)
            P_r = float(rec['P_if']) + dcb_corr + sat_clk - rho_corr
            sv_p_res.setdefault(sv_id, []).append(P_r)
    sv_bias = {}
    for sv, vals in sv_p_res.items():
        if len(vals) >= 3: sv_bias[sv] = float(np.median(vals))
    sv_bias_ref = float(np.mean(list(sv_bias.values()))) if sv_bias else 0.0
    for sv in sv_bias: sv_bias[sv] -= sv_bias_ref

    # Init EKF
    r0_ecef = interpolate_ref(ref_orbit, epochs[0])
    v0_ecef = interpolate_ref(ref_vel, epochs[0])
    if r0_ecef is None or v0_ecef is None: return None
    MJD_J2000 = 51544.5
    mjd_start = MJD_J2000 + epochs[0] / 86400.0
    r0_eci, v0_eci = ecef_to_eci(r0_ecef, v0_ecef, mjd_start)

    ekf_cfg = {
        'Cd': 2.2, 'CR': 1.3, 'area_drag': 0.68, 'area_srp': 3.4, 'mass': 580.0,
        'bodies': ['Sun', 'Moon'], 'Cnm': Cnm, 'Snm': Snm,
        'GM_grav': GM_grav, 'R_grav': R_grav, 'gravity_nmax': Nmax,
        'sigma_acc_process': 1e-3, 'tau_emp': 600.0, 'sigma_emp_ss': 1e-8,
        'sigma_zwd_rw': 1e-9, 'sigma_phase': 0.20, 'sigma_code': 0.30,
        'chi2_threshold': chi2_thresh, 'el_min': 0.087,
        'use_phase_windup': True, 'use_relativity': True, 'use_cycle_slip': False,
        'ar_min_epochs': 6, 'antex_data': antex, 'dcb_data': dcb_pair,
        'elev_exp_phase': 1.0, 'elev_exp_code': 0.70 if arc_hours >= 0.3 else 1.0,
        'clock_rw': 0.0004 if arc_hours < 0.3 else 0.001, 'mw_max_epochs': 200,
    }

    ekf = SequentialEKF(ekf_cfg)
    state = ekf.initialize(r0_eci, v0_eci, mjd_start, epochs[0])
    dr_vals = []; n_rej_total = 0
    pass1_geometry = []

    SEC_PER_DAY = 86400.0; C_LIGHT = 299792458.0; OMEGA_E = 7.2921151467e-5

    for i_ep, gps_sod in enumerate(epochs):
        mjd_utc = MJD_J2000 + gps_sod / 86400.0
        if i_ep > 0:
            mjd_prev = MJD_J2000 + epochs[i_ep-1] / 86400.0
            state = ekf.predict(state, gps_sod, mjd_prev, mjd_prev + 69.184/86400.0)
        rcv_ecef, _ = eci_to_ecef(state.r_eci, state.v_eci, mjd_utc)
        ep_data = compute_epoch_geometry(gps_sod, gps1b, sp3, rcv_ecef, clk_data)
        if not ep_data: continue
        doy = 120
        state, stats = ekf.process_epoch(state, ep_data, sp3, sv_bias, sv_bias_ref,
                                          mjd_utc, mjd_utc + 69.184/86400.0, doy)
        r_ecef, _ = eci_to_ecef(state.r_eci, state.v_eci, mjd_utc)
        r_gnv = interpolate_ref(ref_orbit, gps_sod)
        if r_gnv is not None:
            dr_vals.append(np.linalg.norm(r_ecef - r_gnv))
        n_rej_total += stats['n_rej']

        # Store geometry for batch solver
        from src.troposphere import saastamoinen_zhd
        lat_rad = np.arcsin(rcv_ecef[2] / np.linalg.norm(rcv_ecef))
        h_m = np.linalg.norm(rcv_ecef) - 6378137.0; zhd = saastamoinen_zhd(lat_rad, h_m)
        for d in ep_data:
            sv = d['sv']
            se = np.asarray(d['sat_pos'], dtype=float); sc = float(d.get('sat_clk', 0))
            el = float(d.get('el', 0.5))
            rho = np.linalg.norm(se - rcv_ecef)
            sag = (OMEGA_E/C_LIGHT)*(se[0]*rcv_ecef[1] - se[1]*rcv_ecef[0])
            mf = 1.001/np.sqrt(0.002001+np.sin(el)**2)
            dcb_c = compute_dcb_if_correction(dcb_pair, sv)
            d['_geo_full'] = rho + sag - sc + zhd * mf
            d['_obs_code'] = float(d.get('P_if_raw', 0)) + dcb_c - sv_bias.get(sv, 0.0)
            d['_obs_phase'] = float(d.get('L_if_raw', 0)) - sv_bias.get(sv, 0.0)
        pass1_geometry.append(ep_data)

    rms_3d = np.sqrt(np.mean([d**2 for d in dr_vals])) if dr_vals else 0
    ekf_phase_rms = stats['rms_phase'] if rms_3d > 0 else 0

    # Batch solver on EKF orbit
    from src.batch_solver import BatchLinearSolver
    bls = BatchLinearSolver(pass1_geometry, sigma_phase=0.20, sigma_code=0.30)
    bls_sol = bls.solve()
    batch_phase_rms = bls_sol['rms_phase']
    batch_code_rms = bls_sol['rms_code']

    return (rms_3d, ekf_phase_rms, batch_phase_rms, batch_code_rms,
            state.n_sv, (gps_start + 300 - GPS_SOD_START) / 3600.0)


# ── Run ──
print("\nRunning 0.17h arcs every 2h (12 segments)...")
results = []
for h in range(0, 24, 2):
    gps_start = GPS_SOD_START + h * 3600
    ret = run_one_segment(gps_start, 0.17, chi2_thresh=25)
    if ret:
        rms, epr, bpr, bcr, nsv, hour = ret
        imp = (epr - bpr)/epr*100 if epr > 0 else 0
        results.append({'hour': hour, 'rms_3d': rms, 'ekf_phase': epr,
                        'batch_phase': bpr, 'batch_code': bcr, 'n_sv': nsv,
                        'improvement': imp})
        print(f"  hour {h:2d}: EKF={rms:.3f}m  phase: {epr:.3f}m→{bpr:.3f}m "
              f"({imp:+.0f}%)  SVs={nsv}")

# ── Summary ──
h = [r['hour'] for r in results]
ekf_p = [r['ekf_phase'] for r in results]
bat_p = [r['batch_phase'] for r in results]
imp_p = [r['improvement'] for r in results]
rms_v = [r['rms_3d'] for r in results]

print(f"\n=== Full-Day Batch Solver Assessment ({len(results)} segments) ===")
print(f"  EKF phase RMS:   mean={np.mean(ekf_p):.3f}m  best={np.min(ekf_p):.3f}m  worst={np.max(ekf_p):.3f}m")
print(f"  Batch phase RMS: mean={np.mean(bat_p):.3f}m  best={np.min(bat_p):.3f}m  worst={np.max(bat_p):.3f}m")
print(f"  Improvement:     mean={np.mean(imp_p):+.1f}%  best={np.max(imp_p):+.1f}%  worst={np.min(imp_p):+.1f}%")

# ── Plot ──
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), gridspec_kw={'height_ratios': [2,1]})

ax1.plot(h, ekf_p, 'o-', color='#FF9800', linewidth=1.5, markersize=6, label='EKF sequential')
ax1.plot(h, bat_p, 's-', color='#2196F3', linewidth=1.5, markersize=6, label='Batch solver (Framework v3)')
ax1.set_ylabel('Phase RMS [m]', fontsize=12)
ax1.set_title('GRACE-FO C Phase Residuals — EKF vs Batch Solver\n'
              f'({len(results)} × 0.17h arcs, 2024-04-29)',
              fontsize=12)
ax1.legend(fontsize=10); ax1.grid(True, alpha=0.3)
ax1.set_xlim(-0.5, 24.5)

ax2.bar(np.array(h)-0.25, imp_p, 0.45, color='#4CAF50', label='Phase RMS improvement')
ax2.axhline(y=0, color='gray', linewidth=0.5)
ax2.set_ylabel('Improvement [%]', fontsize=12); ax2.set_xlabel('GPS Time of Day [hours]', fontsize=12)
ax2.legend(fontsize=10); ax2.grid(True, alpha=0.3)
ax2.set_xlim(-0.5, 24.5)

plt.tight_layout()
png_path = OUT_DIR / 'batch_v3_fullday.png'
plt.savefig(str(png_path), dpi=150, bbox_inches='tight')
print(f"\n  Plot: {png_path}")

pkl_path = OUT_DIR / 'batch_v3_fullday.pkl'
pickle.dump(results, open(str(pkl_path), 'wb'))
print(f"  Data: {pkl_path}")
