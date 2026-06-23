"""Full-day accuracy assessment: 10-min segments at 6 representative hours.

Processes 6 × 6 = 36 ten-minute segments spread across the day
(every 4 hours), computes 3D RMS per segment, and saves a PNG plot.
"""
import sys, os, pickle, math
from pathlib import Path
from datetime import datetime, timedelta
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))

DP = Path(r"d:\prj\gnss_pod\data")
DATE_STR = "2024-04-29"
GRACE_ID = "C"
GPS_SOD_START = 767620800
INTERVAL = 30
HOURS_TO_TEST = [0, 4, 8, 12, 16, 20]

# ── Load data once ──
print("Loading...")
gps1b_raw = pickle.load(open(str(DP / f"GPS1B_{DATE_STR}_{GRACE_ID}_04.pkl"), "rb"))
from run_sequential_pod import load_gnv1b
gnv_path = DP / 'gracefo' / '2024' / DATE_STR / f'GNV1B_{DATE_STR}_{GRACE_ID}_04.txt'
ref_orbit, ref_vel = load_gnv1b(str(gnv_path))
sp3 = pickle.load(open(str(DP / "CODE/2024/COD0OPSFIN_20241200000_01D_05M_ORB.pkl"), "rb"))
from src.precision_products import read_rinex_clk, read_antex, load_code_dcb_pair, compute_dcb_if_correction, setup_iers_from_c04
clk_data = read_rinex_clk(str(DP / "CODE/2024/COD0OPSFIN_20241200000_01D_30S_CLK.CLK"))
antex = read_antex(str(DP / "igs14.atx"))
dcb_pair = load_code_dcb_pair(str(DP / "CODE/2024/P1C12404.DCB"),
                               str(DP / "CODE/2024/P1P22404.DCB"))
dcb_if = {prn: compute_dcb_if_correction(dcb_pair, prn) for prn in dcb_pair}
setup_iers_from_c04(str(DP / "IERS/eopc04_IAU2000.txt"))
from src.gravity_model import read_icgem_gfc
gfc_path = DP / 'gravity' / 'GGM05C.gfc'
Cnm, Snm, Nmax, GM_grav, R_grav = read_icgem_gfc(str(gfc_path))
Nmax = min(Nmax, 90)
print(f"Loaded: SP3={len(sp3['ts'])} CLK={len(clk_data)} grav Nmax={Nmax}")

from src.coordinates import ecef_to_eci, eci_to_ecef
from src.sequential_filter import SequentialEKF
from run_sequential_pod import compute_epoch_geometry, interpolate_ref
from run_sequential_pod import get_sat_geometry

J2000 = datetime(2000, 1, 1, 12, 0, 0)
MJD_J2000 = 51544.5

ekf_cfg = {
    'Cd': 2.2, 'CR': 1.3, 'area_drag': 0.68, 'area_srp': 3.4, 'mass': 580.0,
    'dt_integ': 10.0, 'bodies': ['Sun', 'Moon'],
    'Cnm': Cnm, 'Snm': Snm, 'GM_grav': GM_grav, 'R_grav': R_grav,
    'gravity_nmax': Nmax,
    'sigma_acc_process': 1e-3, 'tau_emp': 600.0, 'sigma_emp_ss': 1e-8,
    'sigma_zwd_rw': 1e-9, 'sigma_phase': 0.20, 'sigma_code': 0.30,
    'chi2_threshold': 100.0, 'el_min': 0.087,
    'use_phase_windup': True, 'use_relativity': True, 'use_cycle_slip': False,
    'ar_min_epochs': 6,
    'antex_data': antex, 'dcb_data': dcb_pair,
    'elev_exp_phase': 1.0, 'elev_exp_code': 0.70,
}

results_all = []

for test_hour in HOURS_TO_TEST:
    gps_start = GPS_SOD_START + test_hour * 3600
    # Process 6 consecutive 10-minute segments within this hour
    for seg_offset in range(0, 3600, 600):
        seg_start = gps_start + seg_offset
        seg_end = seg_start + 600

        segment_epochs = []
        for gps_sod in sorted(gps1b_raw.keys()):
            if seg_start <= gps_sod < seg_end:
                dt_ep = gps_sod - seg_start
                nearest = round(dt_ep / INTERVAL) * INTERVAL
                if abs(dt_ep - nearest) <= 2.0:
                    segment_epochs.append(gps_sod)
        segment_epochs = sorted(set(segment_epochs))
        if len(segment_epochs) < 6:
            continue

        # Per-SV code bias
        sv_bias_local = {}
        sv_code_res = {}
        for gps_sod in segment_epochs[:min(30, len(segment_epochs))]:
            ref_pos = interpolate_ref(ref_orbit, gps_sod)
            if ref_pos is None:
                continue
            utc_dt = J2000 + timedelta(seconds=gps_sod)
            recs = gps1b_raw.get(gps_sod, {})
            for sv_id, rec in recs.items():
                if 'L_if' not in rec:
                    continue
                sat_pos, sat_clk, rho_corr = get_sat_geometry(sp3, sv_id, utc_dt, ref_pos, clk_data)
                if sat_pos is None:
                    continue
                rho_tmp = float(np.linalg.norm(sat_pos - ref_pos))
                P_r = float(rec['P_if']) + dcb_if.get(sv_id, 0.0) + sat_clk - rho_tmp
                sv_code_res.setdefault(sv_id, []).append(P_r)
        for sv_id, vals in sv_code_res.items():
            if len(vals) >= 3:
                sv_bias_local[sv_id] = float(np.median(vals))
        if sv_bias_local:
            ref_m = float(np.mean(list(sv_bias_local.values())))
            for sv in sv_bias_local:
                sv_bias_local[sv] -= ref_m
        sv_bias_ref = ref_m if sv_bias_local else 0.0

        # Initialize EKF
        r0_ecef = interpolate_ref(ref_orbit, segment_epochs[0])
        v0_ecef = interpolate_ref(ref_vel, segment_epochs[0])
        if r0_ecef is None or v0_ecef is None:
            continue
        mjd_s = MJD_J2000 + segment_epochs[0] / 86400.0
        r0_eci, v0_eci = ecef_to_eci(r0_ecef, v0_ecef, mjd_s)

        ekf = SequentialEKF(ekf_cfg)
        state = ekf.initialize(r0_eci, v0_eci, mjd_s, segment_epochs[0])
        segment_dr = []

        for i_ep, gps_sod in enumerate(segment_epochs):
            mjd_utc = MJD_J2000 + gps_sod / 86400.0
            mjd_tt = mjd_utc + 69.184 / 86400.0
            if i_ep > 0:
                mjd_up = MJD_J2000 + segment_epochs[i_ep - 1] / 86400.0
                state = ekf.predict(state, gps_sod, mjd_up, mjd_up + 69.184 / 86400.0)
            rcv_ecef, _ = eci_to_ecef(state.r_eci, state.v_eci, mjd_utc)
            ep_data = compute_epoch_geometry(gps_sod, gps1b_raw, sp3, rcv_ecef, clk_data)
            if not ep_data:
                continue
            doy = 120
            state, _ = ekf.process_epoch(state, ep_data, sp3, sv_bias_local, sv_bias_ref, mjd_utc, mjd_tt, doy)
            r_ecef, _ = eci_to_ecef(state.r_eci, state.v_eci, mjd_utc)
            r_gnv = interpolate_ref(ref_orbit, gps_sod)
            if r_gnv is not None:
                segment_dr.append(np.linalg.norm(r_ecef - r_gnv))

        if len(segment_dr) >= 6:
            rms_3d = np.sqrt(np.mean([d**2 for d in segment_dr]))
            hour = (seg_start + 300 - GPS_SOD_START) / 3600.0
            results_all.append({'hour': hour, 'rms_3d': rms_3d, 'n_epochs': len(segment_dr)})

    n_done = len(results_all)
    print(f"  hour {test_hour:2d}: {sum(1 for r in results_all if abs(r['hour'] - test_hour) < 0.5)} segments done, total={n_done}")

# ── Results ──
hours = [r['hour'] for r in results_all]
rms_vals = [r['rms_3d'] for r in results_all]
print(f"\n=== Summary ===")
print(f"  Segments: {len(results_all)}")
print(f"  Mean RMS: {np.mean(rms_vals):.3f}m")
print(f"  Median:   {np.median(rms_vals):.3f}m")
print(f"  Best:     {np.min(rms_vals):.3f}m at {hours[np.argmin(rms_vals)]:.1f}h")
print(f"  Worst:    {np.max(rms_vals):.3f}m at {hours[np.argmax(rms_vals)]:.1f}h")

# ── Plot ──
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), gridspec_kw={'height_ratios': [3, 1]})

# Main scatter
ax1.scatter(hours, rms_vals, c='#2196F3', s=40, alpha=0.7, edgecolors='#1565C0', linewidth=0.5)
ax1.axhline(y=0.50, color='green', linestyle='--', alpha=0.5, label='0.50 m')
ax1.axhline(y=1.00, color='orange', linestyle='--', alpha=0.5, label='1.00 m')
ax1.axhline(y=np.median(rms_vals), color='red', linestyle='-', alpha=0.7, label=f'Median {np.median(rms_vals):.3f}m')
ax1.set_ylabel('3D RMS [m]', fontsize=12)
ax1.set_xlabel('GPS Time of Day [hours]', fontsize=12)
ax1.set_title(f'GRACE-FO C POD Accuracy — 2024-04-29 ({len(results_all)}×10min segments)\n'
              f'Mean={np.mean(rms_vals):.3f}m  Median={np.median(rms_vals):.3f}m  '
              f'Best={np.min(rms_vals):.3f}m  Worst={np.max(rms_vals):.3f}m',
              fontsize=13)
ax1.legend(fontsize=10, loc='upper right')
ax1.set_xlim(-0.5, 24.5)
ax1.grid(True, alpha=0.3)

# Histogram
ax2.hist(rms_vals, bins=20, color='#2196F3', alpha=0.7, edgecolor='#1565C0')
ax2.axvline(x=0.50, color='green', linestyle='--', alpha=0.5)
ax2.axvline(x=np.median(rms_vals), color='red', linestyle='-', alpha=0.7)
ax2.set_xlabel('3D RMS [m]', fontsize=12)
ax2.set_ylabel('Count', fontsize=12)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
out_png = Path('results') / 'accuracy_2024-04-29.png'
out_png.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(str(out_png), dpi=150, bbox_inches='tight')
print(f"\n  Saved: {out_png}")

# Also save data
pickle.dump({'hours': hours, 'rms_vals': rms_vals, 'results': results_all},
            open(str(out_png.with_suffix('.pkl')), 'wb'))
