"""V2.2.4 Multi-hour validation on 2024-04-29.

Processes 24 hours of data into 0.17h and 0.5h arc segments,
computing per-segment 3D RMS vs GNV1B reference.
Generates summary statistics table and PNG plot.

This evaluates algorithm consistency across a full day with
the available CODE precision products.
"""
import sys, os, pickle, math
from pathlib import Path
from datetime import datetime, timedelta
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'src'))

DATE_STR = "2024-04-29"
GRACE_ID = "C"
DATA_ROOT = ROOT / 'data'
GPS_SOD_START = 767620800
INTERVAL = 30
OUT_DIR = ROOT / 'results' / 'multiday'
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Load global data ──
print("Loading...")
gps1b_raw = pickle.load(open(str(DATA_ROOT / f"GPS1B_{DATE_STR}_{GRACE_ID}_04.pkl"), "rb"))
sp3 = pickle.load(open(str(DATA_ROOT / "CODE/2024/COD0OPSFIN_20241200000_01D_05M_ORB.pkl"), "rb"))
from precision_products import (read_rinex_clk, read_antex, load_code_dcb_pair,
                                 compute_dcb_if_correction, setup_iers_from_c04)
clk_data = read_rinex_clk(str(DATA_ROOT / "CODE/2024/COD0OPSFIN_20241200000_01D_30S_CLK.CLK"))
antex = read_antex(str(DATA_ROOT / "igs14.atx"))
dcb_pair = load_code_dcb_pair(str(DATA_ROOT / "CODE/2024/P1C12404.DCB"),
                               str(DATA_ROOT / "CODE/2024/P1P22404.DCB"))
setup_iers_from_c04(str(DATA_ROOT / "IERS/eopc04_IAU2000.txt"))
from gravity_model import read_icgem_gfc
Cnm, Snm, Nmax, GM_grav, R_grav = read_icgem_gfc(str(DATA_ROOT / 'gravity' / 'GGM05C.gfc'))
Nmax = min(Nmax, 90)

print(f"GPS1B={len(gps1b_raw)}ep  SP3={len(sp3['ts'])}ep  grav Nmax={Nmax}")


def run_one_segment(gps_start, arc_hours, chi2_thresh):
    """Run EKF for one segment, return (3D RMS, metrics dict)."""
    from run_sequential_pod import load_gnv1b, compute_epoch_geometry, interpolate_ref
    from sequential_filter import SequentialEKF
    from coordinates import ecef_to_eci, eci_to_ecef

    gps_end = gps_start + int(arc_hours * 3600)
    epochs = sorted(set(g for g in gps1b_raw.keys()
                        if gps_start <= g <= gps_end
                        and abs((g - gps_start) % INTERVAL) <= 2.0))
    if len(epochs) < 6:
        return None, {}

    # Load ref orbit
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
        recs = gps1b_raw.get(int(gps_sod), gps1b_raw.get(gps_sod, {}))
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

    # Init
    r0_ecef = interpolate_ref(ref_orbit, epochs[0])
    v0_ecef = interpolate_ref(ref_vel, epochs[0])
    if r0_ecef is None or v0_ecef is None: return None, {}
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
        'clock_rw': 0.0004 if arc_hours < 0.3 else 0.001,
        'mw_max_epochs': 200,
    }

    ekf = SequentialEKF(ekf_cfg)
    state = ekf.initialize(r0_eci, v0_eci, mjd_start, epochs[0])
    dr_vals = []; sv_counts = []; rej_total = 0

    for i_ep, gps_sod in enumerate(epochs):
        mjd_utc = MJD_J2000 + gps_sod / 86400.0
        if i_ep > 0:
            mjd_prev = MJD_J2000 + epochs[i_ep-1] / 86400.0
            state = ekf.predict(state, gps_sod, mjd_prev, mjd_prev + 69.184 / 86400.0)
        rcv_ecef, _ = eci_to_ecef(state.r_eci, state.v_eci, mjd_utc)
        ep_data = compute_epoch_geometry(gps_sod, gps1b_raw, sp3, rcv_ecef, clk_data)
        if not ep_data: continue
        doy = 120
        state, stats = ekf.process_epoch(state, ep_data, sp3, sv_bias, sv_bias_ref,
                                          mjd_utc, mjd_utc + 69.184 / 86400.0, doy)
        r_ecef, _ = eci_to_ecef(state.r_eci, state.v_eci, mjd_utc)
        r_gnv = interpolate_ref(ref_orbit, gps_sod)
        if r_gnv is not None:
            dr_vals.append(np.linalg.norm(r_ecef - r_gnv))
        sv_counts.append(state.n_sv)
        rej_total += stats['n_rej']

    if len(dr_vals) < 6:
        return None, {}

    rms = np.sqrt(np.mean([d**2 for d in dr_vals]))
    return rms, {
        'rms_3d': rms, 'mean_3d': float(np.mean(dr_vals)),
        'max_3d': float(np.max(dr_vals)), 'min_3d': float(np.min(dr_vals)),
        'n_epochs': len(dr_vals), 'n_sv': int(np.mean(sv_counts)),
        'rej': rej_total,
    }


# ── Run: 24 hours of 0.17h arcs + 6 hours of 0.5h arcs ──
print("\nRunning 0.17h segments (every 30 min, 24h)...")
results_017 = []
for h in range(0, 24):
    for m in range(0, 60, 30):
        gps_start = GPS_SOD_START + h * 3600 + m * 60
        rms, metrics = run_one_segment(gps_start, 0.17, chi2_thresh=25)
        if rms is not None:
            hour = (gps_start + 300 - GPS_SOD_START) / 3600.0
            results_017.append({'hour': hour, 'rms': rms, **metrics})
    if h % 4 == 0:
        n = len([r for r in results_017 if h <= r['hour'] < h+1])
        print(f"  hour {h:2d}: {n} segments")

print(f"\nRunning 0.5h segments (every 2h, 24h)...")
results_05 = []
for h in range(0, 24, 2):
    gps_start = GPS_SOD_START + h * 3600
    rms, metrics = run_one_segment(gps_start, 0.5, chi2_thresh=100)
    if rms is not None:
        hour = (gps_start + 900 - GPS_SOD_START) / 3600.0
        results_05.append({'hour': hour, 'rms': rms, **metrics})
    print(f"  hour {h:2d}: {len([r for r in results_05 if abs(r['hour']-h-0.5)<0.3])} segments")

# ── Statistics ──
def stats_summary(data, label):
    r = [d['rms'] for d in data]
    if not r: return
    print(f"\n{label} ({len(r)} segments):")
    print(f"  Mean={np.mean(r):.3f}m  Median={np.median(r):.3f}m  "
          f"Std={np.std(r):.3f}m")
    print(f"  Best={np.min(r):.3f}m  Worst={np.max(r):.3f}m")
    print(f"  <0.5m: {sum(1 for v in r if v<0.5)}/{len(r)}  "
          f"<0.8m: {sum(1 for v in r if v<0.8)}/{len(r)}")

stats_summary(results_017, "0.17h arcs")
stats_summary(results_05, "0.5h arcs")

# ── Plot ──
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), gridspec_kw={'height_ratios': [3, 1]})

# 0.17h (blue) + 0.5h (orange)
h017 = [r['hour'] for r in results_017]
r017 = [r['rms'] for r in results_017]
h05  = [r['hour'] for r in results_05]
r05  = [r['rms'] for r in results_05]

ax1.scatter(h017, r017, c='#2196F3', s=30, alpha=0.7, marker='o', edgecolors='#1565C0', linewidth=0.4, label=f'0.17h (n={len(r017)})')
ax1.scatter(h05, r05, c='#FF9800', s=50, alpha=0.8, marker='s', edgecolors='#E65100', linewidth=0.6, label=f'0.5h (n={len(r05)})')
ax1.axhline(y=0.30, color='green', linestyle='--', alpha=0.4)
ax1.axhline(y=1.00, color='gray', linestyle='--', alpha=0.3)
ax1.set_ylabel('3D RMS [m]', fontsize=12)
ax1.set_xlabel('GPS Time of Day [hours]', fontsize=12)
ax1.set_title(f'V2.2.4 POD — GRACE-FO C, 2024-04-29, Full-Day Validation\n'
              f'0.17h: mean={np.mean(r017):.3f}m median={np.median(r017):.3f}m best={np.min(r017):.3f}m  '
              f'0.5h: mean={np.mean(r05):.3f}m median={np.median(r05):.3f}m best={np.min(r05):.3f}m',
              fontsize=12)
ax1.legend(fontsize=8)
ax1.set_xlim(-0.5, 24.5)
ax1.grid(True, alpha=0.3)

# Histogram
bins = np.linspace(0, max(max(r017), max(r05)), 25)
ax2.hist(r017, bins=bins, color='#2196F3', alpha=0.5, edgecolor='#1565C0', label=f'0.17h')
ax2.hist(r05, bins=bins, color='#FF9800', alpha=0.5, edgecolor='#E65100', label=f'0.5h')
ax2.axvline(x=0.30, color='green', linestyle='--', alpha=0.4)
ax2.axvline(x=1.0, color='gray', linestyle='--', alpha=0.3)
ax2.set_xlabel('3D RMS [m]', fontsize=12)
ax2.set_ylabel('Count', fontsize=12)
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
png_path = OUT_DIR / 'fullday_v224.png'
plt.savefig(str(png_path), dpi=150, bbox_inches='tight')
print(f"\n  Plot: {png_path}")

pkl_path = OUT_DIR / 'fullday_v224.pkl'
pickle.dump({'results_017': results_017, 'results_05': results_05}, open(str(pkl_path), 'wb'))
print(f"  Data: {pkl_path}")
