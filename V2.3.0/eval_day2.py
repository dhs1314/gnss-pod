"""Full-day assessment: run EKF at 6 representative hours, plot 3D RMS.

Uses the existing run_sequential_pod data loading pipeline but starts
from different GPS seconds of day for each test hour.
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
SEGMENT_HOURS = 0.5  # each test run covers 30 min
HOURS = [0, 4, 8, 12, 16, 20]

# ── Load all shared data ──
print("Loading data...")
gps1b_raw = pickle.load(open(str(DP / f"GPS1B_{DATE_STR}_{GRACE_ID}_04.pkl"), "rb"))
epochs_all = sorted(gps1b_raw.keys())
sp3 = pickle.load(open(str(DP / "CODE/2024/COD0OPSFIN_20241200000_01D_05M_ORB.pkl"), "rb"))

# Run each hour
results = []
for test_hour in HOURS:
    gps_sod_start = GPS_SOD_START + test_hour * 3600
    gps_sod_end = gps_sod_start + SEGMENT_HOURS * 3600

    # Select epochs
    epochs = []
    for gps_sod in epochs_all:
        if not (gps_sod_start <= gps_sod <= gps_sod_end):
            continue
        dt_ep = gps_sod - gps_sod_start
        nearest = round(dt_ep / INTERVAL) * INTERVAL
        if abs(dt_ep - nearest) <= 2.0:
            epochs.append(gps_sod)
    epochs = sorted(set(epochs))

    if len(epochs) < 10:
        print(f"  Hour {test_hour}: only {len(epochs)} epochs, skip")
        continue

    print(f"  Hour {test_hour}: {len(epochs)} epochs, processing...")

    # ── Re-use the main EKF pipeline for this hour ──
    # Load GNV1B ref
    from run_sequential_pod import load_gnv1b, compute_epoch_geometry, interpolate_ref
    gnv_path = DP / 'gracefo' / '2024' / DATE_STR / f'GNV1B_{DATE_STR}_{GRACE_ID}_04.txt'
    ref_orbit, ref_vel = load_gnv1b(str(gnv_path))

    # Load CLK
    from src.precision_products import read_rinex_clk, read_antex, load_code_dcb_pair, compute_dcb_if_correction, setup_iers_from_c04
    clk_data = read_rinex_clk(str(DP / "CODE/2024/COD0OPSFIN_20241200000_01D_30S_CLK.CLK"))
    antex = read_antex(str(DP / "igs14.atx"))
    dcb_pair = load_code_dcb_pair(str(DP / "CODE/2024/P1C12404.DCB"),
                                   str(DP / "CODE/2024/P1P22404.DCB"))
    dcb_if = {prn: compute_dcb_if_correction(dcb_pair, prn) for prn in dcb_pair}
    setup_iers_from_c04(str(DP / "IERS/eopc04_IAU2000.txt"))

    # Load gravity
    from src.gravity_model import read_icgem_gfc
    gfc_path = DP / 'gravity' / 'GGM05C.gfc'
    Cnm, Snm, Nmax, GM_grav, R_grav = read_icgem_gfc(str(gfc_path))
    Nmax = min(Nmax, 90)

    # Estimate per-SV code bias
    from run_sequential_pod import get_sat_geometry
    J2000 = datetime(2000, 1, 1, 12, 0, 0)
    N_BIAS = min(60, len(epochs))
    sv_p_residuals = {}
    for gps_sod in epochs[:N_BIAS]:
        ref_pos = interpolate_ref(ref_orbit, gps_sod)
        if ref_pos is None:
            continue
        utc_dt = J2000 + timedelta(seconds=gps_sod)
        recs = gps1b_raw.get(int(gps_sod), gps1b_raw.get(gps_sod, {}))
        for sv_id, rec in recs.items():
            if 'P_if' not in rec:
                continue
            sat_pos, sat_clk, rho_corr = get_sat_geometry(sp3, sv_id, utc_dt, ref_pos, clk_data)
            if sat_pos is None:
                continue
            P_r = float(rec['P_if']) + dcb_if.get(sv_id, 0) + sat_clk - rho_corr
            sv_p_residuals.setdefault(sv_id, []).append(P_r)
    sv_bias = {}
    for sv, vals in sv_p_residuals.items():
        if len(vals) >= 3:
            sv_bias[sv] = float(np.median(vals))
    if sv_bias:
        sv_bias_ref = float(np.mean(list(sv_bias.values())))
        for sv in sv_bias:
            sv_bias[sv] -= sv_bias_ref
    else:
        sv_bias_ref = 0.0

    # Initialize EKF
    from src.coordinates import ecef_to_eci, eci_to_ecef
    from src.sequential_filter import SequentialEKF
    MJD_J2000 = 51544.5

    r0_ecef = interpolate_ref(ref_orbit, epochs[0])
    v0_ecef = interpolate_ref(ref_vel, epochs[0])
    if r0_ecef is None or v0_ecef is None:
        continue
    mjd_start = MJD_J2000 + epochs[0] / 86400.0
    r0_eci, v0_eci = ecef_to_eci(r0_ecef, v0_ecef, mjd_start)

    ekf_cfg = {
        'Cd': 2.2, 'CR': 1.3, 'area_drag': 0.68, 'area_srp': 3.4, 'mass': 580.0,
        'bodies': ['Sun', 'Moon'],
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

    ekf = SequentialEKF(ekf_cfg)
    state = ekf.initialize(r0_eci, v0_eci, mjd_start, epochs[0])

    # Run EKF
    epoch_rms = []
    for i_ep, gps_sod in enumerate(epochs):
        mjd_utc = MJD_J2000 + gps_sod / 86400.0
        mjd_tt = mjd_utc + 69.184 / 86400.0
        if i_ep > 0:
            mjd_prev = MJD_J2000 + epochs[i_ep - 1] / 86400.0
            state = ekf.predict(state, gps_sod, mjd_prev, mjd_prev + 69.184 / 86400.0)
        rcv_ecef, _ = eci_to_ecef(state.r_eci, state.v_eci, mjd_utc)
        ep_data = compute_epoch_geometry(gps_sod, gps1b_raw, sp3, rcv_ecef, clk_data)
        if not ep_data:
            continue
        doy = 120
        state, stats = ekf.process_epoch(state, ep_data, sp3, sv_bias, sv_bias_ref, mjd_utc, mjd_tt, doy)
        r_ecef, _ = eci_to_ecef(state.r_eci, state.v_eci, mjd_utc)
        r_gnv = interpolate_ref(ref_orbit, gps_sod)
        if r_gnv is not None:
            epoch_rms.append({'hour': (gps_sod - GPS_SOD_START) / 3600.0,
                              'dr': np.linalg.norm(r_ecef - r_gnv)})

    # Compute 10-min window RMS
    for window_start_sec in range(0, int(SEGMENT_HOURS * 3600), 600):
        win_end = window_start_sec + 600
        win_dr = []
        for e in epoch_rms:
            sod = e['hour'] * 3600 - test_hour * 3600
            if window_start_sec <= sod < win_end:
                win_dr.append(e['dr'])
        if len(win_dr) >= 6:
            win_rms = np.sqrt(np.mean([d**2 for d in win_dr]))
            win_hour = test_hour + (window_start_sec + 300) / 3600.0
            results.append({'hour': win_hour, 'rms_3d': win_rms, 'n_epochs': len(win_dr)})

    # Also add overall RMS for this hour
    if epoch_rms:
        hour_rms = np.sqrt(np.mean([e['dr']**2 for e in epoch_rms]))
        print(f"    Overall 0.5h RMS: {hour_rms:.3f}m, {len(results)} total windows")

# ── Results ──
hours = [r['hour'] for r in results]
rms_vals = [r['rms_3d'] for r in results]

if len(results) >= 6:
    print(f"\n === Full-Day Summary ===")
    print(f"  Windows: {len(results)}")
    print(f"  Mean RMS:   {np.mean(rms_vals):.3f}m")
    print(f"  Median RMS: {np.median(rms_vals):.3f}m")
    print(f"  Best:       {np.min(rms_vals):.3f}m at {hours[np.argmin(rms_vals)]:.1f}h")
    print(f"  Worst:      {np.max(rms_vals):.3f}m at {hours[np.argmax(rms_vals)]:.1f}h")
    print(f"  <0.5m:      {sum(1 for v in rms_vals if v < 0.5)}/{len(rms_vals)}")
    print(f"  <0.8m:      {sum(1 for v in rms_vals if v < 0.8)}/{len(rms_vals)}")

    # ── Plot ──
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8),
                                    gridspec_kw={'height_ratios': [3, 1]})
    ax1.scatter(hours, rms_vals, c='#2196F3', s=50, alpha=0.8, edgecolors='#1565C0', linewidth=0.8, zorder=5)
    ax1.axhline(y=0.50, color='green', linestyle='--', alpha=0.4, label='0.50 m')
    ax1.axhline(y=1.00, color='orange', linestyle='--', alpha=0.4, label='1.00 m')
    ax1.axhline(y=np.median(rms_vals), color='red', linestyle='-', alpha=0.6,
                label=f'Median: {np.median(rms_vals):.3f} m')
    ax1.set_ylabel('3D RMS [m]', fontsize=12)
    ax1.set_xlabel('GPS Time of Day [hours]', fontsize=12)
    ax1.set_title(
        f'GRACE-FO C POD Accuracy — 2024-04-29\n'
        f'{len(results)} × 10-min windows | '
        f'Mean={np.mean(rms_vals):.3f}m  Median={np.median(rms_vals):.3f}m  '
        f'Best={np.min(rms_vals):.3f}m  Worst={np.max(rms_vals):.3f}m',
        fontsize=13)
    ax1.legend(fontsize=9, loc='upper right')
    ax1.set_xlim(-1, 25)
    ax1.grid(True, alpha=0.3)

    # Histogram
    ax2.hist(rms_vals, bins=15, color='#2196F3', alpha=0.7, edgecolor='#1565C0')
    ax2.axvline(x=np.median(rms_vals), color='red', linestyle='-', alpha=0.7,
                label=f'Median: {np.median(rms_vals):.3f}m')
    ax2.axvline(x=0.50, color='green', linestyle='--', alpha=0.4)
    ax2.set_xlabel('3D RMS [m]', fontsize=12)
    ax2.set_ylabel('Count', fontsize=12)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    out_png = Path('results') / 'accuracy_2024-04-29.png'
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(out_png), dpi=150, bbox_inches='tight')
    print(f"\n  Saved: {out_png}")

    # Save raw data
    pickle.dump({'hours': hours, 'rms_vals': rms_vals, 'results': results},
                open(str(out_png.with_suffix('.pkl')), 'wb'))
else:
    print(f"\n  ERROR: only {len(results)} valid windows found")
