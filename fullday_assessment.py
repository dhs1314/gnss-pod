"""Full-day accuracy assessment: 10-minute segments across 24 hours.

Processes 10-minute arc segments contiguously for the full day,
computing 3D RMS vs GNV1B per segment. Produces a PNG plot.
"""
import sys, os, pickle, math
from pathlib import Path
from datetime import datetime, timedelta
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

# ── Load data once ──
DP = Path(r"d:\prj\gnss_pod\data")
DATE_STR = "2024-04-29"
GRACE_ID = "C"
y, m, d = 2024, 4, 29

print("Loading data...")
gps1b_raw = pickle.load(open(str(DP / f"GPS1B_{DATE_STR}_{GRACE_ID}_04.pkl"), "rb"))
ref_orbit, ref_vel = None, None
gnv_path = DP / 'gracefo' / str(y) / DATE_STR / f'GNV1B_{DATE_STR}_{GRACE_ID}_04.txt'
from run_sequential_pod import load_gnv1b, get_sat_geometry
ref_orbit, ref_vel = load_gnv1b(str(gnv_path))
print(f"  GNV1B: {len(ref_orbit)} pos, {len(ref_vel)} vel")

# Load SP3 from pickle cache (pre-parsed)
sp3 = pickle.load(open(str(DP / "CODE/2024/COD0OPSFIN_20241200000_01D_05M_ORB.pkl"), "rb"))
print(f"  SP3: {len(sp3['ts'])} epochs")

# Load CLK
from src.precision_products import read_rinex_clk
clk_data = read_rinex_clk(str(DP / "CODE/2024/COD0OPSFIN_20241200000_01D_30S_CLK.CLK"))
print(f"  CLK: {len(clk_data)} SVs")

# Load ANTEX+DCB+IERS
from src.precision_products import read_antex, load_code_dcb_pair, compute_dcb_if_correction
antex = read_antex(str(DP / "igs14.atx"))
dcb_pair = load_code_dcb_pair(str(DP / "CODE/2024/P1C12404.DCB"),
                               str(DP / "CODE/2024/P1P22404.DCB"))
dcb_if = {}
for prn in dcb_pair:
    dcb_if[prn] = compute_dcb_if_correction(dcb_pair, prn)
from src.precision_products import setup_iers_from_c04
setup_iers_from_c04(str(DP / "IERS/eopc04_IAU2000.txt"))
print("  ANTEX+DCB+IERS loaded")

# Load gravity
from src.gravity_model import read_icgem_gfc
gfc_path = DP / 'gravity' / 'GGM05C.gfc'
Cnm, Snm, Nmax, GM_grav, R_grav = read_icgem_gfc(str(gfc_path))
Nmax = min(Nmax, 90)
print(f"  Gravity: Nmax={Nmax}")

# ── Set up ──
from src.coordinates import ecef_to_eci, eci_to_ecef
from src.sequential_filter import SequentialEKF
from run_sequential_pod import compute_epoch_geometry, interpolate_ref

GPS_SOD_START = 767620800  # 2024-04-29 00:00:00 GPS time
GPS_SOD_END = GPS_SOD_START + 86400
SEGMENT_DT = 600  # 10 minutes
INTERVAL = 30  # 30s EKF step

J2000 = datetime(2000, 1, 1, 12, 0, 0)
MJD_J2000 = 51544.5
GPS_UTC_OFFSET = 0  # simplified

MJD_START = 60429 + GPS_SOD_START / 86400.0  # approx

# ── EKF config ──
ekf_config = {
    'Cd': 2.2, 'CR': 1.3, 'area_drag': 0.68, 'area_srp': 3.4, 'mass': 580.0,
    'dt_integ': 10.0, 'bodies': ['Sun', 'Moon'],
    'Cnm': Cnm, 'Snm': Snm, 'GM_grav': GM_grav, 'R_grav': R_grav,
    'gravity_nmax': Nmax,
    'sigma_acc_process': 1e-3, 'tau_emp': 600.0, 'sigma_emp_ss': 1e-8,
    'sigma_zwd_rw': 1e-9, 'sigma_phase': 0.20, 'sigma_code': 0.30,
    'chi2_threshold': 100.0, 'el_min': 0.087,
    'use_phase_windup': True, 'use_relativity': True, 'use_cycle_slip': False,
    'ar_min_epochs': 6,
    'antex_data': antex,
    'dcb_data': dcb_pair,
    'elev_exp_phase': 1.0,
    'elev_exp_code': 0.70,  # 10min is short arc, but safe for all
}

# Pre-compute per-SV code biases (first 60 epochs for reference)
from run_sequential_pod import compute_epoch_geometry, interpolate_ref
sv_bias = {}
sv_bias_ref = 0.0

print("\nProcessing segments...")
results_all = []  # list of (gps_sod_mid, rms_3d, n_epochs, n_sv, n_rej)

# Process in 10-min windows
for seg_start in range(GPS_SOD_START, GPS_SOD_END, SEGMENT_DT):
    seg_end = seg_start + SEGMENT_DT  # exclusive

    # Collect epochs in this segment
    segment_epochs = []
    for gps_sod in sorted(gps1b_raw.keys()):
        if seg_start <= gps_sod < seg_end:
            dt_ep = gps_sod - seg_start
            nearest = round(dt_ep / INTERVAL) * INTERVAL
            if abs(dt_ep - nearest) <= 2.0:
                segment_epochs.append(gps_sod)
    segment_epochs = sorted(set(segment_epochs))

    if len(segment_epochs) < 6:
        continue  # too few epochs

    # Estimate per-SV code biases from first 30 epochs of this segment
    sv_code_res = {}
    for gps_sod in segment_epochs[:30]:
        ref_pos = interpolate_ref(ref_orbit, gps_sod)
        if ref_pos is None:
            continue
        utc_dt = J2000 + timedelta(seconds=gps_sod)
        from run_sequential_pod import get_sat_geometry
        recs = gps1b_raw.get(gps_sod, {})
        for sv_id, rec in recs.items():
            if 'L_if' not in rec:
                continue
            sat_pos, sat_clk, rho_corr = get_sat_geometry(sp3, sv_id, utc_dt, ref_pos, clk_data)
            if sat_pos is None:
                continue
            rho_tmp = float(np.linalg.norm(sat_pos - ref_pos))
            sag_tmp = 0  # simplify
            P_r = float(rec['P_if']) + dcb_if.get(sv_id, 0.0) + sat_clk - rho_tmp
            sv_code_res.setdefault(sv_id, []).append(P_r)

    sv_bias_local = {}
    for sv_id, vals in sv_code_res.items():
        if len(vals) >= 3:
            sv_bias_local[sv_id] = float(np.median(vals))
    if sv_bias_local:
        ref = float(np.mean(list(sv_bias_local.values())))
        for sv in sv_bias_local:
            sv_bias_local[sv] -= ref
        sv_bias_ref_local = ref
    else:
        sv_bias_ref_local = 0.0

    # Initialize EKF
    r0_ecef = interpolate_ref(ref_orbit, segment_epochs[0])
    v0_ecef = interpolate_ref(ref_vel, segment_epochs[0])
    if r0_ecef is None or v0_ecef is None:
        continue
    mjd_start = MJD_J2000 + (segment_epochs[0] - GPS_UTC_OFFSET) / 86400.0
    r0_eci, v0_eci = ecef_to_eci(r0_ecef, v0_ecef, mjd_start)

    ekf = SequentialEKF(ekf_config)
    state = ekf.initialize(r0_eci, v0_eci, mjd_start, segment_epochs[0])
    segment_dr = []
    n_rej_total = 0

    for i_ep, gps_sod in enumerate(segment_epochs):
        mjd_utc = MJD_J2000 + (gps_sod - GPS_UTC_OFFSET) / 86400.0
        mjd_tt = mjd_utc + 69.184 / 86400.0

        if i_ep > 0:
            dt = gps_sod - segment_epochs[i_ep - 1]
            mjd_utc_prev = MJD_J2000 + (segment_epochs[i_ep - 1] - GPS_UTC_OFFSET) / 86400.0
            mjd_tt_prev = mjd_utc_prev + 69.184 / 86400.0
            state = ekf.predict(state, gps_sod, mjd_utc_prev, mjd_tt_prev)

        rcv_ecef, _ = eci_to_ecef(state.r_eci, state.v_eci, mjd_utc)
        ep_data = compute_epoch_geometry(gps_sod, gps1b_raw, sp3, rcv_ecef, clk_data)

        if not ep_data:
            continue

        doy = (datetime(y, m, d) - datetime(y, 1, 1)).days + 1
        state, stats = ekf.process_epoch(state, ep_data, sp3, sv_bias_local,
                                          sv_bias_ref_local, mjd_utc, mjd_tt, doy)

        r_ecef, _ = eci_to_ecef(state.r_eci, state.v_eci, mjd_utc)
        r_gnv = interpolate_ref(ref_orbit, gps_sod)
        if r_gnv is not None:
            segment_dr.append(np.linalg.norm(r_ecef - r_gnv))
        n_rej_total += stats['n_rej']

    if len(segment_dr) >= 6:
        rms_3d = np.sqrt(np.mean([d**2 for d in segment_dr]))
        gps_mid = seg_start + SEGMENT_DT / 2
        hour = (gps_mid - GPS_SOD_START) / 3600.0
        results_all.append({
            'hour': hour,
            'rms_3d': rms_3d,
            'n_epochs': len(segment_dr),
            'n_sv': state.n_sv,
            'n_rej': n_rej_total,
        })

    seg_idx = len(results_all)
    if seg_idx % 15 == 0 or seg_idx <= 3:
        rms_str = f"{rms_3d:.3f}m" if 'rms_3d' in locals() and segment_dr else "skip"
        print(f"  [{seg_idx:3d}] hour={hour:.2f}  RMS={rms_str}  "
              f"n_ep={len(segment_dr)}  n_sv={state.n_sv}  rej={n_rej_total}")

# ── Save results ──
out_dir = Path("results") / "fullday"
out_dir.mkdir(parents=True, exist_ok=True)
pickle.dump(results_all, open(str(out_dir / f"fullday_2024-04-29.pkl"), "wb"))

# ── Print summary ──
rms_vals = [r['rms_3d'] for r in results_all]
hours = [r['hour'] for r in results_all]
print(f"\n=== Full-Day Summary ===")
print(f"  Segments processed: {len(results_all)}")
print(f"  Best:  {min(rms_vals):.3f}m at {hours[np.argmin(rms_vals)]:.1f}h")
print(f"  Worst: {max(rms_vals):.3f}m at {hours[np.argmax(rms_vals)]:.1f}h")
print(f"  Mean:  {np.mean(rms_vals):.3f}m")
print(f"  Median:{np.median(rms_vals):.3f}m")
print(f"  Std:   {np.std(rms_vals):.3f}m")
print(f"  <0.5m segments: {sum(1 for v in rms_vals if v < 0.5)}/{len(rms_vals)}")
print(f"  <0.8m segments: {sum(1 for v in rms_vals if v < 0.8)}/{len(rms_vals)}")
