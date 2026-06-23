"""Test 9-parameter GN loop (BatchOrbitLSQv3) on 2024-04-29 0.17h.

Key insight: GN force model must match EKF force model exactly.
We use the EKF's internal force_model closure for consistency.
"""
import sys, os, pickle, numpy as np
from pathlib import Path
from datetime import datetime, timedelta

ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'src'))

DATE_STR = "2024-04-29"; GRACE_ID = "C"; INTERVAL = 30; ARC_HOURS = 0.17
DATA_ROOT = ROOT / 'data'

print("=" * 60)
print(f"9-Param GN Loop Test: {DATE_STR} {ARC_HOURS}h")
print("=" * 60)

# ---- Load data ----
print("\nLoading data...")
gps1b = pickle.load(open(str(DATA_ROOT / f"GPS1B_{DATE_STR}_{GRACE_ID}_04.pkl"), "rb"))
sp3 = pickle.load(open(str(DATA_ROOT / "CODE/2024/COD0OPSFIN_20241200000_01D_05M_ORB.pkl"), "rb"))
from precision_products import *
clk_data = read_rinex_clk(str(DATA_ROOT / "CODE/2024/COD0OPSFIN_20241200000_01D_30S_CLK.CLK"))
antex = read_antex(str(DATA_ROOT / "igs14.atx"))
dcb_pair = load_code_dcb_pair(str(DATA_ROOT / "CODE/2024/P1C12404.DCB"),
                               str(DATA_ROOT / "CODE/2024/P1P22404.DCB"))
setup_iers_from_c04(str(DATA_ROOT / "IERS/eopc04_IAU2000.txt"))
from gravity_model import read_icgem_gfc
Cnm, Snm, Nmax_grav, GM_grav, R_grav = read_icgem_gfc(str(DATA_ROOT / 'gravity' / 'GGM05C.gfc'))
Nmax_grav = min(Nmax_grav, 90)
print(f"GPS1B={len(gps1b)}ep SP3={len(sp3['ts'])}ep grav Nmax={Nmax_grav}")

# ---- Time setup ----
gps_start = min(gps1b.keys())
gps_end = gps_start + int(ARC_HOURS * 3600)
epochs = sorted(set(g for g in gps1b.keys()
                    if gps_start <= g <= gps_end
                    and abs((g - gps_start) % INTERVAL) <= 2.0))
print(f"Epochs: {len(epochs)} ({epochs[0]} to {epochs[-1]})")

MJD_J2000 = 51544.5
mjd_start = MJD_J2000 + epochs[0] / 86400.0
mjd_tt_start = mjd_start + 69.184 / 86400.0
t_epochs = np.array([(g - epochs[0]) for g in epochs], dtype=float)

# ---- Reference orbit ----
from run_sequential_pod import load_gnv1b, compute_epoch_geometry, interpolate_ref
gnv_path = DATA_ROOT / 'gracefo' / '2024' / DATE_STR / f'GNV1B_{DATE_STR}_{GRACE_ID}_04.txt'
ref_orbit, ref_vel = load_gnv1b(str(gnv_path))

r0_ecef = interpolate_ref(ref_orbit, epochs[0])
v0_ecef = interpolate_ref(ref_vel, epochs[0])
from coordinates import ecef_to_eci, eci_to_ecef
r0_eci, v0_eci = ecef_to_eci(r0_ecef, v0_ecef, mjd_start)
print(f"r0_eci: [{r0_eci[0]:.1f}, {r0_eci[1]:.1f}, {r0_eci[2]:.1f}]")

# ---- Per-SV code bias ----
N_BIAS = min(60, len(epochs))
sv_p_res = {}
J2000 = datetime(2000, 1, 1, 12, 0, 0)
for gps_sod in epochs[:N_BIAS]:
    utc_dt = J2000 + timedelta(seconds=gps_sod)
    recs = gps1b.get(int(gps_sod), gps1b.get(gps_sod, {}))
    from run_sequential_pod import get_sat_geometry
    for sv_id, rec in recs.items():
        if 'P_if' not in rec: continue
        rcv_ref = interpolate_ref(ref_orbit, gps_sod)
        if rcv_ref is None: continue
        sat_pos, sat_clk, rho_corr = get_sat_geometry(sp3, sv_id, utc_dt, rcv_ref, clk_data)
        if sat_pos is None: continue
        dcb_corr = compute_dcb_if_correction(dcb_pair, sv_id)
        P_r = float(rec['P_if']) + dcb_corr + sat_clk - rho_corr
        sv_p_res.setdefault(sv_id, []).append(P_r)
sv_bias = {}
for sv, vals in sv_p_res.items():
    if len(vals) >= 3: sv_bias[sv] = float(np.median(vals))
sv_bias_ref = float(np.mean(list(sv_bias.values()))) if sv_bias else 0.0
for sv in sv_bias: sv_bias[sv] -= sv_bias_ref
print(f"SV biases: {len(sv_bias)} SVs")

# ---- Run EKF pass 1 ----
print("\nRunning EKF pass 1...")
from sequential_filter import SequentialEKF, I_AMB_START

# Build EKF with simplified dynamics (same force model as GN)
ekf_cfg = {
    'dynamics_mode': 'simplified',
    'Cd': 2.2, 'CR': 1.3, 'area_drag': 0.68, 'area_srp': 3.4, 'mass': 580.0,
    'bodies': ['Sun', 'Moon'], 'Cnm': Cnm, 'Snm': Snm,
    'GM_grav': GM_grav, 'R_grav': R_grav, 'gravity_nmax': Nmax_grav,
    'sigma_acc_process': 1e-3, 'tau_emp': 600.0, 'sigma_emp_ss': 1e-8,
    'sigma_zwd_rw': 1e-9, 'sigma_phase': 0.20, 'sigma_code': 0.30,
    'chi2_threshold': 25, 'el_min': 0.087,
    'use_phase_windup': True, 'use_relativity': True, 'use_cycle_slip': False,
    'ar_min_epochs': 6, 'antex_data': antex, 'dcb_data': dcb_pair,
    'elev_exp_phase': 1.0, 'elev_exp_code': 1.0,
    'clock_rw': 0.0004, 'mw_max_epochs': 200,
}

ekf = SequentialEKF(ekf_cfg)
state = ekf.initialize(r0_eci, v0_eci, mjd_start, epochs[0])
dr_vals = []

SEC_PER_DAY = 86400.0; C_LIGHT = 299792458.0; OMEGA_E = 7.2921151467e-5
pass1_geometry = []

from src.troposphere import saastamoinen_zhd

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

    # Store geometry (ECEF-based, same as EKF uses internally)
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

ekf_rms = np.sqrt(np.mean([d**2 for d in dr_vals])) if dr_vals else 0
print(f"EKF 3D RMS: {ekf_rms:.3f}m  phase RMS: {stats.get('rms_phase', 0):.3f}m")
print(f"Collected {sum(len(ep) for ep in pass1_geometry)} obs in {len(pass1_geometry)} epochs")

# ---- Baseline: BatchLinearSolver on EKF orbit ----
print("\nBaseline: BatchLinearSolver on EKF orbit...")
from src.batch_solver import BatchLinearSolver
bls = BatchLinearSolver(pass1_geometry, sigma_phase=0.20, sigma_code=0.30)
bls_sol = bls.solve()
print(f"  Phase RMS: {bls_sol['rms_phase']:.3f}m  Code RMS: {bls_sol['rms_code']:.3f}m")

# ---- Build GN force function matching EKF ----
# The EKF uses SequentialEKF's internal force_model closure with total_acc_eci,
# solid tide corrections, and all gravity parameters.
# We replicate this exactly.
print("\nBuilding GN force function (matching EKF simplified dynamics)...")
from src.orbit_dynamics import total_acc_eci
from src.orbit_integrator import integrate_orbit_eci_with_stm

# Compute tide corrections (same as EKF)
tide_corr = {}
try:
    from src.solid_tides import (compute_solid_tide_corrections,
                                 compute_time_varying_gravity, merge_tide_corrections)
    # Solid tides at the arc start
    tide_corr = compute_solid_tide_corrections(mjd_start, mjd_tt_start)
    # Time-varying gravity
    tvgrav = compute_time_varying_gravity(mjd_tt_start)
    tide_corr = merge_tide_corrections(tide_corr, tvgrav)
    print(f"  Tide corrections: loaded")
except Exception as e:
    print(f"  Tide corrections: unavailable ({e})")

def gn_force_fn(pos_eci, vel_eci, CD=2.2, CR=1.3,
                 area_drag=0.68, area_srp=3.4, mass=580.0,
                 empirical_acc_rtn=None, bodies=None, mjd_utc=None, mjd_tt=None,
                 **kwargs):
    """Force function matching EKF simplified dynamics exactly."""
    return total_acc_eci(
        pos_eci, vel_eci,
        mjd_tt=mjd_tt if mjd_tt is not None else mjd_tt_start,
        mjd_utc=mjd_utc if mjd_utc is not None else mjd_start,
        Cnm=Cnm, Snm=Snm, Nmax=Nmax_grav,
        CD=CD, CR=CR,
        area_drag=area_drag, area_srp=area_srp, mass=mass,
        empirical_acc_rtn=empirical_acc_rtn,
        tide_corrections=tide_corr,
        bodies=bodies if bodies else ['Sun', 'Moon'],
        GM_gravity=GM_grav, R_gravity=R_grav)

# ---- Compare EKF vs GN propagated orbits (sanity check) ----
print("\nComparing EKF vs GN propagated orbits...")
integ_gn = integrate_orbit_eci_with_stm(
    r0_eci, v0_eci, (0, max(t_epochs) + 30), gn_force_fn,
    Cd=2.2, CR=1.3, area_drag=0.68, area_srp=3.4, mass=580.0,
    empirical_acc_rtn=np.zeros(3), param_names=['aR','aT','aN'], dt=10.0,
    mjd_utc=mjd_start, mjd_tt=mjd_tt_start, bodies=['Sun','Moon'])

# Get EKF initial propagation (pre-fit, from GNV1B prior)
# The EKF's first predict() propagates from GNV1B initial state
# We can compare endpoints
for label, r, v in [("GNV1B", r0_eci, v0_eci),
                      ("GN-prop", integ_gn['r'][-1], integ_gn['v'][-1])]:
    print(f"  {label}: r=[{r[0]:.1f},{r[1]:.1f},{r[2]:.1f}] v=[{v[0]:.3f},{v[1]:.3f},{v[2]:.3f}]")

# Compare GN propagated end position vs. EKF state at last epoch
ekf_last = state.r_eci
t_end = max(t_epochs)
idx_end = np.argmin(np.abs(integ_gn['t'] - t_end))
gn_end = integ_gn['r'][idx_end]
dr_ekf_gn = np.linalg.norm(ekf_last - gn_end)
print(f"  EKF end vs GN end: dr={dr_ekf_gn:.3f}m")
if dr_ekf_gn > 5.0:
    print(f"  WARNING: Large EKF-GN mismatch! Force models may differ.")
    print(f"  EKF end: [{ekf_last[0]:.1f},{ekf_last[1]:.1f},{ekf_last[2]:.1f}]")
    print(f"  GN  end: [{gn_end[0]:.1f},{gn_end[1]:.1f},{gn_end[2]:.1f}]")

# ---- Run 9-parameter GN loop ----
print("\n" + "=" * 60)
print("Running 9-parameter GN loop...")
print("=" * 60)

from src.batch_orbit_v3 import BatchOrbitLSQv3

gn_solver = BatchOrbitLSQv3(
    pass1_geometry, gn_force_fn, t_epochs,
    mjd_utc_start=mjd_start, mjd_tt_start=mjd_tt_start,
    sigma_phase=0.20, sigma_code=0.30,
    max_iter=8,
    dx_tol_pos=0.002, dx_tol_vel=0.0002, dx_tol_emp=1e-10,
    damping=0.5, prior_r0=1.0, prior_v0=0.01, prior_emp=1e-7,
)

sol = gn_solver.solve(r0_eci, v0_eci)

# ---- Evaluate GN orbit vs GNV1B ----
print("\n" + "=" * 60)
print("Results")
print("=" * 60)
print(f"Converged: {sol['converged']} in {sol['iterations']} iterations")
print(f"r0 delta: {np.linalg.norm(sol['r0'] - r0_eci):.3f}m")
print(f"v0 delta: {np.linalg.norm(sol['v0'] - v0_eci):.6f}m/s")
a_e = sol['a_emp']
print(f"a_emp: [{a_e[0]:.2e}, {a_e[1]:.2e}, {a_e[2]:.2e}] m/s2")
print(f"Batch Phase RMS: {sol['rms_phase']:.3f}m  Code RMS: {sol['rms_code']:.3f}m")

# 3D RMS vs GNV1B
gn_dr = []
for i_ep, t_ep in enumerate(t_epochs):
    gps_sod = epochs[i_ep]
    r_gnv = interpolate_ref(ref_orbit, gps_sod)
    if r_gnv is not None:
        mjd_utc = MJD_J2000 + gps_sod / 86400.0
        r_ecef, _ = eci_to_ecef(sol['r_eci'][i_ep], sol['v_eci'][i_ep], mjd_utc)
        gn_dr.append(np.linalg.norm(r_ecef - r_gnv))
gn_rms = np.sqrt(np.mean([d**2 for d in gn_dr])) if gn_dr else float('inf')

print(f"\nEKF 3D RMS: {ekf_rms:.3f}m (vs GNV1B)")
print(f"GN  3D RMS: {gn_rms:.3f}m (vs GNV1B)")
if ekf_rms > 0:
    change = (gn_rms - ekf_rms) / ekf_rms * 100
    print(f"Change: {change:+.1f}%")

print(f"\nBaseline Batch Phase RMS (EKF orbit): {bls_sol['rms_phase']:.3f}m")
print(f"GN      Batch Phase RMS (GN orbit):   {sol['rms_phase']:.3f}m")
