"""Test Orekit GN -- fully operational after Orekit v13 API fixes.

Uses Orekit full dynamics (Nmax=150 gravity + tides + drag + SRP)
with FD STM for the GN outer loop Jacobian.
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
print(f"Orekit GN Loop Test: {DATE_STR} {ARC_HOURS}h")
print("=" * 60)

# ---- Check Orekit availability ----
from src.orekit_bridge import is_orekit_available, get_orekit_error
if not is_orekit_available():
    print(f"Orekit NOT available: {get_orekit_error()}")
    sys.exit(1)
print("Orekit available.")

# ---- Create Orekit propagator ----
print("\nSetting up Orekit propagator...")
from src.orekit_bridge import OrekitPropagator
orekit_prop = OrekitPropagator(
    gravity_field=str(DATA_ROOT / 'gravity' / 'GGM05C.gfc'),
    gravity_degree=150,
    solid_tides=True, ocean_tides=True, ocean_tide_degree=50,
    third_body='lunisolar', srp_model='isotropic', relativity=True,
    drag_model='exponential',
    mass=580.0, area_drag=0.68, area_srp=3.4, CR=1.3, CD=2.2,
    stm_perturb=1.0, integrator_tol=1e-12,
)
print("Orekit propagator created.")

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
grav_nmax = min(Nmax_grav, 90)
print(f"GPS1B={len(gps1b)}ep SP3={len(sp3['ts'])}ep grav Nmax={grav_nmax}")

# ---- Time ----
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

# ---- Run EKF pass 1 to collect geometry ----
from src.troposphere import saastamoinen_zhd
from sequential_filter import SequentialEKF

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
        rcv_ref = interpolate_ref(ref_orbit, gps_sod)
        if rcv_ref is None: continue
        sat_pos, sat_clk, rho_corr = get_sat_geometry(sp3, sv_id, utc_dt, rcv_ref, clk_data)
        if sat_pos is None: continue
        dcb_c = compute_dcb_if_correction(dcb_pair, sv_id)
        P_r = float(rec['P_if']) + dcb_c + sat_clk - rho_corr
        sv_p_res.setdefault(sv_id, []).append(P_r)
sv_bias = {}; sv_bias_ref = 0.0
for sv, vals in sv_p_res.items():
    if len(vals) >= 3: sv_bias[sv] = float(np.median(vals))
if sv_bias:
    sv_bias_ref = float(np.mean(list(sv_bias.values())))
    for sv in sv_bias: sv_bias[sv] -= sv_bias_ref

print("Running EKF pass 1...")
ekf_cfg = {
    'dynamics_mode': 'simplified', 'Cd': 2.2, 'CR': 1.3,
    'area_drag': 0.68, 'area_srp': 3.4, 'mass': 580.0,
    'bodies': ['Sun', 'Moon'], 'Cnm': Cnm, 'Snm': Snm,
    'GM_grav': GM_grav, 'R_grav': R_grav, 'gravity_nmax': grav_nmax,
    'sigma_acc_process': 1e-3, 'tau_emp': 600.0, 'sigma_emp_ss': 1e-8,
    'sigma_zwd_rw': 1e-9, 'sigma_phase': 0.20, 'sigma_code': 0.30,
    'chi2_threshold': 100 if ARC_HOURS >= 0.3 else 25, 'el_min': 0.087,
    'use_phase_windup': True, 'use_relativity': True, 'use_cycle_slip': False,
    'ar_min_epochs': 6, 'antex_data': antex, 'dcb_data': dcb_pair,
    'elev_exp_phase': 1.0, 'elev_exp_code': 0.70 if ARC_HOURS >= 0.3 else 1.0,
    'clock_rw': 0.001 if ARC_HOURS >= 0.3 else 0.0004, 'mw_max_epochs': 200,
}
ekf = SequentialEKF(ekf_cfg)
state = ekf.initialize(r0_eci, v0_eci, mjd_start, epochs[0])
SEC_PER_DAY = 86400.0; C_LIGHT = 299792458.0; OMEGA_E = 7.2921151467e-5
pass1_geometry = []

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
    lat_rad = np.arcsin(rcv_ecef[2] / np.linalg.norm(rcv_ecef))
    h_m = np.linalg.norm(rcv_ecef) - 6378137.0; zhd = saastamoinen_zhd(lat_rad, h_m)
    for d in ep_data:
        sv = d['sv']; se = np.asarray(d['sat_pos'], dtype=float)
        sc = float(d.get('sat_clk', 0)); el = float(d.get('el', 0.5))
        rho = np.linalg.norm(se - rcv_ecef)
        sag = (OMEGA_E/C_LIGHT)*(se[0]*rcv_ecef[1] - se[1]*rcv_ecef[0])
        mf = 1.001/np.sqrt(0.002001+np.sin(el)**2)
        dcb_c = compute_dcb_if_correction(dcb_pair, sv)
        d['_geo_full'] = rho + sag - sc + zhd * mf
        d['_obs_code'] = float(d.get('P_if_raw', 0)) + dcb_c - sv_bias.get(sv, 0.0)
        d['_obs_phase'] = float(d.get('L_if_raw', 0)) - sv_bias.get(sv, 0.0)
    pass1_geometry.append(ep_data)

ekf_rms = stats.get('rms_3d', 0.0)
print(f"EKF 3D RMS: {ekf_rms:.3f}m  phase RMS: {stats.get('rms_phase', 0):.3f}m")
print(f"Collected {sum(len(ep) for ep in pass1_geometry)} obs in {len(pass1_geometry)} epochs")

# ---- Build Python force function (for fallback) ----
from src.orbit_dynamics import total_acc_eci
try:
    from src.solid_tides import (compute_solid_tide_corrections,
                                 compute_time_varying_gravity, merge_tide_corrections)
    tide_corr = compute_solid_tide_corrections(mjd_start, mjd_tt_start)
    tvgrav = compute_time_varying_gravity(mjd_tt_start)
    tide_corr = merge_tide_corrections(tide_corr, tvgrav)
except Exception as e:
    tide_corr = {}
    print(f"  Tide corrections unavailable: {e}")

_G = {'Cnm': Cnm, 'Snm': Snm, 'Nmax': grav_nmax, 'GM': GM_grav, 'R': R_grav}

def gn_force_fn(pos_eci, vel_eci, CD=2.2, CR=1.3,
                 area_drag=0.68, area_srp=3.4, mass=580.0,
                 empirical_acc_rtn=None, bodies=None, mjd_utc=None, mjd_tt=None,
                 **kwargs):
    return total_acc_eci(
        pos_eci, vel_eci,
        mjd_tt=mjd_tt if mjd_tt is not None else mjd_tt_start,
        mjd_utc=mjd_utc if mjd_utc is not None else mjd_start,
        Cnm=_G['Cnm'], Snm=_G['Snm'], Nmax=_G['Nmax'],
        CD=CD, CR=CR, area_drag=area_drag, area_srp=area_srp, mass=mass,
        empirical_acc_rtn=empirical_acc_rtn,
        tide_corrections=tide_corr,
        bodies=bodies if bodies else ['Sun', 'Moon'],
        GM_gravity=_G['GM'], R_gravity=_G['R'])

# ---- Run GN loop with Orekit (MSISE00 drag) ----
print("\n" + "=" * 60)
print("GN Loop: Orekit (exponential drag)")
print("=" * 60)
from src.batch_orbit_v3 import BatchOrbitLSQv3

gn_ok = BatchOrbitLSQv3(
    pass1_geometry, gn_force_fn, t_epochs,
    mjd_utc_start=mjd_start, mjd_tt_start=mjd_tt_start,
    sigma_phase=0.20, sigma_code=0.30,
    max_iter=6,
    prior_r0=1.0, prior_v0=0.01, prior_emp=1e-7,
    damping=0.5,
    orekit_prop=orekit_prop,
    estimate_cd_cr=False,
)
import time
t0 = time.time()
sol_ok = gn_ok.solve(r0_eci, v0_eci)
dt_gn = time.time() - t0
print(f"  [TIMING] GN solve: {dt_gn:.1f}s")

# ---- GNV1B 3D RMS ----
def compute_3d_rms(r_eci_arc, v_eci_arc, epochs, ref_orbit):
    dr = []
    for i_ep, gps_sod in enumerate(epochs):
        r_gnv = interpolate_ref(ref_orbit, gps_sod)
        if r_gnv is not None:
            mjd_utc = MJD_J2000 + gps_sod / 86400.0
            r_ecef, _ = eci_to_ecef(r_eci_arc[i_ep], v_eci_arc[i_ep], mjd_utc)
            dr.append(np.linalg.norm(r_ecef - r_gnv))
    return np.sqrt(np.mean([d**2 for d in dr])) if dr else 0

rms_ok = compute_3d_rms(sol_ok['r_eci'], sol_ok['v_eci'], epochs, ref_orbit)

print("\n" + "=" * 60)
print("Results")
print("=" * 60)
print(f"  Orekit GN  3D RMS: {rms_ok:.3f}m  Phase: {sol_ok['rms_phase']:.3f}m")
print(f"  Converged: {sol_ok['converged']} ({sol_ok['iterations']} iter)")
print(f"  r0: [{sol_ok['r0'][0]:.1f}, {sol_ok['r0'][1]:.1f}, {sol_ok['r0'][2]:.1f}]")
print(f"  v0: [{sol_ok['v0'][0]:.3f}, {sol_ok['v0'][1]:.3f}, {sol_ok['v0'][2]:.3f}]")

if sol_ok.get('a_emp') is not None:
    N_segs = sol_ok.get('N_segs', 1)
    if N_segs == 1:
        print(f"  a_emp: [{sol_ok['a_emp'][0]:.2e}, "
              f"{sol_ok['a_emp'][1]:.2e}, {sol_ok['a_emp'][2]:.2e}]")
    else:
        a_e = sol_ok['a_emp']
        for s in range(N_segs):
            b = s * 3
            print(f"  a_emp[{s}]: [{a_e[b]:.2e}, {a_e[b+1]:.2e}, {a_e[b+2]:.2e}]")
if sol_ok.get('Cd') is not None:
    print(f"  Cd: {sol_ok['Cd']:.4f}  CR: {sol_ok['CR']:.4f}")
