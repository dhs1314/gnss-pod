"""Debug script: compare ECI and ECEF force model accuracy."""
import numpy as np
from src.orbit_dynamics import total_acc, total_acc_eci, two_body_acc, j2_acc, drag_acc
from src.coordinates import ecef_to_eci, eci_to_ecef
from src.gravity_model import read_icgem_gfc, compute_gravity_acceleration
from src.orbit_integrator import rk4_step_eci, rk4_step
from src.third_body import compute_total_third_body
from pathlib import Path
from datetime import datetime

SEC_PER_DAY = 86400.0
MJD_J2000 = 51544.5
OMEGA_E = 7.2921151467e-5
date_str = '2024-04-29'

# Load GNV1B
dp = Path('data')
y, m, d = 2024, 4, 29
gnv_path = dp / 'gracefo' / str(y) / date_str / 'GNV1B_2024-04-29_C_04.txt'
ref_orbit = {}
with open(str(gnv_path)) as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        try:
            t = float(parts[0])
            flag = parts[2]
            if flag not in ('C', 'E'):
                continue
            X, Y, Z = float(parts[3]), float(parts[4]), float(parts[5])
            if abs(X) < 1e3:
                continue
            ref_orbit[t] = np.array([X, Y, Z])
        except:
            continue

ts = sorted(ref_orbit.keys())
t0 = ts[0]


def interp(t):
    for i, ti in enumerate(ts):
        if ti >= t:
            if i == 0:
                return ref_orbit[ts[0]]
            a = (t - ts[i - 1]) / (ts[i] - ts[i - 1])
            return ref_orbit[ts[i - 1]] * (1 - a) + ref_orbit[ts[i]] * a
    return ref_orbit[ts[-1]]


r0 = interp(t0)
v_avg = (interp(t0 + 30) - r0) / 30.0
r_mid = 0.5 * (r0 + interp(t0 + 30))
a_mid = total_acc(r_mid, v_avg, 2.2)
v0_ecef = v_avg - 0.5 * a_mid * 30.0
mjd_utc = MJD_J2000 + t0 / SEC_PER_DAY

r0_eci, v0_eci = ecef_to_eci(r0, v0_ecef, mjd_utc)

Cnm, Snm, _, GM_grav, R_grav = read_icgem_gfc('data/gravity/GGM05C.gfc')

# --- Compare gravity models in ECEF ---
print('=== ECEF gravity comparison ===')
a_2b = two_body_acc(r0)
a_j2 = j2_acc(r0)
a_ggm = compute_gravity_acceleration(r0, Cnm, Snm, 90, GM=GM_grav, R=R_grav)
print(f'two-body:  [{a_2b[0]:.4f}, {a_2b[1]:.4f}, {a_2b[2]:.4f}] |a|={np.linalg.norm(a_2b):.4f}')
print(f'J2:        [{a_j2[0]:.4f}, {a_j2[1]:.4f}, {a_j2[2]:.4f}] |a|={np.linalg.norm(a_j2):.4f}')
print(f'GGM05C 90: [{a_ggm[0]:.4f}, {a_ggm[1]:.4f}, {a_ggm[2]:.4f}] |a|={np.linalg.norm(a_ggm):.4f}')
print(f'  GGM - (2b+J2): [{a_ggm[0]-a_2b[0]-a_j2[0]:.6f}, {a_ggm[1]-a_2b[1]-a_j2[1]:.6f}, {a_ggm[2]-a_2b[2]-a_j2[2]:.6f}]')

# --- Compare ECI force model components ---
print()
print('=== ECI force model breakdown ===')
mjd_tt = mjd_utc + 69.184 / SEC_PER_DAY
a_total_eci = total_acc_eci(r0_eci, v0_eci, mjd_tt, mjd_utc, Cnm, Snm, 90,
                            CD=0, CR=0, area_drag=0, area_srp=0, mass=580,
                            GM_gravity=GM_grav, R_gravity=R_grav)
print(f'total_acc_eci (no drag/SRP): |a|={np.linalg.norm(a_total_eci):.4f}')

a_ecef_2b_j2 = a_2b + a_j2
a_ecef_2b_j2_eci, _ = ecef_to_eci(a_ecef_2b_j2, np.zeros(3), mjd_utc)
print(f'ECEF 2b+J2 in ECI: |a|={np.linalg.norm(a_ecef_2b_j2_eci):.4f}')
print(f'  diff GGM-ECI vs 2b+J2-ECI: {np.linalg.norm(a_total_eci - a_ecef_2b_j2_eci):.6f}')

# --- Check third-body ---
a_3b = compute_total_third_body(r0_eci, mjd_tt, ['Sun', 'Moon'])
print(f'Third-body: [{a_3b[0]:.6f}, {a_3b[1]:.6f}, {a_3b[2]:.6f}] |a|={np.linalg.norm(a_3b):.6f}')

# --- Multi-step prediction test ---
print()
print('=== Multi-step prediction test (ECI, GGM05C 90x90) ===')
r_cur, v_cur = r0.copy(), v0_ecef.copy()
dt_s = 30.0
for step in range(5):
    t_cur = t0 + step * dt_s
    t_next = t_cur + dt_s
    mjd_cur = MJD_J2000 + t_cur / SEC_PER_DAY
    mjd_tt_cur = mjd_cur + 69.184 / SEC_PER_DAY

    r_eci, v_eci = ecef_to_eci(r_cur, v_cur, mjd_cur)
    r_pred_eci, v_pred_eci = rk4_step_eci(
        r_eci, v_eci, dt_s, total_acc_eci,
        mjd_tt=mjd_tt_cur, mjd_utc=mjd_cur,
        Cnm=Cnm, Snm=Snm, Nmax=90,
        CD=2.2, CR=1.3,
        area_drag=0.68, area_srp=3.4, mass=580.0,
        empirical_acc_rtn=None,
        bodies=['Sun', 'Moon'],
        GM_gravity=GM_grav, R_gravity=R_grav,
    )
    mjd_next = MJD_J2000 + t_next / SEC_PER_DAY
    r_pred, v_pred = eci_to_ecef(r_pred_eci, v_pred_eci, mjd_next)

    r_gnv = interp(t_next)
    pred_err = np.linalg.norm(r_pred - r_gnv)
    print(f'  step {step}: pred_err={pred_err:.3f}m  |r_pred|={np.linalg.norm(r_pred)/1000:.3f}km  |r_gnv|={np.linalg.norm(r_gnv)/1000:.3f}km')

    # Use GNV1B position + dynamics-corrected velocity
    r_cur = r_gnv
    v_cur = v_pred + (r_gnv - r_pred) / dt_s

# --- Multi-step ECEF comparison ---
print()
print('=== Multi-step prediction test (ECEF, 2b+J2+drag) ===')
r_cur, v_cur = r0.copy(), v0_ecef.copy()
for step in range(5):
    t_cur = t0 + step * dt_s
    t_next = t_cur + dt_s

    r_pred, v_pred = rk4_step(r_cur, v_cur, 2.2, area_to_mass=0.002, dt=dt_s)

    r_gnv = interp(t_next)
    pred_err = np.linalg.norm(r_pred - r_gnv)
    print(f'  step {step}: pred_err={pred_err:.3f}m')

    r_cur = r_gnv
    v_cur = v_pred + (r_gnv - r_pred) / dt_s
