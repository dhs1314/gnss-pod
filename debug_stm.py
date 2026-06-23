"""Verify STM-based observation partials against numerical finite differences."""
import sys, numpy as np
from pathlib import Path
sys.path.insert(0, 'src')
from orbit_dynamics import total_acc_eci
from coordinates import ecef_to_eci, eci_to_ecef
from gravity_model import read_icgem_gfc
from orbit_integrator import integrate_orbit_eci_with_stm

SEC_PER_DAY = 86400.0
MJD_J2000 = 51544.5

# Load GNV1B
dp = Path('data')
date_str = '2024-04-29'; y, m, d = 2024, 4, 29; grace_id = 'C'

ref_orbit = {}
with open(str(dp / 'gracefo' / str(y) / date_str / f'GNV1B_{date_str}_{grace_id}_04.txt'), encoding='utf-8', errors='replace') as f:
    for line in f:
        if line.startswith('#') or not line.strip(): continue
        parts = line.split()
        if len(parts) < 6: continue
        try:
            t, flag = float(parts[0]), parts[2]
            if flag not in ('C', 'E'): continue
            X, Y, Z = float(parts[3]), float(parts[4]), float(parts[5])
            if abs(X) < 1e3: continue
            ref_orbit[t] = np.array([X, Y, Z])
        except: continue

ts = sorted(ref_orbit.keys())
gps0 = ts[0]
r0_ecef = ref_orbit[gps0].copy()
i0 = ts.index(gps0)
v0_ecef = (ref_orbit[ts[i0+1]] - ref_orbit[ts[i0]]) / (ts[i0+1] - ts[i0])
mjd0 = MJD_J2000 + gps0 / SEC_PER_DAY
r0_eci, v0_eci = ecef_to_eci(r0_ecef, v0_ecef, mjd0)
mjd_tt0 = mjd0 + 69.184 / SEC_PER_DAY

# Load gravity
Cnm, Snm, Nmax, GM_grav, R_grav = read_icgem_gfc(str(dp / 'gravity' / 'GGM05C.gfc'))

def force_model(pos_eci, vel_eci, **kwargs):
    return total_acc_eci(pos_eci, vel_eci,
                         Cnm=Cnm, Snm=Snm, Nmax=90,
                         GM_gravity=GM_grav, R_gravity=R_grav,
                         **kwargs)

# Reference integration
t_span = (0.0, 300.0)  # 5 minutes
integ = integrate_orbit_eci_with_stm(
    r0_eci, v0_eci, t_span, force_model,
    Cd=2.2, CR=1.3,
    area_drag=0.68, area_srp=3.4, mass=580.0,
    empirical_acc_rtn=None,
    param_names=['Cd', 'CR'],
    dt=10.0,
    mjd_tt=mjd_tt0, mjd_utc=mjd0,
    bodies=['Sun', 'Moon'],
)

print("STM-based partials vs numerical partials:")
print(f'{"t(s)":>6s}  {"dr0_x":>12s}  {"dr0_y":>12s}  {"dr0_z":>12s}  {"dv0_x":>12s}  {"dCd":>12s}  {"dCR":>12s}')

# Test at a few epochs
for t_test in [30, 60, 120, 180, 300]:
    i_step = int(round(t_test / 10.0))
    if i_step >= len(integ['t']):
        continue

    phi_ep = integ['phi'][i_step]
    S_ep = integ['S'][i_step]
    param_names = integ['param_names']

    phi_rr = phi_ep[0:3, 0:3]
    phi_rv = phi_ep[0:3, 3:6]
    S_r = S_ep[0:3, :]
    col_Cd = param_names.index('Cd') if 'Cd' in param_names else None
    col_CR = param_names.index('CR') if 'CR' in param_names else None

    # Reference ECI position at t_test
    r_ref_eci = integ['r'][i_step]

    # === Numerical partials ===
    eps = 1.0  # 1m perturbation for position, 0.001 m/s for velocity

    # dr0_x
    r0_plus, v0_plus = ecef_to_eci(r0_ecef + np.array([eps, 0, 0]), v0_ecef, mjd0)
    integ_plus = integrate_orbit_eci_with_stm(
        r0_plus, v0_plus, (0.0, t_test), force_model,
        Cd=2.2, CR=1.3,
        area_drag=0.68, area_srp=3.4, mass=580.0,
        empirical_acc_rtn=None,
        param_names=['Cd', 'CR'],
        dt=10.0, mjd_tt=mjd_tt0, mjd_utc=mjd0, bodies=['Sun', 'Moon'],
    )
    r_plus = integ_plus['r'][-1]
    num_drx = (r_plus - r_ref_eci) / eps

    # dr0_y
    r0_plus, v0_plus = ecef_to_eci(r0_ecef + np.array([0, eps, 0]), v0_ecef, mjd0)
    integ_plus = integrate_orbit_eci_with_stm(
        r0_plus, v0_plus, (0.0, t_test), force_model,
        Cd=2.2, CR=1.3,
        area_drag=0.68, area_srp=3.4, mass=580.0,
        empirical_acc_rtn=None,
        param_names=['Cd', 'CR'],
        dt=10.0, mjd_tt=mjd_tt0, mjd_utc=mjd0, bodies=['Sun', 'Moon'],
    )
    r_plus = integ_plus['r'][-1]
    num_dry = (r_plus - r_ref_eci) / eps

    # dr0_z
    r0_plus, v0_plus = ecef_to_eci(r0_ecef + np.array([0, 0, eps]), v0_ecef, mjd0)
    integ_plus = integrate_orbit_eci_with_stm(
        r0_plus, v0_plus, (0.0, t_test), force_model,
        Cd=2.2, CR=1.3,
        area_drag=0.68, area_srp=3.4, mass=580.0,
        empirical_acc_rtn=None,
        param_names=['Cd', 'CR'],
        dt=10.0, mjd_tt=mjd_tt0, mjd_utc=mjd0, bodies=['Sun', 'Moon'],
    )
    r_plus = integ_plus['r'][-1]
    num_drz = (r_plus - r_ref_eci) / eps

    # Compare: STM says ∂r/∂r0 = phi_rr
    # phi_rr[:, 0] should match num_drx, etc.

    stm_drx = phi_rr[:, 0]
    stm_dry = phi_rr[:, 1]
    stm_drz = phi_rr[:, 2]

    # Compare first component (x) of each
    print(f'{t_test:6d}s  '
          f'{stm_drx[0]:12.6f}  {stm_dry[0]:12.6f}  {stm_drz[0]:12.6f}  '
          f'{"...":>12s}  {"...":>12s}  {"...":>12s}')
    print(f'{"NUM":>6s}  '
          f'{num_drx[0]:12.6f}  {num_dry[0]:12.6f}  {num_drz[0]:12.6f}')

    # Check the full norm
    print(f'{"|err|":>6s}  '
          f'{np.linalg.norm(stm_drx-num_drx):12.6f}  {np.linalg.norm(stm_dry-num_dry):12.6f}  {np.linalg.norm(stm_drz-num_drz):12.6f}')
    print()
