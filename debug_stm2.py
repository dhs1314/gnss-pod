"""Verify STM partials — perturb ECI initial state directly."""
import sys, numpy as np
from pathlib import Path
sys.path.insert(0, 'src')
from orbit_dynamics import total_acc_eci
from coordinates import ecef_to_eci, eci_to_ecef
from gravity_model import read_icgem_gfc
from orbit_integrator import integrate_orbit_eci_with_stm

SEC_PER_DAY = 86400.0
MJD_J2000 = 51544.5

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

Cnm, Snm, Nmax, GM_grav, R_grav = read_icgem_gfc(str(dp / 'gravity' / 'GGM05C.gfc'))

def force_model(pos_eci, vel_eci, **kwargs):
    return total_acc_eci(pos_eci, vel_eci,
                         Cnm=Cnm, Snm=Snm, Nmax=90,
                         GM_gravity=GM_grav, R_gravity=R_grav, **kwargs)

def integrate(t_end, r0, v0):
    """Integrate from 0 to t_end with the full force model."""
    if t_end < 1e-6:
        return r0, v0
    integ = integrate_orbit_eci_with_stm(
        r0, v0, (0.0, t_end), force_model,
        Cd=2.2, CR=1.3,
        area_drag=0.68, area_srp=3.4, mass=580.0,
        empirical_acc_rtn=None,
        param_names=['Cd', 'CR'],
        dt=10.0, mjd_tt=mjd_tt0, mjd_utc=mjd0, bodies=['Sun', 'Moon'],
    )
    return integ['r'][-1], integ['v'][-1]

# Reference integration
integ_ref = integrate_orbit_eci_with_stm(
    r0_eci, v0_eci, (0.0, 300.0), force_model,
    Cd=2.2, CR=1.3,
    area_drag=0.68, area_srp=3.4, mass=580.0,
    empirical_acc_rtn=None,
    param_names=['Cd', 'CR'],
    dt=10.0, mjd_tt=mjd_tt0, mjd_utc=mjd0, bodies=['Sun', 'Moon'],
)

eps_pos = 1.0      # 1m position perturbation in ECI
eps_vel = 0.001    # 0.001 m/s velocity perturbation in ECI

print("Comparing STM partials with numerical (ECI perturbations):")
print(f'{"t(s)":>6s}  {"STM_rr[0,:]":>30s}  {"NUM_rr[0,:]":>30s}  {"err_norm":>10s}')

for t_test in [30, 120, 300]:
    i_step = int(round(t_test / 10.0))
    phi_ep = integ_ref['phi'][i_step]
    S_ep = integ_ref['S'][i_step]
    r_ref = integ_ref['r'][i_step]

    phi_rr = phi_ep[0:3, 0:3]
    phi_rv = phi_ep[0:3, 3:6]

    # Numerical: perturb r0_ECI in x, y, z
    num_phi_rr = np.zeros((3, 3))
    num_phi_rv = np.zeros((3, 3))

    for j in range(3):
        pert = np.zeros(3)
        pert[j] = eps_pos
        r_plus, v_plus = integrate(t_test, r0_eci + pert, v0_eci)
        num_phi_rr[:, j] = (r_plus - r_ref) / eps_pos

        pert_v = np.zeros(3)
        pert_v[j] = eps_vel
        r_plus_v, v_plus_v = integrate(t_test, r0_eci, v0_eci + pert_v)
        num_phi_rv[:, j] = (r_plus_v - r_ref) / eps_vel

    # Compare STM vs numerical
    for j in range(3):
        err_rr = np.linalg.norm(phi_rr[:, j] - num_phi_rr[:, j])
        err_rv = np.linalg.norm(phi_rv[:, j] - num_phi_rv[:, j])
        if j == 0:
            print(f'{t_test:6d}s  STM rr[{j}]={phi_rr[0,j]:9.6f} rv[{j}]={phi_rv[0,j]:9.6f}')
            print(f'{"":>6s}  NUM rr[{j}]={num_phi_rr[0,j]:9.6f} rv[{j}]={num_phi_rv[0,j]:9.6f}')
            print(f'{"":>6s}  |err_rr|={err_rr:.6f}  |err_rv|={err_rv:.6f}')

    # Check full matrix error
    err_matrix = np.linalg.norm(phi_rr - num_phi_rr)
    print(f'{"":>6s}  Full |phi_rr - num| = {err_matrix:.6f}')
    print()
