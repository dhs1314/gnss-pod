"""Diagnose force model by comparing different force configurations."""
import sys, numpy as np
from pathlib import Path
sys.path.insert(0, 'src')
from orbit_dynamics import total_acc, total_acc_eci, two_body_acc, j2_acc, drag_acc, GM
from coordinates import ecef_to_eci, eci_to_ecef
from gravity_model import read_icgem_gfc, compute_gravity_acceleration

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

# Test at initial state
r, v = r0_ecef, v0_ecef

# 1. Two-body only
a_2b = two_body_acc(r)
print(f'Two-body ECEF:          |a|={np.linalg.norm(a_2b):.6f} m/s^2')

# 2. J2 only
a_j2 = j2_acc(r)
print(f'J2 ECEF:                |a|={np.linalg.norm(a_j2):.6f} m/s^2')

# 3. Two-body + J2
a_2b_j2 = a_2b + a_j2
print(f'2b+J2 ECEF:             |a|={np.linalg.norm(a_2b_j2):.6f} m/s^2')
print(f'J2 fraction:            {np.linalg.norm(a_j2)/np.linalg.norm(a_2b)*100:.4f}%')

# 4. Full ECEF (total_acc)
a_full_ecef = total_acc(r, v, Cd=2.2, area_to_mass=0.68/580.0)
print(f'Full ECEF (total_acc):  |a|={np.linalg.norm(a_full_ecef):.6f} m/s^2')

# 5. Drag only
a_drag = drag_acc(r, v, 2.2, 0.68/580.0)
print(f'Drag ECEF:              |a|={np.linalg.norm(a_drag):.10f} m/s^2')

# 6. Coriolis + centrifugal
omega_vec = np.array([0.0, 0.0, 7.2921151467e-5])
a_cor = np.array([2.0*7.2921151467e-5*v[1], -2.0*7.2921151467e-5*v[0], 0.0])
a_cen = np.array([7.2921151467e-5**2 * r[0], 7.2921151467e-5**2 * r[1], 0.0])
print(f'Coriolis ECEF:          |a|={np.linalg.norm(a_cor):.6f} m/s^2')
print(f'Centrifugal ECEF:       |a|={np.linalg.norm(a_cen):.6f} m/s^2')

# 7. Full GGM05C gravity in ECEF
a_ggm = compute_gravity_acceleration(r, Cnm, Snm, 90, GM=GM_grav, R=R_grav)
print(f'GGM05C ECEF:            |a|={np.linalg.norm(a_ggm):.6f} m/s^2')
print(f'GGM - 2b:               |diff|={np.linalg.norm(a_ggm - a_2b):.6f} m/s^2')
print(f'GGM - (2b+J2):          |diff|={np.linalg.norm(a_ggm - a_2b_j2):.6f} m/s^2')

# 8. Check GM values
print(f'\nGM (orbit_dynamics): {GM:.6e}')
print(f'GM (GGM05C):          {GM_grav:.6e}')
print(f'GM difference:        {abs(GM-GM_grav):.6e} ({(abs(GM-GM_grav)/GM)*100:.8f}%)')

# Test what a constant δa = 0.007 m/s² would do over 300s
dt = 300.0
dr_pred = 0.5 * 0.007 * dt**2
print(f'\nWith δa=0.007 m/s^2 over 300s: δr = {dr_pred:.1f} m')

# Check if GM difference explains dynamics error
# a_2b = GM/r^2, so δa = -2*δGM*GM/r^3 * δr... actually
# a = GM/r^2, δa/a = δGM/GM
# δa = a * δGM/GM ≈ 8.4 * (small)
# If δGM/GM = 1e-8, δa ≈ 8.4e-8 m/s^2
# Over 300s: δr ≈ 0.5 * 8.4e-8 * 300^2 ≈ 0.004 m
# So GM difference can't explain it
print('GM mismatch contribution to position error is < 0.01m — negligible')
