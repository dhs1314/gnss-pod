"""Minimal ECI two-body integration test — check Keplerian invariants."""
import sys, numpy as np
sys.path.insert(0, 'src')
from coordinates import ecef_to_eci, eci_to_ecef
from orbit_dynamics import GM

SEC_PER_DAY = 86400.0
MJD_J2000 = 51544.5

# Load GNV1B
from pathlib import Path
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
dt_gnv = ts[i0+1] - ts[i0]
v0_ecef = (ref_orbit[ts[i0+1]] - ref_orbit[ts[i0]]) / dt_gnv
mjd0 = MJD_J2000 + gps0 / SEC_PER_DAY

print(f'GNV1B first epoch: gps={gps0}, dt={dt_gnv:.3f}s')
print(f'r0_ecef: {r0_ecef}')
print(f'v0_ecef: {v0_ecef}')
print(f'|r0|: {np.linalg.norm(r0_ecef)/1000:.3f} km')
print(f'|v0|: {np.linalg.norm(v0_ecef):.3f} m/s')

# ECEF -> ECI
r0_eci, v0_eci = ecef_to_eci(r0_ecef, v0_ecef, mjd0)
print(f'r0_eci: {r0_eci}')
print(f'v0_eci: {v0_eci}')

# Keplerian orbital elements in ECI
h_vec = np.cross(r0_eci, v0_eci)
e_vec = np.cross(v0_eci, h_vec) / GM - r0_eci / np.linalg.norm(r0_eci)
a = 1.0 / (2.0 / np.linalg.norm(r0_eci) - np.linalg.norm(v0_eci)**2 / GM)
e = np.linalg.norm(e_vec)
print(f'Semi-major axis: {a/1000:.3f} km')
print(f'Eccentricity: {e:.6f}')
print(f'|h|: {np.linalg.norm(h_vec):.3f} m^2/s')

# Simple RK4 two-body in ECI
def rk4_2b(r, v, dt):
    def f(s):
        ri, vi = s[:3], s[3:6]
        rm = np.linalg.norm(ri)
        a = -GM * ri / rm**3
        return np.concatenate([vi, a])
    s = np.concatenate([r, v])
    k1 = f(s); k2 = f(s + 0.5*dt*k1); k3 = f(s + 0.5*dt*k2); k4 = f(s + dt*k3)
    sn = s + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    return sn[:3], sn[3:6]

# Integrate for 30 min with 10s steps
r, v = r0_eci.copy(), v0_eci.copy()
print(f'\n{"Step":>5s}  {"|r|(km)":>10s}  {"|v|(m/s)":>10s}  {"a(km)":>10s}  {"e":>10s}  {"energy_err":>12s}')
for step in range(181):  # include step 0
    if step % 30 == 0:
        r_mag = np.linalg.norm(r) / 1000
        v_mag = np.linalg.norm(v)
        energy = 0.5 * v_mag**2 - GM / (r_mag * 1000)
        energy0 = 0.5 * np.linalg.norm(v0_eci)**2 - GM / np.linalg.norm(r0_eci)
        h_now = np.cross(r, v)
        e_now = np.cross(v, h_now) / GM - r / np.linalg.norm(r)
        a_now = 1.0 / (2.0 / np.linalg.norm(r) - np.linalg.norm(v)**2 / GM) / 1000
        print(f'{step*10:5d}s  {r_mag:10.3f}  {v_mag:10.3f}  {a_now:10.3f}  {np.linalg.norm(e_now):10.6f}  {abs(energy-energy0)/abs(energy0):12.3e}')
    if step < 180:
        r, v = rk4_2b(r, v, 10.0)

# Also check ECI -> ECEF at various times and compare with GNV1B
print(f'\nECI->ECEF comparison with GNV1B:')
r, v = r0_eci.copy(), v0_eci.copy()
for step in range(181):
    if step % 30 == 0 and step > 0:
        mjd_step = mjd0 + step * 10.0 / SEC_PER_DAY
        r_ecef, _ = eci_to_ecef(r, np.zeros(3), mjd_step)
        gps_sod = gps0 + step * 10.0
        # Interpolate GNV1B
        idx = np.searchsorted(ts, gps_sod)
        if idx < len(ts) and idx > 0:
            t1, t0 = ts[idx], ts[idx-1]
            r_gnv = ref_orbit[t0] + (ref_orbit[t1] - ref_orbit[t0]) * (gps_sod - t0) / (t1 - t0) if t1 != t0 else ref_orbit[t0]
            diff = np.linalg.norm(r_ecef - r_gnv)
            print(f'  t={step*10:5d}s: ECI->ECEF vs GNV1B diff = {diff:.3f} m')
    if step < 180:
        r, v = rk4_2b(r, v, 10.0)
