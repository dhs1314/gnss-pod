"""Compare ECEF vs ECI integration accuracy."""
import sys, numpy as np
from pathlib import Path
from datetime import datetime
sys.path.insert(0, 'src')
from orbit_dynamics import total_acc, total_acc_eci, GM
from coordinates import ecef_to_eci, eci_to_ecef
from gravity_model import read_icgem_gfc

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
print(f'Start GPS SOD: {gps0}, MJD: {mjd0:.6f}')

# Load gravity
Cnm, Snm, Nmax, GM_grav, R_grav = read_icgem_gfc(str(dp / 'gravity' / 'GGM05C.gfc'))
print(f'Gravity loaded: Nmax={Nmax}, GM={GM_grav:.6e}, R={R_grav:.1f}')

# === ECEF Integration (RK4, 10s step) ===
def rk4_step_ecef(r, v, dt):
    def f(s):
        ri, vi = s[:3], s[3:6]
        a = total_acc(ri, vi, Cd=2.2, area_to_mass=0.68/580.0)
        return np.concatenate([vi, a])
    s = np.concatenate([r, v])
    k1 = f(s); k2 = f(s + 0.5*dt*k1); k3 = f(s + 0.5*dt*k2); k4 = f(s + dt*k3)
    sn = s + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    return sn[:3], sn[3:6]

# === ECI Integration with pure two-body ===
def rk4_step_eci_2b(r, v, dt):
    def f(s):
        ri, vi = s[:3], s[3:6]
        r_mag = np.linalg.norm(ri)
        a = -GM * ri / r_mag**3
        return np.concatenate([vi, a])
    s = np.concatenate([r, v])
    k1 = f(s); k2 = f(s + 0.5*dt*k1); k3 = f(s + 0.5*dt*k2); k4 = f(s + dt*k3)
    sn = s + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    return sn[:3], sn[3:6]

# === ECI Integration with full force model (no MJD adjustment per step) ===
def rk4_step_eci_full_no_mjd(r, v, dt, mjd_utc, mjd_tt):
    def f(s):
        ri, vi = s[:3], s[3:6]
        a = total_acc_eci(ri, vi, mjd_tt, mjd_utc,
                         Cnm, Snm, 90,
                         CD=2.2, CR=1.3,
                         area_drag=0.68, area_srp=3.4, mass=580.0,
                         bodies=['Sun', 'Moon'],
                         GM_gravity=GM_grav, R_gravity=R_grav)
        return np.concatenate([vi, a])
    s = np.concatenate([r, v])
    k1 = f(s); k2 = f(s + 0.5*dt*k1); k3 = f(s + 0.5*dt*k2); k4 = f(s + dt*k3)
    sn = s + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    return sn[:3], sn[3:6]

# Integrate ECEF
r_ecef, v_ecef = r0_ecef.copy(), v0_ecef.copy()
pos_ecef = [r_ecef.copy()]
for step in range(180):
    r_ecef, v_ecef = rk4_step_ecef(r_ecef, v_ecef, 10.0)
    if (step + 1) % 3 == 0:
        pos_ecef.append(r_ecef.copy())

# Integrate ECI two-body
r0_eci, v0_eci = ecef_to_eci(r0_ecef, v0_ecef, mjd0)
r_eci_2b, v_eci_2b = r0_eci.copy(), v0_eci.copy()
pos_eci_2b_ecef = []
for step in range(180):
    r_eci_2b, v_eci_2b = rk4_step_eci_2b(r_eci_2b, v_eci_2b, 10.0)
    if (step + 1) % 3 == 0:
        mjd_step = mjd0 + (step + 1) * 10.0 / SEC_PER_DAY
        r_ecef_2b, _ = eci_to_ecef(r_eci_2b, np.zeros(3), mjd_step)
        pos_eci_2b_ecef.append(r_ecef_2b)

# Integrate ECI full force (WITHOUT MJD per-step adjustment)
r_eci_full, v_eci_full = r0_eci.copy(), v0_eci.copy()
mjd_tt0 = mjd0 + 69.184 / SEC_PER_DAY
pos_eci_full_ecef = []
for step in range(180):
    r_eci_full, v_eci_full = rk4_step_eci_full_no_mjd(r_eci_full, v_eci_full, 10.0, mjd0, mjd_tt0)
    if (step + 1) % 3 == 0:
        mjd_step = mjd0 + (step + 1) * 10.0 / SEC_PER_DAY
        r_ecef_full, _ = eci_to_ecef(r_eci_full, np.zeros(3), mjd_step)
        pos_eci_full_ecef.append(r_ecef_full)

print(f'\n{"Step":>5s}  {"ECEF_vs_GNV1B":>13s}  {"ECI-2b_vs_GNV1B":>16s}  {"ECI-full_vs_GNV1B":>17s}')
for i in range(min(len(pos_ecef), len(pos_eci_2b_ecef), len(pos_eci_full_ecef))):
    t_s = (i + 1) * 30
    gps_sod = gps0 + t_s
    idx = np.searchsorted(ts, gps_sod)
    if idx >= len(ts): break
    t1, t0 = ts[idx], ts[idx-1] if idx > 0 else ts[0]
    if t1 == t0: r_gnv = ref_orbit[t0]
    else: r_gnv = ref_orbit[t0] + (ref_orbit[t1] - ref_orbit[t0]) * (gps_sod - t0) / (t1 - t0)

    d_ecef = np.linalg.norm(pos_ecef[i] - r_gnv)
    d_2b = np.linalg.norm(pos_eci_2b_ecef[i] - r_gnv)
    d_full = np.linalg.norm(pos_eci_full_ecef[i] - r_gnv)
    if i < 10 or i % 10 == 0:
        print(f'{t_s:5.0f}s  {d_ecef:13.3f}m  {d_2b:16.3f}m  {d_full:17.3f}m')

# Summary
diffs_ecef = []
diffs_2b = []
diffs_full = []
for i in range(min(len(pos_ecef), len(pos_eci_2b_ecef), len(pos_eci_full_ecef))):
    t_s = (i + 1) * 30
    gps_sod = gps0 + t_s
    idx = np.searchsorted(ts, gps_sod)
    if idx >= len(ts): break
    t1, t0 = ts[idx], ts[idx-1] if idx > 0 else ts[0]
    r_gnv = ref_orbit[t0] + (ref_orbit[t1] - ref_orbit[t0]) * (gps_sod - t0) / (t1 - t0) if t1 != t0 else ref_orbit[t0]
    diffs_ecef.append(np.linalg.norm(pos_ecef[i] - r_gnv))
    diffs_2b.append(np.linalg.norm(pos_eci_2b_ecef[i] - r_gnv))
    diffs_full.append(np.linalg.norm(pos_eci_full_ecef[i] - r_gnv))

print(f'\nSummary (3D RMS vs GNV1B):')
print(f'  ECEF (2b+J2+drag): {np.sqrt(np.mean(np.array(diffs_ecef)**2)):.3f} m')
print(f'  ECI-2b only:       {np.sqrt(np.mean(np.array(diffs_2b)**2)):.3f} m')
print(f'  ECI-full (no MJD): {np.sqrt(np.mean(np.array(diffs_full)**2)):.3f} m')
