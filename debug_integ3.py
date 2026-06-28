"""Compare ECEF vs ECI-full integration accuracy vs GNV1B."""
import sys, numpy as np
from pathlib import Path
from datetime import datetime
sys.path.insert(0, 'src')
from orbit_dynamics import total_acc, total_acc_eci, GM
from coordinates import ecef_to_eci, eci_to_ecef
from gravity_model import read_icgem_gfc
from orbit_integrator import rk4_step_eci, integrate_orbit_eci

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

# ========== ECEF integration (two-body + J2 + drag + Coriolis) ==========
def rk4_step_ecef(r, v, dt):
    def f(s):
        ri, vi = s[:3], s[3:6]
        a = total_acc(ri, vi, Cd=2.2, area_to_mass=0.68/580.0)
        return np.concatenate([vi, a])
    s = np.concatenate([r, v])
    k1 = f(s); k2 = f(s + 0.5*dt*k1); k3 = f(s + 0.5*dt*k2); k4 = f(s + dt*k3)
    sn = s + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    return sn[:3], sn[3:6]

# ========== ECI full model (uses integrator with MJD adjustment) ==========
def force_model(pos_eci, vel_eci, **kwargs):
    return total_acc_eci(pos_eci, vel_eci,
                         Cnm=Cnm, Snm=Snm, Nmax=90,
                         GM_gravity=GM_grav, R_gravity=R_grav,
                         **kwargs)

# Integrate both for 30 min with 10s steps
dt = 10.0
n_steps = 180  # 1800s

# ECEF
r_ec, v_ec = r0_ecef.copy(), v0_ecef.copy()
# ECI full (using rk4_step_eci with proper MJD per step)
r_ei, v_ei = r0_eci.copy(), v0_eci.copy()

print(f'{"Step":>5s}  {"ECEF_err(m)":>12s}  {"ECI-full_err(m)":>16s}')
diffs_ec = []
diffs_ei = []

for step in range(n_steps + 1):
    t_s = step * dt

    if step % 30 == 0 and step > 0:
        # ECEF comparison
        gps_sod = gps0 + t_s
        idx = np.searchsorted(ts, gps_sod)
        if idx < len(ts) and idx > 0:
            t1, t0 = ts[idx], ts[idx-1]
            r_gnv = ref_orbit[t0] + (ref_orbit[t1] - ref_orbit[t0]) * (gps_sod - t0) / (t1 - t0) if t1 != t0 else ref_orbit[t0]
            d_ec = np.linalg.norm(r_ec - r_gnv)
            diffs_ec.append(d_ec)

            # ECI full comparison: convert to ECEF at epoch
            mjd_ep = mjd0 + t_s / SEC_PER_DAY
            r_ei_ecef, _ = eci_to_ecef(r_ei, np.zeros(3), mjd_ep)
            d_ei = np.linalg.norm(r_ei_ecef - r_gnv)
            diffs_ei.append(d_ei)

            if step <= 300 or step % 300 == 0:
                print(f'{t_s:5.0f}s  {d_ec:12.3f}m  {d_ei:16.3f}m')

    if step < n_steps:
        # ECEF step
        r_ec, v_ec = rk4_step_ecef(r_ec, v_ec, dt)
        # ECI step (with MJD advancement for correct Earth rotation)
        step_mjd_utc = mjd0 + step * dt / SEC_PER_DAY
        step_mjd_tt = mjd_tt0 + step * dt / SEC_PER_DAY
        r_ei, v_ei = rk4_step_eci(r_ei, v_ei, dt, force_model,
                                  mjd_utc=step_mjd_utc, mjd_tt=step_mjd_tt,
                                  CD=2.2, CR=1.3,
                                  area_drag=0.68, area_srp=3.4, mass=580.0,
                                  bodies=['Sun', 'Moon'])

if diffs_ec and diffs_ei:
    diffs_ec = np.array(diffs_ec)
    diffs_ei = np.array(diffs_ei)
    print(f'\nSummary (3D RMS vs GNV1B, 30 min):')
    print(f'  ECEF (2b+J2+drag): {np.sqrt(np.mean(diffs_ec**2)):.3f} m  (mean={np.mean(diffs_ec):.1f}, max={np.max(diffs_ec):.1f})')
    print(f'  ECI-full (GGM05C): {np.sqrt(np.mean(diffs_ei**2)):.3f} m  (mean={np.mean(diffs_ei):.1f}, max={np.max(diffs_ei):.1f})')
