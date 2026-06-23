"""Test Orekit propagator with GRACE-FO initial conditions."""
import os, sys
os.environ['JAVA_HOME'] = r"C:\Program Files\JetBrains\PyCharm Community Edition 2024.3.5\jbr"
os.environ['OREKIT_DATA_PATH'] = r"d:\prj\gnss_pod\data\orekit"

import numpy as np

# Load GGM05C from existing pipeline to compare
print("=" * 60)
print("Testing OrekitPropagator for GRACE-FO POD")
print("=" * 60)

from src.orekit_bridge import OrekitPropagator, is_orekit_available
print(f"Orekit available: {is_orekit_available()}")

# GRACE-FO C initial state at 2024-04-29 12:00:00 (approximate)
# From GNV1B reference orbit
r_eci = np.array([5264840.5, -3202983.7, 4031685.2])   # m
v_eci = np.array([-3786.64, -4212.80, -5071.31])       # m/s
a_rtn = np.array([0.0, 0.0, 0.0])  # zero empirical
dt = 30.0  # 30s EKF step

# MJD for 2024-04-29 12:00:00 UTC
# 2024-04-29 = MJD 60429 (approx), 12:00 = 0.5 day
mjd_utc = 60429.5
mjd_tt = mjd_utc + 69.184 / 86400.0

print(f"\nInitial state:")
print(f"  r_eci = [{r_eci[0]/1000:.3f}, {r_eci[1]/1000:.3f}, {r_eci[2]/1000:.3f}] km")
print(f"  v_eci = [{v_eci[0]:.3f}, {v_eci[1]:.3f}, {v_eci[2]:.3f}] m/s")
print(f"  |r| = {np.linalg.norm(r_eci)/1000:.3f} km")
print(f"  |v| = {np.linalg.norm(v_eci):.3f} m/s")
print(f"  dt = {dt} s")

# Create Orekit propagator
print("\nCreating OrekitPropagator...")
prop = OrekitPropagator(
    gravity_field=r'd:\prj\gnss_pod\data\gravity\GGM05C.gfc',
    gravity_degree=90,
    solid_tides=True,
    ocean_tides=False,  # FES2004 data not available
    relativity=True,
    mass=580.0,
    area_drag=0.68,
    area_srp=3.4,
    CR=1.3,
    CD=2.2,
)

print("\nPropagating 30s with Orekit full dynamics...")
r_new, v_new, Phi, S = prop.propagate(r_eci, v_eci, a_rtn, dt, mjd_utc, mjd_tt)

print(f"\nResults:")
print(f"  r_new = [{r_new[0]/1000:.3f}, {r_new[1]/1000:.3f}, {r_new[2]/1000:.3f}] km")
print(f"  v_new = [{v_new[0]:.3f}, {v_new[1]:.3f}, {v_new[2]:.3f}] m/s")
print(f"  |r_new| = {np.linalg.norm(r_new)/1000:.3f} km")
print(f"  |v_new| = {np.linalg.norm(v_new):.3f} m/s")
print(f"  |dr| = {np.linalg.norm(r_new - r_eci):.3f} m")
print(f"  |dv| = {np.linalg.norm(v_new - v_eci):.3f} m/s")
print(f"  Phi (6x6) = {Phi.shape}, det={np.linalg.det(Phi):.6f}")
print(f"  S (6x3) = {S.shape}")

# Compare with simplified dynamics
print("\nComparing with simplified dynamics...")
from src.orbit_dynamics import total_acc_eci
from src.gravity_model import read_icgem_gfc
from src.orbit_integrator import integrate_orbit_eci_with_stm

gfc_path = r'd:\prj\gnss_pod\data\gravity\GGM05C.gfc'
Cnm, Snm, Nmax, GM_grav, R_grav = read_icgem_gfc(gfc_path)
Nmax = min(Nmax, 90)

def force_model(pos, vel, **kw):
    # Extract integrator-provided time args (f_ext adjusts per sub-step)
    tt = kw.pop('mjd_tt')
    utc = kw.pop('mjd_utc')
    return total_acc_eci(pos, vel, tt, utc,
                         Cnm, Snm, Nmax,
                         GM_gravity=GM_grav, R_gravity=R_grav,
                         **kw)

integ = integrate_orbit_eci_with_stm(
    r_eci, v_eci, (0.0, dt), force_model,
    Cd=2.2, CR=1.3,
    area_drag=0.68, area_srp=3.4,
    mass=580.0,
    empirical_acc_rtn=a_rtn,
    param_names=['aR', 'aT', 'aN'],
    dt=10.0, mjd_tt=mjd_tt, mjd_utc=mjd_utc,
    bodies=['Sun', 'Moon'],
)

r_simple = integ['r'][-1]
v_simple = integ['v'][-1]

print(f"  Simplified r_new: [{r_simple[0]/1000:.3f}, {r_simple[1]/1000:.3f}, {r_simple[2]/1000:.3f}] km")
print(f"  Simplified v_new: [{v_simple[0]:.3f}, {v_simple[1]:.3f}, {v_simple[2]:.3f}] m/s")
print(f"  |dr| simple: {np.linalg.norm(r_simple - r_eci):.3f} m")
print(f"  |dv| simple: {np.linalg.norm(v_simple - v_eci):.3f} m/s")

print(f"\n  Position difference (Orekit - simplified): {np.linalg.norm(r_new - r_simple):.3f} m")
print(f"  Velocity difference (Orekit - simplified): {np.linalg.norm(v_new - v_simple):.6f} m/s")
print("\n=== Orekit propagator test complete ===")
