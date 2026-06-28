#!/usr/bin/env python3
"""Deep dive: what's in P_r for each SV?"""
import sys; sys.path.insert(0, '.')
from src.gps1a_loader import download_gps1a, gps_sod_to_utc
from src.ambiguity import compute_if_m
from src.sp3_loader import get_gps_pos_from_sp3
import pickle, numpy as np
from datetime import datetime, timedelta
from collections import defaultdict

C = 299792458.0
OMEGA_E = 7.2921151467e-5
J2000 = datetime(2000, 1, 1, 12, 0, 0)

gps_obs = download_gps1a(2024, 4, 29, grace_filter='C')
sp3 = pickle.load(open('data/2024/120/igs_sp3_FIN.pkl', 'rb'))

# Load reference orbit
ref_orbit = {}
for line in open('data/gracefo/2024/2024-04-29/GNV1B_2024-04-29_C_04.txt'):
    p = line.split()
    if len(p) < 6: continue
    try:
        t = float(p[0]); flag = p[2]
        if flag in ('C', 'E'):
            ref_orbit[gps_sod_to_utc(t)] = np.array([float(p[3]), float(p[4]), float(p[5])])
    except: pass
orbit_ts = sorted(ref_orbit.keys())

def get_rcv_pos(utc_dt):
    t0, t1 = None, None
    for j, ti in enumerate(orbit_ts):
        if ti >= utc_dt:
            t1 = ti; t0 = orbit_ts[j-1] if j>0 else ti; break
        t0 = ti
    if t0 is None: t0 = orbit_ts[0]
    if t1 is None: t1 = orbit_ts[-1]
    dt_frac = (utc_dt - t0).total_seconds()
    dt_tot = (t1 - t0).total_seconds()
    if dt_tot < 0.01:
        return np.array(ref_orbit[t0], dtype=float)
    a = dt_frac / dt_tot
    return np.array(ref_orbit[t0], dtype=float)*(1-a) + np.array(ref_orbit[t1], dtype=float)*a

# Take one epoch with many SVs and decompose P_r per SV
gps_sods = sorted(gps_obs.keys())
sod = gps_sods[0]
utc_dt = gps_sod_to_utc(sod)
rcv_pos = get_rcv_pos(utc_dt)

print(f"Epoch: gps_sod={sod} utc={utc_dt}")
print(f"rcv_pos: {rcv_pos/1000} km")
print(f"SP3 epochs: {len(sp3['ts'])} from {min(sp3['ts'])} to {max(sp3['ts'])}")
print()

# Check each visible SV
for sv in sorted(gps_obs[sod].keys()):
    rec = gps_obs[sod][sv]
    pos, clk, vel = get_gps_pos_from_sp3(sp3, sv, utc_dt)
    if pos is None:
        print(f"{sv}: NO SP3 DATA!")
        continue

    rho = float(np.linalg.norm(pos - rcv_pos))
    sag = (OMEGA_E / C) * (pos[0] * rcv_pos[1] - pos[1] * rcv_pos[0])
    rho_corr = rho + sag

    L_if, P_if, B_if = compute_if_m(rec['L1_cyc'], rec['L2_cyc'], rec['P1'], rec['P2'])

    P_r = P_if + clk - rho_corr

    print(f"{sv}:")
    print(f"  P1={rec['P1']/1000:.3f}km  P2={rec['P2']/1000:.3f}km  P2-P1={rec['P2']-rec['P1']:.1f}m")
    print(f"  P_if={P_if/1000:.3f}km")
    print(f"  sat_pos ECEF: [{pos[0]/1000:.1f}, {pos[1]/1000:.1f}, {pos[2]/1000:.1f}] km")
    print(f"  rho={rho/1000:.3f}km  sag={sag:.2f}m  rho_corr={rho_corr/1000:.3f}km")
    print(f"  clk_sv={clk:.1f}m ({clk/C*1e6:.3f}us)")
    print(f"  P_if - rho_corr = {(P_if - rho_corr)/1000:.3f}km")
    print(f"  P_r = {P_r/1000:.3f}km")

# Now compare: P_r should ≈ receiver clock offset (common across SVs)
print()
print("=== Summary ===")
P_r_vals = []
for sv in sorted(gps_obs[sod].keys()):
    rec = gps_obs[sod][sv]
    pos, clk, vel = get_gps_pos_from_sp3(sp3, sv, utc_dt)
    if pos is None: continue
    rho = float(np.linalg.norm(pos - rcv_pos))
    sag = (OMEGA_E / C) * (pos[0] * rcv_pos[1] - pos[1] * rcv_pos[0])
    rho_corr = rho + sag
    L_if, P_if, B_if = compute_if_m(rec['L1_cyc'], rec['L2_cyc'], rec['P1'], rec['P2'])
    P_r = P_if + clk - rho_corr
    P_r_vals.append((sv, P_r, P_if - rho_corr, clk))

for sv, Pr, diff, clk in P_r_vals:
    print(f"  {sv}: P_r={Pr/1000:.3f}km  P_if-rho_corr={diff/1000:.3f}km  clk_sv={clk/1000:.3f}km")
