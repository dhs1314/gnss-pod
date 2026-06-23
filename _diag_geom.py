#!/usr/bin/env python3
"""Check geometry: LOS vectors, reference position, coordinate frames"""
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

# First GNV1B epoch
print("First 5 GNV1B entries:")
for ts in orbit_ts[:5]:
    r = ref_orbit[ts]
    print(f"  {ts}: [{r[0]/1000:.3f}, {r[1]/1000:.3f}, {r[2]/1000:.3f}] km, r={np.linalg.norm(r)/1000:.3f}km")

# Check first epoch
gps_sods = sorted(gps_obs.keys())
sod0 = gps_sods[0]
utc0 = gps_sod_to_utc(sod0)
print(f"\nFirst GPS epoch: gps_sod={sod0} utc={utc0}")

# Interpolate
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
        return np.array(ref_orbit[t0], dtype=float), t0
    a = dt_frac / dt_tot
    return np.array(ref_orbit[t0], dtype=float)*(1-a) + np.array(ref_orbit[t1], dtype=float)*a, t0

rcv_pos, t0_ref = get_rcv_pos(utc0)
print(f"rcv_pos: [{rcv_pos[0]/1000:.3f}, {rcv_pos[1]/1000:.3f}, {rcv_pos[2]/1000:.3f}] km")
print(f"  r = {np.linalg.norm(rcv_pos)/1000:.3f} km")

# Check a few SVs: LOS vectors
print("\nLOS vectors:")
for sv in sorted(gps_obs[sod0].keys())[:5]:
    rec = gps_obs[sod0][sv]
    pos, clk, vel = get_gps_pos_from_sp3(sp3, sv, utc0)
    if pos is None: continue

    delta = pos - rcv_pos
    rng = np.linalg.norm(delta)
    los = delta / rng

    print(f"  {sv}: sat=[{pos[0]/1000:.1f},{pos[1]/1000:.1f},{pos[2]/1000:.1f}]km"
          f"  rng={rng/1000:.1f}km  los=[{los[0]:.4f},{los[1]:.4f},{los[2]:.4f}]")

# Check if the reference position at UTC key exists
ref_key = utc0
if ref_key in ref_orbit:
    ref_pos = ref_orbit[ref_key]
    print(f"\nReference at {utc0}: [{ref_pos[0]/1000:.3f}, {ref_pos[1]/1000:.3f}, {ref_pos[2]/1000:.3f}] km")
    print(f"  rcv-ref: {np.linalg.norm(rcv_pos - ref_pos):.1f}m")
else:
    # Find closest
    closest = min(orbit_ts, key=lambda t: abs((t - utc0).total_seconds()))
    ref_pos = ref_orbit[closest]
    print(f"\nClosest ref at {closest} (dt={(closest-utc0).total_seconds():.1f}s):")
    print(f"  [{ref_pos[0]/1000:.3f}, {ref_pos[1]/1000:.3f}, {ref_pos[2]/1000:.3f}] km")
    print(f"  rcv-ref: {np.linalg.norm(rcv_pos - ref_pos):.1f}m")

# Simple LSQ: solve for position with known clock
print("\n\n=== Simple LSQ position check ===")
# Use P_r (code measurements) — no ambiguity issue
# State: [dX, dY, dZ, clk_r]
H_rows, y_rows = [], []
for sv in sorted(gps_obs[sod0].keys()):
    rec = gps_obs[sod0][sv]
    pos, clk, vel = get_gps_pos_from_sp3(sp3, sv, utc0)
    if pos is None: continue
    delta = pos - rcv_pos
    rng = np.linalg.norm(delta)
    los = delta / rng

    L_if, P_if, B_if = compute_if_m(rec['L1_cyc'], rec['L2_cyc'], rec['P1'], rec['P2'])
    sag = (OMEGA_E / C) * (pos[0] * rcv_pos[1] - pos[1] * rcv_pos[0])
    rho_corr = rng + sag
    P_r = P_if + clk - rho_corr

    H_rows.append([-los[0], -los[1], -los[2], 1.0])
    y_rows.append(P_r)

H = np.array(H_rows)
y = np.array(y_rows)
print(f"H shape: {H.shape}")
print(f"y: {y/1000} km")

# Solve LSQ
try:
    dx, residuals, rank, svals = np.linalg.lstsq(H, y, rcond=None)
    print(f"dx = [{dx[0]:.1f}, {dx[1]:.1f}, {dx[2]:.1f}]m  clk={dx[3]/1000:.3f}km")
    print(f"Residuals: mean={np.mean(residuals):.1f}m std={np.std(residuals):.1f}m")

    solved = rcv_pos + dx[:3]
    if ref_key in ref_orbit:
        err = solved - ref_orbit[ref_key]
        print(f"3D error: {np.linalg.norm(err):.1f}m")
except Exception as e:
    print(f"LSQ failed: {e}")

# Check: are the residuals per-SV?
print("\nPer-SV P_r values:")
for sv in sorted(gps_obs[sod0].keys())[:8]:
    rec = gps_obs[sod0][sv]
    pos, clk, vel = get_gps_pos_from_sp3(sp3, sv, utc0)
    if pos is None: continue
    delta = pos - rcv_pos
    rng = np.linalg.norm(delta)
    sag = (OMEGA_E / C) * (pos[0] * rcv_pos[1] - pos[1] * rcv_pos[0])
    L_if, P_if, B_if = compute_if_m(rec['L1_cyc'], rec['L2_cyc'], rec['P1'], rec['P2'])
    P_r = P_if + clk - (rng + sag)

    # What should P_r be? = c*dt_r + T + b_Pif
    # If we subtract the median, the receiver clock cancels:
    print(f"  {sv}: P_r={P_r/1000:.3f}km")
