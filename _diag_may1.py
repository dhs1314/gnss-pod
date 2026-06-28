#!/usr/bin/env python3
"""Diagnose May 1, 2024 GPS1B data quality"""
import sys; sys.path.insert(0, 'src')
import pickle
from sp3_loader import get_gps_pos_from_sp3
import numpy as np
from datetime import datetime, timedelta

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
F1_SQ, F2_SQ = F1*F1, F2*F2
ALPHA = F1_SQ / (F1_SQ - F2_SQ)
BETA  = -F2_SQ / (F1_SQ - F2_SQ)
OMEGA_E = 7.2921151467e-5
J2000 = datetime(2000, 1, 1, 12, 0, 0)

# Load data
gps1b = pickle.load(open('data/gracefo/2024/2024-05-01/GPS1B_2024-05-01_C_04.pkl', 'rb'))
sp3 = pickle.load(open('data/2024/122/igs_sp3_FIN.pkl', 'rb'))

# Load GNV1B
ref_orbit = {}
for line in open('data/gracefo/2024/2024-05-01/GNV1B_2024-05-01_C_04.txt'):
    parts = line.split()
    if len(parts) < 6: continue
    try:
        t = float(parts[0]); flag = parts[2]
        if flag in ('C', 'E'):
            ref_orbit[t] = np.array([float(parts[3]), float(parts[4]), float(parts[5])])
    except: pass

def get_ref_pos(gps_sod):
    ts = sorted(ref_orbit.keys())
    t0 = t1 = None
    for i, ti in enumerate(ts):
        if ti >= gps_sod: t1 = ti; t0 = ts[i-1] if i > 0 else None; break
        t0 = ti
    if t1 is None: t0 = t1 = ts[-1]
    if t0 is None: t0 = ts[0]
    if t0 == t1: return ref_orbit[t0]
    a = (gps_sod - t0) / (t1 - t0)
    return ref_orbit[t0] * (1-a) + ref_orbit[t1] * a

# Check SP3 validity
sp3_ts = sorted(sp3['ts'])
print(f"SP3 times: {min(sp3_ts)} to {max(sp3_ts)}")
print(f"SP3 counts: {len(sp3_ts)}")

# Check GPS1B time range
gps1b_keys = sorted(gps1b.keys())
print(f"GPS1B: {len(gps1b_keys)} epochs, {gps1b_keys[0]} to {gps1b_keys[-1]}")
print(f"GPS1B start UTC: {J2000 + timedelta(seconds=gps1b_keys[0])}")
print(f"GPS1B end UTC:   {J2000 + timedelta(seconds=gps1b_keys[-1])}")

# Check GNV1B
gnv_keys = sorted(ref_orbit.keys())
print(f"GNV1B: {len(gnv_keys)} epochs, {gnv_keys[0]} to {gnv_keys[-1]}")

# Check first epoch
sod = gps1b_keys[0]
utc_dt = J2000 + timedelta(seconds=sod)
ref_pos = get_ref_pos(sod)
print(f"\nFirst epoch: sod={sod} utc={utc_dt}")
print(f"Ref pos: [{ref_pos[0]/1000:.3f}, {ref_pos[1]/1000:.3f}, {ref_pos[2]/1000:.3f}] km")

# Check a few SVs
ep = gps1b[sod]
for sv in sorted(ep.keys())[:3]:
    rec = ep[sv]
    pos, clk, vel = get_gps_pos_from_sp3(sp3, sv, utc_dt)
    if pos is None:
        print(f"  {sv}: NO SP3 DATA")
        continue
    rho = float(np.linalg.norm(pos - ref_pos))
    sag = (OMEGA_E / C) * (pos[0] * ref_pos[1] - pos[1] * ref_pos[0])
    rho_corr = rho + sag
    L_r = rec['L_if'] + clk - rho_corr
    P_r = rec['P_if'] + clk - rho_corr
    B_if = rec['L_if'] - rec['P_if']
    print(f"  {sv}: L_r={L_r:.1f}m P_r={P_r:.1f}m B_if={B_if:.2f}m rho={rho/1000:.1f}km")

# Compare L_if values across data sources
print(f"\n--- Comparing with April 29 data ---")
gps1b_0429 = pickle.load(open('data/gracefo/2024/2024-04-29/GPS1B_2024-04-29_C_04.pkl', 'rb'))
sod_0429 = sorted(gps1b_0429.keys())[0]
print(f"April 29 first epoch: {sod_0429}")
ep_0429 = gps1b_0429[sod_0429]
print(f"April 29 SVs: {sorted(ep_0429.keys())}")
print(f"May 1 SVs:    {sorted(ep.keys())}")

# Check if P_if values have different characteristics
for sv in sorted(set(ep_0429.keys()) & set(ep.keys()))[:3]:
    r0429 = ep_0429[sv]
    r0501 = ep[sv]
    print(f"{sv}:")
    print(f"  Apr29: P1={r0429['P1']:.2f} P2={r0429['P2']:.2f} P_if={r0429['P_if']:.2f} L_if={r0429['L_if']:.2f}")
    print(f"  May01: P1={r0501['P1']:.2f} P2={r0501['P2']:.2f} P_if={r0501['P_if']:.2f} L_if={r0501['L_if']:.2f}")
    print(f"  B_if: Apr29={r0429['L_if']-r0429['P_if']:.2f} May01={r0501['L_if']-r0501['P_if']:.2f}")

# Check SP3 quality: does SP3 have data for all SVs?
print(f"\n--- SP3 SV coverage check ---")
sp3_svs = set()
for t in sp3_ts[:5]:
    for k in sp3.get('data', {}).get(t, {}):
        sp3_svs.add(k)
    if hasattr(sp3, 'values'):
        continue
print(f"SP3 has {len(sp3_svs)} unique SVs in first 5 epochs")
print(f"Sample: {sorted(list(sp3_svs))[:10]}")

# Check if SP3 data is valid for May 1
sp3_t0 = sp3_ts[0]
sp3_dt0 = J2000 + timedelta(seconds=sp3_t0) if sp3_t0 > 1e7 else datetime(2000,1,1,12,0,0) + timedelta(seconds=sp3_t0*86400)
print(f"SP3 first epoch: {sp3_t0} -> {sp3_dt0}")
