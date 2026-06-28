#!/usr/bin/env python3
"""Check P_r/L_r values at first RINEX epoch vs ASCII"""
import sys; sys.path.insert(0, 'src')
import pickle
from gps1b_rnx_loader import load_gps1b_rnx
from sp3_loader import get_gps_pos_from_sp3
import numpy as np
from datetime import datetime, timedelta

C = 299792458.0
OMEGA_E = 7.2921151467e-5
J2000 = datetime(2000, 1, 1, 12, 0, 0)

# Load data
rnx_data = load_gps1b_rnx('data/GPS1B_2024-04-29_C_04.rnx')
ascii_data = pickle.load(open('data/gracefo/2024/2024-04-29/GPS1B_2024-04-29_C_04.pkl', 'rb'))
sp3 = pickle.load(open('data/2024/120/igs_sp3_FIN.pkl', 'rb'))

# Load GNV1B for ref position
ref_orbit = {}
for line in open('data/gracefo/2024/2024-04-29/GNV1B_2024-04-29_C_04.txt'):
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

# Pick first epoch
sod = 767620800
utc_dt = J2000 + timedelta(seconds=sod)
ref_pos = get_ref_pos(sod)
print(f"Epoch: {sod} UTC={utc_dt}")
print(f"Ref pos: {ref_pos}")

common_svs = sorted(set(rnx_data[sod].keys()) & set(ascii_data[sod].keys()))
print(f"\n{'SV':<6s} {'P_r(RNX)':>12s} {'P_r(ASC)':>12s} {'L_r(RNX)':>12s} {'L_r(ASC)':>12s} {'dP_r':>10s} {'dL_r':>10s}")
print("-" * 78)

for sv in common_svs:
    r = rnx_data[sod][sv]
    a = ascii_data[sod][sv]

    pos, clk, vel = get_gps_pos_from_sp3(sp3, sv, utc_dt)
    if pos is None: continue

    rho = float(np.linalg.norm(pos - ref_pos))
    sag = (OMEGA_E / C) * (pos[0] * ref_pos[1] - pos[1] * ref_pos[0])
    rho_corr = rho + sag

    P_r_rnx = r['P_if'] + clk - rho_corr
    P_r_asc = a['P_if'] + clk - rho_corr
    L_r_rnx = r['L_if'] + clk - rho_corr
    L_r_asc = a['L_if'] + clk - rho_corr

    print(f"{sv:<6s} {P_r_rnx:>12.1f} {P_r_asc:>12.1f} {L_r_rnx:>12.1f} {L_r_asc:>12.1f} {P_r_rnx-P_r_asc:>10.1f} {L_r_rnx-L_r_asc:>10.1f}")

# Check median per-SV
print(f"\n--- Per-SV median P_r (bias) ---")
print(f"{'SV':<6s} {'RNX median':>12s} {'ASC median':>12s} {'Diff':>10s}")
for sv in common_svs:
    rnx_Prs = []
    asc_Prs = []
    for s in sorted(rnx_data.keys())[:60]:
        if s not in rnx_data or sv not in rnx_data[s]: continue
        pos, clk, vel = get_gps_pos_from_sp3(sp3, sv, J2000 + timedelta(seconds=s))
        if pos is None: continue
        ref = get_ref_pos(s)
        rho = float(np.linalg.norm(pos - ref))
        sag = (OMEGA_E / C) * (pos[0] * ref[1] - pos[1] * ref[0])
        rho_corr = rho + sag
        rnx_Prs.append(rnx_data[s][sv]['P_if'] + clk - rho_corr)
    for s in sorted(ascii_data.keys())[:60]:
        if s not in ascii_data or sv not in ascii_data[s]: continue
        pos, clk, vel = get_gps_pos_from_sp3(sp3, sv, J2000 + timedelta(seconds=s))
        if pos is None: continue
        ref = get_ref_pos(s)
        rho = float(np.linalg.norm(pos - ref))
        sag = (OMEGA_E / C) * (pos[0] * ref[1] - pos[1] * ref[0])
        rho_corr = rho + sag
        asc_Prs.append(ascii_data[s][sv]['P_if'] + clk - rho_corr)

    if rnx_Prs and asc_Prs:
        rnx_med = np.median(rnx_Prs)
        asc_med = np.median(asc_Prs)
        print(f"{sv:<6s} {rnx_med:>12.1f} {asc_med:>12.1f} {rnx_med-asc_med:>10.1f}")
