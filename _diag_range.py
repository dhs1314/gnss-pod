#!/usr/bin/env python3
"""Diagnose: compute per-SV range error budget including light-time, Sagnac, clock"""
import pickle, numpy as np, math
from datetime import datetime, timedelta
from collections import defaultdict
import sys
sys.path.insert(0, '.')
from src.sp3_loader import get_gps_pos_from_sp3

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
OMEGA_E = 7.2921151467e-5
J2000 = datetime(2000, 1, 1, 12, 0, 0)
ALPHA = F1*F1 / (F1*F1 - F2*F2)
BETA = -F2*F2 / (F1*F1 - F2*F2)

sp3 = pickle.load(open("data/2024/120/igs_sp3_FIN.pkl", "rb"))
gps1b = pickle.load(open("data/gracefo/2024/2024-04-29/GPS1B_2024-04-29_C_04.pkl", "rb"))

ref_orbit = {}
for l in open("data/gracefo/2024/2024-04-29/GNV1B_2024-04-29_C_04.txt"):
    p = l.split()
    if len(p) < 6: continue
    try:
        t = float(p[0]); flag = p[2]
        if flag in ("C", "E"): ref_orbit[t] = np.array([float(p[3]), float(p[4]), float(p[5])])
    except: pass

def interpolate_ref(orbit, gps_sod):
    ts = sorted(orbit.keys())
    t0 = t1 = None
    for i, ti in enumerate(ts):
        if ti >= gps_sod: t1 = ti; t0 = ts[i-1] if i > 0 else None; break
        t0 = ti
    if t1 is None: t0 = t1 = ts[-1]
    if t0 is None: t0 = ts[0]
    if t0 == t1: return orbit[t0]
    a = (gps_sod - t0) / (t1 - t0)
    return orbit[t0] * (1-a) + orbit[t1] * a

# Pick first epoch and trace through every SV
gps_sod = sorted(gps1b.keys())[0]
utc_dt = J2000 + timedelta(seconds=gps_sod)
ref_pos = interpolate_ref(ref_orbit, gps_sod)

print(f"Epoch: gps_sod={gps_sod}")
print(f"GNV1B ref_pos: [{ref_pos[0]/1000:.3f}, {ref_pos[1]/1000:.3f}, {ref_pos[2]/1000:.3f}] km")
print(f"GNV1B radius: {np.linalg.norm(ref_pos)/1000:.3f} km")
print()

header = f"{'SV':5s} {'clk_km':>8s} {'rho_rx_km':>9s} {'rho_tx_km':>9s} {'lt_diff_m':>9s} {'sag_m':>7s} {'L-rho_m':>8s} {'P-rho_m':>8s} {'L-P_m':>7s}"
print(header)
print("-" * len(header))

for sv_id, rec in sorted(gps1b[gps_sod].items()):
    # 1. Satellite position at RECEPTION time (no iteration)
    pos_rx, clk_rx, vel_rx = get_gps_pos_from_sp3(sp3, sv_id, utc_dt)
    if pos_rx is None: continue
    if abs(clk_rx) > 0.1 * C: continue

    rho_rx = float(np.linalg.norm(pos_rx - ref_pos))
    if not (1.8e7 < rho_rx < 2.8e7): continue

    # 2. Satellite position at TRANSMISSION time (light-time iteration)
    pos, clk, vel = pos_rx, clk_rx, vel_rx
    rho = rho_rx
    for _ in range(5):
        travel_time = rho / C
        tx_dt = utc_dt - timedelta(seconds=travel_time)
        pos_tx, clk_tx, vel_tx = get_gps_pos_from_sp3(sp3, sv_id, tx_dt)
        if pos_tx is None: break
        rho_new = float(np.linalg.norm(pos_tx - ref_pos))
        if abs(rho_new - rho) < 1e-8:
            pos, clk, vel = pos_tx, clk_tx, vel_tx
            rho = rho_new
            break
        pos, clk, vel = pos_tx, clk_tx, vel_tx
        rho = rho_new

    rho_tx = rho
    lt_diff = rho_rx - rho_tx  # positive = rho_rx is larger
    sag = (OMEGA_E / C) * (pos[0] * ref_pos[1] - pos[1] * ref_pos[0])
    rho_corr = rho_tx + sag

    L_if = float(rec['L_if'])
    P_if = float(rec['P_if'])

    L_rho = L_if + clk - rho_corr
    P_rho = P_if + clk - rho_corr

    print(f"{sv_id:5s} {clk/1000:8.1f} {rho_rx/1000:9.3f} {rho_tx/1000:9.3f} {lt_diff:9.2f} {sag:7.2f} {L_rho:8.2f} {P_rho:8.2f} {L_if-P_if:7.3f}")

# Now check: do the SP3 clock values change significantly between
# reception time and transmission time? (they shouldn't, much)
print(f"\n--- Clock change during light-time (tx - rx) ---")
for sv_id, rec in sorted(gps1b[gps_sod].items()):
    pos_rx, clk_rx, vel_rx = get_gps_pos_from_sp3(sp3, sv_id, utc_dt)
    if pos_rx is None or abs(clk_rx) > 0.1 * C: continue
    rho_rx = float(np.linalg.norm(pos_rx - ref_pos))
    if not (1.8e7 < rho_rx < 2.8e7): continue

    travel_time = rho_rx / C
    tx_dt = utc_dt - timedelta(seconds=travel_time)
    pos_tx, clk_tx, vel_tx = get_gps_pos_from_sp3(sp3, sv_id, tx_dt)
    if pos_tx is None: continue

    print(f"  {sv_id}: clk_rx={clk_rx:.3f}m clk_tx={clk_tx:.3f}m diff={clk_tx-clk_rx:.4f}m")

# Also check: what if we DON'T apply Sagnac?
print(f"\n--- WITHOUT Sagnac correction ---")
for sv_id, rec in sorted(gps1b[gps_sod].items()):
    pos_rx, clk_rx, vel_rx = get_gps_pos_from_sp3(sp3, sv_id, utc_dt)
    if pos_rx is None or abs(clk_rx) > 0.1 * C: continue

    pos, clk, vel = pos_rx, clk_rx, vel_rx
    rho = float(np.linalg.norm(pos - ref_pos))
    for _ in range(5):
        travel_time = rho / C
        tx_dt = utc_dt - timedelta(seconds=travel_time)
        pos_tx, clk_tx, vel_tx = get_gps_pos_from_sp3(sp3, sv_id, tx_dt)
        if pos_tx is None: break
        rho_new = float(np.linalg.norm(pos_tx - ref_pos))
        if abs(rho_new - rho) < 1e-8:
            pos, clk, vel = pos_tx, clk_tx, vel_tx
            rho = rho_new
            break
        pos, clk, vel = pos_tx, clk_tx, vel_tx
        rho = rho_new
    if not (1.8e7 < rho < 2.8e7): continue

    L_if = float(rec['L_if']); P_if = float(rec['P_if'])
    L_rho_no_sag = L_if + clk - rho
    P_rho_no_sag = P_if + clk - rho
    sag = (OMEGA_E / C) * (pos[0] * ref_pos[1] - pos[1] * ref_pos[0])
    print(f"  {sv_id}: L-rho={L_rho_no_sag:.2f}m P-rho={P_rho_no_sag:.2f}m (sag was {sag:.2f}m)")

# Check SP3 clock values directly
print(f"\n--- SP3 clock values at reception time ---")
for sv_id, rec in sorted(gps1b[gps_sod].items()):
    pos_rx, clk_rx, vel_rx = get_gps_pos_from_sp3(sp3, sv_id, utc_dt)
    if pos_rx is None: continue
    print(f"  {sv_id}: clk={clk_rx:.1f}m ({clk_rx/C*1e6:.3f} us)")
