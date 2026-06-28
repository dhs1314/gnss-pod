#!/usr/bin/env python3
"""Verify Lagrange interpolation fixes the range-rate discrepancy"""
import pickle, numpy as np
from datetime import datetime, timedelta

C = 299792458.0
F1 = 1575.42e6
F2 = 1227.60e6
F1_SQ = F1 * F1
F2_SQ = F2 * F2
ALPHA = F1_SQ / (F1_SQ - F2_SQ)
BETA = -F2_SQ / (F1_SQ - F2_SQ)
J2000 = datetime(2000, 1, 1, 12, 0, 0)

sp3 = pickle.load(open("data/2024/120/igs_sp3_FIN.pkl", "rb"))
gps1b = pickle.load(open("data/gracefo/2024/2024-04-29/GPS1B_2024-04-29_C_04.pkl", "rb"))

# GNV1B
ref_orbit = {}
for l in open("data/gracefo/2024/2024-04-29/GNV1B_2024-04-29_C_04.txt"):
    p = l.split()
    if len(p) < 6:
        continue
    try:
        t = float(p[0])
        flag = p[2]
        if flag in ("C", "E"):
            ref_orbit[t] = np.array([float(p[3]), float(p[4]), float(p[5])])
    except:
        pass

# NEW Lagrange interpolation
from src.sp3_loader import get_gps_pos_from_sp3

sv = "G05"
gps_sods = sorted(gps1b.keys())
gps0 = gps_sods[0]
utc0 = J2000 + timedelta(seconds=gps0)
utc1 = J2000 + timedelta(seconds=gps_sods[1])

# GNV1B interpolation helper
ts_ref = sorted(ref_orbit.keys())

def interp_gnv1b(gps_sod):
    t0 = t1 = None
    for j, tj in enumerate(ts_ref):
        if tj >= gps_sod:
            t1 = tj
            t0 = ts_ref[j - 1] if j > 0 else None
            break
        t0 = tj
    if t1 is None:
        t0 = t1 = ts_ref[-1]
    if t0 is None:
        t0 = ts_ref[0]
    if t0 == t1:
        return ref_orbit[t0]
    a = (gps_sod - t0) / (t1 - t0)
    return ref_orbit[t0] * (1 - a) + ref_orbit[t1] * a

# OLD linear interpolation (for comparison)
def old_get_sat_pos(sp3_data, sv, utc_dt):
    ts = sp3_data['ts']
    epochs = sp3_data['epochs']
    t0 = t1 = None
    for i, ti in enumerate(ts):
        if ti >= utc_dt:
            t1 = ti
            t0 = ts[i - 1] if i > 0 else None
            break
        t0 = ti
    if t1 is None:
        t0 = t1 = ts[-1]
    if t0 is None:
        t0 = ts[0]
    p0 = epochs.get(t0, {}).get(sv)
    p1 = epochs.get(t1, {}).get(sv)
    dt_sec = (utc_dt - t0).total_seconds()
    dt_tot = (t1 - t0).total_seconds()
    if dt_tot == 0:
        return np.array(p0[:3]), float(p0[3])
    a = dt_sec / dt_tot
    pos = np.array(p0[:3]) * (1 - a) + np.array(p1[:3]) * a
    clk = p0[3] * (1 - a) + p1[3] * a
    return pos, clk

# ===== EPOCH 0 =====
pos0_new, clk0_new, vel0_new = get_gps_pos_from_sp3(sp3, sv, utc0)
pos0_old, clk0_old = old_get_sat_pos(sp3, sv, utc0)
g0 = interp_gnv1b(gps0)
rng0_new = np.linalg.norm(pos0_new - g0)
rng0_old = np.linalg.norm(pos0_old - g0)

print("=== Epoch 0 ===")
print(f"Lagrange pos: {pos0_new}")
print(f"Linear   pos: {pos0_old}")
print(f"Pos diff: {np.linalg.norm(pos0_new - pos0_old):.3f} m")
print(f"Lagrange vel: {vel0_new}  |v|={np.linalg.norm(vel0_new):.1f} m/s")
print(f"Lagrange rng: {rng0_new/1000:.3f} km")
print(f"Linear   rng: {rng0_old/1000:.3f} km")

# ===== EPOCH 1 (10s later) =====
pos1_new, clk1_new, vel1_new = get_gps_pos_from_sp3(sp3, sv, utc1)
pos1_old, clk1_old = old_get_sat_pos(sp3, sv, utc1)
g1 = interp_gnv1b(gps_sods[1])
rng1_new = np.linalg.norm(pos1_new - g1)
rng1_old = np.linalg.norm(pos1_old - g1)
dt_gps = gps_sods[1] - gps0

print(f"\n=== Epoch 1 (dt={dt_gps}s) ===")
print(f"Lagrange pos: {pos1_new}")
print(f"Linear   pos: {pos1_old}")
print(f"Pos diff: {np.linalg.norm(pos1_new - pos1_old):.3f} m")
print(f"Lagrange vel: {vel1_new}  |v|={np.linalg.norm(vel1_new):.1f} m/s")
print(f"Lagrange rng: {rng1_new/1000:.3f} km")
print(f"Linear   rng: {rng1_old/1000:.3f} km")

# ===== RATE CHECK =====
rng_rate_new = (rng1_new - rng0_new) / dt_gps
rng_rate_old = (rng1_old - rng0_old) / dt_gps

L1_0 = gps1b[gps0][sv]["L1_phase"]
L1_1 = gps1b[gps_sods[1]][sv]["L1_phase"]
L1_rate = (L1_1 - L1_0) / dt_gps

print(f"\n=== Rate Comparison ===")
print(f"L1_phase rate:              {L1_rate:10.1f} m/s")
print(f"Geometric range rate (OLD): {rng_rate_old:10.1f} m/s  diff: {L1_rate - rng_rate_old:.1f} m/s")
print(f"Geometric range rate (NEW): {rng_rate_new:10.1f} m/s  diff: {L1_rate - rng_rate_new:.1f} m/s")

# Apply satellite clock to L1_phase for proper comparison
# L1_phase + clk_s should be close to geometric_range + clk_r + iono
# The rate of (L1_phase + clk_s) should match geometric_range rate closely
L1_clk_rate_new = (L1_1 + clk1_new - (L1_0 + clk0_new)) / dt_gps
L1_clk_rate_old = (L1_1 + clk1_old - (L1_0 + clk0_old)) / dt_gps
print(f"\nL1+clk rate (NEW):          {L1_clk_rate_new:10.1f} m/s  diff: {L1_clk_rate_new - rng_rate_new:.1f} m/s")
print(f"L1+clk rate (OLD):          {L1_clk_rate_old:10.1f} m/s  diff: {L1_clk_rate_old - rng_rate_old:.1f} m/s")

print(f"\n=== Improvement ===")
improvement = abs(L1_rate - rng_rate_old) - abs(L1_rate - rng_rate_new)
print(f"Rate discrepancy reduced by: {improvement:.1f} m/s")
print(f"Lagrange velocity magnitude: {np.linalg.norm(vel0_new):.1f} m/s (expected ~3874 for GPS MEO)")
