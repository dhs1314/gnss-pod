#!/usr/bin/env python3
"""Diagnose: verify SP3 clock values and correction"""
import pickle, numpy as np, math
from datetime import datetime, timedelta
import sys
sys.path.insert(0, '.')
from src.sp3_loader import get_gps_pos_from_sp3

C = 299792458.0
J2000 = datetime(2000, 1, 1, 12, 0, 0)

sp3 = pickle.load(open("data/2024/120/igs_sp3_FIN.pkl", "rb"))
gps1b = pickle.load(open("data/gracefo/2024/2024-04-29/GPS1B_2024-04-29_C_04.pkl", "rb"))

# Check SP3 epoch range
print("SP3 info:")
print(f"  epochs: {len(sp3['ts'])}")
print(f"  time range: {sp3['ts'][0]} to {sp3['ts'][-1]}")
print(f"  dt between epochs: {(sp3['ts'][1] - sp3['ts'][0]).total_seconds()}s")
print()

# Check GPS1B epoch vs SP3 epochs
gps_sod = sorted(gps1b.keys())[0]
utc_dt = J2000 + timedelta(seconds=gps_sod)
print(f"GPS1B first epoch: gps_sod={gps_sod}, utc={utc_dt}")

# Find SP3 epochs near this time
ts = sp3['ts']
for i, t in enumerate(ts):
    if t >= utc_dt:
        print(f"  SP3 epoch before: {ts[i-1]} (dt={(utc_dt-ts[i-1]).total_seconds():.0f}s)")
        print(f"  SP3 epoch after:  {t} (dt={(t-utc_dt).total_seconds():.0f}s)")
        break

# Load GNV1B
ref_orbit = {}
for l in open("data/gracefo/2024/2024-04-29/GNV1B_2024-04-29_C_04.txt"):
    p = l.split()
    if len(p) < 6: continue
    try:
        t = float(p[0]); flag = p[2]
        if flag in ("C", "E"): ref_orbit[t] = np.array([float(p[3]), float(p[4]), float(p[5])])
    except: pass

ref_pos = ref_orbit[gps_sod]  # exact epoch match
print(f"\nGNV1B at {gps_sod}: {ref_pos/1000} km")

# For each SV: compare L_if - rho (no clock) vs -clk
print(f"\n{'SV':5s} {'clk_m':>10s} {'L_if-rho':>10s} {'L-rho_corr':>10s} {'P_if-rho':>10s} {'P-rho_corr':>10s} {'check':>10s}")
print("-" * 70)
for sv_id, rec in sorted(gps1b[gps_sod].items()):
    pos_rx, clk_rx, vel_rx = get_gps_pos_from_sp3(sp3, sv_id, utc_dt)
    if pos_rx is None: continue
    if abs(clk_rx) > 0.1 * C: continue

    # Simple range (no light-time, no Sagnac)
    rho = float(np.linalg.norm(pos_rx - ref_pos))
    if not (1.8e7 < rho < 2.8e7): continue

    L_if = float(rec['L_if'])
    P_if = float(rec['P_if'])

    L_minus_rho = L_if - rho
    P_minus_rho = P_if - rho

    L_corr = L_if + clk_rx - rho  # with clock
    P_corr = P_if + clk_rx - rho

    # If clk is correct, L_corr should be ~constant across SVs
    # And L_minus_rho should ≈ -clk_rx
    check = L_minus_rho + clk_rx  # should be same as L_corr

    print(f"{sv_id:5s} {clk_rx:10.1f} {L_minus_rho:10.1f} {L_corr:10.2f} {P_minus_rho:10.1f} {P_corr:10.2f} {check:10.2f}")

# Check the SP3 clock interpolation accuracy at a non-nodal point
print(f"\n--- Clock interpolation test: compare exact-node vs interpolated ---")
# Find an SP3 epoch that's an exact node
test_sv = 'G05'
for t in sp3['ts'][:3]:
    pos, clk, vel = get_gps_pos_from_sp3(sp3, test_sv, t)
    epoch_data = sp3['epochs'][t]
    if test_sv in epoch_data:
        clk_direct = float(epoch_data[test_sv][3])
        print(f"  {test_sv} at {t}: interpolated clk={clk:.6f}m, direct={clk_direct:.6f}m, diff={clk-clk_direct:.10f}m")

# Check at mid-point between two SP3 epochs
t_mid = sp3['ts'][0] + (sp3['ts'][1] - sp3['ts'][0]) / 2
pos_mid, clk_mid, vel_mid = get_gps_pos_from_sp3(sp3, test_sv, t_mid)
# Compare with linear interpolation of clock
t0 = sp3['ts'][0]; t1 = sp3['ts'][1]
p0 = sp3['epochs'][t0].get(test_sv)
p1 = sp3['epochs'][t1].get(test_sv)
if p0 and p1:
    clk_lin = p0[3] + (p1[3] - p0[3]) * (t_mid - t0).total_seconds() / (t1 - t0).total_seconds()
    print(f"  {test_sv} at mid {t_mid}: barycentric={clk_mid:.6f}m, linear={clk_lin:.6f}m, diff={clk_mid-clk_lin:.6f}m")
