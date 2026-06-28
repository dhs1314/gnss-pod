#!/usr/bin/env python3
"""Test: run PPP WITHOUT satellite clock, light-time, Sagnac — like batch_v12_fixed.py"""
import pickle, numpy as np, math
from datetime import datetime, timedelta
from collections import defaultdict
import sys
sys.path.insert(0, '.')
from src.sp3_loader import get_gps_pos_from_sp3

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
F_IF = (F1*F1 - F2*F2) / (F1 + F2)
LAM_IF = C / F_IF   # ~0.1070 m
LAM1 = C / F1
LAM2 = C / F2
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

# First: check if BATCH pickle has L1/L2 in cycles or meters
# The reference data has L1/L2 in cycles. Let's check our BATCH.
gps_sod_start = min(gps1b.keys())
first_sod = sorted(gps1b.keys())[0]
print("=== Data Format Check ===")
for sv_id, rec in sorted(gps1b[first_sod].items()):
    if 'L1_phase' not in rec: continue

    # Try interpreting L1_phase as meters
    L1_m = float(rec['L1_phase'])
    L2_m = float(rec['L2_phase'])

    # Try interpreting L1_phase as cycles
    L1_cyc = float(rec['L1_phase'])
    L2_cyc = float(rec['L2_phase'])

    # Get satellite range to check which interpretation makes sense
    utc_dt = J2000 + timedelta(seconds=first_sod)
    ref_pos = interpolate_ref(ref_orbit, first_sod)
    pos, clk, vel = get_gps_pos_from_sp3(sp3, sv_id, utc_dt)
    if pos is None: continue

    rho = float(np.linalg.norm(pos - ref_pos))

    # If L1 is in meters: L1_m should be ~ rho + clock_bias
    # If L1 is in cycles: L1_cyc * LAM1 should be ~ rho + clock_bias
    L1_from_cyc = L1_cyc * LAM1  # cycles → meters

    print(f"\n{sv_id}: rho={rho:.1f}m")
    print(f"  L1 as meters: {L1_m:.1f}m  (diff from rho: {L1_m - rho:.1f}m)")
    print(f"  L1 as cycles*lambda: {L1_from_cyc:.1f}m  (diff from rho: {L1_from_cyc - rho:.1f}m)")
    print(f"  L1_m / L1_from_cyc ratio: {L1_m / L1_from_cyc:.6f}")

    # Check L_if interpretation
    L_if_stored = float(rec['L_if'])
    L_if_from_m = ALPHA * L1_m + BETA * L2_m
    L_if_from_cyc = (ALPHA * L1_cyc + BETA * L2_cyc) * LAM_IF

    print(f"  L_if stored: {L_if_stored:.1f}")
    print(f"  L_if from meters: {L_if_from_m:.1f}")
    print(f"  L_if from cycles*LAM_IF: {L_if_from_cyc:.1f}")
    break

# Now test: run WITHOUT satellite clock, like the reference
print(f"\n=== Test WITHOUT satellite clock (like batch_v12_fixed.py) ===")

gps_sod_end = gps_sod_start + 4 * 3600
interval = 30.0

epochs_to_test = []
for gps_sod in sorted(gps1b.keys()):
    if not (gps_sod_start <= gps_sod <= gps_sod_end): continue
    dt_from_start = gps_sod - gps_sod_start
    nearest = round(dt_from_start / interval) * interval
    if abs(dt_from_start - nearest) > 2.0: continue
    if len(epochs_to_test) < 30:
        epochs_to_test.append(gps_sod)

# Pass 1: float ambiguity (no clock needed since phase-code diff)
sv_diff = defaultdict(list)
for gps_sod in epochs_to_test:
    utc_dt = J2000 + timedelta(seconds=gps_sod)
    for sv_id, rec in gps1b[gps_sod].items():
        # Only need L_if and P_if - no satellite position or clock
        sv_diff[sv_id].append(float(rec['L_if']) - float(rec['P_if']))

float_amb = {}
for sv, diffs in sv_diff.items():
    if len(diffs) >= 10:
        arr = np.array(diffs)
        float_amb[sv] = float(np.median(arr))

# Pass 2: Test WITHOUT satellite clock, WITHOUT light-time, WITHOUT Sagnac
print(f"\n{'='*70}")
print(f"epoch  | clk_est | trop_est | phase_res_std | code_res_std")
print(f"{'='*70}")

for gps_sod in epochs_to_test[:5]:
    utc_dt = J2000 + timedelta(seconds=gps_sod)
    ref_pos = interpolate_ref(ref_orbit, gps_sod)

    H_rows, y_rows, w_rows = [], [], []
    sv_info = []

    for sv_id, rec in gps1b[gps_sod].items():
        if sv_id not in float_amb: continue

        # SIMPLE: satellite position at RECEPTION time (no light-time iteration)
        pos, clk, vel = get_gps_pos_from_sp3(sp3, sv_id, utc_dt)
        if pos is None: continue

        # SIMPLE: Euclidean range (no Sagnac, no light-time)
        rho = float(np.linalg.norm(pos - ref_pos))
        if not (1.8e7 < rho < 2.8e7): continue

        el = math.asin(abs(pos[2] - ref_pos[2]) / rho)
        if el < 0.087: continue
        mf = 1.0 / max(math.sin(el), 0.05)

        L_if = float(rec['L_if'])   # meters
        P_if = float(rec['P_if'])   # meters

        # NO satellite clock! (like batch_v12_fixed.py)
        L_corr = L_if - float_amb[sv_id]  # just subtract ambiguity
        P_corr = P_if                      # code as-is

        h = np.array([1.0, mf])

        w_l = 1.0 / (0.003**2 / (math.sin(el)**2 + 0.01))
        H_rows.append(h); y_rows.append(L_corr - rho); w_rows.append(w_l)

        w_p = 1.0 / (0.3**2 / (math.sin(el)**2 + 0.01))
        H_rows.append(h); y_rows.append(P_corr - rho); w_rows.append(w_p)

        sv_info.append((sv_id, L_corr - rho, P_corr - rho))

    if len(H_rows) < 4: continue

    H = np.array(H_rows); y = np.array(y_rows); w_diag = np.array(w_rows)
    HtWH = H.T @ (w_diag[:, None] * H)
    HtWy = H.T @ (w_diag * y)
    dx = np.linalg.solve(HtWH, HtWy)
    clk_est, trop_est = dx[0], dx[1]

    model = H @ dx
    residuals = y - model
    phase_res = residuals[0::2]
    code_res = residuals[1::2]

    print(f"{gps_sod:.0f} | {clk_est:7.1f} | {trop_est:8.3f} | {np.std(phase_res):12.3f} | {np.std(code_res):11.3f}")

# Also test: with satellite clock BUT no light-time/Sagnac
print(f"\n=== Test WITH satellite clock, no light-time/Sagnac ===")
print(f"{'='*70}")

for gps_sod in epochs_to_test[:5]:
    utc_dt = J2000 + timedelta(seconds=gps_sod)
    ref_pos = interpolate_ref(ref_orbit, gps_sod)

    H_rows, y_rows, w_rows = [], [], []

    for sv_id, rec in gps1b[gps_sod].items():
        if sv_id not in float_amb: continue

        pos, clk, vel = get_gps_pos_from_sp3(sp3, sv_id, utc_dt)
        if pos is None: continue
        if abs(clk) > 0.1 * C: continue

        rho = float(np.linalg.norm(pos - ref_pos))
        if not (1.8e7 < rho < 2.8e7): continue

        el = math.asin(abs(pos[2] - ref_pos[2]) / rho)
        if el < 0.087: continue
        mf = 1.0 / max(math.sin(el), 0.05)

        L_if = float(rec['L_if'])
        P_if = float(rec['P_if'])

        # WITH satellite clock
        L_corr = L_if + clk - float_amb[sv_id]
        P_corr = P_if + clk

        h = np.array([1.0, mf])

        w_l = 1.0 / (0.003**2 / (math.sin(el)**2 + 0.01))
        H_rows.append(h); y_rows.append(L_corr - rho); w_rows.append(w_l)

        w_p = 1.0 / (0.3**2 / (math.sin(el)**2 + 0.01))
        H_rows.append(h); y_rows.append(P_corr - rho); w_rows.append(w_p)

    if len(H_rows) < 4: continue

    H = np.array(H_rows); y = np.array(y_rows); w_diag = np.array(w_rows)
    HtWH = H.T @ (w_diag[:, None] * H)
    HtWy = H.T @ (w_diag * y)
    dx = np.linalg.solve(HtWH, HtWy)
    clk_est, trop_est = dx[0], dx[1]

    model = H @ dx
    residuals = y - model
    phase_res = residuals[0::2]
    code_res = residuals[1::2]

    print(f"{gps_sod:.0f} | {clk_est:7.1f} | {trop_est:8.3f} | {np.std(phase_res):12.3f} | {np.std(code_res):11.3f}")
