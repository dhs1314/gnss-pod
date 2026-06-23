#!/usr/bin/env python3
"""Diagnose: compare code-level L_if vs phase-level L_if from raw L1/L2 phase"""
import pickle, numpy as np, math
from datetime import datetime, timedelta
from collections import defaultdict
import sys
sys.path.insert(0, '.')
from src.sp3_loader import get_gps_pos_from_sp3

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
F1_SQ, F2_SQ = F1*F1, F2*F2
ALPHA = F1_SQ / (F1_SQ - F2_SQ)
BETA = -F2_SQ / (F1_SQ - F2_SQ)
LAM_W = C / (F1 - F2)
OMEGA_E = 7.2921151467e-5
J2000 = datetime(2000, 1, 1, 12, 0, 0)

sp3 = pickle.load(open("data/2024/120/igs_sp3_FIN.pkl", "rb"))

# Use gps1b cache (from gps1b_loader, rebuilt)
gps1b = pickle.load(open("data/gracefo/2024/2024-04-29/gps1b_C.pkl", "rb"))

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

gps_sod = sorted(gps1b.keys())[0]
utc_dt = J2000 + timedelta(seconds=gps_sod)
ref_pos = interpolate_ref(ref_orbit, gps_sod)

print("=== Comparing code-level L_if vs phase-level L_if ===")
print(f"Epoch: gps_sod={gps_sod}")
print()
print(f"{'SV':5s} {'L_if_code':>12s} {'L_if_phase':>12s} {'diff_m':>8s} {'P_if':>12s}")
print(f"{'':5s} {'(from L*_range)':>12s} {'(from L*_phase)':>12s} {'':>8s} {'(from P*)':>12s}")
print("-" * 60)

for sv_id, rec in sorted(gps1b[gps_sod].items()):
    if 'L1_phase' not in rec: continue

    L1_phase_m = float(rec['L1_phase'])
    L2_phase_m = float(rec['L2_phase'])
    L1_range_m = float(rec['L1_range'])
    L2_range_m = float(rec['L2_range'])

    # Code-level IF (what the gps1b loader computes: from L*_range)
    L_if_code = ALPHA * L1_range_m + BETA * L2_range_m

    # Phase-level IF (from L*_phase — true carrier phase in meters)
    L_if_phase = ALPHA * L1_phase_m + BETA * L2_phase_m

    P_if_val = float(rec['P_if'])

    print(f"{sv_id:5s} {L_if_code:12.3f} {L_if_phase:12.3f} {L_if_phase-L_if_code:8.4f} {P_if_val:12.3f}")

    if sv_id >= 'G07': break

# Check across epochs: is L_if_phase - L_if_code CONSTANT?
print(f"\n=== Checking phase-code difference STABILITY ===")

# Get 30 epochs
gps_sod_start = min(gps1b.keys())
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

# Track L_if_phase - L_if_code for G05 and G07 across epochs
for sv_track in ['G05', 'G07', 'G16']:
    phase_diffs = []
    code_diffs = []
    for gps_sod in epochs_to_test:
        if sv_track in gps1b[gps_sod]:
            rec = gps1b[gps_sod][sv_track]
            if 'L1_phase' not in rec: continue

            L1_phase_m = float(rec['L1_phase'])
            L2_phase_m = float(rec['L2_phase'])
            L1_range_m = float(rec['L1_range'])
            L2_range_m = float(rec['L2_range'])
            L1 = float(rec['L1'])  # = L1_range (preferred)
            L2 = float(rec['L2'])

            L_if_phase = ALPHA * L1_phase_m + BETA * L2_phase_m
            L_if_code = ALPHA * L1_range_m + BETA * L2_range_m
            P_if_val = float(rec['P_if'])

            phase_diffs.append(L_if_phase - P_if_val)
            code_diffs.append(L_if_code - P_if_val)

    if phase_diffs:
        p_arr = np.array(phase_diffs)
        c_arr = np.array(code_diffs)
        print(f"{sv_track}:")
        print(f"  L_if_phase - P_if: median={np.median(p_arr):.4f}m std={np.std(p_arr):.4f}m n={len(p_arr)}")
        print(f"  L_if_code  - P_if: median={np.median(c_arr):.4f}m std={np.std(c_arr):.4f}m n={len(c_arr)}")
        print(f"  phase-noise advantage: std_ratio = {np.std(c_arr)/np.std(p_arr):.1f}x")

# The KEY test: compute L-rho with phase-level L_if and see if residuals are better
print(f"\n=== Phase-level IF: fixed-position residuals ===")
for gps_sod in epochs_to_test[:3]:
    utc_dt = J2000 + timedelta(seconds=gps_sod)
    ref_pos = interpolate_ref(ref_orbit, gps_sod)

    code_vals = []
    phase_vals = []
    for sv_id, rec in gps1b[gps_sod].items():
        if 'L1_phase' not in rec: continue

        pos, clk, vel = get_gps_pos_from_sp3(sp3, sv_id, utc_dt)
        if pos is None or abs(clk) > 0.1 * C: continue

        # Light-time iteration
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

        sag = (OMEGA_E / C) * (pos[0] * ref_pos[1] - pos[1] * ref_pos[0])
        rho_corr = rho + sag

        L1_phase_m = float(rec['L1_phase'])
        L2_phase_m = float(rec['L2_phase'])
        L1_range_m = float(rec['L1_range'])
        L2_range_m = float(rec['L2_range'])

        L_if_phase = ALPHA * L1_phase_m + BETA * L2_phase_m
        L_if_code = ALPHA * L1_range_m + BETA * L2_range_m
        P_if = float(rec['P_if'])

        # Phase-level L-rho
        L_rho_phase = L_if_phase + clk - rho_corr
        # Code-level L-rho
        L_rho_code = L_if_code + clk - rho_corr
        # Code P-rho
        P_rho = P_if + clk - rho_corr

        code_vals.append(P_rho)
        phase_vals.append((sv_id, L_rho_phase, L_rho_code, P_rho))

    if not code_vals: continue
    simple_clk = np.median(code_vals)

    print(f"\nepoch {gps_sod}: simple_clk={simple_clk:.3f}m")
    for sv_id, L_rp, L_rc, P_r in phase_vals:
        print(f"  {sv_id}: phase-L-rho={L_rp:.3f}m code-L-rho={L_rc:.3f}m P-rho={P_r:.3f}m  "
              f"phase_res={L_rp-simple_clk:.3f}m code_res={L_rc-simple_clk:.3f}m")
