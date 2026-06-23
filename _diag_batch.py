#!/usr/bin/env python3
"""Diagnose: check if BATCH pickle L_if/P_if are correctly formed from raw L1/L2"""
import pickle, numpy as np, math
from datetime import datetime, timedelta
import sys
sys.path.insert(0, '.')
from src.sp3_loader import get_gps_pos_from_sp3

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
F1_SQ, F2_SQ = F1*F1, F2*F2
ALPHA = F1_SQ / (F1_SQ - F2_SQ)
BETA = -F2_SQ / (F1_SQ - F2_SQ)
LAM_IF_STD = C / (F1*F1 - F2*F2) * (F1 + F2)  # standard IF wavelength
F_IF_NL = (F1*F1 - F2*F2) / (F1 + F2)  # = F1 - F2
LAM_WL = C / (F1 - F2)

print(f"ALPHA={ALPHA:.6f}, BETA={BETA:.6f}, ALPHA+BETA={ALPHA+BETA:.10f}")
print(f"LAM_IF (standard narrow-lane) = {LAM_IF_STD:.6f} m")
print(f"F_IF (computed as F1-F2) = {F_IF_NL:.0f} Hz")
print(f"LAM_IF (if using F1-F2) = {C/F_IF_NL:.6f} m  (wide-lane wavelength!)")
print(f"LAM_WL = {C/(F1-F2):.6f} m")
print()

sp3 = pickle.load(open("data/2024/120/igs_sp3_FIN.pkl", "rb"))
gps1b = pickle.load(open("data/gracefo/2024/2024-04-29/GPS1B_2024-04-29_C_04.pkl", "rb"))

gps_sod = sorted(gps1b.keys())[0]
utc_dt = J2000 = datetime(2000, 1, 1, 12, 0, 0) + timedelta(seconds=gps_sod)

# Load GNV1B
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

ref_pos = interpolate_ref(ref_orbit, gps_sod)
C = 299792458.0

# For each SV, compute L_if TWO WAYS:
# Method A: Use the pre-computed L_if from the BATCH pickle
# Method B: Compute L_if from raw L1/L2 phase using standard IF formula
# Method C: Check IF computed as (F1²*L1 - F2²*L2)/(F1²-F2²) — what gps1b_loader does

print(f"{'SV':5s} {'L_if_A':>12s} {'L_if_B':>12s} {'A-B':>8s} {'P_if_A':>12s} {'P_if_B':>12s} {'A-B':>8s}")
print(f"{'':5s} {'(from pickle)':>12s} {'(raw L1/L2)':>12s} {'':>8s} {'(from pickle)':>12s} {'(raw P1/P2)':>12s} {'':>8s}")
print("-" * 80)

for sv_id, rec in sorted(gps1b[gps_sod].items()):
    if 'L1_phase' not in rec: continue

    L1_m = float(rec['L1_phase'])
    L2_m = float(rec['L2_phase'])
    P1_m = float(rec.get('P1', 0))
    P2_m = float(rec.get('P2', 0))

    # Method A: pre-computed
    L_if_A = float(rec['L_if'])
    P_if_A = float(rec['P_if'])

    # Method B: standard IF from raw L1/L2 in meters
    L_if_B = ALPHA * L1_m + BETA * L2_m
    P_if_B = ALPHA * P1_m + BETA * P2_m

    # Method C: what gps1b_loader.py does: (f1²*L1 - f2²*L2)/(f1²-f2²)
    # This is identical to ALPHA*L1 + BETA*L2 since ALPHA = f1²/(f1²-f2²)
    L_if_C = (F1_SQ * L1_m - F2_SQ * L2_m) / (F1_SQ - F2_SQ)

    print(f"{sv_id:5s} {L_if_A:12.3f} {L_if_B:12.3f} {L_if_A-L_if_B:8.4f} {P_if_A:12.3f} {P_if_B:12.3f} {P_if_A-P_if_B:8.4f}")

    if sv_id >= 'G07': break

# Now: check if the L1/L2 phase values are SMOOTHED RANGE (in meters)
# or ACTUAL PHASE (in cycles). Check by looking at L1 and P1:
print(f"\n--- Check L1 vs P1 (phase vs code) ---")
for sv_id, rec in sorted(gps1b[gps_sod].items()):
    if 'L1_phase' not in rec: continue

    L1_m = float(rec['L1_phase'])
    P1_m = float(rec.get('P1', 0))

    # For smoothed range, L1 ≈ P1 (within code noise)
    # For raw phase, L1 ≈ P1 + N*lambda (differ by integer cycles * 0.19m)
    diff = L1_m - P1_m
    diff_cycles = diff / (C/F1)

    print(f"  {sv_id}: L1={L1_m:.3f}m P1={P1_m:.3f}m diff={diff:.3f}m ({diff_cycles:.1f} L1 cycles)")

    if sv_id >= 'G07': break

# Also check: what's the "true" IF wavelength for converting IF cycles to meters?
# IF combination: L_if = α*L1 + β*L2 where L1, L2 in meters
# IF in cycles: L_if_cyc = α*L1/λ1 + β*L2/λ2 = α*L1*f1/c + β*L2*f2/c
# L_if_m / L_if_cyc = ???
print(f"\n--- IF wavelength analysis ---")
for sv_id, rec in sorted(gps1b[gps_sod].items()):
    if 'L1_phase' not in rec: continue

    L1_m = float(rec['L1_phase'])
    L2_m = float(rec['L2_phase'])

    # Treat as meters
    L_if_m = ALPHA * L1_m + BETA * L2_m

    # Treat as cycles
    L1_cyc = L1_m  # just reinterpret the same number as cycles
    L2_cyc = L2_m
    # Convert cycles → meters properly
    L_if_from_cyc_m = ALPHA * L1_cyc * (C/F1) + BETA * L2_cyc * (C/F2)
    # Or: treat L_if_cyc = ALPHA*L1_cyc + BETA*L2_cyc, then × λ_IF
    L_if_cyc = ALPHA * L1_cyc + BETA * L2_cyc

    # What λ would make L_if_cyc * λ = L_if_m ?
    implied_lam = L_if_m / L_if_cyc if abs(L_if_cyc) > 1 else 0

    print(f"  {sv_id}: L_if_m={L_if_m:.3f}  L_if_cyc={L_if_cyc:.3f}  implied_λ={implied_lam:.6f}m")
    print(f"          L_if_from_cyc_m={L_if_from_cyc_m:.3f} (correct conversion)")
    print(f"          ratio L_if_m/L_if_from_cyc_m = {L_if_m/L_if_from_cyc_m:.6f}")
    break
