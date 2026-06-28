#!/usr/bin/env python3
"""Test: per-SV bias estimation to reduce code-level residuals"""
import pickle, numpy as np, math
from datetime import datetime, timedelta
from collections import defaultdict
import sys
sys.path.insert(0, '.')
from src.sp3_loader import get_gps_pos_from_sp3

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
LAM_W = C / (F1 - F2)
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

# Select test epochs
gps_sod_start = min(gps1b.keys())
gps_sod_end = gps_sod_start + 4 * 3600
interval = 30.0

epochs_to_test = []
for gps_sod in sorted(gps1b.keys()):
    if not (gps_sod_start <= gps_sod <= gps_sod_end): continue
    dt_from_start = gps_sod - gps_sod_start
    nearest = round(dt_from_start / interval) * interval
    if abs(dt_from_start - nearest) > 2.0: continue
    if len(epochs_to_test) < 60:
        epochs_to_test.append(gps_sod)

# Step 1: Estimate float ambiguities (phase-code diff) and per-SV biases
# Collect all L-P differences per SV
sv_lp_diff = defaultdict(list)
sv_sagnac = defaultdict(list)
sv_elevation = defaultdict(list)

for gps_sod in epochs_to_test:
    utc_dt = J2000 + timedelta(seconds=gps_sod)
    ref_pos = interpolate_ref(ref_orbit, gps_sod)

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
        el = math.asin(abs(pos[2] - ref_pos[2]) / rho)

        L1_phase_m = float(rec['L1_phase'])
        L2_phase_m = float(rec['L2_phase'])
        L_if = ALPHA * L1_phase_m + BETA * L2_phase_m
        P_if = float(rec['P_if'])

        L_r = L_if + clk - rho_corr
        P_r = P_if + clk - rho_corr

        sv_lp_diff[sv_id].append(L_if - P_if)  # ambiguity
        sv_sagnac[sv_id].append(sag)
        sv_elevation[sv_id].append(el)

# Float ambiguities
float_amb = {}
for sv, diffs in sv_lp_diff.items():
    if len(diffs) >= 10:
        float_amb[sv] = float(np.median(diffs))

print("Float ambiguities (phase-code diff):")
for sv in sorted(float_amb.keys()):
    arr = np.array(sv_lp_diff[sv])
    print(f"  {sv}: B={float_amb[sv]:.4f}m std={np.std(arr):.4f}m n={len(arr)}")

# Step 2: At each epoch, estimate [clk, trop] + per-SV biases
# The per-SV bias is the DCB-like term: observed - model without bias
# We estimate it jointly with clock and tropo

print(f"\n=== Joint estimation: clk + trop + per-SV biases ===")

# First, estimate mean per-SV bias across all epochs
sv_bias_mean = {}
for sv in float_amb:
    biases = []
    for gps_sod in epochs_to_test:
        if gps_sod not in gps1b or sv not in gps1b[gps_sod]: continue
        rec = gps1b[gps_sod][sv]
        if 'L1_phase' not in rec: continue

        utc_dt = J2000 + timedelta(seconds=gps_sod)
        ref_pos = interpolate_ref(ref_orbit, gps_sod)

        pos, clk, vel = get_gps_pos_from_sp3(sp3, sv, utc_dt)
        if pos is None or abs(clk) > 0.1 * C: continue

        rho = float(np.linalg.norm(pos - ref_pos))
        for _ in range(5):
            travel_time = rho / C
            tx_dt = utc_dt - timedelta(seconds=travel_time)
            pos_tx, clk_tx, vel_tx = get_gps_pos_from_sp3(sp3, sv, tx_dt)
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
        L_if = ALPHA * L1_phase_m + BETA * L2_phase_m
        P_if = float(rec['P_if'])

        L_r = L_if + clk - rho_corr

        # Bias = L-rho - float_amb - (median of all code residuals at this epoch)
        # We'll iterate to estimate clock + biases jointly
        biases.append((gps_sod, L_r))

    if len(biases) >= 5:
        # Simple: bias = L-rho - float_amb, then median across epochs
        b_vals = [b[1] - float_amb[sv] for b in biases]
        sv_bias_mean[sv] = float(np.median(b_vals))

print("\nPer-SV mean bias (L-rho - B_float, median across epochs):")
for sv in sorted(sv_bias_mean.keys()):
    print(f"  {sv}: bias={sv_bias_mean[sv]:.3f}m")

# Step 3: Apply per-SV bias correction and check residuals
print(f"\n=== Residuals after bias correction ===")
for gps_sod in epochs_to_test[:5]:
    utc_dt = J2000 + timedelta(seconds=gps_sod)
    ref_pos = interpolate_ref(ref_orbit, gps_sod)

    code_vals = []
    phase_vals = []
    for sv_id, rec in gps1b[gps_sod].items():
        if sv_id not in float_amb or sv_id not in sv_bias_mean: continue
        if 'L1_phase' not in rec: continue

        pos, clk, vel = get_gps_pos_from_sp3(sp3, sv_id, utc_dt)
        if pos is None or abs(clk) > 0.1 * C: continue

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
        L_if = ALPHA * L1_phase_m + BETA * L2_phase_m
        P_if = float(rec['P_if'])

        L_r = L_if + clk - rho_corr
        P_r = P_if + clk - rho_corr

        # Apply float ambiguity + per-SV bias correction
        L_corrected = L_r - float_amb[sv_id] - sv_bias_mean[sv_id]
        P_corrected = P_r  # code has no ambiguity

        code_vals.append(P_corrected)
        phase_vals.append((sv_id, L_corrected, P_corrected))

    if len(code_vals) < 4: continue

    clk_est = np.median(code_vals)
    phase_res = np.array([pv[1] - clk_est for pv in phase_vals])

    print(f"epoch {gps_sod}: clk={clk_est:.3f}m  phase_res_std={np.std(phase_res):.3f}m  "
          f"phase_res_range=[{np.min(phase_res):.3f}, {np.max(phase_res):.3f}]m")

# Step 4: Full LSQ at each epoch, with per-SV bias pre-applied
# Design: [dX, dY, dZ, clk_r, trop_wet]
# Phase measurements: corrected for float_amb + per-SV bias
# Code measurements: used as-is (bias is in phase domain)

print(f"\n=== Full LSQ: pos + clk + trop (bias-corrected phase) ===")
pos_errors = []
for gps_sod in epochs_to_test:
    utc_dt = J2000 + timedelta(seconds=gps_sod)
    ref_pos = interpolate_ref(ref_orbit, gps_sod)

    sv_list = []
    A_rows, y_rows, w_rows = [], [], []

    for sv_id, rec in gps1b[gps_sod].items():
        if sv_id not in float_amb or sv_id not in sv_bias_mean: continue
        if 'L1_phase' not in rec: continue

        pos, clk, vel = get_gps_pos_from_sp3(sp3, sv_id, utc_dt)
        if pos is None or abs(clk) > 0.1 * C: continue

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
        el = math.asin(abs(pos[2] - ref_pos[2]) / rho)
        if el < 0.087: continue
        mf = 1.0 / max(math.sin(el), 0.05)

        L1_phase_m = float(rec['L1_phase'])
        L2_phase_m = float(rec['L2_phase'])
        L_if = ALPHA * L1_phase_m + BETA * L2_phase_m
        P_if = float(rec['P_if'])

        L_r = L_if + clk - rho_corr
        P_r = P_if + clk - rho_corr

        e = (pos - ref_pos) / rho
        h = np.array([-e[0], -e[1], -e[2], 1.0, mf])

        # Phase: corrected for float ambiguity AND per-SV bias
        A_rows.append(h)
        y_rows.append(L_r - float_amb[sv_id] - sv_bias_mean[sv_id])
        w_rows.append(1.0 / 0.01**2)  # optimistic: ~1cm for bias-corrected phase

        # Code
        A_rows.append(h)
        y_rows.append(P_r)
        w_rows.append(1.0 / 0.3**2)

        sv_list.append(sv_id)

    if len(A_rows) < 10: continue

    A = np.array(A_rows)
    y = np.array(y_rows)
    W = np.diag(w_rows)

    AtWA = A.T @ W @ A
    AtWy = A.T @ W @ y
    dx = np.linalg.solve(AtWA, AtWy)

    pos_est = ref_pos + dx[:3]
    clk_est = dx[3]
    trop_est = dx[4]

    model = A @ dx
    residuals = y - model
    phase_res = residuals[0::2]
    code_res = residuals[1::2]
    d3 = np.linalg.norm(pos_est - ref_pos)
    pos_errors.append(d3)

    if gps_sod in epochs_to_test[:3]:
        print(f"epoch {gps_sod}:")
        print(f"  pos={pos_est/1000} km")
        print(f"  3D error={d3:.3f}m clk={clk_est:.3f}m trop={trop_est:.3f}m")
        print(f"  phase post-fit: median={np.median(phase_res):.3f}m std={np.std(phase_res):.3f}m "
              f"range=[{np.min(phase_res):.3f}, {np.max(phase_res):.3f}]m")

if pos_errors:
    print(f"\n  3D RMS over {len(pos_errors)} epochs: {np.sqrt(np.mean(np.array(pos_errors)**2)):.3f}m")

# Step 5: Two-pass approach — estimate biases from code, then phase-only positioning
print(f"\n=== Two-pass: code-based bias estimation, phase-based positioning ===")
print(f"Pass 1: Estimate per-SV biases using code measurements only")
# Re-estimate biases using P_r only (no phase ambiguity involved)
sv_bias_v2 = {}
for sv in float_amb:
    p_vals = []
    for gps_sod in epochs_to_test:
        if gps_sod not in gps1b or sv not in gps1b[gps_sod]: continue
        rec = gps1b[gps_sod][sv]
        if 'P_if' not in rec: continue

        utc_dt = J2000 + timedelta(seconds=gps_sod)
        ref_pos = interpolate_ref(ref_orbit, gps_sod)
        pos, clk, vel = get_gps_pos_from_sp3(sp3, sv, utc_dt)
        if pos is None or abs(clk) > 0.1 * C: continue

        rho = float(np.linalg.norm(pos - ref_pos))
        for _ in range(5):
            travel_time = rho / C
            tx_dt = utc_dt - timedelta(seconds=travel_time)
            pos_tx, clk_tx, vel_tx = get_gps_pos_from_sp3(sp3, sv, tx_dt)
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
        P_r = float(rec['P_if']) + clk - rho_corr
        p_vals.append(P_r)

    if len(p_vals) >= 5:
        sv_bias_v2[sv] = float(np.median(p_vals))

print("Per-SV code bias (P_r median):")
for sv in sorted(sv_bias_v2.keys()):
    print(f"  {sv}: {sv_bias_v2[sv]:.3f}m")

print(f"\nPass 2: Phase-only positioning with bias-corrected phase")
pos_errors_v2 = []
for gps_sod in epochs_to_test:
    utc_dt = J2000 + timedelta(seconds=gps_sod)
    ref_pos = interpolate_ref(ref_orbit, gps_sod)

    A_rows, y_rows, w_rows = [], [], []
    for sv_id, rec in gps1b[gps_sod].items():
        if sv_id not in float_amb or sv_id not in sv_bias_v2: continue
        if 'L1_phase' not in rec: continue

        pos, clk, vel = get_gps_pos_from_sp3(sp3, sv_id, utc_dt)
        if pos is None or abs(clk) > 0.1 * C: continue

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
        el = math.asin(abs(pos[2] - ref_pos[2]) / rho)
        if el < 0.087: continue
        mf = 1.0 / max(math.sin(el), 0.05)

        L1_phase_m = float(rec['L1_phase'])
        L2_phase_m = float(rec['L2_phase'])
        L_if = ALPHA * L1_phase_m + BETA * L2_phase_m

        L_r = L_if + clk - rho_corr
        e = (pos - ref_pos) / rho
        h = np.array([-e[0], -e[1], -e[2], 1.0, mf])

        # Phase corrected for float amb + code-derived bias
        A_rows.append(h)
        y_rows.append(L_r - float_amb[sv_id] - sv_bias_v2[sv_id])
        w_rows.append(1.0 / 0.01**2)

    if len(A_rows) < 4: continue

    A = np.array(A_rows)
    y = np.array(y_rows)
    W = np.diag(w_rows)

    AtWA = A.T @ W @ A
    AtWy = A.T @ W @ y
    dx = np.linalg.solve(AtWA, AtWy)

    pos_est = ref_pos + dx[:3]
    d3 = np.linalg.norm(pos_est - ref_pos)
    pos_errors_v2.append(d3)

    model = A @ dx
    phase_res = y - model

if pos_errors_v2:
    print(f"  3D RMS over {len(pos_errors_v2)} epochs: {np.sqrt(np.mean(np.array(pos_errors_v2)**2)):.3f}m")
    print(f"  3D mean: {np.mean(pos_errors_v2):.3f}m")
    print(f"  3D max:  {np.max(pos_errors_v2):.3f}m")

# Step 6: Kalman Filter with per-SV bias correction
# State: [X, Y, Z, clk_r, trop_wet]
# Uses bias-corrected phase + code measurements
# Process noise: position random walk, clock white noise, tropo random walk

print(f"\n=== Kalman Filter: 5-state [pos, clk, trop] ===")

std_phase = 0.2   # post-fit phase noise level
std_code = 0.3    # code noise
dt = 30.0         # epoch interval

# Sort epochs
epochs_sorted = sorted(epochs_to_test)

# Find which SVs have bias estimates
sv_with_bias = set(sv_bias_v2.keys()) & set(float_amb.keys())

# Initial state and covariance
x = np.zeros(5)  # [dX, dY, dZ, clk_r, trop_wet]
P = np.diag([1e4**2, 1e4**2, 1e4**2, 1e4**2, 0.1**2])

# Process noise per second
Q_per_s = np.diag([1e-4, 1e-4, 1e-4, 1e8, 1e-16])  # m²/s for pos, m²/s for clk, m²/s for trop
Q_dt = Q_per_s * dt

# Transition matrix (position is static in ECI? no — but for short arcs, identity is ok)
# Actually GRACE-FO moves at ~7.6 km/s. With dt=30s, that's ~230 km displacement.
# The KF state is dX from reference, so identity transition is wrong for long arcs.
# But if we re-linearize at each epoch (IEKF style), it works.
F = np.eye(5)

kf_errors = []
prev_gps_sod = None
for gps_sod in epochs_sorted:
    utc_dt = J2000 + timedelta(seconds=gps_sod)
    ref_pos = interpolate_ref(ref_orbit, gps_sod)

    # Predict step
    if prev_gps_sod is not None:
        actual_dt = gps_sod - prev_gps_sod
        Q = Q_per_s * actual_dt
        x = F @ x  # identity, so no change
        P = F @ P @ F.T + Q
    prev_gps_sod = gps_sod

    # Build measurements at this epoch
    H_rows, z_rows, R_diag = [], [], []
    for sv_id, rec in gps1b[gps_sod].items():
        if sv_id not in sv_with_bias: continue
        if 'L1_phase' not in rec: continue

        pos, clk, vel = get_gps_pos_from_sp3(sp3, sv_id, utc_dt)
        if pos is None or abs(clk) > 0.1 * C: continue

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
        el = math.asin(abs(pos[2] - ref_pos[2]) / rho)
        if el < 0.087: continue
        mf = 1.0 / max(math.sin(el), 0.05)
        sin_el = max(math.sin(el), 0.1)

        L1_phase_m = float(rec['L1_phase'])
        L2_phase_m = float(rec['L2_phase'])
        L_if = ALPHA * L1_phase_m + BETA * L2_phase_m
        P_if = float(rec['P_if'])

        L_r = L_if + clk - rho_corr
        P_r = P_if + clk - rho_corr

        e = (pos - ref_pos) / rho
        h = np.array([-e[0], -e[1], -e[2], 1.0, mf])

        # Phase: corrected for float amb + code-derived bias
        H_rows.append(h)
        z_rows.append(L_r - float_amb[sv_id] - sv_bias_v2[sv_id])
        R_diag.append((std_phase / sin_el)**2)

        # Code
        H_rows.append(h)
        z_rows.append(P_r)
        R_diag.append((std_code / sin_el)**2)

    if len(H_rows) < 8: continue

    H_mat = np.array(H_rows)
    z = np.array(z_rows)
    R = np.diag(R_diag)

    # KF update
    S = H_mat @ P @ H_mat.T + R
    K = P @ H_mat.T @ np.linalg.inv(S)

    # Innovation
    innov = z - H_mat @ x
    dx = K @ innov

    # Joseph-form covariance update for numerical stability
    I_KH = np.eye(5) - K @ H_mat
    P = I_KH @ P @ I_KH.T + K @ R @ K.T

    x = x + dx

    pos_est = ref_pos + x[:3]
    d3 = np.linalg.norm(pos_est - ref_pos)
    kf_errors.append(d3)

if kf_errors:
    errs = np.array(kf_errors)
    print(f"  KF 3D RMS over {len(errs)} epochs: {np.sqrt(np.mean(errs**2)):.3f}m")
    print(f"  KF 3D mean: {np.mean(errs):.3f}m")
    print(f"  KF 3D median: {np.median(errs):.3f}m")
    print(f"  KF 3D max: {np.max(errs):.3f}m")
    print(f"  Fraction < 1m: {np.sum(errs < 1.0) / len(errs) * 100:.1f}%")
    print(f"  Fraction < 2m: {np.sum(errs < 2.0) / len(errs) * 100:.1f}%")
    print(f"  Fraction < 5m: {np.sum(errs < 5.0) / len(errs) * 100:.1f}%")

# Step 7: Try with realistic phase noise (matched to post-fit residuals)
print(f"\n=== KF with tuned noise: phase=0.2m, code=0.3m, Q scaled ===")

x2 = np.zeros(5)
P2 = np.diag([1e4**2, 1e4**2, 1e4**2, 1e4**2, 0.1**2])
Q_per_s2 = np.diag([1e-2, 1e-2, 1e-2, 1e6, 1e-16])

kf_errors2 = []
prev_gps_sod = None
for gps_sod in epochs_sorted:
    utc_dt = J2000 + timedelta(seconds=gps_sod)
    ref_pos = interpolate_ref(ref_orbit, gps_sod)

    if prev_gps_sod is not None:
        actual_dt = gps_sod - prev_gps_sod
        Q = Q_per_s2 * actual_dt
        P2 = P2 + Q  # F=I
    prev_gps_sod = gps_sod

    H_rows, z_rows, R_diag = [], [], []
    for sv_id, rec in gps1b[gps_sod].items():
        if sv_id not in sv_with_bias: continue
        if 'L1_phase' not in rec: continue

        pos, clk, vel = get_gps_pos_from_sp3(sp3, sv_id, utc_dt)
        if pos is None or abs(clk) > 0.1 * C: continue

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
        el = math.asin(abs(pos[2] - ref_pos[2]) / rho)
        if el < 0.087: continue
        mf = 1.0 / max(math.sin(el), 0.05)
        sin_el = max(math.sin(el), 0.1)

        L1_phase_m = float(rec['L1_phase'])
        L2_phase_m = float(rec['L2_phase'])
        L_if = ALPHA * L1_phase_m + BETA * L2_phase_m
        P_if = float(rec['P_if'])

        L_r = L_if + clk - rho_corr
        P_r = P_if + clk - rho_corr

        e = (pos - ref_pos) / rho
        h = np.array([-e[0], -e[1], -e[2], 1.0, mf])

        H_rows.append(h)
        z_rows.append(L_r - float_amb[sv_id] - sv_bias_v2[sv_id])
        R_diag.append((0.2 / sin_el)**2)

        H_rows.append(h)
        z_rows.append(P_r)
        R_diag.append((0.3 / sin_el)**2)

    if len(H_rows) < 8: continue

    H_mat = np.array(H_rows)
    z = np.array(z_rows)
    R = np.diag(R_diag)

    S = H_mat @ P2 @ H_mat.T + R
    K = P2 @ H_mat.T @ np.linalg.inv(S)
    innov = z - H_mat @ x2
    dx = K @ innov
    I_KH = np.eye(5) - K @ H_mat
    P2 = I_KH @ P2 @ I_KH.T + K @ R @ K.T
    x2 = x2 + dx

    pos_est = ref_pos + x2[:3]
    d3 = np.linalg.norm(pos_est - ref_pos)
    kf_errors2.append(d3)

if kf_errors2:
    errs2 = np.array(kf_errors2)
    print(f"  KF-v2 3D RMS: {np.sqrt(np.mean(errs2**2)):.3f}m")
    print(f"  KF-v2 3D mean: {np.mean(errs2):.3f}m")
    print(f"  KF-v2 3D median: {np.median(errs2):.3f}m")
    print(f"  Fraction < 1m: {np.sum(errs2 < 1.0) / len(errs2) * 100:.1f}%")
    print(f"  Fraction < 2m: {np.sum(errs2 < 2.0) / len(errs2) * 100:.1f}%")

# Step 8: KF with proper re-linearization — reset position state each epoch
# Since we re-linearize at the TRUE position, dX should be zero-mean each epoch.
# We must NOT propagate position errors across epochs with changing reference.
print(f"\n=== KF-v3: Reset position state per epoch (correct re-linearization) ===")

x3_clk_trop = np.zeros(2)  # [clk_r, trop_wet] — persistent states
P3 = np.diag([1e4**2, 0.1**2])  # initial uncertainty

Q_clk_per_s = 1e6     # m²/s — large but not infinite (allows some smoothing)
Q_trop_per_s = 1e-16   # m²/s — tiny random walk

kf_errors3 = []
prev_gps_sod = None
for gps_sod in epochs_sorted:
    utc_dt = J2000 + timedelta(seconds=gps_sod)
    ref_pos = interpolate_ref(ref_orbit, gps_sod)

    # Predict clock + tropo states
    if prev_gps_sod is not None:
        actual_dt = gps_sod - prev_gps_sod
        P3[0, 0] += Q_clk_per_s * actual_dt  # clock: nearly white noise
        P3[1, 1] += Q_trop_per_s * actual_dt  # tropo: slow random walk
    prev_gps_sod = gps_sod

    # Position state: [dX, dY, dZ] — reset each epoch (re-linearized at truth)
    x_pos = np.zeros(3)
    P_pos = np.diag([1e4**2, 1e4**2, 1e4**2])  # reset position covariance each epoch

    H_rows, z_rows, R_diag = [], [], []
    for sv_id, rec in gps1b[gps_sod].items():
        if sv_id not in sv_with_bias: continue
        if 'L1_phase' not in rec: continue

        pos, clk, vel = get_gps_pos_from_sp3(sp3, sv_id, utc_dt)
        if pos is None or abs(clk) > 0.1 * C: continue

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
        el = math.asin(abs(pos[2] - ref_pos[2]) / rho)
        if el < 0.087: continue
        mf = 1.0 / max(math.sin(el), 0.05)
        sin_el = max(math.sin(el), 0.1)

        L1_phase_m = float(rec['L1_phase'])
        L2_phase_m = float(rec['L2_phase'])
        L_if = ALPHA * L1_phase_m + BETA * L2_phase_m
        P_if = float(rec['P_if'])

        L_r = L_if + clk - rho_corr
        P_r = P_if + clk - rho_corr

        e = (pos - ref_pos) / rho
        # H = [dX, dY, dZ | clk_r, trop_wet]
        h_full = np.array([-e[0], -e[1], -e[2], 1.0, mf])

        H_rows.append(h_full)
        z_rows.append(L_r - float_amb[sv_id] - sv_bias_v2[sv_id])
        R_diag.append((0.2 / sin_el)**2)

        H_rows.append(h_full)
        z_rows.append(P_r)
        R_diag.append((0.3 / sin_el)**2)

    if len(H_rows) < 8: continue

    H_mat = np.array(H_rows)
    z = np.array(z_rows)
    R = np.diag(R_diag)

    # Full state: [pos(3), clk(1), trop(1)]
    x_full = np.concatenate([x_pos, x3_clk_trop])
    P_full = np.zeros((5, 5))
    P_full[:3, :3] = P_pos
    P_full[3:, 3:] = P3

    S = H_mat @ P_full @ H_mat.T + R
    K = P_full @ H_mat.T @ np.linalg.inv(S)
    innov = z - H_mat @ x_full
    dx = K @ innov

    # Joseph-form
    I_KH = np.eye(5) - K @ H_mat
    P_full = I_KH @ P_full @ I_KH.T + K @ R @ K.T

    x_full = x_full + dx

    # Extract states
    pos_est = ref_pos + x_full[:3]
    x3_clk_trop = x_full[3:]
    P3 = P_full[3:, 3:]

    d3 = np.linalg.norm(pos_est - ref_pos)
    kf_errors3.append(d3)

if kf_errors3:
    errs3 = np.array(kf_errors3)
    print(f"  KF-v3 3D RMS: {np.sqrt(np.mean(errs3**2)):.3f}m")
    print(f"  KF-v3 3D mean: {np.mean(errs3):.3f}m")
    print(f"  KF-v3 3D median: {np.median(errs3):.3f}m")
    print(f"  KF-v3 3D max: {np.max(errs3):.3f}m")
    print(f"  Fraction < 1m: {np.sum(errs3 < 1.0) / len(errs3) * 100:.1f}%")
    print(f"  Fraction < 2m: {np.sum(errs3 < 2.0) / len(errs3) * 100:.1f}%")

# Step 9: KF with matched LSQ weights (phase σ=0.01m, code σ=0.3m)
print(f"\n=== KF-v4: Matched LSQ weights (phase sigma=0.01m) ===")
std_p9, std_c9 = 0.01, 0.3

x9 = np.zeros(2)
P9 = np.diag([1e4**2, 0.1**2])
Q_clk9, Q_trop9 = 1e6, 1e-16

kf_errors9 = []
prev_gps_sod = None
for gps_sod in epochs_sorted:
    utc_dt = J2000 + timedelta(seconds=gps_sod)
    ref_pos = interpolate_ref(ref_orbit, gps_sod)

    if prev_gps_sod is not None:
        actual_dt = gps_sod - prev_gps_sod
        P9[0, 0] += Q_clk9 * actual_dt
        P9[1, 1] += Q_trop9 * actual_dt
    prev_gps_sod = gps_sod

    x_pos = np.zeros(3)
    P_pos = np.diag([1e4**2, 1e4**2, 1e4**2])

    H_rows, z_rows, R_diag = [], [], []
    for sv_id, rec in gps1b[gps_sod].items():
        if sv_id not in sv_with_bias: continue
        if 'L1_phase' not in rec: continue

        pos, clk, vel = get_gps_pos_from_sp3(sp3, sv_id, utc_dt)
        if pos is None or abs(clk) > 0.1 * C: continue

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
        el = math.asin(abs(pos[2] - ref_pos[2]) / rho)
        if el < 0.087: continue
        mf = 1.0 / max(math.sin(el), 0.05)
        sin_el = max(math.sin(el), 0.1)

        L1_phase_m = float(rec['L1_phase'])
        L2_phase_m = float(rec['L2_phase'])
        L_if = ALPHA * L1_phase_m + BETA * L2_phase_m
        P_if = float(rec['P_if'])

        L_r = L_if + clk - rho_corr
        P_r = P_if + clk - rho_corr

        e = (pos - ref_pos) / rho
        h_full = np.array([-e[0], -e[1], -e[2], 1.0, mf])

        H_rows.append(h_full)
        z_rows.append(L_r - float_amb[sv_id] - sv_bias_v2[sv_id])
        R_diag.append((std_p9 / sin_el)**2)

        H_rows.append(h_full)
        z_rows.append(P_r)
        R_diag.append((std_c9 / sin_el)**2)

    if len(H_rows) < 8: continue

    H_mat = np.array(H_rows)
    z = np.array(z_rows)
    R = np.diag(R_diag)

    x_full = np.concatenate([x_pos, x9])
    P_full = np.zeros((5, 5))
    P_full[:3, :3] = P_pos
    P_full[3:, 3:] = P9

    S = H_mat @ P_full @ H_mat.T + R
    K = P_full @ H_mat.T @ np.linalg.inv(S)
    innov = z - H_mat @ x_full
    dx = K @ innov
    I_KH = np.eye(5) - K @ H_mat
    P_full = I_KH @ P_full @ I_KH.T + K @ R @ K.T
    x_full = x_full + dx

    pos_est = ref_pos + x_full[:3]
    x9 = x_full[3:]
    P9 = P_full[3:, 3:]

    d3 = np.linalg.norm(pos_est - ref_pos)
    kf_errors9.append(d3)

if kf_errors9:
    errs9 = np.array(kf_errors9)
    print(f"  KF-v4 3D RMS: {np.sqrt(np.mean(errs9**2)):.3f}m")
    print(f"  KF-v4 3D mean: {np.mean(errs9):.3f}m")
    print(f"  KF-v4 3D median: {np.median(errs9):.3f}m")
    print(f"  Fraction < 1m: {np.sum(errs9 < 1.0) / len(errs9) * 100:.1f}%")
    print(f"  Fraction < 2m: {np.sum(errs9 < 2.0) / len(errs9) * 100:.1f}%")

# Step 10: LSQ with realistic phase noise (σ=0.2m) to see true performance
print(f"\n=== LSQ with realistic phase noise (σ=0.2m) ===")
pos_errors_realistic = []
for gps_sod in epochs_to_test:
    utc_dt = J2000 + timedelta(seconds=gps_sod)
    ref_pos = interpolate_ref(ref_orbit, gps_sod)

    A_rows, y_rows, w_rows = [], [], []
    for sv_id, rec in gps1b[gps_sod].items():
        if sv_id not in sv_with_bias: continue
        if 'L1_phase' not in rec: continue

        pos, clk, vel = get_gps_pos_from_sp3(sp3, sv_id, utc_dt)
        if pos is None or abs(clk) > 0.1 * C: continue

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
        el = math.asin(abs(pos[2] - ref_pos[2]) / rho)
        if el < 0.087: continue
        mf = 1.0 / max(math.sin(el), 0.05)

        L1_phase_m = float(rec['L1_phase'])
        L2_phase_m = float(rec['L2_phase'])
        L_if = ALPHA * L1_phase_m + BETA * L2_phase_m
        P_if = float(rec['P_if'])

        L_r = L_if + clk - rho_corr
        P_r = P_if + clk - rho_corr

        e = (pos - ref_pos) / rho
        h = np.array([-e[0], -e[1], -e[2], 1.0, mf])

        A_rows.append(h)
        y_rows.append(L_r - float_amb[sv_id] - sv_bias_v2[sv_id])
        w_rows.append(1.0 / 0.2**2)  # realistic phase noise

        A_rows.append(h)
        y_rows.append(P_r)
        w_rows.append(1.0 / 0.3**2)  # realistic code noise

    if len(A_rows) < 8: continue

    A = np.array(A_rows)
    y = np.array(y_rows)
    W = np.diag(w_rows)
    AtWA = A.T @ W @ A
    AtWy = A.T @ W @ y
    dx = np.linalg.solve(AtWA, AtWy)
    pos_est = ref_pos + dx[:3]
    d3 = np.linalg.norm(pos_est - ref_pos)
    pos_errors_realistic.append(d3)

if pos_errors_realistic:
    errs_r = np.array(pos_errors_realistic)
    print(f"  LSQ-realistic 3D RMS: {np.sqrt(np.mean(errs_r**2)):.3f}m")
    print(f"  LSQ-realistic 3D mean: {np.mean(errs_r):.3f}m")
