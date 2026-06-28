#!/usr/bin/env python3
"""Critical test: fix position at GNV1B, estimate clock+tropo only, check phase residuals"""
import pickle, numpy as np, math
from datetime import datetime, timedelta
from collections import defaultdict
import sys
sys.path.insert(0, '.')
from src.sp3_loader import get_gps_pos_from_sp3

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
LAM_W = C / (F1 - F2)
F_IF = (F1*F1 - F2*F2) / (F1 + F2)
LAM_IF = C / F_IF
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

# First pass: estimate float ambiguities (phase-code difference, median across epochs)
gps_sod_start = min(gps1b.keys())
gps_sod_end = gps_sod_start + 4 * 3600
interval = 30.0

sv_diff = defaultdict(list)
epochs_to_test = []

for gps_sod in sorted(gps1b.keys()):
    if not (gps_sod_start <= gps_sod <= gps_sod_end): continue
    dt_from_start = gps_sod - gps_sod_start
    nearest = round(dt_from_start / interval) * interval
    if abs(dt_from_start - nearest) > 2.0: continue
    if len(epochs_to_test) < 30:
        epochs_to_test.append(gps_sod)

    utc_dt = J2000 + timedelta(seconds=gps_sod)
    ref_pos = interpolate_ref(ref_orbit, gps_sod)

    for sv_id, rec in gps1b[gps_sod].items():
        pos, clk, vel = get_gps_pos_from_sp3(sp3, sv_id, utc_dt)
        if pos is None: continue
        if abs(clk) > 0.1 * C: continue
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

        L_if = float(rec['L_if'])
        P_if = float(rec['P_if'])
        sv_diff[sv_id].append(L_if - P_if)

# Compute float ambiguities
float_amb = {}
for sv, diffs in sv_diff.items():
    if len(diffs) >= 10:
        arr = np.array(diffs)
        float_amb[sv] = float(np.median(arr))
        print(f"{sv}: B_if={float_amb[sv]:.3f}m n={len(diffs)} std={np.std(arr):.3f}m")

# Second pass: fix position at GNV1B, estimate only clock+tropo
print(f"\n── Test: position fixed at GNV1B, estimate clock+tropo only ──")
for gps_sod in epochs_to_test[:5]:
    utc_dt = J2000 + timedelta(seconds=gps_sod)
    ref_pos = interpolate_ref(ref_orbit, gps_sod)

    H_rows, y_rows, w_rows = [], [], []
    sv_info = []
    for sv_id, rec in gps1b[gps_sod].items():
        if sv_id not in float_amb: continue
        pos, clk, vel = get_gps_pos_from_sp3(sp3, sv_id, utc_dt)
        if pos is None: continue
        if abs(clk) > 0.1 * C: continue
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

        L_if = float(rec['L_if'])
        P_if = float(rec['P_if'])
        el = math.asin(abs(pos[2] - ref_pos[2]) / rho)
        if el < 0.087: continue
        mf = 1.0 / max(math.sin(el), 0.05)

        # Phase measurement (ambiguity-corrected)
        L_corr = L_if + clk - float_amb[sv_id]
        # Code measurement
        P_corr = P_if + clk

        # Observation: [clk, trop]
        h = np.array([1.0, mf])

        # Phase
        w_l = 1.0 / (0.003**2 / (math.sin(el)**2 + 0.01))
        H_rows.append(h); y_rows.append(L_corr - rho_corr); w_rows.append(w_l)

        # Code
        w_p = 1.0 / (0.3**2 / (math.sin(el)**2 + 0.01))
        H_rows.append(h); y_rows.append(P_corr - rho_corr); w_rows.append(w_p)

        sv_info.append((sv_id, L_corr - rho_corr, P_corr - rho_corr))

    if len(H_rows) < 4: continue

    H = np.array(H_rows); y = np.array(y_rows); w_diag = np.array(w_rows)
    HtWH = H.T @ (w_diag[:, None] * H)
    HtWy = H.T @ (w_diag * y)
    dx = np.linalg.solve(HtWH, HtWy)
    clk_est, trop_est = dx[0], dx[1]

    # Post-fit residuals
    model = H @ dx
    residuals = y - model
    phase_res = residuals[0::2]  # every other residual is phase
    code_res = residuals[1::2]

    print(f"\nepoch {gps_sod}: clk={clk_est:.3f}m trop={trop_est:.3f}m")
    print(f"  Phase post-fit: median={np.median(phase_res):.3f}m std={np.std(phase_res):.3f}m")
    print(f"  Code post-fit:  median={np.median(code_res):.3f}m std={np.std(code_res):.3f}m")

    # Check: what if we use the SIMPLE median clock (as in pass 1 debug)?
    code_prefits = [si[2] for si in sv_info]
    simple_clk = np.median(code_prefits)
    simple_phase_res = [si[1] - simple_clk for si in sv_info]
    print(f"  Simple clk (median code): {simple_clk:.3f}m")
    print(f"  Simple phase res: median={np.median(simple_phase_res):.3f}m std={np.std(simple_phase_res):.3f}m")
    for i, si in enumerate(sv_info):
        print(f"    {si[0]}: L-rho={si[1]:.3f}m P-rho={si[2]:.3f}m simple_phase_res={simple_phase_res[i]:.3f}m")
