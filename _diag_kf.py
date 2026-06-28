#!/usr/bin/env python3
"""Diagnose the KF measurements after bias correction"""
import sys; sys.path.insert(0, '.')
from src.gps1a_loader import download_gps1a, gps_sod_to_utc
from src.ambiguity import compute_if_m
from src.sp3_loader import get_gps_pos_from_sp3
import pickle, numpy as np
from datetime import datetime, timedelta
from collections import defaultdict

C = 299792458.0
OMEGA_E = 7.2921151467e-5
J2000 = datetime(2000, 1, 1, 12, 0, 0)

gps_obs = download_gps1a(2024, 4, 29, grace_filter='C')
sp3 = pickle.load(open('data/2024/120/igs_sp3_FIN.pkl', 'rb'))

# Load reference orbit
ref_orbit = {}
for line in open('data/gracefo/2024/2024-04-29/GNV1B_2024-04-29_C_04.txt'):
    p = line.split()
    if len(p) < 6: continue
    try:
        t = float(p[0]); flag = p[2]
        if flag in ('C', 'E'):
            ref_orbit[gps_sod_to_utc(t)] = np.array([float(p[3]), float(p[4]), float(p[5])])
    except: pass
orbit_ts = sorted(ref_orbit.keys())

def get_rcv_pos(utc_dt):
    t0, t1 = None, None
    for j, ti in enumerate(orbit_ts):
        if ti >= utc_dt:
            t1 = ti; t0 = orbit_ts[j-1] if j>0 else ti; break
        t0 = ti
    if t0 is None: t0 = orbit_ts[0]
    if t1 is None: t1 = orbit_ts[-1]
    dt_frac = (utc_dt - t0).total_seconds()
    dt_tot = (t1 - t0).total_seconds()
    if dt_tot < 0.01:
        return np.array(ref_orbit[t0], dtype=float)
    a = dt_frac / dt_tot
    return np.array(ref_orbit[t0], dtype=float)*(1-a) + np.array(ref_orbit[t1], dtype=float)*a

def get_sat_geometry(sv, utc_dt, rcv_pos):
    from src.sp3_loader import get_gps_pos_from_sp3
    pos, clk, vel = get_gps_pos_from_sp3(sp3, sv, utc_dt)
    if pos is None or abs(clk) > 0.1 * C:
        return None, None, None
    rho = float(np.linalg.norm(pos - rcv_pos))
    for _ in range(5):
        travel_time = rho / C
        tx_dt = utc_dt - timedelta(seconds=travel_time)
        pos_tx, clk_tx, vel_tx = get_gps_pos_from_sp3(sp3, sv, tx_dt)
        if pos_tx is None: break
        rho_new = float(np.linalg.norm(pos_tx - rcv_pos))
        if abs(rho_new - rho) < 1e-8: pos, clk = pos_tx, clk_tx; rho = rho_new; break
        pos, clk = pos_tx, clk_tx; rho = rho_new
    if not (1.8e7 < rho < 2.8e7): return None, None, None
    sag = (OMEGA_E / C) * (pos[0] * rcv_pos[1] - pos[1] * rcv_pos[0])
    return pos, clk, rho + sag

# Build records (simplified, first 100 epochs)
gps_sods = sorted(gps_obs.keys())
t_start = min(ref_orbit.keys())

records = []
for gps_sod in gps_sods[:600]:  # first ~100 min of 10s data
    utc_dt = gps_sod_to_utc(gps_sod)
    if utc_dt < t_start: continue
    rcv_pos = get_rcv_pos(utc_dt)

    epoch_recs = []
    for sv, rec in gps_obs[gps_sod].items():
        sat_pos, clk_sv, rho_corr = get_sat_geometry(sv, utc_dt, rcv_pos)
        if sat_pos is None: continue
        delta = sat_pos - rcv_pos
        rng = float(np.linalg.norm(delta))
        if not (2e7 < rng < 5e7): continue

        L_if, P_if, B_if = compute_if_m(rec['L1_cyc'], rec['L2_cyc'], rec['P1'], rec['P2'])
        L_r = L_if + clk_sv - rho_corr
        P_r = P_if + clk_sv - rho_corr
        epoch_recs.append({'sv': sv, 'P_r': P_r, 'L_r': L_r, 'B_if': B_if,
                          'el': 45.0, 'gps_sod': gps_sod})

    if epoch_recs:
        records.append({'utc': utc_dt, 'rcv_pos': rcv_pos, 'sv_data': epoch_recs})

print(f"Built {len(records)} records")

# Estimate per-SV biases from first 30 records
sv_p_vals = defaultdict(list)
for ep in records[:30]:
    for d in ep['sv_data']:
        sv_p_vals[d['sv']].append(d['P_r'])

sv_bias = {}
for sv, p_vals in sv_p_vals.items():
    if len(p_vals) >= 5:
        sv_bias[sv] = float(np.median(p_vals))

print(f"SV biases: {len(sv_bias)} SVs")

# Check bias-corrected P_r at a few epochs
for idx in [0, 15, 30, 60]:
    if idx >= len(records): break
    ep = records[idx]
    print(f"\nEpoch {idx}: utc={ep['utc']}")
    corrected = []
    for d in ep['sv_data']:
        sv = d['sv']
        if sv not in sv_bias: continue
        P_c = d['P_r'] - sv_bias[sv]
        corrected.append((sv, P_c))

    if corrected:
        vals = [c[1] for c in corrected]
        print(f"  P_r-bias: mean={np.mean(vals):.1f}m std={np.std(vals):.1f}m "
              f"range=[{np.min(vals):.1f},{np.max(vals):.1f}]m")
        for sv, Pc in corrected[:6]:
            print(f"    {sv}: P_r={d['P_r']:.1f} bias={sv_bias[sv]:.1f} corrected={Pc:.1f}m")

# Now check: what's the typical position error caused by a measurement bias?
# If one SV has a -4000m residual and another has +4000m, what position error does that cause?
print("\n\n=== Geometry analysis ===")
ep = records[0]
rcv_pos = ep['rcv_pos']
print(f"rcv_pos = {rcv_pos/1000} km")

# Simple: solve for position using only code measurements with bias correction
# H rows: [-los_x, -los_y, -los_z, 1, mf], y = P_r - bias
H_rows, y_rows = [], []
for d in ep['sv_data']:
    sv = d['sv']
    if sv not in sv_bias: continue
    sat_pos, clk_sv, rho_corr = get_sat_geometry(sv, ep['utc'], rcv_pos)
    if sat_pos is None: continue
    los = sat_pos - rcv_pos
    los = los / np.linalg.norm(los)
    mf = 1.0 / max(np.sin(np.radians(45.0)), 0.1)
    P_c = d['P_r'] - sv_bias[sv]
    H_rows.append([-los[0], -los[1], -los[2], 1.0, mf])
    y_rows.append(P_c)

H = np.array(H_rows)
y = np.array(y_rows)
print(f"H shape: {H.shape}, y: {y}")
print(f"y mean={np.mean(y):.1f}m std={np.std(y):.1f}m")

# Solve
try:
    P0 = np.eye(5) * 1e6
    W = np.eye(len(y_rows)) / 0.3**2
    S = H @ P0 @ H.T + np.linalg.inv(W)
    K = P0 @ H.T @ np.linalg.inv(S)
    dx = K @ y
    print(f"dx = {dx} (dX={dx[0]:.1f} dY={dx[1]:.1f} dZ={dx[2]:.1f} clk={dx[3]:.1f} trop={dx[4]:.4f})")

    # Compute position error
    solved = rcv_pos + dx[:3]
    # No reference for this epoch, just show the correction
    print(f"Position correction norm: {np.linalg.norm(dx[:3]):.1f}m")

    # Post-fit residuals
    v = y - H @ dx
    print(f"Post-fit residuals: mean={np.mean(v):.1f}m std={np.std(v):.1f}m")

    # Is the solution geometry good?
    # Check DOP
    Q = np.linalg.inv(H.T @ W @ H)
    gdop = np.sqrt(np.trace(Q[:3,:3]) + Q[4,4])
    pdop = np.sqrt(np.trace(Q[:3,:3]))
    print(f"PDOP={pdop:.2f} GDOP={gdop:.2f}")
except Exception as e:
    print(f"LSQ failed: {e}")
