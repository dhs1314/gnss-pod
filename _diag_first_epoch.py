#!/usr/bin/env python3
"""Compare RINEX vs ASCII pipeline at first epoch in detail"""
import sys; sys.path.insert(0, 'src')
import pickle, math
from gps1b_rnx_loader import load_gps1b_rnx
from sp3_loader import get_gps_pos_from_sp3 as _sp3_get
import numpy as np
from datetime import datetime, timedelta

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
F1_SQ, F2_SQ = F1*F1, F2*F2
ALPHA = F1_SQ / (F1_SQ - F2_SQ)
BETA  = -F2_SQ / (F1_SQ - F2_SQ)
OMEGA_E = 7.2921151467e-5
J2000 = datetime(2000, 1, 1, 12, 0, 0)

# Load data
rnx_data = load_gps1b_rnx('data/GPS1B_2024-04-29_C_04.rnx')
ascii_data = pickle.load(open('data/gracefo/2024/2024-04-29/GPS1B_2024-04-29_C_04.pkl', 'rb'))
sp3 = pickle.load(open('data/2024/120/igs_sp3_FIN.pkl', 'rb'))

# Load GNV1B
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

def get_sat_geometry(sp3, sv, utc_dt, rcv_pos):
    pos, clk, vel = _sp3_get(sp3, sv, utc_dt)
    if pos is None or abs(clk) > 0.1 * C:
        return None, None, None
    rho = float(np.linalg.norm(pos - rcv_pos))
    for _ in range(5):
        travel_time = rho / C
        tx_dt = utc_dt - timedelta(seconds=travel_time)
        pos_tx, clk_tx, vel_tx = _sp3_get(sp3, sv, tx_dt)
        if pos_tx is None: break
        rho_new = float(np.linalg.norm(pos_tx - rcv_pos))
        if abs(rho_new - rho) < 1e-8:
            pos, clk = pos_tx, clk_tx; rho = rho_new; break
        pos, clk = pos_tx, clk_tx; rho = rho_new
    if not (1.8e7 < rho < 2.8e7):
        return None, None, None
    sag = (OMEGA_E / C) * (pos[0] * rcv_pos[1] - pos[1] * rcv_pos[0])
    return pos, clk, rho + sag

# First epoch
sod = 767620800
utc_dt = J2000 + timedelta(seconds=sod)
ref_pos = get_ref_pos(sod)

print(f"Epoch: {sod}, UTC={utc_dt}")
print(f"Ref pos (ECEF km): [{ref_pos[0]/1000:.3f}, {ref_pos[1]/1000:.3f}, {ref_pos[2]/1000:.3f}]")
print(f"r = {np.linalg.norm(ref_pos)/1000:.3f} km")

common_svs = sorted(set(rnx_data[sod].keys()) & set(ascii_data[sod].keys()))

N_BIAS = 60
# Compute sv_bias for both (using first 60 epochs)
sods_for_bias = sorted(rnx_data.keys())[:N_BIAS]

rnx_sv_bias = {}
ascii_sv_bias = {}
for sv in common_svs:
    rnx_Prs = []; asc_Prs = []
    for s in sods_for_bias:
        if s not in rnx_data or sv not in rnx_data[s]: continue
        ref = get_ref_pos(s)
        # RINEX
        r = rnx_data[s][sv]
        pos, clk, rho_corr = get_sat_geometry(sp3, sv, J2000 + timedelta(seconds=s), ref)
        if pos is not None:
            rnx_Prs.append(r['P_if'] + clk - rho_corr)
        # ASCII
        if s in ascii_data and sv in ascii_data[s]:
            a = ascii_data[s][sv]
            # Recompute using same geometry
            asc_Prs.append(a['P_if'] + clk - rho_corr)
    if len(rnx_Prs) >= 10:
        rnx_sv_bias[sv] = float(np.median(rnx_Prs))
    if len(asc_Prs) >= 10:
        ascii_sv_bias[sv] = float(np.median(asc_Prs))

print(f"\n{'SV':<6s} {'P_r(RNX)':>10s} {'sv_bias':>10s} {'L_r-sv_bias':>14s} {'P_r(ASC)':>10s} {'sv_b_asc':>10s} {'L_r-sv_b_asc':>14s}")
print("-" * 76)

for sv in common_svs:
    r = rnx_data[sod][sv]
    a = ascii_data[sod][sv]
    pos, clk, rho_corr = get_sat_geometry(sp3, sv, utc_dt, ref_pos)
    if pos is None: continue

    P_r_rnx = r['P_if'] + clk - rho_corr
    P_r_asc = a['P_if'] + clk - rho_corr
    L_r_rnx = r['L_if'] + clk - rho_corr
    L_r_asc = a['L_if'] + clk - rho_corr

    rb = rnx_sv_bias.get(sv, float('nan'))
    ab = ascii_sv_bias.get(sv, float('nan'))

    print(f"{sv:<6s} {P_r_rnx:>10.1f} {rb:>10.1f} {L_r_rnx - rb:>14.1f} "
          f"{P_r_asc:>10.1f} {ab:>10.1f} {L_r_asc - ab:>14.1f}")

# Now run the full KF for the first epoch manually
print(f"\n--- Manual KF first epoch ---")

# Initialize B states
b_rnx = {}; b_asc = {}
for sv in common_svs:
    r = rnx_data[sod][sv]; a = ascii_data[sod][sv]
    b_rnx[sv] = r['L_if'] - r['P_if']  # B_if ≈ 0
    b_asc[sv] = a['L_if'] - a['P_if']  # B_if ≈ -0.7

    if sv not in rnx_sv_bias or sv not in ascii_sv_bias:
        continue

# Build and solve for RINEX
sv_list = sorted([sv for sv in common_svs if sv in rnx_sv_bias])
N_B = len(sv_list)
N_STATE = 5 + N_B

x_prior = np.zeros(N_STATE)
P_prior = np.eye(N_STATE) * 1e6
for i, sv in enumerate(sv_list):
    x_prior[5 + i] = b_rnx[sv]
    P_prior[5 + i, 5 + i] = 0.5
P_prior[4, 4] = 1.0  # trop

H_rows, y_rows, w_rows = [], [], []
for i, d_sv in enumerate(sv_list):
    r = rnx_data[sod][d_sv]
    pos, clk, rho_corr = get_sat_geometry(sp3, d_sv, utc_dt, ref_pos)
    if pos is None: continue

    los = (pos - ref_pos) / rho_corr
    el = math.asin(abs(pos[2] - ref_pos[2]) / rho_corr)
    mf = 1.0 / max(math.sin(el), 0.1)

    L_r = r['L_if'] + clk - rho_corr
    P_r = r['P_if'] + clk - rho_corr

    # Phase
    h_p = np.zeros(N_STATE)
    h_p[0:3] = -los; h_p[3] = 1.0; h_p[4] = mf; h_p[5 + i] = 1.0
    H_rows.append(h_p); y_rows.append(L_r - rnx_sv_bias[d_sv]); w_rows.append(1.0/0.01**2)

    # Code
    h_c = np.zeros(N_STATE)
    h_c[0:3] = -los; h_c[3] = 1.0; h_c[4] = mf
    H_rows.append(h_c); y_rows.append(P_r - rnx_sv_bias[d_sv]); w_rows.append(1.0/0.30**2)

H = np.array(H_rows); y = np.array(y_rows); W = np.diag(w_rows)

try:
    S = H @ P_prior @ H.T + np.linalg.inv(W)
    K = P_prior @ H.T @ np.linalg.inv(S)
    innov = y - H @ x_prior
    print(f"RINEX KF: {N_B} SVs: {sv_list}")
    print(f"  x_prior B states: {[f'{x_prior[5+i]:.1f}' for i in range(N_B)]}")
    print(f"  y values: {[f'{y[j]:.1f}' for j in range(len(y))]}")
    print(f"  innovation: mean={np.mean(innov):.1f}m std={np.std(innov):.1f}m")
    dx = K @ innov
    print(f"  dx pos: [{dx[0]:.1f}, {dx[1]:.1f}, {dx[2]:.1f}]m clk={dx[3]:.1f}m")
    print(f"  3D error: {np.linalg.norm(dx[:3]):.1f}m")
except np.linalg.LinAlgError as e:
    print(f"RINEX KF LIN ALG ERROR: {e}")

# Build and solve for ASCII
x_prior_a = np.zeros(N_STATE)
P_prior_a = np.eye(N_STATE) * 1e6
for i, sv in enumerate(sv_list):
    x_prior_a[5 + i] = b_asc[sv]
    P_prior_a[5 + i, 5 + i] = 0.5
P_prior_a[4, 4] = 1.0

H_rows_a, y_rows_a, w_rows_a = [], [], []
for i, d_sv in enumerate(sv_list):
    a = ascii_data[sod][d_sv]
    pos, clk, rho_corr = get_sat_geometry(sp3, d_sv, utc_dt, ref_pos)
    if pos is None: continue

    los = (pos - ref_pos) / rho_corr
    el = math.asin(abs(pos[2] - ref_pos[2]) / rho_corr)
    mf = 1.0 / max(math.sin(el), 0.1)

    L_r = a['L_if'] + clk - rho_corr
    P_r = a['P_if'] + clk - rho_corr

    h_p = np.zeros(N_STATE)
    h_p[0:3] = -los; h_p[3] = 1.0; h_p[4] = mf; h_p[5 + i] = 1.0
    H_rows_a.append(h_p); y_rows_a.append(L_r - ascii_sv_bias[d_sv]); w_rows_a.append(1.0/0.01**2)

    h_c = np.zeros(N_STATE)
    h_c[0:3] = -los; h_c[3] = 1.0; h_c[4] = mf
    H_rows_a.append(h_c); y_rows_a.append(P_r - ascii_sv_bias[d_sv]); w_rows_a.append(1.0/0.30**2)

H_a = np.array(H_rows_a); y_a = np.array(y_rows_a); W_a = np.diag(w_rows_a)

try:
    S_a = H_a @ P_prior_a @ H_a.T + np.linalg.inv(W_a)
    K_a = P_prior_a @ H_a.T @ np.linalg.inv(S_a)
    innov_a = y_a - H_a @ x_prior_a
    print(f"\nASCII KF: {N_B} SVs: {sv_list}")
    print(f"  x_prior B states: {[f'{x_prior_a[5+i]:.1f}' for i in range(N_B)]}")
    print(f"  y values: {[f'{y_a[j]:.1f}' for j in range(len(y_a))]}")
    print(f"  innovation: mean={np.mean(innov_a):.1f}m std={np.std(innov_a):.1f}m")
    dx_a = K_a @ innov_a
    print(f"  dx pos: [{dx_a[0]:.1f}, {dx_a[1]:.1f}, {dx_a[2]:.1f}]m clk={dx_a[3]:.1f}m")
    print(f"  3D error: {np.linalg.norm(dx_a[:3]):.1f}m")
except np.linalg.LinAlgError as e:
    print(f"ASCII KF LIN ALG ERROR: {e}")

# Compare y vectors directly
print(f"\n--- Direct y comparison ---")
print(f"{'idx':<5s} {'y_phase_vals':>30s} {'y_code_vals':>30s}")
for i, sv in enumerate(sv_list):
    r = rnx_data[sod][sv]; a = ascii_data[sod][sv]
    pos, clk, rho_corr = get_sat_geometry(sp3, sv, utc_dt, ref_pos)
    L_r_rnx = r['L_if'] + clk - rho_corr
    P_r_rnx = r['P_if'] + clk - rho_corr
    L_r_asc = a['L_if'] + clk - rho_corr
    P_r_asc = a['P_if'] + clk - rho_corr
    yp_rnx = L_r_rnx - rnx_sv_bias[sv]
    yp_asc = L_r_asc - ascii_sv_bias[sv]
    yc_rnx = P_r_rnx - rnx_sv_bias[sv]
    yc_asc = P_r_asc - ascii_sv_bias[sv]
    print(f"{i:<5d} phase: R={yp_rnx:>8.1f} A={yp_asc:>8.1f}   code: R={yc_rnx:>8.1f} A={yc_asc:>8.1f}")
