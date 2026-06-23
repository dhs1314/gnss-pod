#!/usr/bin/env python3
"""Debug the float-ambiguity KF measurement model"""
import sys; sys.path.insert(0, '.')
from src.gps1a_loader import download_gps1a, gps_sod_to_utc
from src.ambiguity import compute_if_m
from src.sp3_loader import get_gps_pos_from_sp3
import pickle, numpy as np
from datetime import datetime, timedelta
from collections import defaultdict
from pathlib import Path

C = 299792458.0
OMEGA_E = 7.2921151467e-5
J2000 = datetime(2000, 1, 1, 12, 0, 0)

# Use the actual pipeline to build records
import run_gps1a_ppp as ppp
gps_obs = download_gps1a(2024, 4, 29, grace_filter='C')
sp3 = pickle.load(open('data/2024/120/igs_sp3_FIN.pkl', 'rb'))

gnv_path = Path("data/gracefo/2024/2024-04-29/GNV1B_2024-04-29_C_04.txt")
ref_orbit = ppp.load_gnv1b(str(gnv_path), 4)
ref_ts = sorted(ref_orbit.keys())
t_start = ref_ts[0]

records = ppp.build_gps1a_records(gps_obs, ref_orbit, t_start, 1, 30, sp3)
print(f"Records: {len(records)}")

# Now run the KF manually, printing diagnostics
b_states = {}
b_vars = {}
P_trop = 100.0
last_epoch = None

for ep_idx in range(min(10, len(records))):
    ep = records[ep_idx]
    utc_dt = ep['utc']
    rcv_pos = ep['rcv_pos']

    # Initialize B states for new SVs
    for d in ep['sv_data']:
        sv = d['sv']
        if sv not in b_states:
            b_states[sv] = d['L_r'] - d['P_r']
            b_vars[sv] = 0.5**2
            print(f"  Init {sv}: B_if = {b_states[sv]:.1f}m = {b_states[sv]/1000:.1f}km")

    active_svs = [d['sv'] for d in ep['sv_data']]
    N_B = len(active_svs)
    N_STATE = 5 + N_B

    x_prior = np.zeros(N_STATE)
    P_prior = np.eye(N_STATE) * 1e6
    for i, sv in enumerate(active_svs):
        x_prior[5 + i] = b_states[sv]
        P_prior[5 + i, 5 + i] = b_vars[sv]
    P_prior[4, 4] = P_trop

    H_rows, y_rows, w_rows = [], [], []
    for i, d in enumerate(ep['sv_data']):
        sv = d['sv']
        los = d['sat_pos'] - rcv_pos
        los = los / np.linalg.norm(los)
        mf = 1.0 / max(np.sin(np.radians(d['el'])), 0.1)

        # Phase
        row_p = np.zeros(N_STATE)
        row_p[0:3] = -los
        row_p[3] = 1.0
        row_p[4] = mf
        row_p[5 + i] = 1.0
        H_rows.append(row_p)
        y_rows.append(d['L_r'])
        w_rows.append(1.0 / 0.005**2)

        # Code
        row_c = np.zeros(N_STATE)
        row_c[0:3] = -los
        row_c[3] = 1.0
        row_c[4] = mf
        H_rows.append(row_c)
        y_rows.append(d['P_r'])
        w_rows.append(1.0 / 0.30**2)

    H = np.array(H_rows)
    y = np.array(y_rows)
    W = np.diag(w_rows)

    try:
        S = H @ P_prior @ H.T + np.linalg.inv(W)
        K = P_prior @ H.T @ np.linalg.inv(S)
        innovation = y - H @ x_prior

        print(f"\nEpoch {ep_idx} ({utc_dt}): {N_B} SVs: {active_svs}")
        print(f"  Innovation: mean={np.mean(innovation):.1f}m std={np.std(innovation):.1f}m")
        print(f"  H@x_prior: mean={np.mean(H @ x_prior):.1f}m std={np.std(H @ x_prior):.1f}m")
        print(f"  y: mean={np.mean(y):.1f}m std={np.std(y)/1000:.1f}km")

        dx = K @ innovation
        x = x_prior + dx

        print(f"  dx: pos=[{x[0]:.1f}, {x[1]:.1f}, {x[2]:.1f}]m clk={x[3]/1000:.2f}km trop={x[4]:.4f}m")
        print(f"  B states: first={x[5]:.1f}m last={x[-1]:.1f}m")

        P_post = (np.eye(N_STATE) - K @ H) @ P_prior
        for i, sv in enumerate(active_svs):
            b_states[sv] = float(x[5 + i])
            b_vars[sv] = float(P_post[5 + i, 5 + i])
        P_trop = float(P_post[4, 4])

        solved_pos = rcv_pos + x[:3]
        print(f"  |dx| = {np.linalg.norm(x[:3]):.1f}m")

        # Check reference
        ref = ref_orbit.get(utc_dt, rcv_pos)
        if utc_dt in ref_orbit:
            err = solved_pos - ref
            print(f"  3D error = {np.linalg.norm(err):.1f}m")
        elif ep_idx > 0:
            # Try interpolation
            pass

    except np.linalg.LinAlgError as e:
        print(f"  LIN ALG ERROR: {e}")

print(f"\nFinal B states:")
for sv, val in sorted(b_states.items()):
    print(f"  {sv}: B_if = {val:.1f}m")
