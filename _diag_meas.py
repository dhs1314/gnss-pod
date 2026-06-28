#!/usr/bin/env python3
"""Diagnose P_r and L_r measurement values"""
import sys; sys.path.insert(0, '.')
from src.gps1a_loader import download_gps1a, gps_sod_to_utc
from src.ambiguity import compute_if_m
from src.sp3_loader import get_gps_pos_from_sp3
import pickle, numpy as np
from datetime import datetime, timedelta

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

gps_sods = sorted(gps_obs.keys())
sod0 = gps_sods[0]
utc0 = gps_sod_to_utc(sod0)

# Interpolate GNV1B to get rcv_pos
orbit_ts = sorted(ref_orbit.keys())
t0_orb = t1_orb = None
for j, ti in enumerate(orbit_ts):
    if ti >= utc0:
        t1_orb = ti
        t0_orb = orbit_ts[j-1] if j > 0 else ti
        break
    t0_orb = ti
if t0_orb is None:
    t0_orb = orbit_ts[0]
if t1_orb is None:
    t1_orb = orbit_ts[-1]

dt_frac = (utc0 - t0_orb).total_seconds()
dt_tot = (t1_orb - t0_orb).total_seconds()
if dt_tot < 0.01:
    rcv_pos = np.array(ref_orbit[t0_orb], dtype=float)
else:
    a = dt_frac / dt_tot
    rcv_pos = np.array(ref_orbit[t0_orb], dtype=float) * (1-a) + np.array(ref_orbit[t1_orb], dtype=float) * a

print(f'utc0={utc0} rcv_r={np.linalg.norm(rcv_pos)/1000:.1f}km')
print()

# Check first 5 SVs
for sv in sorted(gps_obs[sod0].keys())[:6]:
    rec = gps_obs[sod0][sv]
    L_if, P_if, B_if = compute_if_m(rec['L1_cyc'], rec['L2_cyc'], rec['P1'], rec['P2'])

    pos, clk, vel = get_gps_pos_from_sp3(sp3, sv, utc0)
    if pos is None: continue
    rho = float(np.linalg.norm(pos - rcv_pos))
    sag = (OMEGA_E / C) * (pos[0] * rcv_pos[1] - pos[1] * rcv_pos[0])
    rho_corr = rho + sag

    L_r = L_if + clk - rho_corr
    P_r = P_if + clk - rho_corr

    print(f'{sv}: P_if={P_if/1000:.3f}km L_if={L_if/1000:.3f}km')
    print(f'  rho={rho/1000:.3f}km clk_m={clk:.1f}m sag={sag:.1f}m')
    print(f'  P_r={P_r:.1f}m  L_r={L_r:.1f}m  B_if={B_if:.1f}m')
    print(f'  P_if-rho={P_if-rho:.1f}m')

# Aggregate across all SVs in first epoch
print(f'\n--- All SVs at epoch 0 ---')
P_r_all, L_r_all, names = [], [], []
for sv in sorted(gps_obs[sod0].keys()):
    rec = gps_obs[sod0][sv]
    pos, clk, vel = get_gps_pos_from_sp3(sp3, sv, utc0)
    if pos is None: continue
    L_if, P_if, B_if = compute_if_m(rec['L1_cyc'], rec['L2_cyc'], rec['P1'], rec['P2'])
    rho = float(np.linalg.norm(pos - rcv_pos))
    sag = (OMEGA_E / C) * (pos[0] * rcv_pos[1] - pos[1] * rcv_pos[0])
    rho_corr = rho + sag
    P_r = P_if + clk - rho_corr
    L_r = L_if + clk - rho_corr
    P_r_all.append(P_r)
    L_r_all.append(L_r)
    names.append(sv)

P_arr = np.array(P_r_all)
L_arr = np.array(L_r_all)
print(f'{len(names)} SVs: {names}')
print(f'P_r: mean={np.mean(P_arr):.1f}m  std={np.std(P_arr):.1f}m  min={np.min(P_arr):.1f}m  max={np.max(P_arr):.1f}m')
print(f'L_r: mean={np.mean(L_arr):.1f}m  std={np.std(L_arr):.1f}m')

# Now check: does clk_sv vary a lot between SVs?
print(f'\n--- SV Clock values ---')
for sv in sorted(gps_obs[sod0].keys())[:8]:
    pos, clk, vel = get_gps_pos_from_sp3(sp3, sv, utc0)
    if pos is not None:
        print(f'  {sv}: clk={clk:.1f}m = {clk/C*1e6:.3f}us')

# Check if P_if and rho_corr are well-matched
print(f'\n--- P_if vs rho ---')
for sv in sorted(gps_obs[sod0].keys())[:8]:
    rec = gps_obs[sod0][sv]
    pos, clk, vel = get_gps_pos_from_sp3(sp3, sv, utc0)
    if pos is None: continue
    L_if, P_if, B_if = compute_if_m(rec['L1_cyc'], rec['L2_cyc'], rec['P1'], rec['P2'])
    rho = float(np.linalg.norm(pos - rcv_pos))
    sag = (OMEGA_E / C) * (pos[0] * rcv_pos[1] - pos[1] * rcv_pos[0])
    print(f'  {sv}: P_if={P_if:.1f}m rho={rho:.1f}m diff={P_if-rho:.1f}m sag={sag:.1f}m')
