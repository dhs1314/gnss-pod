#!/usr/bin/env python3
"""Diagnose P_if-rho variation across SVs — is it per-SV DCB or something else?"""
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

gps_sods = sorted(gps_obs.keys())
orbit_ts = sorted(ref_orbit.keys())

# Interpolate GNV1B position
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

# Check 5 epochs spread over 30 minutes
step = max(1, len(gps_sods) // 10)
for idx in [0, 30, 60, 90, 120]:
    if idx >= len(gps_sods): break
    sod = gps_sods[idx]
    utc_dt = gps_sod_to_utc(sod)
    rcv_pos = get_rcv_pos(utc_dt)

    print(f'\nEpoch {idx}: gps_sod={sod} utc={utc_dt}')
    diffs = []
    for sv in sorted(gps_obs[sod].keys()):
        rec = gps_obs[sod][sv]
        pos, clk, vel = get_gps_pos_from_sp3(sp3, sv, utc_dt)
        if pos is None: continue
        L_if, P_if, B_if = compute_if_m(rec['L1_cyc'], rec['L2_cyc'], rec['P1'], rec['P2'])
        rho = float(np.linalg.norm(pos - rcv_pos))
        diffs.append((sv, P_if - rho, clk, rec['P1'], rec['P2']))

    d_arr = np.array([d[1] for d in diffs])
    print(f'  P_if-rho: mean={np.mean(d_arr)/1000:.1f}km std={np.std(d_arr)/1000:.1f}km '
          f'range=[{np.min(d_arr)/1000:.1f},{np.max(d_arr)/1000:.1f}]km')

    # Show per-SV values
    for sv, diff, clk, P1, P2 in diffs[:5]:
        print(f'  {sv}: P_if-rho={diff/1000:.1f}km  P1={P1/1e6:.3f}M  P2={P2/1e6:.3f}M  P2-P1={P2-P1:.1f}m')

# Now check: what if we compare against a simple rho from satellite position (no receiver clock)?
print('\n\n=== Checking if P1/P2 include receiver clock ===')
# For one SV across multiple epochs
sv_test = 'G05'
sv_diffs = []
for idx in range(0, min(200, len(gps_sods))):
    sod = gps_sods[idx]
    utc_dt = gps_sod_to_utc(sod)
    rcv_pos = get_rcv_pos(utc_dt)
    if sv_test not in gps_obs[sod]: continue
    rec = gps_obs[sod][sv_test]
    pos, clk, vel = get_gps_pos_from_sp3(sp3, sv_test, utc_dt)
    if pos is None: continue
    rho = float(np.linalg.norm(pos - rcv_pos))
    P_if, _, _ = compute_if_m(rec['L1_cyc'], rec['L2_cyc'], rec['P1'], rec['P2'])
    sv_diffs.append((idx, P_if - rho, rec['P1'], rec['P2']))

if sv_diffs:
    d_arr = np.array([d[1] for d in sv_diffs])
    print(f'{sv_test}: P_if-rho across {len(sv_diffs)} epochs: '
          f'mean={np.mean(d_arr)/1000:.1f}km std={np.std(d_arr)/1000:.1f}km')

# Compare P_if vs geometric range for a different SV
print('\n=== Check raw P1 values ===')
sod0 = gps_sods[0]
for sv in sorted(gps_obs[sod0].keys())[:5]:
    rec = gps_obs[sod0][sv]
    print(f'{sv}: P1={rec["P1"]:.3f}m  P2={rec["P2"]:.3f}m  '
          f'P1/c={rec["P1"]/C*1000:.3f}ms')

# For reference: what range do we expect?
print(f'\nGRACE pos norm: {np.linalg.norm(get_rcv_pos(gps_sod_to_utc(gps_sods[0])))/1000:.1f}km')
print(f'GPS sat pos norm: typically ~26600km')
print(f'Expected range: ~19000-27000km')
