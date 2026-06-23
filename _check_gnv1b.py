#!/usr/bin/env python3
"""Check GNV1B data quality"""
import numpy as np
import sys

orbit = {}
all_vals = []
for l in open('data/gracefo/2024/2024-04-29/GNV1B_2024-04-29_C_04.txt'):
    p = l.split()
    if len(p) < 6: continue
    try:
        t = float(p[0]); flag = p[2]
        all_vals.append((t, flag, float(p[3]), float(p[4]), float(p[5])))
        if flag in ('C', 'E'): orbit[t] = np.array([float(p[3]), float(p[4]), float(p[5])])
    except: pass

print(f'Total records: {len(all_vals)}')
flags = {}
for v in all_vals:
    flags[v[1]] = flags.get(v[1], 0) + 1
for k, v in sorted(flags.items()):
    print(f'  Flag {k}: {v}')

# Check specific epoch
gps_sod = 767620800.0
ts = sorted(orbit.keys())
print(f'\nOrbit epoch range: {ts[0]} to {ts[-1]}')
print(f'At gps_sod={gps_sod}:')
if gps_sod in orbit:
    print(f'  Exact match: {orbit[gps_sod]/1000} km')
else:
    for i, ti in enumerate(ts):
        if ti >= gps_sod:
            print(f'  t0={ts[i-1]}, t1={ti}')
            p0 = orbit[ts[i-1]]; p1 = orbit[ti]
            print(f'  p0={p0/1000} km')
            print(f'  p1={p1/1000} km')
            print(f'  distance: {np.linalg.norm(p1-p0):.3f} m')
            break

# Orbit velocity
diffs = []
for i in range(1, min(100, len(ts))):
    dt = ts[i] - ts[i-1]
    dp = np.linalg.norm(orbit[ts[i]] - orbit[ts[i-1]])
    diffs.append(dp / dt)
if diffs:
    print(f'\nOrbit velocity (first 100 epochs): mean={np.mean(diffs):.1f} m/s std={np.std(diffs):.3f} m/s')

# Check a few epochs to verify geocentric radius
for t in ts[:5]:
    p = orbit[t]
    r = np.linalg.norm(p)
    print(f'  t={t}: r={r/1000:.3f} km (alt={r/1000-6371:.1f} km)')
