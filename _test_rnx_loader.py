#!/usr/bin/env python3
"""Quick test of GPS1B RINEX loader"""
import sys; sys.path.insert(0, 'src')
from gps1b_rnx_loader import load_gps1b_rnx

data = load_gps1b_rnx('data/GPS1B_2024-04-29_C_04.rnx')
if data:
    sods = sorted(data.keys())
    print(f'Epochs: {len(sods)}')
    print(f'First SOD: {sods[0]}, Last SOD: {sods[-1]}')
    first = data[sods[0]]
    print(f'SVs in first epoch: {sorted(first.keys())}')
    for sv in sorted(first.keys())[:3]:
        rec = first[sv]
        L1_m = rec['L1_cyc'] * 299792458.0 / 1575.42e6
        print(f'  {sv}: L1={rec["L1_cyc"]:.0f}cyc L2={rec["L2_cyc"]:.0f}cyc P1={rec["P1"]:.1f}m P2={rec["P2"]:.1f}m')
        print(f'       L1_m={L1_m:.3f}m')
        print(f'       P_if={rec["P_if"]:.3f}m L_if={rec["L_if"]:.3f}m B_if={rec["L_if"]-rec["P_if"]:.2f}m')
else:
    print('FAILED to load data')
