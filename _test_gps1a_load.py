#!/usr/bin/env python3
"""Quick test of GPS1A loader"""
import sys; sys.path.insert(0, '.')
from src.gps1a_loader import download_gps1a

gps_obs = download_gps1a(2024, 4, 29, grace_filter='C')
if gps_obs:
    sods = sorted(gps_obs.keys())
    print(f"Epochs: {len(sods)}")
    print(f"First: {sods[0]}, Last: {sods[-1]}")
    for sv, rec in sorted(gps_obs[sods[0]].items())[:3]:
        print(f"  {sv}: L1_cyc={rec['L1_cyc']:.3f} L2_cyc={rec['L2_cyc']:.3f} "
              f"P1={rec['P1']:.3f}m slip=({rec['slip_L1']},{rec['slip_L2']})")
else:
    print("FAILED to load GPS1A")
