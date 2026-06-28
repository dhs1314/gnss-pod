#!/usr/bin/env python3
"""Check D satellite position errors"""
import csv
from pathlib import Path

# Read C and D results
for sat in ['C', 'D']:
    tag = '20240429'
    csv_path = Path(f'output/ppp_rnx_{tag}.csv')
    if sat == 'C':
        # Re-run needed with different grace-id output name...
        pass

# Just check the last run (which was D)
import pandas as pd
try:
    df = pd.read_csv('output/ppp_rnx_20240429.csv')
    print("D satellite errors:")
    print(df[['dE_m', 'dN_m', 'dU_m', 'd3_m']].describe())
    print(f"\nEpochs with d3 > 5m: {(df['d3_m'] > 5).sum()}")
    print(f"Epochs with d3 > 10m: {(df['d3_m'] > 10).sum()}")
    print(f"\nFirst 20 epochs:")
    print(df[['dE_m', 'dN_m', 'dU_m', 'd3_m', 'n_sat']].head(20))
except ImportError:
    # Manual parsing
    d3_vals = []
    with open('output/ppp_rnx_20240429.csv') as f:
        reader = csv.DictReader(f)
        for row in reader:
            d3_vals.append(float(row['d3_m']))

    import numpy as np
    d3 = np.array(d3_vals)
    print(f"D satellite: {len(d3)} epochs")
    print(f"  d3 mean={np.mean(d3):.2f}m std={np.std(d3):.2f}m")
    print(f"  d3 min={np.min(d3):.2f}m max={np.max(d3):.2f}m")
    print(f"  d3 < 2m: {np.sum(d3 < 2)}/{len(d3)}")
    print(f"  d3 < 5m: {np.sum(d3 < 5)}/{len(d3)}")
    print(f"  d3 > 10m: {np.sum(d3 > 10)}/{len(d3)}")
    print(f"  First 5: {d3[:5]}")
    print(f"  Last 5: {d3[-5:]}")
