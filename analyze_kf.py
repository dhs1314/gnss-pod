"""Analyze KF accuracy from CSV output."""
import csv, numpy as np

rows = []
with open('output/ppp_rd_20240429_C.csv') as f:
    for row in csv.DictReader(f):
        try:
            rows.append({k: float(v) for k, v in row.items()})
        except ValueError:
            rows.append({k: float(v) for k, v in row.items() if k != 'time'})

print(f'Total epochs: {len(rows)}')

d3 = np.array([r['d3_m'] for r in rows])
dE = np.array([r['dE_m'] for r in rows])
dN = np.array([r['dN_m'] for r in rows])
dU = np.array([r['dU_m'] for r in rows])

print(f'\n=== 3D Error Distribution ===')
print(f'  Mean:   {np.mean(d3):.3f} m')
print(f'  Median: {np.median(d3):.3f} m')
print(f'  Std:    {np.std(d3):.3f} m')
print(f'  RMS:    {np.sqrt(np.mean(d3**2)):.3f} m')
print(f'  <0.1m:  {np.sum(d3<0.1)} epochs ({100*np.sum(d3<0.1)/len(d3):.1f}%)')
print(f'  <0.2m:  {np.sum(d3<0.2)} epochs ({100*np.sum(d3<0.2)/len(d3):.1f}%)')
print(f'  <0.5m:  {np.sum(d3<0.5)} epochs ({100*np.sum(d3<0.5)/len(d3):.1f}%)')
print(f'  <1.0m:  {np.sum(d3<1.0)} epochs ({100*np.sum(d3<1.0)/len(d3):.1f}%)')

print(f'\n=== Component RMS ===')
print(f'  East:   {np.sqrt(np.mean(dE**2)):.3f} m')
print(f'  North:  {np.sqrt(np.mean(dN**2)):.3f} m')
print(f'  Up:     {np.sqrt(np.mean(dU**2)):.3f} m')

# Per-hour breakdown
print(f'\n=== Per-Hour RMS ===')
for h in range(2):
    mask = np.arange(len(d3)) // 120 == h  # 120 epochs per 30s = 1 hour
    if np.sum(mask) > 0:
        rms_h = np.sqrt(np.mean(d3[mask]**2))
        print(f'  Hour {h}: {rms_h:.3f} m ({np.sum(mask)} epochs)')

# Phase residual analysis
ph = np.array([r['ph_res_std_m'] for r in rows])
print(f'\n=== Phase Residual (std) ===')
print(f'  Mean:   {np.mean(ph):.3f} m')
print(f'  Median: {np.median(ph):.3f} m')
print(f'  Max:    {np.max(ph):.3f} m')

# Worst 5 epochs
worst_idx = np.argsort(d3)[-5:]
print(f'\n=== Worst 5 Epochs ===')
for idx in worst_idx:
    r = rows[idx]
    print(f'  ep {idx+1:3d}: 3D={r["d3_m"]:.3f}m E={r["dE_m"]:.3f} N={r["dN_m"]:.3f} U={r["dU_m"]:.3f} ph_res={r["ph_res_std_m"]:.3f} n_sat={int(r["n_sat"])}')

# Best 5 epochs
best_idx = np.argsort(d3)[:5]
print(f'\n=== Best 5 Epochs ===')
for idx in best_idx:
    r = rows[idx]
    print(f'  ep {idx+1:3d}: 3D={r["d3_m"]:.3f}m E={r["dE_m"]:.3f} N={r["dN_m"]:.3f} U={r["dU_m"]:.3f} ph_res={r["ph_res_std_m"]:.3f} n_sat={int(r["n_sat"])}')

# Correlation: 3D error vs phase residual
corr = np.corrcoef(d3, ph)[0, 1]
print(f'\n=== Correlation (3D err vs phase residual): {corr:.3f} ===')

# Number of satellites vs 3D error
n_sat = np.array([int(r['n_sat']) for r in rows])
for ns in sorted(set(n_sat)):
    mask = n_sat == ns
    if np.sum(mask) > 0:
        rms_ns = np.sqrt(np.mean(d3[mask]**2))
        print(f'  n_sat={ns}: {np.sum(mask)} epochs, 3D RMS={rms_ns:.3f}m')
