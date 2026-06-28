#!/usr/bin/env python3
"""
验证GPS1B L1/L2的正确单位转换
ICD: L1_phase/L2_phase 是 Accumulated Phase in microseconds (μs)
GPS1B文件: L1_range/L2_range 单位是 meters

转换: L1_phase(μs) × C/1e6 = L1_phase(meters)
     L2_phase(μs) × C/1e6 = L2_phase(meters)
"""
import sys, os, math
sys.path.insert(0, '/workspace/gnss_pod/src')

C = 299792458.0
F1 = 1575.42e6
F2 = 1227.60e6
LAM1 = C/F1
LAM2 = C/F2
LAM_IF = C / (F1**2/(F1**2 - F2**2))  # ~0.1070 m

fname = '/workspace/gnss_pod/data/GPS1B_2024-04-29_C_04.txt'
lines = open(fname, 'r').readlines()

# 找数据行
data_lines = [l for l in lines if not l.startswith('#') and not l.startswith('header:') and not l.startswith('  dimensions:') and not l.startswith('    - ')]
# 找第一个GPS卫星
from gps1b_loader import parse_gps1b_record

print("=== 原始GPS1B L1/L2字段分析 ===")
count = 0
for line in lines:
    parts = line.strip().split()
    if len(parts) < 12: continue
    rec = parse_gps1b_record(parts, 'C')
    if not rec or not rec.get('sv_prn','').startswith('G'): continue
    if rec.get('L1_phase') is None or rec.get('L2_phase') is None: continue
    
    sv = rec['sv_prn']
    L1_raw = rec['L1_phase']  # μs
    L2_raw = rec['L2_phase']  # μs
    L1_range = rec.get('L1_range', 0)
    
    # 转换 μs → meters
    L1_m = L1_raw * C/1e6
    L2_m = L2_raw * C/1e6
    
    # L_if (无电离层组合, meters)
    F1_SQ, F2_SQ = F1**2, F2**2
    L_if_m = (F1_SQ*L1_m - F2_SQ*L2_m)/(F1_SQ - F2_SQ)
    
    print(f"{sv}: L1_raw={L1_raw:.6f} μs")
    print(f"     L1(m)={L1_m:.3f} ({L1_m/1000:.6f} km)")
    print(f"     L_if(m)={L_if_m:.3f} ({L_if_m/1000:.6f} km)")
    print(f"     期望GPS范围≈20265 km")
    print(f"     GPS范围误差={abs(L_if_m-20265000)/20265000*100:.2f}%")
    print(f"     L1_range={L1_range:.3f} m")
    print()
    count += 1
    if count >= 2: break

print("\n=== 结论 ===")
print("GPS1B L1_phase/L2_phase 单位是 μs (microseconds)")
print("转换: L1_meters = L1_μs × C/1e6 = L1_μs × 299792458/1e6")
print("或: L1_cycles = L1_μs × 1e6 / LAM1 (→ 再转meters)")
print()
print("修复方案:")
print("  src/gps1b_loader.py 中: L1 = L1_phase(μs), L2 = L2_phase(μs)")
print("  src/run_ppp.py 中: L1(m) = L1_μs × 299.792458")
print("  (或者在loader中: L1 = L1_phase × 299792458/1e6)")