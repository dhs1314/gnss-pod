#!/usr/bin/env python3
"""
Apr29诊断: 完整测试GPS1B loader的L1/L2单位
"""
import sys, os
sys.path.insert(0, '/workspace/gnss_pod/src')

C = 299792458.0
F1 = 1575.42e6
F2 = 1227.60e6
LAM1 = C/F1
LAM2 = C/F2
LAM_IF = 0.10697  # ~0.1070 m

# 读GPS1B数据
fname = '/workspace/gnss_pod/data/GPS1B_2024-04-29_C_04.txt'
lines = open(fname, 'r').readlines()

# 找数据开始行
data_start = 0
for i, line in enumerate(lines):
    parts = line.strip().split()
    if len(parts) == 37 and parts[1] in ('US', 'CN'):
        data_start = i
        break

print(f"数据从第{data_start}行开始, 共{len(lines)-data_start}行")

# 分析前10个GPS卫星的L1/L2
from gps1b_loader import parse_gps1b_record

print("\n=== GPS卫星 L1/L2 分析 ===")
for line in lines[data_start:data_start+500]:
    parts = line.strip().split()
    if len(parts) != 37: continue
    rec = parse_gps1b_record(parts, 'C')
    if not rec or not rec.get('sv_prn','').startswith('G'): continue
    
    sv = rec['sv_prn']
    L1_phase = rec.get('L1_phase', 0)
    L2_phase = rec.get('L2_phase', 0)
    L1_range = rec.get('L1_range', 0)
    L2_range = rec.get('L2_range', 0)
    
    # 原始loader: L1 = L1_raw (cycles), 直接传入PPP
    # L_if = (F1²*L1 - F2²*L2)/(F1²-F2²) in cycles
    F1_SQ, F2_SQ = F1**2, F2**2
    L_if_cycles = (F1_SQ*L1_phase - F2_SQ*L2_phase)/(F1_SQ - F2_SQ)
    L_if_m = L_if_cycles * LAM_IF  # batch_v12中 L_if * LAM_IF
    
    print(f"{sv}: L1_phase={L1_phase:.3f} cycles, L2_phase={L2_phase:.3f} cycles")
    print(f"     L_if(cycles)={L_if_cycles:.3f}, L_if(m)={L_if_m:.1f} ({L_if_m/1000:.3f} km)")
    print(f"     期望GPS范围≈20265 km")
    break

# 测试: 如果L1_phase是μs而不是cycles,会怎样?
print("\n=== 假设L1_phase是μs (ICD标准) ===")
for line in lines[data_start:data_start+500]:
    parts = line.strip().split()
    if len(parts) != 37: continue
    rec = parse_gps1b_record(parts, 'C')
    if not rec or not rec.get('sv_prn','').startswith('G'): continue
    
    sv = rec['sv_prn']
    L1_raw = rec.get('L1_phase', 0)  # μs
    L2_raw = rec.get('L2_phase', 0)  # μs
    
    # 如果L1是μs,则转换为米: L1(μs) × C/1e6 = L1(m)
    L1_m = L1_raw * C/1e6
    L2_m = L2_raw * C/1e6
    
    # L_if(米) = (F1²*L1 - F2²*L2)/(F1²-F2²)  -- 这里L1,L2已经是米
    L_if_m2 = (F1_SQ*L1_m - F2_SQ*L2_m)/(F1_SQ - F2_SQ)
    
    print(f"{sv}: L1_phase(μs)={L1_raw:.3f} → L1(m)={L1_m:.1f}")
    print(f"     L_if(m)={L_if_m2:.1f} ({L_if_m2/1000:.3f} km)")
    print(f"     期望GPS范围=20265 km, 误差={abs(L_if_m2-20265000)/20265000*100:.1f}%")
    break

print("\n=== 结论: 正确的单位转换 ===")
print("GPS1B L1/L2字段是 μs (ICD: seconds)")
print("GPS1B L1/L2 在 loader 中不需要转换 (保持μs)")
print("run_ppp.py 的 ionospheric_free 需要 L1/L2 in meters")
print("→ 在 run_ppp 中: L1(m) = L1(μs) × C/1e6 = L1 × 299.792458")