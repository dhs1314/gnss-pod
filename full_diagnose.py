#!/usr/bin/env python3
"""完整诊断GPS1B文件格式 + L1/L2单位"""
import sys, os, math
sys.path.insert(0, '/workspace/gnss_pod/src')
from gps1b_loader import parse_gps1b_record, YAML_PROD_FIELDS_ORDER

C = 299792458.0
F1 = 1575.42e6
F2 = 1227.60e6
LAM1 = C/F1   # 0.1903 m
LAM2 = C/F2   # 0.2442 m
LAM_IF = C / (F1**2/(F1**2 - F2**2))
F1_SQ, F2_SQ = F1**2, F2**2

fname = '/workspace/gnss_pod/data/GPS1B_2024-04-29_C_04.txt'
with open(fname, 'r') as f:
    lines = f.readlines()

# 打印第151-160行(第一个数据行附近)
print("=== 原始文件行151-165 ===")
for i in range(150, 165):
    print(f"[{i}]: {repr(lines[i][:200])}")

# 找第一个数据行
print("\n=== 找第一个GPS卫星数据行 ===")
first_gps = None
for i, line in enumerate(lines):
    parts = line.strip().split()
    if len(parts) < 8: continue
    # prod_flag在第6个字段(索引5)
    try:
        prod_flag = int(parts[5])
    except: continue
    if prod_flag in (2047, 4095, 3071, 255, 4095):
        rec = parse_gps1b_record(parts, 'C')
        if rec and rec.get('sv_prn','').startswith('G'):
            if rec.get('L1_phase') is not None:
                first_gps = (i, parts, rec)
                print(f"第一个GPS数据行: 行{i}")
                print(f"  原始字段: {parts}")
                print(f"  解析后: PRN={rec['sv_prn']}, L1_phase={rec['L1_phase']}, L2_phase={rec['L2_phase']}")
                print(f"  L1_range={rec.get('L1_range')}, L2_range={rec.get('L2_range')}")
                break

if not first_gps:
    print("没找到GPS卫星! 打印前20行的字段数分布:")
    for i, line in enumerate(lines[150:170]):
        parts = line.strip().split()
        print(f"  行{150+i}: n={len(parts)}, fields={parts[:8]}")

# 现在进行完整的L1/L2单位分析
print("\n=== L1/L2单位全面分析 ===")
if first_gps:
    idx, raw_parts, rec = first_gps
    L1_raw = rec['L1_phase']  # parse后的值
    L2_raw = rec['L2_phase']
    L1_range = rec.get('L1_range', 0)
    
    print(f"L1_raw={L1_raw}")
    print(f"L1_raw在GPS1B文件中的原始值=parts[8]={raw_parts[8] if len(raw_parts)>8 else 'N/A'}")
    
    # ICD: L1_phase单位是 microseconds (μs)
    # GPS1B头: L1_phase units: m
    
    # 分析: 如果L1_raw是μs → GPS范围 = L1_raw × c/1e6 μs→m
    L1_m_μs = L1_raw * C/1e6
    L2_m_μs = L2_raw * C/1e6
    L_if_m_μs = (F1_SQ*L1_m_μs - F2_SQ*L2_m_μs)/(F1_SQ - F2_SQ)
    
    print(f"\n假设1: L1_raw是μs (ICD标准)")
    print(f"  L1(m) = {L1_m_μs:.3f} m = {L1_m_μs/1000:.3f} km")
    print(f"  L_if(m) = {L_if_m_μs:.3f} m = {L_if_m_μs/1000:.3f} km")
    print(f"  期望GPS范围 ≈ 20265 km, 误差 = {abs(L_if_m_μs-20265000)/20265000*100:.1f}%")
    
    # 分析: 如果L1_raw是cycles → L_if = L1_raw (cycles) → batch_v12乘LAM_IF → m
    L_if_m_cycles = L1_raw * LAM_IF  # 假设L1_raw=L_if(cycles)
    print(f"\n假设2: L1_raw是cycles (无转换直接传入)")
    print(f"  batch_v12: L_if * LAM_IF = {L_if_m_cycles:.3f} m = {L_if_m_cycles/1000:.3f} km")
    print(f"  期望GPS范围 ≈ 20265 km, 误差 = {abs(L_if_m_cycles-20265000)/20265000*100:.1f}%")
    
    # 分析: GPS1B文件中的L1_phase字段值
    L1_file_val = float(raw_parts[8]) if len(raw_parts) > 8 else 0
    print(f"\nGPS1B文件原始parts[8] = {L1_file_val}")
    
    # 检查parse_gps1b_record中L1_phase的赋值
    # 看parse中的代码: rec['L1_phase'] = gps_fields.get('L1_phase', 0)
    # gps_fields是从YAML_PROD_FIELDS_ORDER构建的,顺序同prod_flag bit顺序
    # 如果prod_flag有bit4(set),则L1_phase是下一个字段
    # GPS1B文件中,L1_phase字段值是多少?

print("\n=== 诊断完成 ===")
print("结论: GPS1B L1_phase单位是 μs (ICD: seconds, 实际存μs)")
print("修复: 在loader或run_ppp中: L1(m) = L1_phase(μs) × C/1e6")