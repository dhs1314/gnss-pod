#!/usr/bin/env python3
"""GPS1B完整诊断: 分析channel→PRN映射 + L1/L2单位"""
import sys, os, math
sys.path.insert(0, '/workspace/gnss_pod/src')
from gps1b_loader import parse_gps1b_record

C = 299792458.0
F1 = 1575.42e6
F2 = 1227.60e6
LAM1 = C/F1
LAM2 = C/F2
LAM_IF = C / (F1**2/(F1**2 - F2**2))
F1_SQ, F2_SQ = F1**2, F2**2

fname = '/workspace/gnss_pod/data/GPS1B_2024-04-29_C_04.txt'
with open(fname, 'r') as f:
    lines = f.readlines()

# 数据从行196开始
data_start = 196
print(f"数据行数: {len(lines)-data_start}")

# 直接分析前20行的原始字段
print("\n=== 前20个数据行(原始字段) ===")
for i in range(data_start, data_start+20):
    parts = lines[i].strip().split()
    if len(parts) < 10:
        continue
    chan = parts[3]
    L1_raw = float(parts[8]) if len(parts) > 8 else 0
    L1_ph = float(parts[10]) if len(parts) > 10 else 0
    print(f"行{i}: chan={chan}, L1_range={parts[7] if len(parts)>7 else '?'}, L1_phase(字段10)={L1_ph:.3f}")

# 解析前5个卫星并打印
print("\n=== 解析前5个卫星 ===")
count = 0
for i in range(data_start, min(data_start+5000, len(lines))):
    parts = lines[i].strip().split()
    if len(parts) < 10:
        continue
    rec = parse_gps1b_record(parts, 'C')
    if not rec:
        continue
    sv = rec.get('sv_prn', '?')
    L1_ph = rec.get('L1_phase')
    L2_ph = rec.get('L2_phase')
    L1_rng = rec.get('L1_range', 0)
    L2_rng = rec.get('L2_range', 0)
    
    print(f"\n行{i}: sv={sv}, L1_phase={L1_ph}, L2_phase={L2_ph}")
    if L1_ph and L2_ph:
        # 计算L_if在各种假设下
        L_if_m_μs = (F1_SQ*(L1_ph*C/1e6) - F2_SQ*(L2_ph*C/1e6))/(F1_SQ - F2_SQ)
        L_if_m_cyc = L1_ph * LAM_IF
        print(f"  L1_range(字段7)={L1_rng:.3f}")
        print(f"  假设μs: L_if(m)={L_if_m_μs:.1f} ({L_if_m_μs/1000:.3f} km)")
        print(f"  假设cycles×LAM_IF: L_if(m)={L_if_m_cyc:.1f} ({L_if_m_cyc/1000:.3f} km)")
    
    count += 1
    if count >= 5:
        break

print(f"\n解析了{count}个卫星")