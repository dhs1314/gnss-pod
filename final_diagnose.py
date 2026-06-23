#!/usr/bin/env python3
"""GPS1B完整诊断: 找到真实数据行并分析L1/L2单位"""
import sys, os, math
sys.path.insert(0, '/workspace/gnss_pod/src')
from gps1b_loader import parse_gps1b_record, N_HEADER_FIELDS, YAML_PROD_FIELDS_ORDER

C = 299792458.0
F1 = 1575.42e6
F2 = 1227.60e6
LAM1 = C/F1
LAM2 = C/F2
LAM_IF = C / (F1**2/(F1**2 - F2**2))
F1_SQ, F2_SQ = F1**2, F2**2

fname = '/workspace/gnss_pod/data/GPS1B_2024-04-29_C_04.txt'
with open(fname, 'r') as f:
    raw = f.read()

lines = raw.split('\n')
print(f"总行数: {len(lines)}")

# 找数据开始行
data_start = 0
for i, line in enumerate(lines):
    parts = line.strip().split()
    if len(parts) > 7:
        try:
            int(parts[0])
            data_start = i
            break
        except:
            pass

print(f"数据从行{data_start}开始")
print(f"第{data_start}行: {lines[data_start].rstrip()[:200]}")

# 检查前10个数据行的字段
print("\n=== 前10个数据行 ===")
for i in range(data_start, data_start+10):
    parts = lines[i].strip().split()
    if len(parts) > 7:
        try:
            pf = int(parts[5])
            sv_type = '?' 
            for j in range(data_start, i+1):
                rec = parse_gps1b_record(lines[j].strip().split(), 'C')
                if rec:
                    sv_type = rec.get('sv_prn', '?')
                    break
            print(f"行{i}: nfields={len(parts)}, prod={pf}, sv={sv_type}, fields={parts[:8]}")
        except Exception as e:
            print(f"行{i}: nfields={len(parts)}, error={e}")

# 找第一个GPS卫星(含L1_phase)
print("\n=== 找第一个GPS卫星(含L1_phase) ===")
first_gps = None
for i in range(data_start, min(data_start+20000, len(lines))):
    parts = lines[i].strip().split()
    if len(parts) < 10:
        continue
    try:
        prod_flag = int(parts[5])
        bit4 = (prod_flag >> 4) & 1
        if not bit4:
            continue
        rec = parse_gps1b_record(parts, 'C')
        if not rec:
            continue
        sv = rec.get('sv_prn', '')
        if not sv.startswith('G'):
            continue
        L1 = rec.get('L1_phase')
        L2 = rec.get('L2_phase')
        if L1 is None or L2 is None:
            continue
        first_gps = (i, parts, rec)
        print(f"行{i}: PRN={sv}, L1_phase={L1}, L2_phase={L2}")
        print(f"  L1_range={rec.get('L1_range',0)}, L2_range={rec.get('L2_range',0)}")
        
        # 各种单位假设下的L_if计算
        L1_m_μs = L1 * C/1e6
        L2_m_μs = L2 * C/1e6
        L_if_m_μs = (F1_SQ*L1_m_μs - F2_SQ*L2_m_μs)/(F1_SQ - F2_SQ)
        L_if_m_cycles = L1 * LAM_IF
        
        print(f"\n--- L1_phase单位分析 ---")
        print(f"L1_phase={L1}")
        print(f"假设μs: L1(m)={L1_m_μs:.1f}, L_if(m)={L_if_m_μs:.1f}")
        print(f"假设cycles×LAM_IF: L_if(m)={L_if_m_cycles:.1f}")
        print(f"期望GPS范围≈20265000 m")
        
        break
    except Exception as e:
        pass

if not first_gps:
    print("未找到带L1_phase的GPS卫星")
    # 打印prod_flag分布
    pfs = {}
    for i in range(data_start, min(data_start+1000, len(lines))):
        parts = lines[i].strip().split()
        if len(parts) > 5:
            try:
                pf = int(parts[5])
                pfs[pf] = pfs.get(pf, 0) + 1
            except:
                pass
    print(f"prod_flag分布: {dict(sorted(pfs.items())[:10])}")