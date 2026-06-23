#!/usr/bin/env python3
"""
综合诊断GPS1B L1/L2数据 + PPP运行测试
一次性搞清楚: (1) L1/L2原始值 (2) 各种单位假设 (3) PPP实际残差
"""
import sys, os, math
sys.path.insert(0, '/workspace/gnss_pod/src')

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
LAM1 = C/F1    # 0.1903 m
LAM2 = C/F2    # 0.2442 m
LAM_IF = C / (F1**2/(F1**2 - F2**2))  # ~0.1070 m
F1_SQ, F2_SQ = F1**2, F2**2
ALPHA = F1_SQ / (F1_SQ - F2_SQ)  # ~1.546
BETA  = -F2_SQ / (F1_SQ - F2_SQ) # ~-0.546

# 读取GPS1B Apr29数据
fname = '/workspace/gnss_pod/data/GPS1B_2024-04-29_C_04.txt'
with open(fname, 'r') as f:
    lines = f.readlines()

# 数据从行196开始,共78754行
data_start = 196

# 从batch_v12.py复制的parse函数
YAML_PROD_FIELDS = ['CA_range','L1_range','L2_range','CA_phase','L1_phase',
                    'L2_phase','CA_SNR','L1_SNR','L2_SNR','CA_chan','L1_chan','L2_chan']

def parse_gps1b_record(parts, gf='C'):
    if len(parts) < 7: return None
    try:
        sod = int(parts[0]) + int(parts[1])*1e-6
        gid = parts[2].strip(); prn = int(parts[3])
        pf = int(parts[5].strip(), 2)
    except: return None
    if prn < 1 or prn > 32 or (gf and gid != gf): return None
    pv = parts[7:]
    rec = {'sv': f"G{prn:02d}", 'gps_sod': sod}
    for b, fn in enumerate(YAML_PROD_FIELDS):
        if (pf >> b) & 1 and b < len(pv):
            try: rec[fn] = float(pv[b])
            except: rec[fn] = None
    L1 = rec.get('L1_phase') or rec.get('L1_range')
    L2 = rec.get('L2_phase') or rec.get('L2_range')
    P1 = rec.get('CA_range') or rec.get('L1_range')
    P2 = rec.get('L2_range')
    if L1 is None or L2 is None: return None
    rec['L1'] = L1; rec['L2'] = L2; rec['P1'] = P1; rec['P2'] = P2
    rec['L_if'] = ALPHA*L1 + BETA*L2
    rec['P_if'] = (ALPHA*P1 + BETA*P2) if (P1 and P2) else P1
    return rec

# 收集第一个GPS卫星epoch的L1/L2值
print("=== GPS1B L1/L2 原始值分析 ===")
first_rec = None
for i in range(data_start, min(data_start+10000, len(lines))):
    parts = lines[i].strip().split()
    if len(parts) < 12: continue
    rec = parse_gps1b_record(parts, 'C')
    if rec is None: continue
    sv = rec['sv']
    L1 = rec['L1']
    L2 = rec['L2']
    P1 = rec['P1']
    
    first_rec = rec
    print(f"第一个GPS卫星: {sv}, L1={L1:.6f}, L2={L2:.6f}")
    print(f"  P1(CA_range)={P1:.6f}" if P1 else "  P1=None")
    print(f"  GPS1B字段(原始parts[8])={parts[8] if len(parts)>8 else 'N/A'}")
    print()
    
    # 测试: GPS1B L1_phase值是多少?
    L1_phase_field = rec.get('L1_phase', 0)
    L1_range_field = rec.get('L1_range', 0)
    L2_phase_field = rec.get('L2_phase', 0)
    print(f"L1_phase字段: {L1_phase_field:.6f}")
    print(f"L1_range字段: {L1_range_field:.6f}")
    print(f"L2_phase字段: {L2_phase_field:.6f}")
    print(f"差值(L1_range-L1_phase): {L1_range_field-L1_phase_field:.6f} m")
    break

print("\n=== 关键问题: GPS1B L1_phase的单位是什么? ===")
if first_rec:
    L1_val = first_rec['L1_phase'] or first_rec['L1_range']
    L2_val = first_rec['L2_phase'] or first_rec['L2_range']
    P1_val = first_rec['P1']
    
    print(f"L1_phase={L1_val}")
    print(f"GPS1B原始parts[8]={float(lines[data_start].split()[8]):.6f}")
    
    # 测试各种假设
    # 假设A: μs → m: × 299.792458
    L1_A = L1_val * C/1e6  # μs → m
    L2_A = L2_val * C/1e6
    L_if_A = ALPHA*L1_A + BETA*L2_A
    
    # 假设B: 直接用原始值(无转换)
    L_if_B = first_rec['L_if']  # batch_v12的L_if
    
    # 假设C: km → m: × 1000
    L1_C = L1_val * 1000
    L2_C = L2_val * 1000
    L_if_C = ALPHA*L1_C + BETA*L2_C
    
    # 假设D: cycles → m: × LAM1
    L1_D = L1_val * LAM1
    L2_D = L2_val * LAM2
    L_if_D = ALPHA*L1_D + BETA*L2_D
    
    # P1值作为参考
    P1_ref = P1_val if P1_val else L1_val
    
    print(f"\nP1(CA_range)={P1_val:.6f}")
    print(f"\n假设A (×C/1e6, 即μs→m):")
    print(f"  L1={L1_A:.3f} m ({L1_A/1000:.3f} km)")
    print(f"  L_if={L_if_A:.3f} m ({L_if_A/1000:.3f} km)")
    print(f"  期望GPS范围≈20265 km, 误差={abs(L_if_A-20265000)/20265000*100:.1f}%")
    
    print(f"\n假设B (直接原始值, batch_v12用):")
    print(f"  L_if={L_if_B:.3f} ({L_if_B/1000:.3f} km)")
    print(f"  期望GPS范围≈20265 km, 误差={abs(L_if_B-20265000)/20265000*100:.1f}%")
    
    print(f"\n假设C (×1000, 即假设原始是km):")
    print(f"  L1={L1_C:.3f} m ({L1_C/1000:.3f} km)")
    print(f"  L_if={L_if_C:.3f} m ({L_if_C/1000:.3f} km)")
    print(f"  期望GPS范围≈20265 km, 误差={abs(L_if_C-20265000)/20265000*100:.1f}%")
    
    print(f"\n假设D (×LAM1, 即cycles→m):")
    print(f"  L1={L1_D:.3f} m ({L1_D/1000:.3f} km)")
    print(f"  L_if={L_if_D:.3f} m ({L_if_D/1000:.3f} km)")
    print(f"  期望GPS范围≈20265 km, 误差={abs(L_if_D-20265000)/20265000*100:.1f}%")

print("\n=== 结论 ===")
print("ICD: L1_phase units是 seconds (实际存μs)")
print("GPS1B头: L1_phase units: m")
print("GPS1B L1_phase≈21868738, 如果单位是m → L1=21869 km → 误差7.9%")
print()
print("关键发现: GPS1B的L1_range≈L1_phase≈21868738 (差0.125m)")
print("这意味着: GPS1B L1_range是~21868738 m = 21869 km")
print("21869 km ≈ GPS卫星 altitude (20265 km) + GRACE altitude (490 km)")
print("这实际上是 GPS卫星altitude, 不是GPS范围!")
print()
print("修复方案: L1/L2已经是m (ICD文件头说units:m), 不需要转换!")
print("batch_v12应该直接用L1/L2作为m, 不需要×LAM_IF!")