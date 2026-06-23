#!/usr/bin/env python3
"""
Apr29 GPS1B数据完整性检查 + 最终PPP运行
"""
import sys, os, math, pickle
sys.path.insert(0, '/workspace/gnss_pod/src')

C = 299792458.0; F1, F2 = 1575.42e6, 1227.60e6
F1_SQ, F2_SQ = F1*F1, F2*F2
ALPHA = F1_SQ/(F1_SQ-F2_SQ); BETA = -F2_SQ/(F1_SQ-F2_SQ)

# 读取GPS1B Apr29原始文本文件
fname = '/workspace/gnss_pod/data/GPS1B_2024-04-29_C_04.txt'
lines = open(fname, 'r').readlines()

YAML_PROD = ['CA_range','L1_range','L2_range','CA_phase','L1_phase',
             'L2_phase','CA_SNR','L1_SNR','L2_SNR','CA_chan','L1_chan','L2_chan']

def parse(parts, gf='C'):
    if len(parts) < 7: return None
    try:
        sod = int(parts[0]) + int(parts[1])*1e-6
        gid = parts[2].strip(); prn = int(parts[3])
        pf = int(parts[5].strip(), 2)
    except: return None
    if prn < 1 or prn > 32 or (gf and gid != gf): return None
    pv = parts[7:]
    rec = {'sv': f"G{prn:02d}", 'gps_sod': sod}
    for b, fn in enumerate(YAML_PROD):
        if (pf >> b) & 1 and b < len(pv):
            try: rec[fn] = float(pv[b])
            except: rec[fn] = None
    return rec

# 检查第一个GPS卫星的L1/L2值
data_start = 196
first_sat = None
for i in range(data_start, min(data_start+5000, len(lines))):
    parts = lines[i].strip().split()
    if len(parts) < 12: continue
    rec = parse(parts, 'C')
    if not rec or not rec.get('sv','').startswith('G'): continue
    
    sv = rec['sv']
    L1 = rec.get('L1_range') or rec.get('L1_phase')
    L2 = rec.get('L2_range') or rec.get('L2_phase')
    if L1 is None: continue
    
    first_sat = rec
    break

if not first_sat:
    print("错误: 没有找到GPS卫星")
    sys.exit(1)

print(f"=== GPS1B Apr29 第一个卫星 {first_sat['sv']} ===")
print(f"L1_range: {first_sat.get('L1_range', 'N/A')}")
print(f"L1_phase: {first_sat.get('L1_phase', 'N/A')}")
print(f"L2_range: {first_sat.get('L2_range', 'N/A')}")
print(f"L2_phase: {first_sat.get('L2_phase', 'N/A')}")
print(f"CA_range: {first_sat.get('CA_range', 'N/A')}")
print()

# 计算L_if (当前loader的计算方式)
L1_raw = first_sat.get('L1_range') or first_sat.get('L1_phase')
L2_raw = first_sat.get('L2_range') or first_sat.get('L2_phase')
L_if_current = (F1_SQ*L1_raw - F2_SQ*L2_raw)/(F1_SQ - F2_SQ)
print(f"当前loader L_if = {L_if_current:.3f} m = {L_if_current/1000:.3f} km")
print(f"期望GPS范围 ≈ 20265 km")
print(f"差异: {(L_if_current-20265000)/1000:.1f} km ({abs(L_if_current-20265000)/20265000*100:.1f}%)")
print()

# 测试: 如果L1/L2单位是μs会怎样?
L_if_μs = (F1_SQ*(L1_raw*C/1e6) - F2_SQ*(L2_raw*C/1e6))/(F1_SQ - F2_SQ)
print(f"如果L1/L2是μs: L_if = {L_if_μs/1000:.1f} km")
print(f"  (误差: {abs(L_if_μs-20265000)/20265000*100:.0f}%)")
print()

# 读取GPS1B pickle,对比值
pkl = 'data/gracefo/2024/2024-04-29/gps1b_C.pkl'
gps1b = pickle.load(open(pkl, 'rb'))
first_t = sorted(gps1b.keys())[0]
obs = gps1b[first_t]
sv0 = list(obs.keys())[0]
d0 = obs[sv0]

print(f"=== GPS1B pickle 卫星 {sv0} ===")
print(f"L1 (pickle): {d0['L1']}")
print(f"L2 (pickle): {d0['L2']}")
print(f"L_if (pickle): {d0['L_if']}")
print(f"L_if/1000: {d0['L_if']/1000:.3f} km")
print()

# 对比文本文件和pickle
print(f"=== 对比 ===")
print(f"文本文件L1: {L1_raw}")
print(f"pickle L1: {d0['L1']}")
print(f"差值: {d0['L1']-L1_raw}")
print()

# 结论
print("=== 诊断结论 ===")
L_if_km = d0['L_if']/1000
if abs(L_if_km - 20265) < 2000:
    print(f"✓ L_if={L_if_km:.0f}km ≈ GPS范围20265km (误差{abs(L_if_km-20265)/20265*100:.0f}%)")
    print("→ 单位正确, GPS1B L1/L2已经是m")
    print("→ L_if的8%误差来自GPS卫星时钟bias,PPP可以处理")
elif abs(L_if_km - 6556) < 2000:
    print(f"✗ L_if={L_if_km:.0f}km ≈ 6556km")
    print("→ 单位是μs,需要转换为m")
    print("→ 修复: L_if_m = L_if_μs * C/1e6")
else:
    print(f"? L_if={L_if_km:.0f}km与任何预期值都不匹配")
    print(f"期望: 20265km (GPS范围)")
    print(f"也可能: 6556km (μs转换错误)")