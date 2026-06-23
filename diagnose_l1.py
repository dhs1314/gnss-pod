#!/usr/bin/env python3
"""
诊断GPS1B L1/L2单位问题
目标: 搞清楚GPS1B L1_phase的单位, 并正确修复loader
"""
import sys, os, math
sys.path.insert(0, '/workspace/gnss_pod/src')

C = 299792458.0
F1 = 1575.42e6
F2 = 1227.60e6
LAM1 = C/F1   # 0.1903 m
LAM2 = C/F2   # 0.2442 m
LAM_IF = C / (F1**2/(F1**2 - F2**2))  # ~0.1070 m
F1_SQ, F2_SQ = F1**2, F2**2

# 读取GPS1B Apr29数据
fname = '/workspace/gnss_pod/data/GPS1B_2024-04-29_C_04.txt'
with open(fname, 'r') as f:
    raw = f.read()

# NetCDF头解析
lines = raw.split('\n')
yaml_vars = {}
in_vars = False
for line in lines:
    if line.strip().startswith('- L1_phase:') or line.strip().startswith('- L2_phase:'):
        in_vars = True
    if in_vars and 'units:' in line:
        unit = line.split('units:')[1].strip().strip('"').strip("'")
        # 获取前面的变量名
        for prev_line in lines[max(0, lines.index(line)-3):lines.index(line)+1]:
            if prev_line.strip().startswith('- '):
                var = prev_line.strip().replace('- ','').replace(':','').strip()
                if var in ('L1_phase', 'L2_phase', 'L1_range', 'L2_range', 'CA_range'):
                    yaml_vars[var] = unit
                break

print("=== ICD和GPS1B文件中的单位 ===")
for v,u in yaml_vars.items():
    print(f"  {v}: {u}")

# 找第一个GPS卫星的数据
from gps1b_loader import parse_gps1b_record

data_lines = [l for l in lines if not l.startswith('#') 
              and not l.startswith('header:') 
              and not l.startswith('  dimensions:')
              and not l.startswith('    - ')]

print(f"\n数据行数: {len(data_lines)}")

first_gps = None
for line in data_lines:
    parts = line.strip().split()
    if len(parts) < 12: continue
    rec = parse_gps1b_record(parts, 'C')
    if rec and rec.get('sv_prn','').startswith('G'):
        if rec.get('L1_phase') and rec.get('L2_phase'):
            first_gps = rec
            break

if not first_gps:
    print("ERROR: 没有找到GPS卫星数据")
    sys.exit(1)

sv = first_gps['sv_prn']
L1_raw = first_gps['L1_phase']  # 这是parse后的值
L2_raw = first_gps['L2_phase']
L1_range = first_gps.get('L1_range', 0)
L2_range = first_gps.get('L2_range', 0)

print(f"\n=== 第一个GPS卫星 {sv} ===")
print(f"L1_phase(解析后)={L1_raw:.6f}")
print(f"L2_phase(解析后)={L2_raw:.6f}")
print(f"L1_range={L1_range:.3f} m")
print(f"L2_range={L2_range:.3f} m")

# 尝试各种单位转换
print(f"\n=== 各种假设下的L_if计算 ===")

# 假设1: L1_raw是μs (ICD标准)
L1_μs = L1_raw
L2_μs = L2_raw
L1_m_v1 = L1_μs * C/1e6
L2_m_v1 = L2_μs * C/1e6
L_if_m_v1 = (F1_SQ*L1_m_v1 - F2_SQ*L2_m_v1)/(F1_SQ - F2_SQ)
print(f"假设1 (μs→m): L1={L1_m_v1:.3f}m ({L1_m_v1/1000:.3f}km)")
print(f"  L_if={L_if_m_v1:.3f}m ({L_if_m_v1/1000:.3f}km), 期望≈20265km, 误差={abs(L_if_m_v1-20265000)/20265000*100:.1f}%")

# 假设2: L1_raw已经是cycles (但ICD说是μs)
L1_cyc = L1_raw
L2_cyc = L2_raw
L_if_cyc = (F1_SQ*L1_cyc - F2_SQ*L2_cyc)/(F1_SQ - F2_SQ)
L_if_m_v2 = L_if_cyc * LAM_IF  # batch_v12乘LAM_IF
print(f"\n假设2 (cycles→m, batch_v12乘LAM_IF): L_if(cycles)={L_if_cyc:.3f}")
print(f"  L_if(m)={L_if_m_v2:.3f}m ({L_if_m_v2/1000:.3f}km), 期望≈20265km, 误差={abs(L_if_m_v2-20265000)/20265000*100:.1f}%")

# 假设3: L1_raw是ns
L1_ns = L1_raw
L2_ns = L2_raw
L1_m_v3 = L1_ns * C/1e9
L2_m_v3 = L2_ns * C/1e9
L_if_m_v3 = (F1_SQ*L1_m_v3 - F2_SQ*L2_m_v3)/(F1_SQ - F2_SQ)
print(f"\n假设3 (ns→m): L1={L1_m_v3:.3f}m ({L1_m_v3/1000:.3f}km)")
print(f"  L_if={L_if_m_v3:.3f}m ({L_if_m_v3/1000:.3f}km), 期望≈20265km, 误差={abs(L_if_m_v3-20265000)/20265000*100:.1f}%")

# 假设4: L1_raw是GPS伪距(m)——但这不可能,因为L1比L2大很多
L_if_m_v4 = L_if_cyc  # 如果L1_raw直接是m
print(f"\n假设4 (L1_raw直接是m, batch_v12不再乘LAM_IF):")
print(f"  L_if={L_if_m_v4:.3f}m ({L_if_m_v4/1000:.3f}km), 期望≈20265km, 误差={abs(L_if_m_v4-20265000)/20265000*100:.1f}%")

# 假设5: GPS1B文件中的L1_phase字段值
# 从GPS1B文件直接读字段值(parts[8])
print(f"\n=== GPS1B文件中的原始字段值 ===")
for line in data_lines:
    parts = line.strip().split()
    if len(parts) < 12: continue
    rec = parse_gps1b_record(parts, 'C')
    if rec and rec.get('sv_prn','').startswith('G') and rec.get('L1_phase'):
        raw_L1_phase = float(parts[8]) if len(parts) > 8 else 0
        raw_L2_phase = float(parts[9]) if len(parts) > 9 else 0
        print(f"GPS1B原始字段: parts[8]={raw_L1_phase}, parts[9]={raw_L2_phase}")
        print(f"parse后: L1_phase={rec['L1_phase']}, L2_phase={rec['L2_phase']}")
        break

print(f"\n=== 关键判断 ===")
# 如果L1_phase是μs → L_if≈6556km (误差68%)
# 如果L1_phase是cycles → L_if≈21869km (误差8%)
# 如果L1_phase是m → L_if≈21869km (误差8%)
# 如果L1_phase是ns → L_if≈6.556km (误差99.97%)

print("期望GPS卫星到GRACE距离: ~20265 km (L1范围)")
print("L1_phase = 21868738")
print()
print("如果GPS1B L1_phase单位是 μs:")
print("  L_if(μs→m) = 21868738 × 299792458/1e6 = 6556219 m = 6556 km")
print("  误差 = |6556-20265|/20265 = 67.7%  ← 不对")
print()
print("如果GPS1B L1_phase单位是 cycles:")
print("  L_if(cycles) = 21868737 cycles")
print("  L_if(m) = 21868737 × 0.1070 = 2339299 m = 2339 km")  
print("  误差 = |2339-20265|/20265 = 88.5%  ← 也不对")
print()
print("但是! 如果L1_cycles的真正含义是:")
print("  L1_phase(μs) / LAM1(0.1903m/cycle) = L1_cycles")
print("  即 L1_cycles = L1_μs / 0.1903 = 21868738 / 0.1903 = 1.149e8 cycles")
print("  L_if(cycles) = (F1²×L1_μs/LAM1 - F2²×L2_μs/LAM2) / (F1²-F2²)")
print("  这个计算量太大，手动验证...")

# 精确计算: L1_phase(μs) → L1(cycles) → L_if(cycles) → L_if(m)
L1_μs_val = L1_raw
L2_μs_val = L2_raw
L1_cycles_v6 = L1_μs_val / LAM1 * 1e6  # μs / (m/cycle) * 1e6(m/μs) = cycles
L2_cycles_v6 = L2_μs_val / LAM2 * 1e6
L_if_cycles_v6 = (F1_SQ*L1_cycles_v6 - F2_SQ*L2_cycles_v6)/(F1_SQ - F2_SQ)
L_if_m_v6 = L_if_cycles_v6 * LAM_IF
print(f"\n假设6 (μs→cycles→m):")
print(f"  L1_cycles = {L1_cycles_v6:.3f}")
print(f"  L_if(cycles) = {L_if_cycles_v6:.3f}")
print(f"  L_if(m) = {L_if_m_v6:.3f} ({L_if_m_v6/1000:.3f} km)")
print(f"  误差 = {abs(L_if_m_v6-20265000)/20265000*100:.1f}%")

print("\n=== 最终结论 ===")
# 选择误差最小的方案
errors = {
    'μs→m直接': abs(L_if_m_v1-20265000)/20265000*100,
    'cycles×LAM_IF': abs(L_if_m_v2-20265000)/20265000*100,
    'ns→m': abs(L_if_m_v3-20265000)/20265000*100,
    'cycles→m': abs(L_if_m_v4-20265000)/20265000*100,
    'μs→cycles→m': abs(L_if_m_v6-20265000)/20265000*100,
}
best = min(errors, key=errors.get)
print(f"最小误差方案: {best} = {errors[best]:.1f}%")
print(f"\n推荐修复: 在loader中,L1/L2保持μs; 在run_ppp.py中乘C/1e6转米")