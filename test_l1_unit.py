#!/usr/bin/env python3
"""
GPS1B L1/L2符号和单位全面测试
用Apr29数据验证所有假设
"""
import sys, os, math, pickle
sys.path.insert(0, '/workspace/gnss_pod/src')

C = 299792458.0; F1, F2 = 1575.42e6, 1227.60e6
LAM1 = C/F1; LAM2 = C/F2
F1_SQ, F2_SQ = F1*F1, F2*F2
ALPHA = F1_SQ/(F1_SQ-F2_SQ); BETA = -F2_SQ/(F1_SQ-F2_SQ)

# 读取Apr29 GPS1B数据
pkl = 'data/gracefo/2024/2024-04-29/gps1b_C.pkl'
gps1b = pickle.load(open(pkl, 'rb'))
first_t = sorted(gps1b.keys())[0]
obs = gps1b[first_t]
sv0 = list(obs.keys())[0]
d0 = obs[sv0]

L1 = d0['L1']; L2 = d0['L2']
P1 = d0['P1']; P2 = d0['P2']

print(f"=== GPS1B Apr29 卫星 {sv0} ===")
print(f"L1 = {L1:.6f} ({L1/1000:.3f} km)")
print(f"L2 = {L2:.6f} ({L2/1000:.3f} km)")
print(f"P1 = {P1:.6f} ({P1/1000:.3f} km)")
print(f"P2 = {P2:.6f} ({P2/1000:.3f} km)")
print(f"L_if = {d0['L_if']:.6f} ({d0['L_if']/1000:.3f} km)")
print(f"P_if = {d0['P_if']:.6f} ({d0['P_if']/1000:.3f} km)")
print()

# 测试所有假设
tests = [
    # (名称, L1计算方式, L2计算方式)
    ('原始(m)', L1, L2),
    ('符号反(m)', -L1, -L2),
    ('μs→m', L1*C/1e6, L2*C/1e6),
    ('-μs→m', -L1*C/1e6, -L2*C/1e6),
    ('ns→m', L1*C/1e9, L2*C/1e9),
    ('-ns→m', -L1*C/1e9, -L2*C/1e9),
    ('cycles(LAM1)', L1*LAM1, L2*LAM2),
    ('-cycles(LAM1)', -L1*LAM1, -L2*LAM2),
]

# GPS理论范围 ~20265 km
R_GPS = 20265000  # m
print(f"期望GPS范围: {R_GPS/1000:.0f} km")
print()
print(f"{'假设':20s} {'L_if(km)':12s} {'差异(km)':12s} {'误差%':8s}")
print("-" * 55)
for name, L1_test, L2_test in tests:
    L_if = ALPHA*L1_test + BETA*L2_test
    diff_km = (L_if - R_GPS)/1000
    pct = abs(diff_km) / (R_GPS/1000) * 100
    status = "✓" if pct < 10 else "✗"
    print(f"{name:20s} {L_if/1000:12.3f} {diff_km:12.1f} {pct:7.1f}% {status}")

# 关键测试: P1/P2是否是有效的pseudo-range?
print()
print("=== P1/P2 分析 ===")
P_if_tests = [
    ('原始P_if', d0['P_if']),
    ('P1/2原始(m)', ALPHA*P1 + BETA*P2),
    ('-P1/2原始(m)', ALPHA*(-P1) + BETA*(-P2)),
]
for name, pif in P_if_tests:
    diff_km = (pif - R_GPS)/1000
    print(f"  {name}: P_if={pif/1000:.3f}km, 与GPS范围差={diff_km:.1f}km")

print()
print("=== 最终诊断 ===")
# 检查L1和P1的关系
print(f"L1 - P1 = {L1-P1:.3f} m = {L1-P1:.1f} cm")
print(f"L2 - P2 = {L2-P2:.3f} m = {L2-P2:.1f} cm")
print(f"L_if - P_if = {d0['L_if']-d0['P_if']:.3f} m = {d0['L_if']-d0['P_if']:.1f} cm")
print()

# 关键: L1_phase/L1_range差值
L1_phase = d0.get('L1_phase', 0)
L1_range = d0.get('L1_range', 0)
L2_phase = d0.get('L2_phase', 0)
L2_range = d0.get('L2_range', 0)
print(f"L1_phase - L1_range = {L1_phase-L1_range:.3f} m")
print(f"L2_phase - L2_range = {L2_phase-L2_range:.3f} m")

# 这个差值是否等于 iono delay?
iono_L1 = -BETA*(P2-P1) if P1 and P2 else 0
print(f"iono_L1 ≈ {iono_L1:.3f} m (BETA*(P2-P1))")
print(f"L1_phase ≈ L1_range - iono_L1? = {L1_range - iono_L1:.3f}m vs L1_phase={L1_phase:.3f}m")
if abs(L1_phase - (L1_range - iono_L1)) < 1:
    print("  ✓ 物理上合理: carrier phase ≈ range - iono")