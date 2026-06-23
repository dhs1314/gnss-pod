#!/usr/bin/env python3
"""
GPS1B L1/L2全面测试: 尝试所有可能的单位转换,找出物理上正确的那个
"""
import sys, os, math
sys.path.insert(0, '/workspace/gnss_pod/src')

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
LAM1 = C/F1   # 0.1903 m
LAM2 = C/F2   # 0.2442 m
F1_SQ, F2_SQ = F1*F1, F2*F2
ALPHA = F1_SQ/(F1_SQ - F2_SQ)  # 1.546
BETA  = -F2_SQ/(F1_SQ - F2_SQ)  # -0.546
LAM_IF = C / (F1_SQ/(F1_SQ-F2_SQ))  # 0.1070 m

# 真实GPS1B数据
L1_RANGE = 21868738.004  # meters (per ICD)
L1_PHASE = 21868738.163  # meters (per ICD header)
L2_RANGE = 21868738.441  # meters
L2_PHASE = 21868738.266  # meters

# GPS参数
GPS_SMA = 26560e3       # GPS semi-major axis: 26,560 km
R_EARTH  = 6371e3       # Earth radius: 6,371 km
R_GRACE  = 6871e3       # GRACE orbital radius: ~6,871 km (490 km altitude)
R_GPS_TYPICAL = 20265e3  # 典型GPS范围: ~20,265 km

print("=== GPS1B 原始值 ===")
print(f"L1_RANGE={L1_RANGE}  ({L1_RANGE/1000:.3f} km)")
print(f"L1_PHASE={L1_PHASE}  ({L1_PHASE/1000:.3f} km)")
print(f"L2_RANGE={L2_RANGE}  ({L2_RANGE/1000:.3f} km)")
print(f"L2_PHASE={L2_PHASE}  ({L2_PHASE/1000:.3f} km)")
print()

# 关键: 如果L1/L2已经是m(ICD), L_if应该约为~21868 km
# 但GPS范围应该是~20265 km
# 差: 21868-20265=1603 km = 8%
# 这个差是否是一个常数bias?

# 测试1: 如果L1_RANGE=GRACE半径(不是GPS范围!)
print("=== 测试1: L1_RANGE = GRACE orbital radius? ===")
if abs(L1_RANGE/1000 - R_GRACE/1000) < 1000:
    print(f"  MATCH! L1_RANGE = {L1_RANGE/1000:.0f} km ≈ GRACE半径 {R_GRACE/1000:.0f} km")
    print(f"  → L1_RANGE 不是GPS范围,而是GRACE轨道半径!")
    print(f"  → 但L1_RANGE(GRACE半径=6871 km) vs 实际GRACE: {R_GRACE/1000:.0f} km")
    print(f"  → 差: {L1_RANGE/1000 - R_GRACE/1000:.0f} km")
else:
    print(f"  NO MATCH. L1_RANGE={L1_RANGE/1000:.0f} km, GRACE={R_GRACE/1000:.0f} km")

print()
print("=== 所有可能的单位转换测试 ===")

# 转换1: 直接用原始值(m)
L_if_1 = ALPHA*L1_RANGE + BETA*L2_RANGE
print(f"1. 直接原始值(m): L_if={L_if_1/1000:.3f} km, 期望={R_GPS_TYPICAL/1000:.0f} km, 差={abs(L_if_1-R_GPS_TYPICAL)/R_GPS_TYPICAL*100:.1f}%")

# 转换2: L1/L2是μs(ICD标准)→m
L_if_2 = ALPHA*(L1_RANGE*C/1e6) + BETA*(L2_RANGE*C/1e6)
print(f"2. μs→m (×C/1e6): L_if={L_if_2/1000:.3f} km, 期望={R_GPS_TYPICAL/1000:.0f} km, 差={abs(L_if_2-R_GPS_TYPICAL)/R_GPS_TYPICAL*100:.1f}%")

# 转换3: L1/L2是cycles (1 cycle = LAM1 m) → m
L_if_3 = ALPHA*(L1_RANGE*LAM1) + BETA*(L2_RANGE*LAM2)
print(f"3. cycles→m (×LAM): L_if={L_if_3/1000:.3f} km, 期望={R_GPS_TYPICAL/1000:.0f} km, 差={abs(L_if_3-R_GPS_TYPICAL)/R_GPS_TYPICAL*100:.1f}%")

# 转换4: L1/L2是ns→m
L_if_4 = ALPHA*(L1_RANGE*C/1e9) + BETA*(L2_RANGE*C/1e9)
print(f"4. ns→m (×C/1e9): L_if={L_if_4/1000:.3f} km, 期望={R_GPS_TYPICAL/1000:.0f} km, 差={abs(L_if_4-R_GPS_TYPICAL)/R_GPS_TYPICAL*100:.1f}%")

# 转换5: 如果L1/L2原始值需要÷1000 (假设km→m)
L_if_5 = ALPHA*(L1_RANGE*1000) + BETA*(L2_RANGE*1000)
print(f"5. ×1000 (km→m): L_if={L_if_5/1000:.3f} km, 期望={R_GPS_TYPICAL/1000:.0f} km, 差={abs(L_if_5-R_GPS_TYPICAL)/R_GPS_TYPICAL*100:.1f}%")

# 转换6: 假设L1/L2是GPS时间(秒)→m
L_if_6 = ALPHA*(L1_RANGE*C) + BETA*(L2_RANGE*C)
print(f"6. 秒→m (×C): L_if={L_if_6/1000:.3f} km, 期望={R_GPS_TYPICAL/1000:.0f} km, 差={abs(L_if_6-R_GPS_TYPICAL)/R_GPS_TYPICAL*100:.1f}%")

# 转换7: 如果L1_RANGE=GRACE半径+GPSaltitude (即L1_RANGE=32867 km)
# 不, 21868 km ≠ 32867 km

# 关键发现测试
print()
print("=== 关键: GPS1B值代表什么? ===")
print(f"L1_RANGE={L1_RANGE/1000:.3f} km")
print(f"L1_PHASE={L1_PHASE/1000:.3f} km")
print(f"差(L1_RANGE-L1_PHASE)={abs(L1_RANGE-L1_PHASE):.3f} m = {abs(L1_RANGE-L1_PHASE)/LAM1:.3f} cycles")
print()

# 如果L1_RANGE是GPS light travel time in μs → m
L1_μs_to_m = L1_RANGE * C/1e6
L2_μs_to_m = L2_RANGE * C/1e6
L_if_μs = ALPHA*L1_μs_to_m + BETA*L2_μs_to_m
print(f"如果L1_RANGE单位是μs: L1={L1_μs_to_m/1000:.3f} km")
print(f"  L_if={L_if_μs/1000:.3f} km")
print(f"  期望GPS范围≈{R_GPS_TYPICAL/1000:.0f} km, 差={abs(L_if_μs-R_GPS_TYPICAL)/R_GPS_TYPICAL*100:.1f}%")
print()

# GPS1B值 ≈ GRACE orbital radius in km?
grace_km = R_GRACE/1000  # 6871 km
gps_sma_km = GPS_SMA/1000  # 26560 km
gps_alt_km = gps_sma_km - R_EARTH/1000  # 20189 km

print(f"GRACE orbital radius: {grace_km:.0f} km")
print(f"GPS orbital radius: {gps_sma_km:.0f} km") 
print(f"GPS altitude: {gps_alt_km:.0f} km")
print()

# L1_RANGE=21868 km接近GPS altitude (20189 km)!
# 差: 21868-20189=1679 km = 8.3%
print(f"!!! L1_RANGE({L1_RANGE/1000:.0f} km) ≈ GPS altitude ({gps_alt_km:.0f} km)")
print(f"    差: {L1_RANGE/1000-gps_alt_km:.0f} km = {abs(L1_RANGE/1000-gps_alt_km)/gps_alt_km*100:.1f}%")

# 或者 L1_RANGE ≈ GPS orbital radius in different units?
print(f"\n    L1_RANGE({L1_RANGE/1000:.0f} km) vs GPS orbital radius: {gps_sma_km:.0f} km")
print(f"    差: {L1_RANGE/1000-gps_sma_km:.0f} km = {abs(L1_RANGE/1000-gps_sma_km)/gps_sma_km*100:.1f}%")

# 最终判断
print()
print("=== 最终结论 ===")
print("GPS1B L1_range和L1_phase值约为21,868 km")
print("这个值接近但不等于:")
print(f"  - GPS altitude:  {gps_alt_km:.0f} km (差{abs(L1_RANGE/1000-gps_alt_km):.0f} km)")
print(f"  - GPS SMA:       {gps_sma_km:.0f} km (差{abs(L1_RANGE/1000-gps_sma_km):.0f} km)")
print(f"  - GPS range:     ~{R_GPS_TYPICAL/1000:.0f} km (差{abs(L1_RANGE/1000-R_GPS_TYPICAL/1000):.0f} km)")
print()
print("→ GPS1B的L1/L2不是直接的GPS范围,而是某种派生量")
print("→ 最可能: L1/L2已经是正确单位(m), 不需要转换")
print("→ L_if ≈ 21868 km, 误差8% (可能被时钟bias吸收)")
print()
print("修复方案: 在loader中,确保L1/L2以m为单位:")
print("  L1(μs) = GPS1B_L1_phase × C/1e6 → m")
print("  L2(μs) = GPS1B_L2_phase × C/1e6 → m")
print("  L_if = α×L1 + β×L2 → m")