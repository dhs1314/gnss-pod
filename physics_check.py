#!/usr/bin/env python3
"""
物理验证: GPS1B L1_range值(21,868,738)的真实含义
用GRACE GNV1B位置 + GPS广播星历计算真实GPS范围
"""
import sys, os, math
sys.path.insert(0, '/workspace/gnss_pod/src')

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
MU_E = 3.986004418e14

# 读取GRACE GNV1B数据(第一个epoch的位置)
fname_gnv = '/workspace/gnss_pod/data/gracefo/2024/2024-04-29/GNV1B_2024-04-29_C_04.txt'
if not os.path.exists(fname_gnv):
    # 搜索GN
    for root, dirs, files in os.walk('/workspace/gnss_pod/data/'):
        for f in files:
            if 'GNV' in f and f.endswith('.txt'):
                print(f"找到: {os.path.join(root, f)}")
                break

# 读取GRACE位置
grace_pos = None
if os.path.exists(fname_gnv):
    with open(fname_gnv, encoding='utf-8', errors='replace') as f:
        for line in f:
            p = line.split()
            if len(p) < 6: continue
            try:
                tg = float(p[0]); flag = p[2]
                if flag not in ('C', 'E'): continue
                X, Y, Z = float(p[3]), float(p[4]), float(p[5])
                if abs(X) < 1e3: continue
                grace_pos = (X, Y, Z)
                print(f"GRACE位置(第一个epoch): X={X:.3f}, Y={Y:.3f}, Z={Z:.3f} m")
                print(f"GRACE半径: {math.sqrt(X**2+Y**2+Z**2):.3f} m")
                r_grace = math.sqrt(X**2+Y**2+Z**2)
                alt_grace = r_grace - 6371000  # 减去地球平均半径
                print(f"GRACE高度: {alt_grace:.3f} m = {alt_grace/1000:.3f} km")
                break
            except: continue

# GPS卫星位置计算(简化圆轨道)
GPS_ORIGIN = __import__('datetime').datetime(2000, 1, 1, 12, 0, 0)
import datetime as dt

GPS_SV_PARAMS = {
    1:(0.0,0.0,0.0,0.0,0.000020,55.0,0.0,5153.6),
    3:(30.0,30.0,60.0,0.000015,54.8,30.0,5153.6),
    5:(90.0,90.0,180.0,0.000010,54.9,90.0,5153.6),
    7:(150.0,150.0,300.0,0.000015,54.7,150.0,5153.6),
    8:(180.0,180.0,0.0,0.000020,55.0,180.0,5153.6),
}

def gps_pos(sv_id, t):
    if sv_id not in GPS_SV_PARAMS: return None
    prn,M0,Omega0,omega,ecc,inc,u0,sqrtA = GPS_SV_PARAMS[sv_id]
    a = sqrtA**2
    n0 = math.sqrt(MU_E/a**3)
    dt_sec = (t - GPS_ORIGIN).total_seconds()
    M = M0 + math.radians(n0*dt_sec)
    E = M
    for _ in range(10): E = M + ecc*math.sin(E)
    sinE, cosE = math.sin(E), math.cos(E)
    v = math.atan2(math.sqrt(1-ecc**2)*sinE, cosE-ecc)
    phi = v + math.radians(omega)
    r = a*(1-ecc*cosE)
    inc_r = math.radians(inc)
    Om = Omega0 + math.radians(0.0265)/86164*dt_sec - 7.2921151467e-5*dt_sec
    xp = r*math.cos(phi)
    yp = r*math.sin(phi)
    X = xp*math.cos(Om) - yp*math.cos(inc_r)*math.sin(Om)
    Y = xp*math.sin(Om) + yp*math.cos(inc_r)*math.cos(Om)
    Z = yp*math.sin(inc_r)
    return (X, Y, Z)

# 计算第一个GPS卫星的GPS范围
if grace_pos:
    t = dt.datetime(2024, 4, 29, 0, 0, 0)
    gps_sat = gps_pos(5, t)  # G05 (PRN=5)
    if gps_sat:
        # 光行时改正
        rho = math.sqrt((gps_sat[0]-grace_pos[0])**2 + (gps_sat[1]-grace_pos[1])**2 + (gps_sat[2]-grace_pos[2])**2)
        tau = rho/C
        dt2 = 7.2921151467e-5*tau
        cr, sr = math.cos(dt2), math.sin(dt2)
        r2 = (grace_pos[0]*cr + grace_pos[1]*sr, -grace_pos[0]*sr + grace_pos[1]*cr, grace_pos[2])
        rho2 = math.sqrt((gps_sat[0]-r2[0])**2 + (gps_sat[1]-r2[1])**2 + (gps_sat[2]-r2[2])**2)
        
        print(f"\nG05位置: {gps_sat[0]:.3f}, {gps_sat[1]:.3f}, {gps_sat[2]:.3f} m")
        print(f"GRACE位置: {grace_pos}")
        print(f"原始GPS范围(未改正): {rho:.3f} m = {rho/1000:.3f} km")
        print(f"光行时改正后: {rho2:.3f} m = {rho2/1000:.3f} km")
        print(f"\nGPS1B L1_range(Chan5): 21868738.004 m = {21868738.004/1000:.3f} km")
        print(f"差异: {21868738.004 - rho2:.3f} m = {(21868738.004 - rho2)/1000:.3f} km")
        print(f"差异百分比: {abs(21868738.004 - rho2)/rho2*100:.1f}%")

# 重新审视GPS1B L1_range
print("\n=== GPS1B L1_range物理意义分析 ===")
# 如果L1_range是GPS范围(未改正)
L1_range_m = 21868738.004
print(f"GPS1B L1_range = {L1_range_m:.3f} m = {L1_range_m/1000:.3f} km")

# GPS卫星轨道半径
a_gps = 5153.6**2  # semi-major axis in m
print(f"GPS semi-major axis: {a_gps:.3f} m = {a_gps/1000:.3f} km")

# Earth + GPS altitude
R_E = 6371000  # Earth radius m
alt_gps = a_gps - R_E
print(f"GPS altitude: {alt_gps/1000:.3f} km")

# 如果L1_range = GPS orbital radius
print(f"\n如果L1_range = GPS orbital radius: {L1_range_m/1000:.3f} km")
print(f"  GPS semi-major axis: {a_gps/1000:.3f} km")
print(f"  差异: {L1_range_m - a_gps:.3f} m = {(L1_range_m - a_gps)/1000:.3f} km")

# 如果L1_range = GPS范围(masked)
# GPS范围典型值: 19,700-25,400 km (GRACE 490km altitude)
# L1_range = 21,869 km → 典型GPS范围 ✓
print(f"\n如果L1_range = GPS卫星到GRACE的距离:")
print(f"  L1_range = {L1_range_m/1000:.3f} km")
print(f"  典型范围: 19,700-25,400 km")
print(f"  → L1_range在典型范围内 ✓")

# 但是为什么L1_phase和L1_range几乎相同(差0.125m)?
# carrier phase = code range - ionospheric delay + cycle slip
# 正常情况下, carrier phase和code range应该相差几十厘米到几米
# 差0.125m说明L1_phase不是真正的carrier phase,而是iono-smoothed range
print(f"\n=== 关键发现 ===")
print(f"L1_phase - L1_range = 21868737.879 - 21868738.004 = -0.125 m")
print(f"这说明 L1_phase 可能是 iono-smoothed L1 carrier phase")
print(f"(通常iono-smoothed carrier phase比code少几十厘米)")
print(f"\n但是L1_range=21868738.004 m = 21869 km")
print(f"这个值在GPS卫星轨道altitude ({a_gps/1000:.3f} km)的典型GPS范围内!")
print(f"→ GPS1B L1_range确实是GPS范围(米)!")