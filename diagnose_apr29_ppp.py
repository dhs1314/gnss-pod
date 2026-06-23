#!/usr/bin/env python3
"""
Apr29 PPP诊断: 直接计算理论测量值vs实际测量值,找出RMS大的原因
"""
import sys, os, math
sys.path.insert(0, '/workspace/gnss_pod/src')

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
LAM1 = C/F1; LAM2 = C/F2
F1_SQ, F2_SQ = F1*F1, F2*F2
ALPHA = F1_SQ/(F1_SQ-F2_SQ)
BETA  = -F2_SQ/(F1_SQ-F2_SQ)

# 读取GPS1B Apr29数据
fname = '/workspace/gnss_pod/data/GPS1B_2024-04-29_C_04.txt'
lines = open(fname, 'r').readlines()

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
    return rec

# 读取GRACE位置(来自GNV1B,用于计算理论GPS范围)
# 读取第一个epoch的GRACE位置
grace_pos = None
gnv_file = '/workspace/gnss_pod/data/gracefo/2024/2024-04-29/GNV1B_2024-04-29_C_04.txt'
if os.path.exists(gnv_file):
    with open(gnv_file, encoding='utf-8', errors='replace') as f:
        for line in f:
            p = line.split()
            if len(p) < 6: continue
            try:
                flag = p[2]
                if flag not in ('C',): continue
                X, Y, Z = float(p[3]), float(p[4]), float(p[5])
                if abs(X) < 1e3: continue
                grace_pos = (X, Y, Z)
                break
            except: continue

if not grace_pos:
    print("找不到GRACE位置数据")
    sys.exit(1)

print(f"GRACE位置: X={grace_pos[0]:.3f}, Y={grace_pos[1]:.3f}, Z={grace_pos[2]:.3f} m")
print(f"GRACE半径: {math.sqrt(sum(x**2 for x in grace_pos)):.3f} m")
print()

# 读取GPS广播星历(解析)
# GPS卫星位置(简化计算)
OMEGA_E = 7.2921151467e-5
MU_E = 3.986004418e14

GPS_PARAMS = {
    'G05': {'M0': 0.0, 'Omega0': 0.0, 'omega': 0.0, 'ecc': 0.000020, 'inc': 55.0, 'sqrtA': 5153.6, 'dn': 0.0, 'Omdot': 0.0},
}

def gps_position(sv, t_sod):
    """计算GPS卫星ECEF位置 (简化,无摄动)"""
    if sv not in GPS_PARAMS:
        return None
    p = GPS_PARAMS[sv]
    a = p['sqrtA']**2
    n0 = math.sqrt(MU_E/a**3)
    t = t_sod  # GPS second of day
    M = p['M0'] + math.radians(n0 * t)
    E = M
    for _ in range(10): E = M + p['ecc']*math.sin(E)
    v = math.atan2(math.sqrt(1-p['ecc']**2)*math.sin(E), math.cos(E)-p['ecc'])
    phi = v + math.radians(p['omega'])
    r = a*(1-p['ecc']*math.cos(E))
    Om = p['Omega0'] + math.radians(p['Omdot'])*t - OMEGA_E*t
    inc = math.radians(p['inc'])
    X = r*(math.cos(phi)*math.cos(Om) - math.sin(phi)*math.cos(inc)*math.sin(Om))
    Y = r*(math.cos(phi)*math.sin(Om) + math.sin(phi)*math.cos(inc)*math.cos(Om))
    Z = r*math.sin(phi)*math.sin(inc)
    return (X, Y, Z)

# 读取Apr29第一个epoch的GPS卫星和测量值
first_epoch = None
obs_list = []
data_start = 196

current_epoch = None
for i in range(data_start, min(data_start+2000, len(lines))):
    parts = lines[i].strip().split()
    if len(parts) < 12: continue
    rec = parse_gps1b_record(parts, 'C')
    if not rec: continue
    
    sv = rec['sv']
    t = rec['gps_sod']
    
    if current_epoch is None:
        current_epoch = t
    
    if t != current_epoch:
        break  # 处理完第一个epoch
    
    L1 = rec.get('L1_range') or rec.get('L1_phase')
    L2 = rec.get('L2_range') or rec.get('L2_phase')
    P1 = rec.get('CA_range') or rec.get('L1_range')
    P2 = rec.get('L2_range')
    if L1 is None: continue
    
    # 计算L_if (无电离层组合, meters)
    L_if = ALPHA*L1 + BETA*L2
    P_if = ALPHA*P1 + BETA*P2
    
    # 计算GPS卫星位置
    gps_pos = gps_position(sv, t)
    if not gps_pos: continue
    
    # 计算几何范围 (考虑光行时)
    rho_est = math.sqrt(sum((gps_pos[j]-grace_pos[j])**2 for j in range(3)))
    tau = rho_est/C
    dtheta = OMEGA_E * tau
    cr, sr = math.cos(dtheta), math.sin(dtheta)
    rcv_rot = (grace_pos[0]*cr + grace_pos[1]*sr, -grace_pos[0]*sr + grace_pos[1]*cr, grace_pos[2])
    rho = math.sqrt(sum((gps_pos[j]-rcv_rot[j])**2 for j in range(3)))
    
    # 残差 (测量值 - 理论值)
    resid_phase = L_if - rho  # 载波相位残差(m)
    resid_code = P_if - rho   # 码残差(m)
    
    obs_list.append({
        'sv': sv, 't': t,
        'L_if': L_if, 'P_if': P_if,
        'rho': rho,
        'resid_phase': resid_phase,
        'resid_code': resid_code,
        'L1': L1, 'L2': L2, 'P1': P1, 'P2': P2,
    })

print(f"=== 第一个epoch ({current_epoch:.0f}s) 的观测 ===")
print(f"卫星数: {len(obs_list)}")
print()

# 分析残差
print(f"{'SV':5s} {'L_if(km)':12s} {'P_if(km)':12s} {'rho(km)':10s} {'res_phase(km)':14s} {'res_code(km)':12s}")
print("-" * 75)
for obs in sorted(obs_list, key=lambda x: x['sv']):
    print(f"{obs['sv']:5s} {obs['L_if']/1000:12.3f} {obs['P_if']/1000:12.3f} {obs['rho']/1000:10.3f} {obs['resid_phase']/1000:14.3f} {obs['resid_code']/1000:12.3f}")

print()
# 统计
res_phase = [o['resid_phase'] for o in obs_list]
res_code = [o['resid_code'] for o in obs_list]
print(f"载波相位残差: mean={sum(res_phase)/len(res_phase)/1000:.1f}km, std={math.sqrt(sum((x-sum(res_phase)/len(res_phase))**2 for x in res_phase)/len(res_phase))/1000:.1f}km")
print(f"码残差:       mean={sum(res_code)/len(res_code)/1000:.1f}km, std={math.sqrt(sum((x-sum(res_code)/len(res_code))**2 for x in res_code)/len(res_code))/1000:.1f}km")

print()
print("=== 关键发现 ===")
# 检查残差的系统性
phase_means = [o['resid_phase']/1000 for o in obs_list]
print(f"所有卫星载波相位残差: {phase_means[0]:.1f} km (应该相似,如果系统误差)")
if max(phase_means) - min(phase_means) < 10:
    print("→ 所有卫星残差相似 → 系统性误差(接收机时钟bias) ✓")
    print(f"→ 接收机时钟偏置 ≈ {phase_means[0]:.1f} km")
else:
    print("→ 不同卫星残差差异大 → 非系统性误差")

# 检查GPS范围是否合理
print(f"\nGPS范围(rho): {obs_list[0]['rho']/1000:.1f} km")
print(f"L_if: {obs_list[0]['L_if']/1000:.1f} km")
print(f"P_if: {obs_list[0]['P_if']/1000:.1f} km")
print(f"期望GPS范围: ~20265 km")
print(f"L_if误差: {abs(obs_list[0]['L_if']/1000-20265):.1f} km ({abs(obs_list[0]['L_if']/1000-20265)/20265*100:.1f}%)")