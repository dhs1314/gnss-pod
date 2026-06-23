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

# GPS1B数据
gps1b_file = '/workspace/gnss_pod/data/gracefo/2024/2024-04-29/GPS1B_2024-04-29_C_04.txt'
# GRACE位置
gnv_file = '/workspace/gnss_pod/data/gracefo/2024/2024-04-29/GNV1B_2024-04-29_C_04.txt'

if not os.path.exists(gps1b_file):
    print(f"GPS1B文件不存在: {gps1b_file}")
    sys.exit(1)

print(f"GPS1B文件: {gps1b_file}")
print(f"GNV1B文件: {gnv_file} (存在: {os.path.exists(gnv_file)})")

# 读取GRACE位置
grace_pos = None
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
                print(f"GRACE位置: X={X:.3f}, Y={Y:.3f}, Z={Z:.3f} m, R={math.sqrt(X**2+Y**2+Z**2):.3f} m")
                break
            except: continue

if not grace_pos:
    print("警告: 找不到GRACE位置,用近似值(GRACE altitude 490km)")
    # GRACEorbital radius ≈ 6871 km
    grace_pos = (0, 0, 6871000)

print()
OMEGA_E = 7.2921151467e-5
MU_E = 3.986004418e14

# GPS广播星历 (简化参数)
GPS_SV = {
    'G05': {'M0': 0.0, 'Omega0': 0.0, 'omega': 0.0, 'ecc': 0.000020, 'inc': 55.0, 'sqrtA': 5153.6},
}

def gps_pos(sv, t_sod):
    if sv not in GPS_SV: return None
    p = GPS_SV[sv]
    a = p['sqrtA']**2
    n0 = math.sqrt(MU_E/a**3)
    M = math.radians(p['M0'] + n0 * t_sod)
    E = M
    for _ in range(8): E = M + p['ecc']*math.sin(E)
    v = math.atan2(math.sqrt(1-p['ecc']**2)*math.sin(E), math.cos(E)-p['ecc'])
    phi = v + math.radians(p['omega'])
    r = a*(1-p['ecc']*math.cos(E))
    Om = math.radians(p['Omega0']) - OMEGA_E*t_sod
    inc = math.radians(p['inc'])
    X = r*(math.cos(phi)*math.cos(Om) - math.sin(phi)*math.cos(inc)*math.sin(Om))
    Y = r*(math.cos(phi)*math.sin(Om) + math.sin(phi)*math.cos(inc)*math.cos(Om))
    Z = r*math.sin(phi)*math.sin(inc)
    return (X, Y, Z)

# 读取GPS1B数据
lines = open(gps1b_file, 'r').readlines()
data_start = 196

# 处理第一个epoch
current_t = None
obs_list = []
for i in range(data_start, min(data_start+5000, len(lines))):
    parts = lines[i].strip().split()
    if len(parts) < 12: continue
    rec = parse_gps1b_record(parts, 'C')
    if not rec: continue
    sv = rec['sv']
    t = rec['gps_sod']
    if current_t is None: current_t = t
    if t != current_t: break
    
    L1 = rec.get('L1_range') or rec.get('L1_phase')
    L2 = rec.get('L2_range') or rec.get('L2_phase')
    P1 = rec.get('CA_range') or rec.get('L1_range')
    P2 = rec.get('L2_range')
    if L1 is None: continue
    
    L_if = ALPHA*L1 + BETA*L2
    P_if = ALPHA*P1 + BETA*P2
    
    gp = gps_pos(sv, t)
    if not gp: continue
    
    rho_est = math.sqrt(sum((gp[j]-grace_pos[j])**2 for j in range(3)))
    tau = rho_est/C
    dtheta = OMEGA_E*tau
    cr, sr = math.cos(dtheta), math.sin(dtheta)
    rcv_rot = (grace_pos[0]*cr+grace_pos[1]*sr, -grace_pos[0]*sr+grace_pos[1]*cr, grace_pos[2])
    rho = math.sqrt(sum((gp[j]-rcv_rot[j])**2 for j in range(3)))
    
    obs_list.append({
        'sv': sv, 't': t, 'L_if': L_if, 'P_if': P_if, 'rho': rho,
        'resid_phase': L_if - rho, 'resid_code': P_if - rho,
        'L1': L1, 'L2': L2, 'P1': P1, 'P2': P2,
    })

print(f"=== Epoch {current_t:.0f}s, {len(obs_list)}个卫星 ===")
print(f"{'SV':5s} {'L_if(km)':12s} {'rho(km)':10s} {'res_phase(km)':14s} {'res_code(km)':12s}")
print("-" * 60)
for obs in sorted(obs_list, key=lambda x: x['sv']):
    print(f"{obs['sv']:5s} {obs['L_if']/1000:12.3f} {obs['rho']/1000:10.3f} {obs['resid_phase']/1000:14.3f} {obs['resid_code']/1000:12.3f}")

print()
print("=== 分析 ===")
if obs_list:
    phase_biases = [o['resid_phase']/1000 for o in obs_list]
    code_biases = [o['resid_code']/1000 for o in obs_list]
    print(f"载波相位残差: mean={sum(phase_biases)/len(phase_biases):.3f}km, std={math.sqrt(sum((x-sum(phase_biases)/len(phase_biases))**2 for x in phase_biases)/len(phase_biases)):.3f}km")
    print(f"码残差:       mean={sum(code_biases)/len(code_biases):.3f}km, std={math.sqrt(sum((x-sum(code_biases)/len(code_biases))**2 for x in code_biases)/len(code_biases)):.3f}km")
    
    # 检查是否有系统性时钟bias
    if max(phase_biases) - min(phase_biases) < 5:
        print("→ 残差一致性好 → 系统性误差(时钟bias)")
        print(f"→ 时钟偏置 ≈ {phase_biases[0]:.3f} km")
    else:
        print(f"→ 残差差异大(最大{ max(phase_biases)-min(phase_biases):.1f}km) → 非系统性")
    
    # 关键检查:L_if是否在GPS范围附近?
    L_if_mean = sum(o['L_if']/1000 for o in obs_list)/len(obs_list)
    rho_mean = sum(o['rho']/1000 for o in obs_list)/len(obs_list)
    print(f"\nL_if平均: {L_if_mean:.1f} km, GPS理论范围: {rho_mean:.1f} km")
    print(f"L_if vs GPS范围差异: {L_if_mean - rho_mean:.1f} km = {abs(L_if_mean-rho_mean)/rho_mean*100:.1f}%")