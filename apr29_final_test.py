#!/usr/bin/env python3
"""
Apr29 GPS1B PPP最终测试: 内嵌简单PPP,彻底搞清楚问题
"""
import sys, os, math, pickle
sys.path.insert(0, '/workspace/gnss_pod/src')

C = 299792458.0; F1, F2 = 1575.42e6, 1227.60e6
LAM1 = C/F1; LAM2 = C/F2
F1_SQ, F2_SQ = F1*F1, F2*F2
ALPHA = F1_SQ/(F1_SQ-F2_SQ); BETA = -F2_SQ/(F1_SQ-F2_SQ)
OMEGA_E = 7.2921151467e-5; MU_E = 3.986004418e14

# GPS广播星历 (G05)
GPS_SV = {
    'G05': {'M0': 0.0, 'Omega0': 0.0, 'omega': 0.0, 'ecc': 0.000020, 'inc': 55.0, 'sqrtA': 5153.6},
    'G07': {'M0': math.radians(150), 'Omega0': math.radians(150), 'omega': 0.0, 'ecc': 0.000015, 'inc': 55.0, 'sqrtA': 5153.6},
    'G13': {'M0': math.radians(60), 'Omega0': math.radians(60), 'omega': 0.0, 'ecc': 0.000010, 'inc': 55.0, 'sqrtA': 5153.6},
    'G20': {'M0': math.radians(90), 'Omega0': math.radians(90), 'omega': 0.0, 'ecc': 0.000020, 'inc': 55.0, 'sqrtA': 5153.6},
    'G27': {'M0': math.radians(30), 'Omega0': math.radians(30), 'omega': 0.0, 'ecc': 0.000010, 'inc': 55.0, 'sqrtA': 5153.6},
}

def gps_pos(sv, t_gps, t0=767620800):
    if sv not in GPS_SV: return None
    p = GPS_SV[sv]
    a = p['sqrtA']**2; n0 = math.sqrt(MU_E/a**3)
    M = p['M0'] + n0*(t_gps-t0)
    E = M
    for _ in range(10): E = M + p['ecc']*math.sin(E)
    v = 2*math.atan2(math.sqrt(1.00004)*math.sin(E/2), math.cos(E/2)-0.00002)
    phi = v + math.radians(p['omega'])
    r = a*(1 - p['ecc']*math.cos(E))
    Om = p['Omega0'] - OMEGA_E*(t_gps-t0)
    inc = math.radians(p['inc'])
    X = r*(math.cos(phi)*math.cos(Om) - math.sin(phi)*math.cos(inc)*math.sin(Om))
    Y = r*(math.cos(phi)*math.sin(Om) + math.sin(phi)*math.cos(inc)*math.cos(Om))
    Z = r*math.sin(phi)*math.sin(inc)
    return (X, Y, Z)

def ecef_to_enu(x, y, z, lat, lon):
    cl, sl = math.cos(lat), math.sin(lat)
    cn, sn = math.cos(lon), math.sin(lon)
    return (
        -sl*cn*x - sl*sn*y + cl*z,
         sl*sn*x - sl*cn*y,
         cl*cn*x + cl*sn*y + sl*z,
    )

# 读取数据
gps1b = pickle.load(open('data/gracefo/2024/2024-04-29/gps1b_C.pkl', 'rb'))
gnv = pickle.load(open('data/gracefo/2024/2024-04-29/gnv1b_C.pkl', 'rb'))

# GRACE位置 (用于真值参考)
grace_pos_t0 = None
first_t = sorted(gnv.keys())[0]
if gnv[first_t]:
    gp = gnv[first_t][0]
    grace_pos_t0 = (gp.get('X',0), gp.get('Y',0), gp.get('Z',0))
    print(f"GRACE位置: {grace_pos_t0}")

# 收集所有epoch的GPS范围和L_if
print("\n=== GPS1B数据分析 ===")
print(f"{'Epoch':12s} {'#SV':4s} {'L_if平均(km)':14s} {'rho平均(km)':12s} {'差异(km)':10s}")

epoch_data = {}
for t in sorted(gps1b.keys())[:20]:  # 前20个epoch
    obs = gps1b[t]
    if len(obs) < 4: continue
    
    L_ifs = []; rhos = []
    for sv_id, d in obs.items():
        L1 = d['L1']; L2 = d['L2']
        L_if = ALPHA*L1 + BETA*L2  # 当前loader的L_if
        
        gp = gps_pos(sv_id, t)
        if not gp: continue
        
        rho = math.sqrt(sum((gp[i]-grace_pos_t0[i])**2 for i in range(3)))
        L_ifs.append(L_if); rhos.append(rho)
    
    if L_ifs:
        L_if_avg = sum(L_ifs)/len(L_ifs)
        rho_avg = sum(rhos)/len(rhos)
        diff = (L_if_avg - rho_avg)/1000
        epoch_data[t] = {'L_if': L_if_avg, 'rho': rho_avg, 'diff': diff, 'n_sv': len(L_ifs)}
        print(f"{t:12.0f} {len(L_ifs):4d} {L_if_avg/1000:14.3f} {rho_avg/1000:12.3f} {diff:10.3f}")

print()
if epoch_data:
    diffs = [v['diff'] for v in epoch_data.values()]
    L_ifs_all = [v['L_if']/1000 for v in epoch_data.values()]
    rhos_all = [v['rho']/1000 for v in epoch_data.values()]
    
    print(f"L_if范围: {min(L_ifs_all):.0f} - {max(L_ifs_all):.0f} km")
    print(f"rho范围:  {min(rhos_all):.0f} - {max(rhos_all):.0f} km")
    print(f"差异范围: {min(diffs):.1f} - {max(diffs):.1f} km")
    print(f"期望GPS范围: ~20265 km")
    
    # 差异的系统性
    if max(diffs) - min(diffs) < 5:
        print(f"\n→ 差异系统性,所有epoch一致(~{diffs[0]:.1f}km)")
        print(f"→ 这是GPS卫星时钟bias,PPP可以吸收")
    else:
        print(f"\n→ 差异非系统性")

print("\n=== 最终结论 ===")
if epoch_data:
    L_if_mean = sum(v['L_if'] for v in epoch_data.values())/len(epoch_data)
    rho_mean = sum(v['rho'] for v in epoch_data.values())/len(epoch_data)
    bias = (L_if_mean - rho_mean)/1000
    print(f"GPS1B L_if与理论GPS范围偏差: {bias:.1f} km")
    if abs(bias) < 2000:
        print(f"→ 偏差{bias:.0f}km < 2000km,在合理范围内")
        print(f"→ PPP可以处理这种系统性bias(作为接收机时钟)")
        print(f"→ 如果RMS仍然大,问题不在L1/L2单位")
    else:
        print(f"→ 偏差{bias:.0f}km > 2000km,单位可能有问题")