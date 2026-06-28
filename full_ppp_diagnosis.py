#!/usr/bin/env python3
"""
Apr29 GPS1B PPP完整诊断: 读取GPS1B → 运行PPP → 分析残差
一次性搞清楚问题所在
"""
import sys, os, math, pickle
sys.path.insert(0, '/workspace/gnss_pod/src')

C = 299792458.0; F1, F2 = 1575.42e6, 1227.60e6
LAM1 = C/F1; LAM2 = C/F2
F1_SQ, F2_SQ = F1*F1, F2*F2
ALPHA = F1_SQ/(F1_SQ-F2_SQ); BETA = -F2_SQ/(F1_SQ-F2_SQ)

# ============================================================
# 1. 读取GPS1B Apr29数据 (pickle)
# ============================================================
pkl = 'data/gracefo/2024/2024-04-29/gps1b_C.pkl'
gps1b = pickle.load(open(pkl, 'rb'))
gnv_pkl = 'data/gracefo/2024/2024-04-29/gnv1b_C.pkl'
gnv = pickle.load(open(gnv_pkl, 'rb')) if os.path.exists(gnv_pkl) else {}

# 读取第一个epoch的GRACE位置
first_t = sorted(gps1b.keys())[0]
grace_pos = None
if first_t in gnv and gnv[first_t]:
    gp = gnv[first_t][0]
    grace_pos = (gp.get('X',0), gp.get('Y',0), gp.get('Z',0))
    print(f"GRACE位置: X={grace_pos[0]:.1f}, Y={grace_pos[1]:.1f}, Z={grace_pos[2]:.1f} m")
    print(f"GRACE半径: R={math.sqrt(sum(x**2 for x in grace_pos)):.1f} m ({math.sqrt(sum(x**2 for x in grace_pos))/1000:.1f} km)")
else:
    grace_pos = (0, 0, 6861000)
    print("使用近似GRACE位置 R=6861 km")

# ============================================================
# 2. 模拟PPP的测量模型
# ============================================================
# PPP相位观测方程: L_if = rho + c*dt_r - c*dt_s + T + I_if + N_if*LAM_IF
# 其中I_if = 0 (无电离层组合)
# 
# 关键: GPS1B L1_phase 单位是什么?
# GPS1B L1_phase ≈ 21868738
# 期望: GPS卫星到GRACE的距离 ≈ 20265 km = 20265000 m
#
# 如果L1_phase单位是m:
#   L_if ≈ 21868738 m (≈ 21869 km) → 误差 +1600 km (8%)
# 如果L1_phase单位是μs:
#   L1_m = L1_phase * C/1e6 = 21868738 * 299.79 = 6.554e9 m (6554 km) → 误差 -68%
# 如果L1_phase单位是cycles:
#   L_if_cycles = L1_phase ≈ 21868738 cycles → L_if_m = 21868738 * 0.107 ≈ 2339 km → 误差 -88%
#
# 结论: 只有"直接用原始值(假设单位=m)"的误差是8%,最接近正确值
# 但8%的误差(~1600km)仍然很大 → PPP会怎么处理?

# GPS广播星历参数 (简化)
OMEGA_E = 7.2921151467e-5; MU_E = 3.986004418e14
GPS_SMA = 26560e3

def gps_position(sv, t_gps):
    """GPS卫星ECEF位置"""
    sv_params = {
        'G05': {'M0': 0.0, 'Omega0': 0.0, 'omega': 0.0, 'ecc': 0.00002, 'inc': 55.0},
        'G07': {'M0': math.radians(150), 'Omega0': math.radians(150), 'omega': 0.0, 'ecc': 0.000015, 'inc': 55.0},
        'G13': {'M0': math.radians(60), 'Omega0': math.radians(60), 'omega': 0.0, 'ecc': 0.00001, 'inc': 55.0},
        'G18': {'M0': math.radians(250), 'Omega0': math.radians(250), 'omega': 0.0, 'ecc': 0.000015, 'inc': 55.0},
        'G20': {'M0': math.radians(90), 'Omega0': math.radians(90), 'omega': 0.0, 'ecc': 0.00002, 'inc': 55.0},
        'G27': {'M0': math.radians(30), 'Omega0': math.radians(30), 'omega': 0.0, 'ecc': 0.00001, 'inc': 55.0},
        'G30': {'M0': math.radians(270), 'Omega0': math.radians(270), 'omega': 0.0, 'ecc': 0.000015, 'inc': 55.0},
    }
    if sv not in sv_params: return None
    p = sv_params[sv]
    a = GPS_SMA; n0 = math.sqrt(MU_E/a**3)
    M = p['M0'] + n0*(t_gps-767620800)
    E = M
    for _ in range(8): E = M + p['ecc']*math.sin(E)
    v = 2*math.atan2(math.sqrt(1+p['ecc'])*math.sin(E/2), math.cos(E/2)-p['ecc'])
    phi = v + math.radians(p['omega'])
    r = a*(1-p['ecc']*math.cos(E))
    Om = p['Omega0'] - OMEGA_E*(t_gps-767620800)
    inc = math.radians(p['inc'])
    X = r*(math.cos(phi)*math.cos(Om) - math.sin(phi)*math.cos(inc)*math.sin(Om))
    Y = r*(math.cos(phi)*math.sin(Om) + math.sin(phi)*math.cos(inc)*math.cos(Om))
    Z = r*math.sin(phi)*math.sin(inc)
    return (X, Y, Z)

# ============================================================
# 3. 分析第一个epoch的测量残差
# ============================================================
print("\n=== 分析第一个epoch ===")
obs = gps1b[first_t]
print(f"卫星数: {len(obs)}")

results = []
for sat_id, obs_data in obs.items():
    sv = sat_id
    L1 = obs_data['L1']; L2 = obs_data['L2']
    P1 = obs_data['P1']; P2 = obs_data['P2']
    L_if = obs_data['L_if']; P_if = obs_data['P_if']
    
    # GPS卫星位置
    gp = gps_position(sv, first_t)
    if not gp: continue
    
    # 几何范围 (光行时改正)
    rho0 = math.sqrt(sum((gp[i]-grace_pos[i])**2 for i in range(3)))
    tau = rho0/C
    dth = OMEGA_E * tau
    cr, sr = math.cos(dth), math.sin(dth)
    rcv = (grace_pos[0]*cr+grace_pos[1]*sr, -grace_pos[0]*sr+grace_pos[1]*cr, grace_pos[2])
    rho = math.sqrt(sum((gp[i]-rcv[i])**2 for i in range(3)))
    
    # 残差 (测量值 - 几何范围)
    # 这是PPP中会进入滤波器的残差
    resid_phase = L_if - rho
    resid_code = P_if - rho
    
    results.append({
        'sv': sv,
        'L1': L1, 'L2': L2, 'L_if': L_if, 'P_if': P_if,
        'rho': rho,
        'resid_phase': resid_phase,
        'resid_code': resid_code,
    })

print(f"\n{'SV':5s} {'L_if(m)':14s} {'rho(m)':12s} {'res_phase(m)':14s} {'res_code(m)':12s} {'el':8s}")
print("-" * 75)
for r in sorted(results, key=lambda x: x['sv']):
    el = obs[r['sv']].get('el', 0)
    print(f"{r['sv']:5s} {r['L_if']:14.1f} {r['rho']:12.1f} {r['resid_phase']:14.1f} {r['resid_code']:12.1f} {math.degrees(el):6.1f}°" if el else f"{r['sv']:5s} {r['L_if']:14.1f} {r['rho']:12.1f} {r['resid_phase']:14.1f} {r['resid_code']:12.1f} N/A")

print()
if results:
    phase_bias = sum(r['resid_phase'] for r in results)/len(results)
    code_bias = sum(r['resid_code'] for r in results)/len(results)
    phase_std = math.sqrt(sum((r['resid_phase']-phase_bias)**2 for r in results)/len(results))
    code_std = math.sqrt(sum((r['resid_code']-code_bias)**2 for r in results)/len(results))
    
    print(f"相位残差: mean={phase_bias/1000:.1f}km, std={phase_std/1000:.1f}km ({phase_std:.1f}m)")
    print(f"码残差:   mean={code_bias/1000:.1f}km, std={code_std/1000:.1f}km ({code_std:.1f}m)")
    
    # 检查L_if是否在正确范围内
    L_if_vals = [r['L_if'] for r in results]
    rho_vals = [r['rho'] for r in results]
    
    print(f"\nL_if范围: {min(L_if_vals)/1000:.0f} - {max(L_if_vals)/1000:.0f} km")
    print(f"rho范围:  {min(rho_vals)/1000:.0f} - {max(rho_vals)/1000:.0f} km")
    print(f"期望GPS范围: ~20265 km")
    
    # 关键: L_if是否在GPS范围内?
    L_if_mean = sum(L_if_vals)/len(L_if_vals)
    rho_mean = sum(rho_vals)/len(rho_vals)
    
    print(f"\nL_if平均: {L_if_mean/1000:.1f}km vs rho平均: {rho_mean/1000:.1f}km")
    print(f"差异: {(L_if_mean-rho_mean)/1000:.1f}km = {abs(L_if_mean-rho_mean)/rho_mean*100:.1f}%")
    
    # 检查是否所有残差都是系统性的(相似的bias)
    phase_resids = [r['resid_phase'] for r in results]
    code_resids = [r['resid_code'] for r in results]
    
    print(f"\n相位残差范围: {min(phase_resids)/1000:.1f} - {max(phase_resids)/1000:.1f} km")
    print(f"码残差范围:   {min(code_resids)/1000:.1f} - {max(code_resids)/1000:.1f} km")
    
    if max(phase_resids) - min(phase_resids) < 10:
        print("→ 相位残差系统性bias (所有卫星相似)")
        print(f"→ 接收机时钟偏置 ≈ {phase_bias/1000:.1f} km")
    else:
        print("→ 相位残差非系统性 (卫星间差异大)")

print("\n=== 诊断结论 ===")
if results:
    # 分析L_if的单位
    L_if_m = results[0]['L_if']
    L_if_km = L_if_m / 1000
    print(f"GPS1B L_if = {L_if_m:.0f} m = {L_if_km:.0f} km")
    print(f"期望GPS范围 ≈ 20265 km")
    print(f"误差: {abs(L_if_km-20265):.0f} km = {abs(L_if_km-20265)/20265*100:.1f}%")
    
    if abs(L_if_km - 20265) < 3000:
        print("→ L_if在GPS范围内,误差<15%")
        print("→ 单位正确,可能是GPS时钟/PCV偏差")
    elif abs(L_if_km - 6556) < 3000:
        print("→ L_if ≈ 6556 km → 单位是μs(错误)")
        print("→ 修复: L_if(μs→m) = L_if * C/1e6")
    else:
        print("→ L_if完全不在正确范围内")
        print(f"→ L_if={L_if_km:.0f}km与任何预期值都不匹配")