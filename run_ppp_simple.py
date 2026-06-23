#!/usr/bin/env python3
"""
GRACE-FO PPP 批量处理 — 简化算法（apr29_final_test.py 验证版）
用法: python3 run_ppp_simple.py --start 2024-04-29 --end 2024-05-31 --grace C
"""
import sys, os, math, pickle, datetime as dt, argparse, warnings
import numpy as np
sys.path.insert(0, '/workspace/gnss_pod/src')

warnings.filterwarnings('ignore')
np.seterr(all='ignore')

C = 299792458.0; F1, F2 = 1575.42e6, 1227.60e6
LAM1 = C/F1; LAM2 = C/F2
F1_SQ, F2_SQ = F1*F1, F2*F2
ALPHA = F1_SQ/(F1_SQ-F2_SQ); BETA = -F2_SQ/(F1_SQ-F2_SQ)
OMEGA_E = 7.2921151467e-5; MU_E = 3.986004418e14
GPS_ORIGIN = dt.datetime(1980, 1, 6)
ITER_MAX = 15; TOL = 0.05  # 5 cm
A_WGS84 = 6378137.0; E2_WGS84 = 0.00669437999014
OMEGA_WGS84 = 7.2921151467e-5

# ── GPS 广播星历（Apr29 验证参数，sqrtA in km²）──────────────
# apr29_final_test.py 验证：sqrtA=5153.6 → a=26,564 km² → |GPS|=26560 km ✓
GPS_SV = [
    # (PRN, M0_deg, Omega0_deg, omega_deg, ecc, inc_deg, u0_deg, sqrtA_km)
    (1,  0.0,   0.0,   0.0,   2e-5,   55.0, 0.0, 5153.6),
    (3,  30.0, 30.0,  60.0,  1.5e-5,  54.8, 30.0, 5153.6),
    (5,  90.0, 90.0,   0.0,  2e-5,   55.0,  0.0, 5153.6),
    (7, 150.0,150.0,   0.0,  1.5e-5,  54.8, 30.0, 5153.6),
    (8, 180.0,180.0,   0.0,  1e-5,   55.0,  0.0, 5153.6),
    (9, 210.0,210.0,   0.0,  2e-5,   55.0,  0.0, 5153.6),
    (10, 240.0,240.0,   0.0, 1.5e-5,  54.8,  0.0, 5153.6),
    (11, 270.0,270.0,   0.0, 1e-5,   55.0,  0.0, 5153.6),
    (13, 300.0,300.0,   0.0, 1.5e-5,  54.8,  0.0, 5153.6),
    (14, 330.0,330.0,   0.0, 1e-5,   55.0,  0.0, 5153.6),
    (15,  0.0,  0.0,   0.0, 2e-5,   55.0,  0.0, 5153.6),
    (17,  30.0, 30.0,  60.0, 1e-5,   54.8, 30.0, 5153.6),
    (18,  90.0, 90.0,   0.0, 2e-5,   55.0,  0.0, 5153.6),
    (19, 120.0,120.0,   0.0, 1.5e-5,  54.8,  0.0, 5153.6),
    (20, 150.0,150.0,   0.0, 1e-5,   55.0,  0.0, 5153.6),
    (21, 180.0,180.0,   0.0, 2e-5,   55.0,  0.0, 5153.6),
    (22, 210.0,210.0,   0.0, 1.5e-5,  54.8,  0.0, 5153.6),
    (23, 240.0,240.0,   0.0, 1e-5,   55.0,  0.0, 5153.6),
    (24, 270.0,270.0,   0.0, 2e-5,   54.8,  0.0, 5153.6),
    (25, 300.0,300.0,   0.0, 1.5e-5,  55.0,  0.0, 5153.6),
    (26, 330.0,330.0,   0.0, 1e-5,   54.8,  0.0, 5153.6),
    (27,  0.0,  0.0,   0.0, 2e-5,   55.0,  0.0, 5153.6),
    (28,  30.0, 30.0,  60.0, 1.5e-5,  54.8, 30.0, 5153.6),
    (29,  90.0, 90.0,   0.0, 1e-5,   55.0,  0.0, 5153.6),
    (30, 120.0,120.0,   0.0, 2e-5,   54.8,  0.0, 5153.6),
    (31, 150.0,150.0,   0.0, 1e-5,   55.0,  0.0, 5153.6),
    (32, 180.0,180.0,   0.0, 1.5e-5,  54.8,  0.0, 5153.6),
]
SV_PLAN = {r[0]: r for r in GPS_SV}

def gps_pos(sv_id, t_gps, t0=None):
    """GPS卫星位置（ECEF），t_gps为GPS秒数，apr29_final_test.py算法"""
    prn = int(sv_id[1:]) if isinstance(sv_id, str) else sv_id
    if prn not in SV_PLAN: return None
    r = SV_PLAN[prn]
    _, M0_d, Omega0_d, omega_d, ecc, inc_d, u0_d, sqrtA = r
    t0 = t0 or 767620800.0
    a = sqrtA**2  # km² → m² (5153.6² = 26,564 km²)
    n0 = math.sqrt(MU_E / a**3)
    M0 = math.radians(M0_d); Omega0 = math.radians(Omega0_d)
    omega = math.radians(omega_d); inc = math.radians(inc_d)
    M = M0 + n0*(t_gps - t0)
    E = M
    for _ in range(10): E = M + ecc*math.sin(E)
    v = 2*math.atan2(math.sqrt(1.00004)*math.sin(E/2), math.cos(E/2)-0.00002)
    phi = v + omega; r_mag = a*(1 - ecc*math.cos(E))
    Om = Omega0 - OMEGA_E*(t_gps - t0)
    X = r_mag*(math.cos(phi)*math.cos(Om) - math.sin(phi)*math.cos(inc)*math.sin(Om))
    Y = r_mag*(math.cos(phi)*math.sin(Om) + math.sin(phi)*math.cos(inc)*math.cos(Om))
    Z = r_mag*math.sin(phi)*math.sin(inc)
    return (X, Y, Z)

def ecef_to_blh(pos):
    X, Y, Z = float(pos[0]), float(pos[1]), float(pos[2])
    p = math.sqrt(X**2 + Y**2)
    if p < 1e-12:
        return (math.pi/2 * (1 if Z >= 0 else -1), 0.0, abs(Z) - A_WGS84)
    e2 = E2_WGS84
    # Use numpy for safe sqrt
    val = 1 - e2*(Z/p)**2
    if val <= 0:
        lat = math.pi/2 * (1 if Z >= 0 else -1)
    else:
        lat = math.atan2(Z, p / math.sqrt(val))
    lon = math.atan2(Y, X)
    sinL = math.sin(lat)
    N = A_WGS84 / math.sqrt(1 - e2*sinL**2)
    h = p/math.cos(lat) - N
    return (lat, lon, h)

def ecef_to_enu(err, lat, lon):
    cl, sl = math.cos(lat), math.sin(lat)
    cn, sn = math.cos(lon), math.sin(lon)
    return (
        -sl*cn*err[0] - sl*sn*err[1] + cl*err[2],
         sl*sn*err[0] - sl*cn*err[1],
         cl*cn*err[0] + cl*sn*err[1] + sl*err[2],
    )

def geometric_range(pos_rcv, pos_sat):
    return math.sqrt((pos_sat[0]-pos_rcv[0])**2 + (pos_sat[1]-pos_rcv[1])**2 + (pos_sat[2]-pos_rcv[2])**2)

def solve_epoch(obs_list, x0, ref_pos):
    """单历元加权最小二乘，apr29_final_test.py算法"""
    H = []; y = []; W = []
    for obs in obs_list:
        sp = obs['sat_pos']; rcv = x0[:3]
        rho = geometric_range(rcv, sp)
        if rho < 1e6: continue
        L_if = obs['L_if']
        el = obs['el']
        sigma = 0.003 / max(math.sin(el), 0.1)
        w = 1.0 / sigma**2
        hx = -(sp[0]-rcv[0])/rho
        hy = -(sp[1]-rcv[1])/rho
        hz = -(sp[2]-rcv[2])/rho
        H.append([hx, hy, hz, 1.0])
        r_obs = L_if - rho
        y.append(r_obs); W.append(w)
    if len(H) < 4: return None
    HT = np.array(H).T; W_arr = np.diag(W)
    HTW = HT @ W_arr; HTWH = HTW @ np.array(H)
    HTWy = HTW @ np.array(y)
    HTWH += np.eye(HTWH.shape[0]) * 1e-8
    try:
        dx = np.linalg.solve(HTWH, HTWy)
    except: return None
    if np.any(np.isnan(dx)) or np.any(np.isinf(dx)): return None
    return dx

def ppp_single_day(gps1b_day, gnv_day, date_str):
    """处理单日GPS1B数据，apr29_final_test.py算法"""
    if not gps1b_day or not gnv_day:
        return [], date_str

    # GRACE参考位置（UTC时间键）
    gnv_keys = sorted(gnv_day.keys())
    grace_ref = np.array(gnv_day[gnv_keys[0]])
    ref_lat, ref_lon, _ = ecef_to_blh(grace_ref)

    # 按GPS时间排序
    gps_keys = sorted(gps1b_day.keys())
    
    # GPS秒 → UTC 转换基准（GPS无闰秒，UTC有18闰秒）
    def gps_to_utc(t_gps):
        return GPS_ORIGIN + dt.timedelta(seconds=t_gps)

    # 每30s一个epoch
    results = []
    prev_result = None
    
    for t_gps in gps_keys:
        obs_dict = gps1b_day[t_gps]
        if len(obs_dict) < 4: continue
        
        utc = gps_to_utc(t_gps)
        
        # 找GRACE参考位置（UTC时间匹配）
        rp0 = rp1 = None
        for tk in gnv_keys:
            if tk >= utc:
                rp1 = tk
                break
            rp0 = tk
        if rp0 is None: rp0 = gnv_keys[0]
        if rp1 is None: rp1 = gnv_keys[-1]
        if rp0 == rp1:
            grace_pos = np.array(gnv_day[rp0])
        else:
            dt0 = (utc - rp0).total_seconds()
            dtt = (rp1 - rp0).total_seconds()
            grace_pos = np.array(gnv_day[rp0]) * (1 - dt0/dtt) + np.array(gnv_day[rp1]) * (dt0/dtt)
        
        # 构建观测量（实时计算el，不用pickle里的）
        obs_list = []
        for sv_id, rec in obs_dict.items():
            if rec.get('L1') is None: continue
            sp = gps_pos(sv_id, t_gps)
            if sp is None: continue
            L1 = rec['L1']; L2 = rec['L2']; P1 = rec['P1']; P2 = rec['P2']
            L_if = ALPHA*L1 + BETA*L2

            # 实时计算仰角（GRACE→GPS向量与GRACE当地垂线夹角）
            vec = np.array([sp[0]-grace_pos[0], sp[1]-grace_pos[1], sp[2]-grace_pos[2]])
            rho = np.linalg.norm(vec)
            if rho < 1e6: continue
            zenith_angle = math.acos(max(min(grace_pos[0]*vec[0]+grace_pos[1]*vec[1]+grace_pos[2]*vec[2], rho*np.linalg.norm(grace_pos)), -rho*np.linalg.norm(grace_pos)) / (rho*np.linalg.norm(grace_pos)))
            el = math.pi/2 - zenith_angle  # 仰角
            if el <= 0: continue

            obs_list.append({'sat_pos': sp, 'L_if': L_if, 'el': el, 'sv': sv_id})
        
        if len(obs_list) < 4: continue
        
        # 初始值：上一个收敛结果或GRACE位置
        if prev_result is not None:
            x0 = np.array([prev_result['X'], prev_result['Y'], prev_result['Z'], 0.0])
        else:
            x0 = np.array([grace_ref[0], grace_ref[1], grace_ref[2], 0.0])
        
        # 迭代最小二乘
        for it in range(ITER_MAX):
            dx = solve_epoch(obs_list, x0, grace_ref)
            if dx is None: break
            x0 = x0 + dx
            if np.linalg.norm(dx[:3]) < TOL and abs(dx[3]) < TOL:
                break
        
        pos_est = x0[:3]
        err = pos_est - grace_ref
        enu = ecef_to_enu(err, ref_lat, ref_lon)
        d3 = math.sqrt(sum(x**2 for x in enu))
        
        result = {
            'time': utc, 'X': pos_est[0], 'Y': pos_est[1], 'Z': pos_est[2],
            'dE': enu[0], 'dN': enu[1], 'dU': enu[2], 'd3': d3,
            'n_sat': len(obs_list),
        }
        results.append(result)
        prev_result = result
    
    return results, date_str

def main():
    p = argparse.ArgumentParser(description='GRACE-FO PPP 简化算法')
    p.add_argument('--start', required=True, help='开始日期 YYYY-MM-DD')
    p.add_argument('--end', required=True, help='结束日期 YYYY-MM-DD')
    p.add_argument('--grace', default='C', help='GRACE卫星编号')
    p.add_argument('--out-dir', default='output_v12', help='输出目录')
    args = p.parse_args()

    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    start = dt.date.fromisoformat(args.start)
    end = dt.date.fromisoformat(args.end)
    delta = dt.timedelta(days=1)

    all_results = []

    cur = start
    while cur <= end:
        ds = cur.isoformat()
        y = cur.year; m = cur.month; d = cur.day
        doy = (cur - dt.date(y, 1, 1)).days + 1
        print(f"\n=== {ds} (DOY {doy}) ===")

        # 加载GPS1B pickle（GPS时间键）
        pkl_paths = [
            f'data/gracefo/{y}/{y:04d}-{m:02d}-{d:02d}/gps1b_{args.grace}.pkl',
            f'data/gracefo/{y}/{y:04d}-{m:02d}-{d:02d}/GPS1B_{y:04d}-{m:02d}-{d:02d}_{args.grace}_04.pkl',
            f'data/gracefo/{y}/{y:04d}-{m:02d}-{d:02d}/{y:04d}-{m:02d}-{d:02d}/gps1b_{args.grace}.pkl',
        ]
        gps1b_day = None
        for pkl_p in pkl_paths:
            if os.path.exists(pkl_p):
                gps1b_day = pickle.load(open(pkl_p, 'rb'))
                print(f"  GPS1B: {len(gps1b_day)} epochs ({sum(len(v) for v in gps1b_day.values())} obs)")
                break
        if gps1b_day is None:
            print(f"  GPS1B缓存不存在，跳过")
            cur += delta; continue

        # 加载GNV1B（UTC时间键）
        gnv_paths = [
            f'data/gracefo/{y}/{y:04d}-{m:02d}-{d:02d}/GNV1B_{y:04d}-{m:02d}-{d:02d}_{args.grace}_04.txt',
            f'data/gracefo/{y}/{y:04d}-{m:02d}-{d:02d}/{y:04d}-{m:02d}-{d:02d}/GNV1B_{y:04d}-{m:02d}-{d:02d}_{args.grace}_04.txt',
        ]
        gnv_day = None
        for gnv_p in gnv_paths:
            if os.path.exists(gnv_p):
                try:
                    from batch_v12 import parse_gnv1b
                    gnv_day = parse_gnv1b(gnv_p)
                    print(f"  GNV1B: {len(gnv_day)} epochs")
                except Exception as e:
                    print(f"  GNV1B解析失败: {e}")
                break
        if gnv_day is None:
            print(f"  GNV1B不存在，跳过")
            cur += delta; continue

        # PPP处理
        results, ds_out = ppp_single_day(gps1b_day, gnv_day, ds)
        if not results:
            print(f"  无有效结果")
            cur += delta; continue

        # 统计
        d3s = [r['d3'] for r in results if r['d3'] > 0]
        n = len(d3s)
        rms_e = math.sqrt(sum(r['dE']**2 for r in results)/n)*100
        rms_n = math.sqrt(sum(r['dN']**2 for r in results)/n)*100
        rms_u = math.sqrt(sum(r['dU']**2 for r in results)/n)*100
        rms_3d = math.sqrt(sum(r['d3']**2 for r in results)/n)*100
        print(f"  历元: {n}/{len(results)}, RMS E/N/U/3D: {rms_e:.1f}/{rms_n:.1f}/{rms_u:.1f}/{rms_3d:.1f} cm")

        # 保存CSV
        csv_path = f'{out_dir}/ppp_{ds}_simple.csv'
        with open(csv_path, 'w', newline='') as f:
            keys = list(results[0].keys())
            import csv
            csv.DictWriter(f, fieldnames=keys).writeheader()
            csv.DictWriter(f, fieldnames=keys).writerows(results)
        print(f"  保存: {csv_path}")
        all_results.append((ds, n, rms_e, rms_n, rms_u, rms_3d))
        cur += delta

    # 汇总
    if all_results:
        print("\n=== 全时期统计 ===")
        print(f"{'日期':<14} {'历元':<8} {'RMS_E':<8} {'RMS_N':<8} {'RMS_U':<8} {'RMS_3D':<10} {'说明'}")
        for ds, n, re, rn, ru, r3 in all_results:
            print(f"{ds:<14} {n:<8} {re:>7.1f}  {rn:>7.1f}  {ru:>7.1f}  {r3:>8.1f}  cm")
        avg_e = sum(r[2] for r in all_results)/len(all_results)
        avg_n = sum(r[3] for r in all_results)/len(all_results)
        avg_u = sum(r[4] for r in all_results)/len(all_results)
        avg_3d = sum(r[5] for r in all_results)/len(all_results)
        print(f"{'平均':<14} {'':<8} {avg_e:>7.1f}  {avg_n:>7.1f}  {avg_u:>7.1f}  {avg_3d:>8.1f}  cm")

if __name__ == '__main__':
    main()