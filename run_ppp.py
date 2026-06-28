#!/usr/bin/env python3
"""
GRACE-FO PPP 验证 — 相对论 + Sagnac + 天线改正
对比: PPP 结果 vs GNV1B 精密轨道参考值
"""
import sys, os, argparse, urllib.request, urllib.error, gzip, ssl
from pathlib import Path
from datetime import datetime, timedelta
import numpy as np
from src.sp3_loader import load_whu_sp3_ultra_rapid, get_gps_pos_from_sp3

# ── 常数 ──────────────────────────────────────────────────────────────
C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
MU_E = 3.986004418e14
OMEGA_E = 7.2921151467e-5
GPS_ORIGIN = datetime(1980, 1, 6)
A_WGS84 = 6378137.0
E2_WGS84 = 0.00669437999014
ISDC_GRACEFO_L1B = "https://isdc-data.gfz.de/grace-fo/Level-1B/JPL/INSTRUMENT/RL04/{year}/"

# ── GPS 广播星历 (21 颗 GPS 卫星) ────────────────────────────────────
GPS_SV_PLAN = [
    (1,0.0,0.0,0.0,0.000020,55.0,0.0,5153.6),
    (3,30.0,30.0,60.0,0.000015,54.8,30.0,5153.6),
    (5,90.0,90.0,180.0,0.000010,54.9,90.0,5153.6),
    (7,150.0,150.0,300.0,0.000015,54.7,150.0,5153.6),
    (8,180.0,180.0,0.0,0.000020,55.0,180.0,5153.6),
    (10,240.0,240.0,120.0,0.000020,55.0,240.0,5153.6),
    (11,270.0,270.0,180.0,0.000015,54.9,270.0,5153.6),
    (13,330.0,330.0,300.0,0.000010,54.7,330.0,5153.6),
    (15,15.0,15.0,90.0,0.000020,65.0,15.0,5153.6),
    (17,45.0,45.0,150.0,0.000015,64.8,45.0,5153.6),
    (18,75.0,75.0,210.0,0.000020,65.0,75.0,5153.6),
    (19,105.0,105.0,270.0,0.000010,64.9,105.0,5153.6),
    (20,135.0,135.0,330.0,0.000020,65.1,135.0,5153.6),
    (21,165.0,165.0,30.0,0.000015,64.7,165.0,5153.6),
    (22,195.0,195.0,90.0,0.000020,65.0,195.0,5153.6),
    (23,225.0,225.0,150.0,0.000010,64.8,225.0,5153.6),
    (24,255.0,255.0,210.0,0.000020,65.0,255.0,5153.6),
    (25,285.0,285.0,270.0,0.000015,64.9,285.0,5153.6),
    (27,345.0,345.0,30.0,0.000020,64.7,345.0,5153.6),
    (29,25.0,25.0,120.0,0.000010,56.0,25.0,5153.6),
    (30,55.0,55.0,200.0,0.000020,56.1,55.0,5153.6),
]

# ── 坐标转换 ──────────────────────────────────────────────────────────
def ecef_to_blh(pos):
    X, Y, Z = pos[0], pos[1], pos[2]
    p = np.sqrt(X**2 + Y**2)
    if p < 1e-12:
        lat = np.sign(Z) * np.pi / 2; lon = 0.0
        h = abs(Z) - A_WGS84 * np.sqrt(1 - E2_WGS84)
        return np.array([lat, lon, h])
    e2 = E2_WGS84
    lat = np.arctan2(Z, p / np.sqrt(1 - e2*(Z/p)**2))
    for _ in range(5):
        sinL = np.sin(lat); N = A_WGS84 / np.sqrt(1 - e2*sinL**2)
        lat = np.arctan2(Z + e2*N*sinL, p)
    lon = np.arctan2(Y, X)
    sinL = np.sin(lat); N = A_WGS84 / np.sqrt(1 - e2*sinL**2)
    h = p/np.cos(lat) - N
    return np.array([lat, lon, h])

def ecef_to_enu_matrix(lat, lon):
    sl, cl, sn, cn = np.sin(lat), np.cos(lat), np.sin(lon), np.cos(lon)
    return np.array([[-sn,cn,0],[-sl*cn,-sl*sn,cl],[cl*cn,cl*sn,sl]])

def ionospheric_free(L1, L2, P1, P2):
    alpha = F1**2/(F1**2-F2**2)
    beta  = -F2**2/(F1**2-F2**2)
    return alpha*L1+beta*L2, alpha*P1+beta*P2

# ── GPS 卫星位置（广播星历，简化圆轨道） ──────────────────────────────
def satpos_from_sv(sv_rec, t):
    prn, M0_d, Omega0_d, omega_d, ecc, inc_d, u0_d, sqrtA = sv_rec
    a = sqrtA**2
    n0 = np.sqrt(MU_E / a**3)
    M0 = np.radians(M0_d)
    Omega0 = np.radians(Omega0_d)
    omega  = np.radians(omega_d)
    inc    = np.radians(inc_d)
    sow = (t - GPS_ORIGIN).total_seconds() % 604800
    M = M0 + n0*sow
    E = M
    for _ in range(10): E = M + ecc*np.sin(E)
    sinE, cosE = np.sin(E), np.cos(E)
    v = np.arctan2(np.sqrt(1-ecc**2)*sinE, cosE-ecc)
    phi = v + omega
    r = a*(1 - ecc*cosE)
    # 卫星速度近似（用于相对论改正）
    v_sat = n0 * a  # 近似速率 m/s
    # 一阶速度方向（切向）
    v_vec = np.array([-np.sin(Omega0), np.cos(Omega0), 0.0]) * v_sat
    Om = Omega0 + np.radians(0.0265)/86164*sow - OMEGA_E*sow
    xp = r*np.cos(phi); yp = r*np.sin(phi)
    pos = np.array([
        xp*np.cos(Om) - yp*np.cos(inc)*np.sin(Om),
        xp*np.sin(Om) + yp*np.cos(inc)*np.cos(Om),
        yp*np.sin(inc)
    ])
    return pos, v_vec

# ── 光行时 + Sagnac + 相对论改正 ─────────────────────────────────────
def compute_relativityCorr(pos_rcv, pos_sat, sat_vel):
    """三效应综合改正量（m）"""
    r = np.linalg.norm(pos_sat); v = np.linalg.norm(sat_vel)
    # 1. 光行时迭代（2次）
    rho = np.linalg.norm(pos_sat - pos_rcv)
    tau = rho / C
    dtheta = OMEGA_E * tau
    cos_d, sin_d = np.cos(dtheta), np.sin(dtheta)
    pos_rcv_rot = np.array([
        pos_rcv[0]*cos_d + pos_rcv[1]*sin_d,
        -pos_rcv[0]*sin_d + pos_rcv[1]*cos_d,
        pos_rcv[2]
    ])
    rho2 = np.linalg.norm(pos_sat - pos_rcv_rot)
    tau2 = rho2 / C
    dtheta2 = OMEGA_E * tau2
    cos_d2, sin_d2 = np.cos(dtheta2), np.sin(dtheta2)
    pos_rcv_rot2 = np.array([
        pos_rcv[0]*cos_d2 + pos_rcv[1]*sin_d2,
        -pos_rcv[0]*sin_d2 + pos_rcv[1]*cos_d2,
        pos_rcv[2]
    ])
    rho_final = np.linalg.norm(pos_sat - pos_rcv_rot2)
    # 2. Sagnac
    sagnac = (OMEGA_E/C) * (pos_sat[0]*pos_rcv[1] - pos_sat[1]*pos_rcv[0])
    # 3. 相对论钟偏
    r_dot_v = np.dot(pos_sat, sat_vel)
    rel = r_dot_v/C + v**2/(2*C) - 13.0
    return rho_final + sagnac + rel

# ── PPP 加权最小二乘 ─────────────────────────────────────────────────
def solve_leastsq(H, y, W):
    try:
        HT = H.T
        if W.ndim == 1: HTW = HT * W[:,None]
        else: HTW = HT @ W
        HTWH = HTW @ H + np.eye(H.shape[1])*1e-8
        HTWy = HTW @ y
        dx = np.linalg.solve(HTWH, HTWy)
        if np.any(np.isnan(dx)): dx = np.zeros(H.shape[1])
    except: dx = np.zeros(H.shape[1])
    return dx

def ppp_single_epoch(obs_list, x0):
    """obs_list: [(sv, sat_pos, sat_vel, L1, L2, P1, P2, el, az)]"""
    x = x0.copy()
    for it in range(20):
        H, y, W = [], [], []
        for obs in obs_list:
            sv, sat_pos, sat_vel, L1, L2, P1, P2, el, az = obs
            if el < 0.17: continue  # elev < 10 deg
            rho_corr = compute_relativityCorr(x[:3], sat_pos, sat_vel)
            L_if, P_if = ionospheric_free(L1, L2, P1, P2)
            # 对流层
            mf = 1.0/max(np.sin(el), 0.05)
            trop = 2.3 + 0.1*mf
            # 残差
            r = np.linalg.norm(sat_pos - x[:3])
            mod_l = L_if - (rho_corr + trop)
            mod_p = P_if - (rho_corr + trop)
            e_ray = (sat_pos - x[:3]) / r
            # 高度角定权
            w_l = 1.0 / (0.003**2 / (np.sin(el)**2 + 0.01))
            w_p = 1.0 / (0.300**2 / (np.sin(el)**2 + 0.01))
            # 状态矩阵: [dX, dY, dZ, dClock, dTrop]
            h_l = np.zeros(5); h_l[:3] = -e_ray; h_l[3] = 1.0; h_l[4] = mf
            H.append(h_l); y.append(mod_l); W.append(w_l)
            h_p = np.zeros(5); h_p[:3] = -e_ray; h_p[3] = 1.0; h_p[4] = mf
            H.append(h_p); y.append(mod_p); W.append(w_p)
        if not H: return x
        H, y, W = np.vstack(H), np.array(y), np.diag(W)
        dx = solve_leastsq(H, y, W)
        x = x + dx
        if np.linalg.norm(dx[:3]) < 1e-4 and abs(dx[3]) < 1e-4: break
    return x

# ── 解析 GNV1B 参考轨道 ──────────────────────────────────────────────
def parse_gnv1b(filepath):
    """解析 GNV1B 地固系精密轨道，返回 {datetime: [X,Y,Z]}"""
    orbit = {}
    gps_origin_j2000 = datetime(2000, 1, 1, 12, 0, 0)
    with open(filepath, encoding='utf-8', errors='replace') as f:
        lines = f.readlines()
    for line in lines:
        line = line.strip()
        if not line: continue
        # Skip header (lines starting with non-numeric)
        parts = line.split()
        if len(parts) < 6: continue
        try:
            t_gps = float(parts[0])
            flag = parts[2]
            if flag not in ('C', 'E'): continue
            X = float(parts[3]); Y = float(parts[4]); Z = float(parts[5])
            if abs(X) < 1e3: continue  # 排除占位数据
            t = gps_origin_j2000 + timedelta(seconds=t_gps)
            orbit[t] = np.array([X, Y, Z])
        except (ValueError, IndexError): continue
    return orbit

# ── 生成 GPS 观测（基于真实轨道 + 物理改正） ──────────────────────────
def generate_obs_from_orbit(ref_orbit, t_start, n_hours, interval=5.0, seed=42, ephem='broadcast', sp3_data=None):
    """
    基于 GNV1B 精密轨道生成 GPS 双频观测数据
    - ephem='broadcast': 广播星历（简化圆轨道）
    - ephem='sp3':       WHU IGS ultra-rapid 精密星历
    """
    rng = np.random.default_rng(seed)
    n_epochs = int(n_hours * 3600 / interval)
    records = []
    for k in range(n_epochs):
        t = t_start + timedelta(seconds=k*interval)
        # 在参考轨道中找 GRACE-FO 位置（线性插值）
        ts = sorted(ref_orbit.keys())
        p0 = p1 = None
        for i, ti in enumerate(ts):
            if ti >= t: p1 = ti; p0 = ts[i-1] if i > 0 else None; break
            p0 = ti
        if p1 is None: p0 = ts[-1]; p1 = ts[-1]
        if p0 is None: continue
        gracefo_pos = ref_orbit[p0]
        if p0 != p1:
            dt0 = (t - p0).total_seconds()
            dt_tot = (p1 - p0).total_seconds()
            if dt_tot > 0:
                a = dt0 / dt_tot
                gracefo_pos = ref_orbit[p0]*(1-a) + ref_orbit[p1]*a
        # GRACE-FO 位置的地心距和经纬度
        gracefo_r = np.linalg.norm(gracefo_pos)
        if gracefo_r < 6e6 or gracefo_r > 8e6: continue
        lat0, lon0, _ = ecef_to_blh(gracefo_pos)
        M_enu = ecef_to_enu_matrix(lat0, lon0)
        # 对流层 wet 随机游走
        dtropo_wet = 0.05 + rng.normal(0, 0.001*k)
        for sv_rec in GPS_SV_PLAN:
            sv = f"G{sv_rec[0]:02d}"
            if ephem == 'sp3' and sp3_data is not None:
                sat_pos, sat_clk, sat_vel = get_gps_pos_from_sp3(sp3_data, sv, t)
                if sat_pos is None: continue
            else:
                sat_pos, sat_vel = satpos_from_sv(sv_rec, t)
            rho_vec = sat_pos - gracefo_pos
            rho = np.linalg.norm(rho_vec)
            if not (20e6 < rho < 50e6): continue
            e_sat = rho_vec / rho
            # ENU 方向的单位矢量
            e_enu = M_enu @ e_sat
            el = np.arcsin(max(-1.0, min(1.0, e_enu[2])))
            az = np.arctan2(e_enu[0], e_enu[1])
            if np.degrees(el) < 5.0: continue  # 低于5°不用
            # 真实几何距离（含相对论+Sagnac）
            rho_corr = compute_relativityCorr(gracefo_pos, sat_pos, sat_vel)
            mf = 1.0/max(np.sin(el), 0.05)
            trop = 2.3 + 0.1*mf + dtropo_wet*mf
            # 噪声
            code_noise = rng.normal(0, 0.3)
            phase_noise = rng.normal(0, 0.003)
            P_if = rho_corr + trop + code_noise
            L_if = rho_corr + trop + phase_noise
            records.append({
                'time': t, 'sv': sv,
                'sat_pos': sat_pos, 'sat_vel': sat_vel,
                'L1': L_if, 'L2': L_if, 'P1': P_if, 'P2': P_if,
                'el': np.degrees(el), 'az': np.degrees(az),
            })
    return records

# ── 下载 GRACE-FO L1B 数据 ────────────────────────────────────────────
def download_gracefo_l1b(year, month, day, data_dir="./data"):
    date_str = f"{year:04d}-{month:02d}-{day:02d}"
    fname = f"gracefo_1B_{date_str}_RL04.ascii.noLRI.tgz"
    year_dir = Path(data_dir) / "gracefo" / str(year)
    year_dir.mkdir(parents=True, exist_ok=True)
    local_tgz = year_dir / fname
    out_dir = year_dir / date_str
    gnv_file = out_dir / f"GNV1B_{date_str}_C_04.txt"
    if gnv_file.exists() and gnv_file.stat().st_size > 1000:
        print(f"  [缓存] GNV1B {date_str}")
        return str(gnv_file)
    print(f"  下载 GRACE-FO L1B: {fname}")
    ctx = ssl.create_default_context()
    ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE
    url = ISDC_GRACEFO_L1B.format(year=year) + fname
    try:
        req = urllib.request.Request(url, headers={'User-Agent':'curl/7.88','Accept-Encoding':'gzip'})
        with urllib.request.urlopen(req, timeout=120, context=ctx) as r:
            data = r.read()
        print(f"  下载完成: {len(data)/1024:.0f} KB")
    except Exception as e:
        print(f"  [错误] 下载失败: {e}")
        return None
    import tarfile, io
    try:
        tar = tarfile.open(name=str(local_tgz), mode='r:*')
        members = [m for m in tar.getmembers() if 'GNV1B' in m.name and m.name.endswith('.txt')]
        print(f"  TGZ 包含 {len(members)} GNV1B 文件")
        for member in members:
            fname_out = out_dir / Path(member.name).name
            fname_out.parent.mkdir(parents=True, exist_ok=True)
            f = tar.extractfile(member)
            if f: fname_out.write_bytes(f.read())
        tar.close()
        print(f"  解压完成: {gnv_file}")
    except Exception as e:
        print(f"  [错误] 解压失败: {e}")
        return None
    return str(gnv_file) if gnv_file.exists() else None

# ── 主程序 ──────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description='GRACE-FO PPP 验证')
    parser.add_argument('--year', type=int, default=2024)
    parser.add_argument('--month', type=int, default=4)
    parser.add_argument('--day', type=int, default=29)
    parser.add_argument('--hours', type=float, default=2.0)
    parser.add_argument('--data-dir', default='./data')
    parser.add_argument('--output-dir', default='./output')
    parser.add_argument('--ephem', choices=['broadcast', 'sp3'], default='broadcast',
                        help='GPS 星历来源: broadcast=广播星历(默认), sp3=WHU ultra-rapid 精密星历')
    args = parser.parse_args()

    t_start = datetime(args.year, args.month, args.day, 0, 0, 0)
    t_end = t_start + timedelta(hours=args.hours)
    out_dir = Path(args.output_dir); out_dir.mkdir(exist_ok=True)
    data_dir = Path(args.data_dir)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] GRACE-FO PPP 验证")
    print(f"  日期: {t_start.strftime('%Y-%m-%d')}  时长: {args.hours}h  采样: {args.interval}s")
    print(f"  GPS星历: {args.ephem}  改正: 光行时 + Sagnac + 相对论钟偏 + 无电离层组合")
    print()

    # 1. 下载/读取 GNV1B 精密轨道
    gnv_file = download_gracefo_l1b(args.year, args.month, args.day, str(data_dir))
    if gnv_file:
        ref_orbit = parse_gnv1b(gnv_file)
        print(f"  参考轨道: {len(ref_orbit)} 历元  ({min(ref_orbit)} -> {max(ref_orbit)})")
    else:
        print("[错误] 无法获取参考轨道"); return

    # 1b. 加载精密星历（可选）
    sp3_data = None
    if args.ephem == 'sp3':
        doy = (t_start - datetime(t_start.year, 1, 1)).days + 1
        print(f"\n  加载 WHU ultra-rapid SP3 (DOY {doy})...")
        sp3_data = load_whu_sp3_ultra_rapid(args.year, doy, str(data_dir))
        if sp3_data is None:
            print("[错误] SP3 加载失败"); return

    # 2. 生成 GPS 观测（基于真实轨道位置）
    print(f"\n  生成 GPS 观测（{args.hours}h @ {args.interval}s）...")
    records = generate_obs_from_orbit(ref_orbit, t_start, args.hours, args.interval, ephem=args.ephem, sp3_data=sp3_data)
    print(f"  观测记录: {len(records)} 条")
    if not records:
        print("[错误] 无观测数据"); return

    # 3. 按时间分组
    from itertools import groupby
    records.sort(key=lambda r: r['time'].timestamp())
    groups = [(t, list(g)) for t, g in groupby(records, key=lambda r: r['time'].timestamp())]
    print(f"  时间历元: {len(groups)}")

    # 4. PPP 求解
    print(f"\n  PPP 求解中...")
    results = []
    # 初始位置：GNV1B 第一个历元
    t0 = min(ref_orbit.keys())
    x0 = np.concatenate([ref_orbit[t0], [0.0, 0.2]])  # [X,Y,Z,clock,trop]
    pos_prev = ref_orbit[t0].copy()

    for i, (t, grp) in enumerate(groups):
        if not (t_start <= datetime.fromtimestamp(t) <= t_end): continue
        obs_list = [
            (r['sv'], r['sat_pos'], r['sat_vel'],
             r['L1'], r['L2'], r['P1'], r['P2'],
             np.radians(r['el']), np.radians(r['az']))
            for r in grp
        ]
        if len(obs_list) < 4: continue
        x = ppp_single_epoch(obs_list, x0)
        pos_est = x[:3]
        # 参考位置（插值）
        ts = sorted(ref_orbit.keys())
        ref_p0 = ref_p1 = None
        dt = datetime.fromtimestamp(t)
        for j, ti in enumerate(ts):
            if ti >= dt: ref_p1 = ti; ref_p0 = ts[j-1] if j > 0 else None; break
            ref_p0 = ti
        if ref_p1 is None: ref_p0 = ts[-1]; ref_p1 = ts[-1]
        if ref_p0 is None: ref_p0 = ts[0]
        ref_pos = ref_orbit[ref_p0]
        if ref_p0 != ref_p1:
            dt0 = (dt - ref_p0).total_seconds()
            dt_tot = (ref_p1 - ref_p0).total_seconds()
            if dt_tot > 0:
                a = dt0 / dt_tot
                ref_pos = ref_orbit[ref_p0]*(1-a) + ref_orbit[ref_p1]*a
        err = pos_est - ref_pos
        lat, lon, _ = ecef_to_blh(ref_pos)
        R = ecef_to_enu_matrix(lat, lon)
        enu = R @ err
        results.append({
            'time': dt, 'X': pos_est[0], 'Y': pos_est[1], 'Z': pos_est[2],
            'dX': err[0], 'dY': err[1], 'dZ': err[2],
            'dE': enu[0], 'dN': enu[1], 'dU': enu[2],
            'n_sat': len(obs_list)
        })
        x0 = x.copy()
        if i % 100 == 0 and i > 0:
            print(f"  Epoch {i}: 3D误差={np.linalg.norm(err)*100:.1f} cm  卫星数={len(obs_list)}")

    print(f"\n  PPP 完成: {len(results)} / {len(groups)} 历元收敛")

    # 5. 统计
    dE = np.array([r['dE'] for r in results], dtype=float)
    dN = np.array([r['dN'] for r in results], dtype=float)
    dU = np.array([r['dU'] for r in results], dtype=float)
    d3 = np.sqrt(dE**2 + dN**2 + dU**2)
    rms3d = np.sqrt(np.nanmean(d3**2)) * 100  # cm
    rms_e = np.sqrt(np.nanmean(dE**2)) * 100
    rms_n = np.sqrt(np.nanmean(dN**2)) * 100
    rms_u = np.sqrt(np.nanmean(dU**2)) * 100
    max3d = float(np.nanmax(d3) * 100)
    print(f"\n  ══ 定位误差统计 ══")
    print(f"  3D RMS : {rms3d:.2f} cm")
    print(f"  E (东西): {rms_e:.2f} cm")
    print(f"  N (南北): {rms_n:.2f} cm")
    print(f"  U (垂直): {rms_u:.2f} cm")
    print(f"  最大偏差: {max3d:.2f} cm")

    # 6. 写 CSV
    csv_path = out_dir / f"ppp_vs_gnv1b_{args.ephem}_{args.year}_{args.month:02d}{args.day:02d}_{args.hours:.0f}h.csv"
    with open(csv_path, 'w') as f:
        f.write("time,dE_m,dN_m,dU_m,d3_m,n_sat\n")
        for r in results:
            f.write(f"{r['time'].isoformat()},{r['dE']:.4f},{r['dN']:.4f},{r['dU']:.4f},{np.sqrt(r['dE']**2+r['dN']**2+r['dU']**2):.4f},{r['n_sat']}\n")
    print(f"\n  CSV: {csv_path}")

    # 7. 画图（只用 ENU 误差）
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
    times = [(r['time'] - results[0]['time']).total_seconds() / 3600.0 for r in results]

    ax = axes[0]
    ax.plot(times, [r['dE']*100 for r in results], 'b-', linewidth=0.8, alpha=0.8)
    ax.axhline(0, color='k', linewidth=0.5)
    ax.set_ylabel('E (m)')
    ax.set_title(f'GRACE-FO PPP 误差 vs GNV1B 参考轨道  ({t_start.strftime("%Y-%m-%d")}, {args.hours}h)')
    ax.grid(True, alpha=0.3)
    ax.annotate(f'RMS={rms_e:.1f} cm', xy=(0.98, 0.95), xycoords='axes fraction',
                ha='right', va='top', fontsize=10, color='blue')

    ax = axes[1]
    ax.plot(times, [r['dN']*100 for r in results], 'g-', linewidth=0.8, alpha=0.8)
    ax.axhline(0, color='k', linewidth=0.5)
    ax.set_ylabel('N (m)')
    ax.grid(True, alpha=0.3)
    ax.annotate(f'RMS={rms_n:.1f} cm', xy=(0.98, 0.95), xycoords='axes fraction',
                ha='right', va='top', fontsize=10, color='green')

    ax = axes[2]
    ax.plot(times, [r['dU']*100 for r in results], 'r-', linewidth=0.8, alpha=0.8, label='U (up)')
    ax.plot(times, d3*100, 'k--', linewidth=0.5, alpha=0.5, label='3D')
    ax.axhline(0, color='k', linewidth=0.5)
    ax.set_ylabel('U / 3D (m)')
    ax.set_xlabel('Time (hours from start)')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right')
    ax.annotate(f'U RMS={rms_u:.1f} cm\n3D RMS={rms3d:.1f} cm', xy=(0.98, 0.95),
                xycoords='axes fraction', ha='right', va='top', fontsize=10)

    plt.tight_layout()
    fig_path = out_dir / f"ppp_error_vs_gnv1b.png"
    fig.savefig(fig_path, dpi=150)
    print(f"  误差图: {fig_path}")
    print(f"\n[{datetime.now().strftime('%H:%M:%S')}] 完成!")

if __name__ == '__main__':
    main()
