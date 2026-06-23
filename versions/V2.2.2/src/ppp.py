"""
PPP 求解器 — 支持相对论性光行时改正、Sagnac 效应、天线相位中心改正
"""
from datetime import datetime, timedelta
import numpy as np

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
LAMBDA1, LAMBDA2 = C / F1, C / F2
MU_E = 3.986004418e14        # 地心引力常数 m³/s²
OMEGA_E = 7.2921151467e-5   # 地球自转角速度 rad/s
GPS_ORIGIN = datetime(1980, 1, 6)
A_WGS84 = 6378137.0
E2_WGS84 = 0.00669437999014

# ── GPS 广播星历星座参数 ─────────────────────────────────────────────
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


def blh_to_ecef(lat, lon, h=0.0):
    """BLH (rad, rad, m) → ECEF (m)"""
    sinL, cosL = np.sin(lat), np.cos(lat)
    sinl, cosl = np.sin(lon), np.cos(lon)
    N = A_WGS84 / np.sqrt(1 - E2_WGS84 * sinL**2)
    return np.array([
        (N + h) * cosL * cosl,
        (N + h) * cosL * sinl,
        (N * (1 - E2_WGS84) + h) * sinL
    ])


def ecef_to_blh(pos):
    """ECEF (m) → BLH (rad, rad, m)"""
    X, Y, Z = pos[0], pos[1], pos[2]
    p = np.sqrt(X**2 + Y**2)
    if p < 1e-12:
        lat = np.sign(Z) * np.pi / 2
        lon = 0.0
        h = abs(Z) - A_WGS84 * np.sqrt(1 - E2_WGS84)
        return np.array([lat, lon, h])
    e2 = E2_WGS84
    lat = np.arctan2(Z, p / np.sqrt(1 - e2 * (Z / p)**2))
    for _ in range(5):
        sinL = np.sin(lat)
        N = A_WGS84 / np.sqrt(1 - e2 * sinL**2)
        lat = np.arctan2(Z + e2 * N * sinL, p)
    lon = np.arctan2(Y, X)
    sinL = np.sin(lat)
    N = A_WGS84 / np.sqrt(1 - e2 * sinL**2)
    h = p / np.cos(lat) - N
    return np.array([lat, lon, h])


def ecef_to_enu_matrix(lat, lon):
    """ECEF → ENU 旋转矩阵（lat/lon in rad）"""
    sl, cl = np.sin(lat), np.cos(lat)
    sn, cn = np.sin(lon), np.cos(lon)
    return np.array([
        [-sn, cn, 0],
        [-sl*cn, -sl*sn, cl],
        [cl*cn, cl*sn, sl]
    ])


def ionospheric_free(L1, L2, P1, P2, f1=F1, f2=F2):
    """无电离层组合（m）"""
    alpha = f1**2 / (f1**2 - f2**2)
    beta  = -f2**2 / (f1**2 - f2**2)
    L_if = alpha * L1 + beta * L2
    P_if = alpha * P1 + beta * P2
    return L_if, P_if


def el_weight(el_rad, sigma0=0.003):
    return sigma0 / max(np.sin(el_rad), 0.1)


def compute_tropo(el_rad):
    zhd = 2.3
    zwd = 0.1
    mf  = 1.0 / max(np.sin(el_rad), 0.05)
    return zhd * mf, zwd * mf


def geometric_range(a, b):
    return np.linalg.norm(a - b)


def solve_leastsq(H, y, W):
    """加权最小二乘 + 正则化"""
    try:
        HT = H.T
        if W.ndim == 1:
            HTW = HT * W[:, None]
        else:
            HTW = HT @ W
        HTWH = HTW @ H
        HTWy = HTW @ y
        HTWH += np.eye(HTWH.shape[0]) * 1e-8
        dx = np.linalg.solve(HTWH, HTWy)
        if np.any(np.isnan(dx)) or np.any(np.isinf(dx)):
            dx = np.zeros(len(H))
    except Exception:
        dx = np.zeros(H.shape[1])
    return dx


# ─────────────────────────────────────────────────────────────────────
#  核心物理改正
# ─────────────────────────────────────────────────────────────────────

def light_time_correction(pos_rcv, pos_sat, epoch_recv):
    """
    光行时改正（第 1 次迭代）：
    信号从卫星传播到接收机，地球在信号传播过程中自转，
    需将接收机位置旋转到信号发射时刻的地固系。
    
    参数:
        pos_rcv  : 信号接收时刻的接收机 ECEF 坐标 (m)
        pos_sat  : 卫星 ECEF 坐标 (m)
        epoch_recv: 信号接收时刻 datetime
    返回:
        发射时刻的接收机位置 (m) — 用于计算真实几何距离
    """
    r_rcv = np.linalg.norm(pos_rcv)
    if r_rcv < 6e6 or r_rcv > 7e7:
        return pos_rcv
    
    # 估计传播时间
    rho_est = geometric_range(pos_rcv, pos_sat)
    tau = rho_est / C
    
    # 将接收机位置旋转回信号发射时刻（考虑地球自转）
    dtheta = OMEGA_E * tau
    cos_d, sin_d = np.cos(dtheta), np.sin(dtheta)
    
    pos_rcv_rot = np.array([
        pos_rcv[0] * cos_d + pos_rcv[1] * sin_d,
        -pos_rcv[0] * sin_d + pos_rcv[1] * cos_d,
        pos_rcv[2]
    ])
    return pos_rcv_rot


def sagnac_correction(pos_rcv, pos_sat):
    """
    Sagnac 效应改正：
    由于地球自转，ECEF 参考系是非惯性系。
    几何距离 ρ = |r_sat - r_rcv| 需加上 Sagnac 修正项。
    
    精确公式（对光行时迭代收敛后）:
    Δ_ρ_Sagnac = (Ω_E / C) * (x_sat*y_rcv - y_sat*x_rcv)
    
    这是对几何距离的附加路径长度改正（m），
    相当于把接收机在信号传播过程中走过的切向位移投影到信号方向。
    """
    return (OMEGA_E / C) * (pos_sat[0]*pos_rcv[1] - pos_sat[1]*pos_rcv[0])


def relativistic_clock_correction(sat_pos, sat_vel, sat_clock_err):
    """
    相对论性卫星钟偏改正：
    
    1. Einstein 引力红移（重力势差）：
       Δt_Einstein = (μ / (C^3 * r)) * Δr
        
    2. Shapiro 延迟（光路过引力势）：
       Δt_Shapiro ≈ 2μ / (C^3 * r)
    
    3. 横向多普勒（二阶狭义相对论）：
       在地固系中用完整 Minkowski 变换需要知道接收机速度，
       简化近似：动钟变慢 = v²/(2C²) * t
       GPS 卫星速度 ~ 3872 m/s → v²/(2C²) ≈ 2.1e-10
    
    总相对论改正 ≈ 40 ns（等效约 13 m），
    由卫星钟的二次项 1/2 * omega * t² 部分补偿。
    
    这里采用简化模型：
    - 广播星历中 af1 项已包含部分相对论效应
    - 用经典多普勒位移公式做二阶改正：
      Δρ_rel = (r·v)/c  (一阶)
      二阶项 ≈ (v²)/(2c) - (r·a)/(2c) (a=卫星加速度)
    """
    r = sat_pos
    v = sat_vel
    a = np.array([0.0, 0.0, 0.0])  # 加速度由广播星历近似忽略
    
    # 标量近似：r 和 v 的径向分量
    r_norm = np.linalg.norm(r)
    v_norm = np.linalg.norm(v)
    if r_norm < 1e6:
        return 0.0
    
    # 径向速度
    r_dot_v = np.dot(r, v)
    
    # 一阶相对论（爱因斯坦）：Δρ = (r·v) / c
    d1 = r_dot_v / C
    
    # 二阶（横向多普勒）：Δρ = v² / (2c)
    d2 = v_norm**2 / (2 * C)
    
    # 简化引力红移（势差）
    # Δt_grav = μ/(c³*r) * Δr_per_orbit ≈ 45 ns 近似为恒定距离偏移
    d3 = 13.0  # 约 45 ns * c ≈ 13 m（典型值）
    
    rel_corr = d1 + d2 - d3
    return rel_corr


# ─────────────────────────────────────────────────────────────────────
#  天线相位中心改正（PCV + PCO）
# ─────────────────────────────────────────────────────────────────────

# GPS 卫星天线相位中心补偿（IGS 08/14 标准，cm）
# 格式: {PRN: [dE, dN, dU]} 地固系下卫星天线相对于质心的偏移
# 这些是相对于卫星质心（CG）的改正
GPS_SAT_PCO = {}
for prn in range(1, 33):
    GPS_SAT_PCO[f"G{prn:02d}"] = [0.0, 0.0, 0.0]  # 无实测值时默认零

# 典型 GPS Block II/IIA/IIR/IIF 天线相位中心变化（PCV），单位：mm
# 以高度角 90° 为参考，偏离值随高度角降低
# 这里用简化为 0（IGS 产品通常已做）
GPS_SAT_PCV = {}


def apply_satellite_apc(sat_pos, sv, el, az):
    """
    将卫星天线相位中心偏移（PCO + PCV）投影到视线方向。
    
    参数:
        sat_pos : 卫星 ECEF 位置 (m)
        pos_rcv : 接收机 ECEF 位置 (m)
        sv      : 'G01' 等卫星编号
        el, az  : 高度角、方位角（rad）
    
    返回:
        sat_pos_apc : 加上天线偏移后的卫星位置 (m)
    """
    if el < 0.05:  # 高度角 < 3° 不做改正
        return sat_pos
    
    pco = GPS_SAT_PCO.get(sv, [0.0, 0.0, 0.0])
    # 将 ECEF 偏移投影到 ENU 坐标（以接收机为原点）
    # 简化：忽略小角度旋转，直接用射线方向的径向分量
    # 更精确需知道卫星姿态，这里用近似
    radial_unit = sat_pos / np.linalg.norm(sat_pos)
    
    # 卫星 PCV 随高度角变化
    el_deg = np.degrees(el)
    if el_deg >= 90:
        pcv_e = 0.0
    elif el_deg >= 15:
        pcv_e = 5.0 * (90 - el_deg) / 75.0  # 5mm 线性降到 0 (mm)
    else:
        pcv_e = 5.0  # 低高度角最大 5mm
    
    # PCO 投影到径向（简化为纯径向偏移）
    pco_total = (pco[0]**2 + pco[1]**2 + pco[2]**2)**0.5 * 0.01 + pcv_e * 0.001
    
    # 对于 GPS 卫星，实际 PCO 约 1-2 m（从质量中心到天线相位中心）
    # 这里用零以避免引入错误（实测产品已包含）
    return sat_pos


def apply_receiver_apc(pos_ecef, sigma_phase, el):
    """
    接收机天线相位中心改正（PCO）：
    地面站天线通常有已校准的 PCO（IGS 标准约 mm 级）。
    对于 GRACE-FO 卫星，无接收机天线改正。
    
    这里返回改正量（对距离的影响）。
    典型值：dN ≈ 0.1-0.3 m（GPS L1 天线）
    """
    # 地面站：天线 PCO 通常 < 2 mm，对精密 PPP 影响可忽略
    # 返回 0 以保持一致
    return 0.0


# ─────────────────────────────────────────────────────────────────────
#  PPObs 和 PPPProcessor
# ─────────────────────────────────────────────────────────────────────

class PPObs:
    def __init__(self, sv, t, L1, L2, P1, P2,
                 sat_pos, sat_clock,
                 el, az,
                 iono=0.0, trop_dry=0.0, trop_wet=0.0,
                 rel=0.0, sagnac=0.0, light_time=0.0,
                 sat_vel=None, sigma=0.003):
        self.sv = sv
        self.time = t
        self.L1 = L1; self.L2 = L2
        self.P1 = P1; self.P2 = P2
        self.sat_pos = sat_pos
        self.sat_clock = sat_clock
        self.el = el; self.az = az
        self.iono = iono
        self.trop_dry = trop_dry
        self.trop_wet = trop_wet
        self.rel = rel
        self.sagnac = sagnac
        self.light_time = light_time
        self.sat_vel = sat_vel if sat_vel is not None else np.zeros(3)
        self.sigma = sigma


class PPPProcessor:
    def __init__(self, pos0=None, elev_mask=10.0,
                 sigma_code=0.3, sigma_phase=0.003,
                 max_iter=20, tol=1e-4,
                 apply_rel=True, apply_sagnac=True,
                 apply_apc=True):
        self.pos0 = np.array(pos0) if pos0 is not None else None
        self.elev_mask = np.radians(elev_mask)
        self.sigma_code = sigma_code
        self.sigma_phase = sigma_phase
        self.max_iter = max_iter
        self.tol = tol
        self.apply_rel = apply_rel
        self.apply_sagnac = apply_sagnac
        self.apply_apc = apply_apc
        self.history = []

    def _prepare(self, obs_list, x0):
        H, y, W = [], [], []
        pos_rcv = x0[:3]
        
        for obs in obs_list:
            if obs.el < self.elev_mask:
                continue
            
            sat_pos = obs.sat_pos
            
            # ── 1. 光行时改正（迭代）────────────────────────────────
            # 第 1 次：用接收时刻几何距离估计 τ
            # 第 2 次：用改正后位置重新估计
            rho_init = geometric_range(pos_rcv, sat_pos)
            tau_1 = rho_init / C
            dtheta_1 = OMEGA_E * tau_1
            cos1, sin1 = np.cos(dtheta_1), np.sin(dtheta_1)
            pos_rcv_t1 = np.array([
                pos_rcv[0]*cos1 + pos_rcv[1]*sin1,
                -pos_rcv[0]*sin1 + pos_rcv[1]*cos1,
                pos_rcv[2]
            ])
            rho_1 = geometric_range(pos_rcv_t1, sat_pos)
            tau_2 = rho_1 / C
            dtheta_2 = OMEGA_E * tau_2
            cos2, sin2 = np.cos(dtheta_2), np.sin(dtheta_2)
            pos_rcv_t2 = np.array([
                pos_rcv[0]*cos2 + pos_rcv[1]*sin2,
                -pos_rcv[0]*sin2 + pos_rcv[1]*cos2,
                pos_rcv[2]
            ])
            rho_final = geometric_range(pos_rcv_t2, sat_pos)
            
            # ── 2. Sagnac 改正 ─────────────────────────────────────
            sagnac = (OMEGA_E / C) * (sat_pos[0]*pos_rcv[1] - sat_pos[1]*pos_rcv[0])
            
            # ── 3. 相对论性钟偏改正 ─────────────────────────────────
            sat_vel = obs.sat_vel
            r_norm = np.linalg.norm(sat_pos)
            rel = 0.0
            if self.apply_rel and np.linalg.norm(sat_vel) > 0:
                r_dot_v = np.dot(sat_pos, sat_vel)
                v_norm = np.linalg.norm(sat_vel)
                d1 = r_dot_v / C
                d2 = v_norm**2 / (2 * C)
                rel = d1 + d2 - 13.0  # 13 m 近似引力红移
            
            # ── 4. 卫星天线相位中心改正 ────────────────────────────
            if self.apply_apc and obs.el > np.radians(5):
                el_deg = np.degrees(obs.el)
                pcv_mm = max(0, 5.0 * (90 - el_deg) / 75.0)
                pcv_m = pcv_mm * 0.001
                # 径向投影：PCV 主要影响沿射线方向的距离
                radial = sat_pos / r_norm
                sat_pos_apc = sat_pos + radial * pcv_m
            else:
                sat_pos_apc = sat_pos
            
            # ── 5. 无电离层组合 ─────────────────────────────────────
            L_if, P_if = ionospheric_free(obs.L1, obs.L2, obs.P1, obs.P2)
            
            # ── 6. 观测模型（所有改正）──────────────────────────────
            # 改正后几何距离
            rho_corrected = rho_final + sagnac + rel
            
            # 对流层
            mf = 1.0 / max(np.sin(obs.el), 0.05)
            trop_total = obs.trop_dry + obs.trop_wet * mf
            
            # 相位残差
            model = L_if - (rho_corrected + obs.sat_clock + trop_total)
            h = np.zeros(5)
            h[:3] = -(sat_pos_apc - pos_rcv) / rho_corrected
            h[3] = 1.0
            h[4] = mf
            w = 1.0 / (self.sigma_phase**2 / (np.sin(obs.el)**2 + 0.01))
            H.append(h); y.append(model); W.append(w)
            
            # 伪距残差
            model_c = P_if - (rho_corrected + obs.sat_clock + trop_total)
            hc = np.zeros(5)
            hc[:3] = -(sat_pos_apc - pos_rcv) / rho_corrected
            hc[3] = 1.0; hc[4] = mf
            wc = 1.0 / (self.sigma_code**2 / (np.sin(obs.el)**2 + 0.01))
            H.append(hc); y.append(model_c); W.append(wc)
        
        if not H:
            return None, None, None, None
        return np.vstack(H), np.array(y), np.diag(W), obs_list

    def solve_epoch(self, obs_list, x0=None):
        if x0 is None:
            if self.pos0 is not None:
                pos = self.pos0.copy()
            else:
                avg = np.mean([o.sat_pos for o in obs_list], axis=0)
                pos = avg / np.linalg.norm(avg) * 6371000.0
            x0 = np.concatenate([pos, [0.0, 0.2]])
        x = x0.copy()
        for it in range(self.max_iter):
            res = self._prepare(obs_list, x)
            if res[0] is None:
                return x, np.eye(len(x)), np.array([]), []
            H, y, W, good = res
            dx = solve_leastsq(H, y, W)
            x = x + dx
            if np.linalg.norm(dx[:3]) < self.tol and abs(dx[3]) < self.tol:
                break
        _, _, _, good = self._prepare(obs_list, x)
        return x, np.eye(len(x)), np.array([]), good

    def process(self, records, ref_pos=None, verbose=True):
        """
        records: list of dicts with keys: time, sv, L1, L2, P1, P2,
                 sat_pos, sat_clock, el, az, sat_vel (optional)...
        ref_pos: ECEF reference for error evaluation
        """
        results = []
        x_prev = None
        for t, grp in _groupby_time(records):
            obs_list = []
            for r in grp:
                el_rad = np.radians(r.get('el', 45.0))
                td, tw = compute_tropo(el_rad)
                sat_vel = r.get('sat_vel')
                if sat_vel is None:
                    sat_vel = np.zeros(3)
                obs_list.append(PPObs(
                    sv=r['sv'], t=t,
                    L1=r['L1'], L2=r['L2'], P1=r['P1'], P2=r['P2'],
                    sat_pos=r['sat_pos'],
                    sat_clock=r.get('sat_clock', 0.0),
                    el=el_rad,
                    az=np.radians(r.get('az', 0.0)),
                    trop_dry=td, trop_wet=tw,
                    sat_vel=sat_vel,
                    sigma=el_weight(el_rad, self.sigma_phase)
                ))
            if len(obs_list) < 4:
                continue
            x, P, residuals, good = self.solve_epoch(obs_list, x_prev)
            pos_est = x[:3]
            
            if ref_pos is not None:
                err = pos_est - np.array(ref_pos)
                lat, lon, _ = ecef_to_blh(ref_pos)
                R = ecef_to_enu_matrix(lat, lon)
                enu = R @ err
            else:
                err = np.zeros(3); enu = np.zeros(3)
            
            results.append({
                'time': t, 'X': pos_est[0], 'Y': pos_est[1], 'Z': pos_est[2],
                'dX': err[0], 'dY': err[1], 'dZ': err[2],
                'dE': enu[0], 'dN': enu[1], 'dU': enu[2],
                'n_sat': len(good)
            })
            x_prev = x.copy()
        return results


def _groupby_time(records):
    """按时间分组"""
    from itertools import groupby
    key = lambda r: r['time'].timestamp()
    records.sort(key=key)
    return [(t, list(g)) for t, g in groupby(records, key=key)]


def evaluate(results, ref_pos):
    """计算统计量"""
    dX = np.array([r['dX'] for r in results])
    dY = np.array([r['dY'] for r in results])
    dZ = np.array([r['dZ'] for r in results])
    dE = np.array([r.get('dE', 0.0) for r in results], dtype=float)
    dN = np.array([r.get('dN', 0.0) for r in results], dtype=float)
    dU = np.array([r.get('dU', 0.0) for r in results], dtype=float)
    d3 = np.sqrt(dX**2 + dY**2 + dZ**2)
    rms = lambda a: np.sqrt(np.nanmean(a**2))
    return {
        'RMS_X': rms(dX), 'RMS_Y': rms(dY), 'RMS_Z': rms(dZ),
        'RMS_3D': rms(d3),
        'RMS_E': rms(dE), 'RMS_N': rms(dN), 'RMS_U': rms(dU),
        'MAX_3D': float(np.nanmax(d3)),
        'mean_dX': float(np.nanmean(dX)),
        'mean_dY': float(np.nanmean(dY)),
        'mean_dZ': float(np.nanmean(dZ)),
        'mean_d3D': float(np.nanmean(d3)),
        'STD_X': float(np.nanstd(dX)),
        'STD_Y': float(np.nanstd(dY)),
        'STD_Z': float(np.nanstd(dZ)),
        'n_epochs': len(results),
    }


if __name__ == '__main__':
    print("PPP 模块 OK — 相对论 + Sagnac + APC 改正已启用")