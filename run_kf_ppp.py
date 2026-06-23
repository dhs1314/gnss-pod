#!/usr/bin/env python3
"""
GRACE-FO Kalman Filter PPP — 使用真实 GPS1B 载波相位 + IGS SP3 精密星历
关键修正:
  1. GPS→UTC: 直接使用 J2000 + gps_sod (GPS1B 时间已是 UTC 参照)
  2. 使用载波相位 (L1_phase, L2_phase) 计算 IF 组合
  3. 预估计 code-phase bias: bias_sv = median(P_if - L_if_phase)
  4. 对载波相位施加 bias 改正, 使其成为绝对测距 (载波精度 + 伪距绝对基准)
  5. SP3 卫星钟差 + 光行时 + Sagnac 改正
  6. KF: [X, Y, Z, clk_r]

用法:
  py -3.12 run_kf_ppp.py                        # 默认: 2024-04-29, 4h
  py -3.12 run_kf_ppp.py --hours 8               # 处理8小时
"""
import sys, os, pickle, csv, argparse, json
from pathlib import Path
from datetime import datetime, timedelta
from collections import defaultdict
import numpy as np

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
F1_SQ, F2_SQ = F1*F1, F2*F2
ALPHA = F1_SQ / (F1_SQ - F2_SQ)
BETA  = -F2_SQ / (F1_SQ - F2_SQ)
OMEGA_E = 7.2921151467e-5
A_WGS84 = 6378137.0
E2_WGS84 = 0.00669437999014
J2000 = datetime(2000, 1, 1, 12, 0, 0)

# ── 坐标转换 ──────────────────────────────────────────────────────────
def ecef_to_blh(pos):
    X, Y, Z = float(pos[0]), float(pos[1]), float(pos[2])
    p = np.sqrt(X**2 + Y**2)
    if p < 1e-6:
        lat = np.pi/2 if Z >= 0 else -np.pi/2
        return np.array([lat, 0.0, 0.0])
    lat = np.arctan2(Z, p)
    for _ in range(10):
        sinL = np.sin(lat)
        N = A_WGS84 / np.sqrt(1 - E2_WGS84 * sinL**2)
        lat_new = np.arctan2(Z + E2_WGS84 * N * sinL, p)
        if abs(lat_new - lat) < 1e-15: lat = lat_new; break
        lat = lat_new
    lon = np.arctan2(Y, X)
    return np.array([lat, lon, 0.0])

def ecef_to_enu_matrix(lat, lon):
    sl, cl = np.sin(lat), np.cos(lat)
    sn, cn = np.sin(lon), np.cos(lon)
    return np.array([[-sn, cn, 0],
                     [-sl*cn, -sl*sn, cl],
                     [cl*cn, cl*sn, sl]])

# ── 几何距离改正 ──────────────────────────────────────────────────────
def geometric_range(rcv_pos, sat_pos, sat_vel):
    rho = np.linalg.norm(sat_pos - rcv_pos)
    for _ in range(2):
        tau = rho / C
        dtheta = OMEGA_E * tau
        cos_d, sin_d = np.cos(dtheta), np.sin(dtheta)
        rcv_rot = np.array([
            rcv_pos[0]*cos_d + rcv_pos[1]*sin_d,
            -rcv_pos[0]*sin_d + rcv_pos[1]*cos_d,
            rcv_pos[2]
        ])
        rho = np.linalg.norm(sat_pos - rcv_rot)
    sagnac = (OMEGA_E / C) * (sat_pos[0]*rcv_pos[1] - sat_pos[1]*rcv_pos[0])
    rel = np.dot(sat_pos, sat_vel) / C if sat_vel is not None and np.linalg.norm(sat_vel) > 0 else 0.0
    return rho + sagnac + rel

# ── SP3 精密星历 ──────────────────────────────────────────────────────
def load_sp3(pkl_path):
    return pickle.load(open(pkl_path, 'rb'))

def get_sat_pos_sp3(sp3_data, sv, utc_dt):
    ts = sp3_data['ts']; epochs = sp3_data['epochs']
    t0 = t1 = None
    for i, ti in enumerate(ts):
        if ti >= utc_dt: t1 = ti; t0 = ts[i-1] if i > 0 else None; break
        t0 = ti
    if t1 is None: t0 = t1 = ts[-1]
    if t0 is None: t0 = ts[0]
    p0 = epochs.get(t0, {}).get(sv); p1 = epochs.get(t1, {}).get(sv)
    if p0 is None and p1 is None: return None, np.zeros(3), 0.0
    if p0 is None: return np.array(p1[:3]), np.zeros(3), float(p1[3]) if len(p1)>3 else 0.0
    if p1 is None: return np.array(p0[:3]), np.zeros(3), float(p0[3]) if len(p0)>3 else 0.0
    dt_sec = (utc_dt - t0).total_seconds(); dt_tot = (t1 - t0).total_seconds()
    if dt_tot == 0: return np.array(p0[:3]), np.zeros(3), float(p0[3]) if len(p0)>3 else 0.0
    a = dt_sec / dt_tot
    pos = np.array(p0[:3]) * (1-a) + np.array(p1[:3]) * a
    vel = (np.array(p1[:3]) - np.array(p0[:3])) / dt_tot if dt_tot > 0 else np.zeros(3)
    clk = p0[3]*(1-a) + p1[3]*a if len(p0)>3 and len(p1)>3 else 0.0
    return pos, vel, clk

# ── GNV1B 参考轨道 ────────────────────────────────────────────────────
def load_gnv1b(filepath):
    orbit = {}
    with open(filepath, encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            parts = line.split()
            if len(parts) < 6: continue
            try:
                t_gps = float(parts[0]); flag = parts[2]
                if flag not in ('C', 'E'): continue
                X, Y, Z = float(parts[3]), float(parts[4]), float(parts[5])
                if abs(X) < 1e3: continue
                orbit[t_gps] = np.array([X, Y, Z])
            except (ValueError, IndexError): continue
    return orbit

def interpolate_ref(orbit, gps_sod):
    ts = sorted(orbit.keys())
    t0 = t1 = None
    for i, ti in enumerate(ts):
        if ti >= gps_sod: t1 = ti; t0 = ts[i-1] if i > 0 else None; break
        t0 = ti
    if t1 is None: t0 = t1 = ts[-1]
    if t0 is None: t0 = ts[0]
    if t0 == t1: return orbit[t0]
    a = (gps_sod - t0) / (t1 - t0)
    return orbit[t0] * (1-a) + orbit[t1] * a

# ── 预估计 code-phase bias ────────────────────────────────────────────
def estimate_code_phase_bias(epoch_obs):
    """对每颗卫星, 计算 P_if - L_if_phase 的中位数作为 code-phase bias

    L_if_phase = ALPHA * L1_phase + BETA * L2_phase  (载波 IF, 含整周模糊度)
    P_if       = ALPHA * P1 + BETA * P2                (伪距 IF, 无模糊度)

    bias = P_if - L_if_phase ≈ 伪距偏差 - 整周模糊度 = 常数
    """
    sv_diffs = defaultdict(list)
    for gps_sod, ep_data in epoch_obs.items():
        for d in ep_data:
            diff = d['P_if'] - d['L_if_phase']
            # 过滤异常值
            if abs(diff) < 1000:  # 合理范围 < 1000m
                sv_diffs[d['sv']].append(diff)

    bias = {}
    for sv, diffs in sv_diffs.items():
        if len(diffs) >= 5:
            bias[sv] = float(np.median(diffs))
            print(f"  Bias {sv}: {bias[sv]:.3f} m (std={np.std(diffs):.3f}m, n={len(diffs)})")
    return bias

# ── Kalman Filter ─────────────────────────────────────────────────────
class KFPPP:
    """KF: [X, Y, Z, clk_r]"""
    def __init__(self, pos0):
        self.x = np.array([pos0[0], pos0[1], pos0[2], 0.0], dtype=float)
        self.P = np.diag([1e3**2, 1e3**2, 1e3**2, 1e5**2])
        self.Q_pos = 500.0**2
        self.Q_clk = 100.0**2

    def predict(self, dt):
        if dt <= 0: return
        Q = np.diag([self.Q_pos*dt, self.Q_pos*dt, self.Q_pos*dt, self.Q_clk*dt])
        self.P = self.P + Q

    def update(self, obs_list):
        """obs_list: [(sat_pos, sat_vel, L_corr, P_corr, el_rad)]"""
        H_list, y_list, w_list = [], [], []
        for sat_pos, sat_vel, L_corr, P_corr, el in obs_list:
            if el < 0.087: continue

            rho = geometric_range(self.x[:3], sat_pos, sat_vel)
            e_vec = (sat_pos - self.x[:3]) / max(rho, 1e6)
            h = np.array([-e_vec[0], -e_vec[1], -e_vec[2], 1.0])

            # 载波相位 (bias-corrected, 绝对测距)
            res_l = L_corr - (rho + self.x[3])
            w_l = 1.0 / (0.005**2 / (np.sin(el)**2 + 0.01))  # 5mm sigma
            H_list.append(h); y_list.append(res_l); w_list.append(w_l)

            # 伪距 (绝对测距, 有噪声)
            res_p = P_corr - (rho + self.x[3])
            w_p = 1.0 / (1.0**2 / (np.sin(el)**2 + 0.01))  # 1m sigma
            H_list.append(h); y_list.append(res_p); w_list.append(w_p)

        if len(H_list) < 8: return False

        H = np.array(H_list); y = np.array(y_list); W = np.diag(w_list)
        try:
            PHT = self.P @ H.T
            S = H @ PHT + np.linalg.inv(W)
            K = PHT @ np.linalg.inv(S)
            self.x = self.x + K @ y
            self.P = (np.eye(4) - K @ H) @ self.P
            return True
        except Exception:
            return False

    def get_position(self):
        return self.x[:3]

# ── 主处理 ────────────────────────────────────────────────────────────
def process_day(date_str, nhours=24, data_dir='./data', grace_id='C', interval=30.0):
    y, m, d = [int(x) for x in date_str.split('-')]
    doy = (datetime(y, m, d) - datetime(y, 1, 1)).days + 1
    data_path = Path(data_dir)

    # 1. 加载 GNV1B 参考轨道
    gnv_path = data_path / 'gracefo' / str(y) / date_str / f'GNV1B_{date_str}_{grace_id}_04.txt'
    if not gnv_path.exists():
        raise FileNotFoundError(f"GNV1B not found: {gnv_path}")
    ref_orbit = load_gnv1b(str(gnv_path))
    print(f"[GNV1B] {len(ref_orbit)} epochs")

    # 2. 加载 GPS1B (优先 batch_v12 格式 = L1_phase)
    gps1b_pkl = data_path / 'gracefo' / str(y) / date_str / f'GPS1B_{date_str}_{grace_id}_04.pkl'
    if not gps1b_pkl.exists():
        gps1b_pkl = data_path / 'gracefo' / str(y) / date_str / f'gps1b_{grace_id}.pkl'
    if not gps1b_pkl.exists():
        raise FileNotFoundError(f"GPS1B not found")
    gps1b_obs = pickle.load(open(str(gps1b_pkl), 'rb'))
    print(f"[GPS1B] {len(gps1b_obs)} epochs")

    # 确定数据格式
    test_sod = sorted(gps1b_obs.keys())[0]
    test_sv = sorted(gps1b_obs[test_sod].keys())[0]
    test_rec = gps1b_obs[test_sod][test_sv]
    use_phase = 'L1_phase' in test_rec and 'L2_phase' in test_rec
    print(f"[DATA] {'L1_phase/L2_phase' if use_phase else 'L1_range/L2_range'}")

    # 3. 加载 SP3
    sp3_pkl = data_path / str(y) / f'{doy:03d}' / 'igs_sp3_FIN.pkl'
    if not sp3_pkl.exists():
        sp3_pkl = data_path / str(y) / f'{doy:03d}' / 'igs_sp3.pkl'
    if not sp3_pkl.exists():
        raise FileNotFoundError(f"SP3 not found for DOY {doy}")
    sp3_data = load_sp3(str(sp3_pkl))
    print(f"[SP3] {len(sp3_data['ts'])} epochs")

    # 4. 构建观测序列
    gps_sod_start = min(gps1b_obs.keys())
    gps_sod_end = gps_sod_start + nhours * 3600

    epoch_obs = {}
    n_skip = {'sp3': 0, 'clk': 0, 'el': 0, 'rng': 0}

    for gps_sod, sv_obs in sorted(gps1b_obs.items()):
        if not (gps_sod_start <= gps_sod <= gps_sod_end): continue
        dt_from_start = gps_sod - gps_sod_start
        nearest = round(dt_from_start / interval) * interval
        if abs(dt_from_start - nearest) > max(2.0, interval * 0.1): continue

        utc_dt = J2000 + timedelta(seconds=gps_sod)  # 不减去闰秒!

        grace_pos = interpolate_ref(ref_orbit, gps_sod)
        lat, lon, _ = ecef_to_blh(grace_pos)
        R_enu = ecef_to_enu_matrix(lat, lon)

        epoch_data = []
        for sv_id, rec in sv_obs.items():
            sat_pos, sat_vel, sat_clk = get_sat_pos_sp3(sp3_data, sv_id, utc_dt)
            if sat_pos is None: n_skip['sp3'] += 1; continue
            if abs(sat_clk) > 0.1 * C: n_skip['clk'] += 1; continue

            rng_raw = np.linalg.norm(sat_pos - grace_pos)
            if not (1.8e7 < rng_raw < 2.8e7): n_skip['rng'] += 1; continue

            los = sat_pos - grace_pos
            e_enu = R_enu @ (los / rng_raw)
            el = float(np.arcsin(np.clip(e_enu[2], -1.0, 1.0)))
            if el < np.radians(5.0): n_skip['el'] += 1; continue

            if use_phase:
                L1_m = float(rec['L1_phase']); L2_m = float(rec['L2_phase'])
            else:
                L1_m = float(rec['L1_range']); L2_m = float(rec['L2_range'])

            L_if_phase = ALPHA * L1_m + BETA * L2_m
            P_if = float(rec['P_if'])

            epoch_data.append({
                'sv': sv_id, 'sat_pos': sat_pos.astype(float),
                'sat_vel': sat_vel.astype(float), 'sat_clk': float(sat_clk),
                'L_if_phase': L_if_phase, 'P_if': P_if, 'el': float(el),
            })

        if len(epoch_data) >= 5:
            epoch_obs[gps_sod] = epoch_data

    print(f"[OBS] {len(epoch_obs)} epochs (skip: {n_skip})")

    if len(epoch_obs) < 10: return None

    # 5. 预估计 code-phase bias
    print("[BIAS] Estimating code-phase bias per satellite...")
    bias = estimate_code_phase_bias(epoch_obs)
    print(f"[BIAS] Estimated bias for {len(bias)} satellites")

    # 6. KF PPP
    sorted_sods = sorted(epoch_obs.keys())
    ref_pos0 = interpolate_ref(ref_orbit, sorted_sods[0])
    kf = KFPPP(ref_pos0)

    results = []
    prev_gps_sod = None

    print(f"[KF] Processing {len(sorted_sods)} epochs...")
    for i, gps_sod in enumerate(sorted_sods):
        ep_data = epoch_obs[gps_sod]

        if prev_gps_sod is not None:
            kf.predict(gps_sod - prev_gps_sod)

        obs_list = []
        for d in ep_data:
            sv = d['sv']
            if sv not in bias: continue  # 跳无bias的卫星

            # 改正: 卫星钟差(SP3) + code-phase bias
            L_corr = d['L_if_phase'] + d['sat_clk'] + bias[sv]
            P_corr = d['P_if'] + d['sat_clk']
            obs_list.append((d['sat_pos'], d['sat_vel'], L_corr, P_corr, d['el']))

        kf.update(obs_list)

        pos_est = kf.get_position()
        ref_pos = interpolate_ref(ref_orbit, gps_sod)
        err = pos_est - ref_pos

        lat_r, lon_r, _ = ecef_to_blh(ref_pos)
        R = ecef_to_enu_matrix(lat_r, lon_r)
        enu = R @ err
        d3 = float(np.linalg.norm(err))

        results.append({
            'time': J2000 + timedelta(seconds=gps_sod),
            'gps_sod': gps_sod,
            'dE': float(enu[0]), 'dN': float(enu[1]), 'dU': float(enu[2]),
            'd3': d3, 'clk': float(kf.x[3]), 'n_sat': len(obs_list),
        })
        prev_gps_sod = gps_sod

        if (i+1) % 100 == 0:
            recent = [r['d3'] for r in results[-100:]]
            print(f"  epoch {i+1}/{len(sorted_sods)}: 3D={np.mean(recent)*100:.1f}cm, "
                  f"clk={kf.x[3]:.1f}m, n_sat={np.mean([r['n_sat'] for r in results[-100:]]):.1f}")

    print(f"[KF] {len(results)} epochs solved")

    # 7. 统计
    dE = np.array([r['dE'] for r in results])
    dN = np.array([r['dN'] for r in results])
    dU = np.array([r['dU'] for r in results])
    d3 = np.array([r['d3'] for r in results])

    def rms(a): return float(np.sqrt(np.nanmean(a**2)))

    stats = {}
    for thresh, label in [(1e9, 'all'), (100.0, '<100m'), (50.0, '<50m'),
                           (10.0, '<10m'), (5.0, '<5m'), (1.0, '<1m')]:
        mask = d3 < thresh if thresh < 1e8 else np.ones(len(d3), dtype=bool)
        n_v = int(mask.sum())
        if n_v > 5:
            stats[label] = {
                'n': n_v,
                'rms_e': rms(dE[mask])*100, 'rms_n': rms(dN[mask])*100,
                'rms_u': rms(dU[mask])*100, 'rms_3d': rms(d3[mask])*100,
                'mean_3d': float(np.nanmean(d3[mask]))*100,
                'max_3d': float(np.nanmax(d3[mask]))*100,
            }

    print("\n" + "="*65)
    print(f"  GRACE-FO KF PPP Results — {date_str} ({nhours}h, SP3, Corrected Phase)")
    print("="*65)
    header = f"  {'Filter':<12s} {'N':>6s} {'E(cm)':>8s} {'N(cm)':>8s} {'U(cm)':>8s} {'3D(cm)':>8s} {'Mean':>8s} {'Max':>8s}"
    print(header); print("  " + "-"*60)
    for label in ['all', '<100m', '<50m', '<10m', '<5m', '<1m']:
        if label in stats:
            s = stats[label]
            print(f"  {label:<12s} {s['n']:>6d} {s['rms_e']:>8.2f} {s['rms_n']:>8.2f} "
                  f"{s['rms_u']:>8.2f} {s['rms_3d']:>8.2f} {s['mean_3d']:>8.2f} {s['max_3d']:>8.2f}")
    print("="*65)

    return {'stats': stats, 'results': results, 'date': date_str, 'bias': bias}

# ── 输出 ──────────────────────────────────────────────────────────────
def save_output(stats, results, date_str, bias, output_dir='./output'):
    out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)
    tag = date_str.replace('-', '')

    csv_path = out / f'kf_ppp_{tag}.csv'
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['time', 'gps_sod', 'dE_m', 'dN_m', 'dU_m', 'd3_m', 'clock_m', 'n_sat'])
        for r in results:
            writer.writerow([
                r['time'].isoformat(), f"{r['gps_sod']:.6f}",
                f"{r['dE']:.4f}", f"{r['dN']:.4f}", f"{r['dU']:.4f}", f"{r['d3']:.4f}",
                f"{r['clk']:.4f}", r['n_sat'],
            ])
    print(f"[CSV] {csv_path}")

    json_path = out / f'kf_ppp_{tag}_stats.json'
    json.dump(stats, open(json_path, 'w'), indent=2)
    print(f"[JSON] {json_path}")

    return csv_path

def generate_plots(results, stats, date_str, output_dir='./output'):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("[PLOT] matplotlib not available")
        return

    out = Path(output_dir); tag = date_str.replace('-', '')
    times_h = [(r['time'] - results[0]['time']).total_seconds()/3600 for r in results]
    dE_cm = [r['dE']*100 for r in results]
    dN_cm = [r['dN']*100 for r in results]
    dU_cm = [r['dU']*100 for r in results]
    d3_cm = [r['d3']*100 for r in results]

    fig, axes = plt.subplots(4, 1, figsize=(14, 10), sharex=True)
    for ax, data, color, ylabel in [
        (axes[0], dE_cm, '#1565C0', 'E (cm)'),
        (axes[1], dN_cm, '#2E7D32', 'N (cm)'),
        (axes[2], dU_cm, '#C62828', 'U (cm)'),
        (axes[3], d3_cm, '#333333', '3D (cm)'),
    ]:
        ax.plot(times_h, data, color=color, linewidth=0.5, alpha=0.7)
        ax.axhline(0, color='k', linewidth=0.5)
        ax.set_ylabel(ylabel); ax.grid(True, alpha=0.3)
    axes[0].set_title(f'GRACE-FO KF PPP Error vs GNV1B — {date_str} (SP3 + Bias-Corrected Phase)')
    axes[3].set_xlabel('Time (hours)')
    if 'all' in stats:
        s = stats['all']
        fig.suptitle(f'3D RMS={s["rms_3d"]:.1f}cm | E={s["rms_e"]:.1f} N={s["rms_n"]:.1f} U={s["rms_u"]:.1f} cm', fontsize=9)
    plt.tight_layout()
    png_path = out / f'kf_ppp_{tag}.png'
    fig.savefig(str(png_path), dpi=150)
    plt.close(fig)
    print(f"[PLOT] {png_path}")

# ═══════════════════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(description='GRACE-FO KF PPP')
    parser.add_argument('--date', default='2024-04-29')
    parser.add_argument('--hours', type=float, default=4.0)
    parser.add_argument('--data-dir', default='./data')
    parser.add_argument('--output-dir', default='./output_kf')
    parser.add_argument('--grace-id', default='C')
    parser.add_argument('--interval', type=float, default=30.0)
    args = parser.parse_args()

    print(f"\n{'='*60}")
    print(f"  GRACE-FO KF PPP (Bias-Corrected Carrier Phase)")
    print(f"  Date: {args.date}  Hours: {args.hours}h  dt: {args.interval}s")
    print(f"{'='*60}")

    result = process_day(date_str=args.date, nhours=args.hours,
                         data_dir=args.data_dir, grace_id=args.grace_id,
                         interval=args.interval)
    if result is None:
        print("[FATAL] Processing failed"); sys.exit(1)

    save_output(result['stats'], result['results'], result['date'],
                result.get('bias', {}), output_dir=args.output_dir)
    generate_plots(result['results'], result['stats'], result['date'],
                   output_dir=args.output_dir)
    print(f"\n[DONE] {args.date}")

if __name__ == '__main__':
    main()
