#!/usr/bin/env python3
"""
GRACE-FO PPP 批量处理 v1.2.0
真实数据链路: GPS1B (ISDC) + IGS Final SP3 (BKG) + MW模糊度固定 + Kalman Filter
处理范围: 2024年5月1日 - 6月30日 (62天)
"""
import sys, os, math, io, json, pickle, ssl
import tarfile, urllib.request
import numpy as np
from datetime import datetime, timedelta
from pathlib import Path
from itertools import groupby
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).parent))
from src.gps1b_loader import gps_sod_to_utc, download_gps1b
from src.sp3_loader import load_igs_sp3, get_gps_pos_from_sp3
from src.ppp import ecef_to_blh, ecef_to_enu_matrix

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
F1_SQ, F2_SQ = F1*F1, F2*F2
ALPHA = F1_SQ / (F1_SQ - F2_SQ)
BETA  = -F2_SQ / (F1_SQ - F2_SQ)
LAM_W = C / (F1 - F2)
MU_E = 3.986004418e14
OMEGA_E = 7.2921151467e-5
GPS_ORIGIN_J2000 = datetime(2000, 1, 1, 12, 0, 0)
ISDC_TGZ = "https://isdc-data.gfz.de/grace-fo/Level-1B/JPL/INSTRUMENT/RL04/{year}/"

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
    L1 = rec.get('L1_phase') or rec.get('L1_range')
    L2 = rec.get('L2_phase') or rec.get('L2_range')
    P1 = rec.get('CA_range') or rec.get('L1_range')
    P2 = rec.get('L2_range')
    if L1 is None or L2 is None: return None
    snr = rec.get('L1_SNR')
    if snr is not None and snr < 5: return None
    rec['L1'] = L1; rec['L2'] = L2; rec['P1'] = P1; rec['P2'] = P2
    rec['L_if'] = ALPHA*L1 + BETA*L2
    rec['P_if'] = (ALPHA*P1 + BETA*P2) if (P1 and P2) else P1
    return rec

def load_gps1b_text(text, gf='C'):
    lines = text.split('\n')
    h_end = next((i for i,l in enumerate(lines) if 'End of YAML header' in l), 0)
    obs = {}
    for line in lines[h_end+1:]:
        if not line.strip() or line.startswith('#'): continue
        parts = line.split()
        rec = parse_gps1b_record(parts, gf)
        if rec is None: continue
        sod = rec['gps_sod']
        obs.setdefault(sod, {})[rec['sv']] = rec
    return obs

def parse_gnv1b(path):
    orbit = {}
    GPS0 = datetime(2000, 1, 1, 12, 0, 0)
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            p = line.split()
            if len(p) < 6: continue
            try:
                tg = float(p[0]); flag = p[2]
                if flag not in ('C','E'): continue
                X,Y,Z = float(p[3]),float(p[4]),float(p[5])
                if abs(X) < 1e3: continue
                orbit[GPS0 + timedelta(seconds=tg)] = np.array([X,Y,Z])
            except: continue
    return orbit

def download_gracefo(year, month, day, data_dir='./data', gf='C'):
    ds = f"{year:04d}-{month:02d}-{day:02d}"
    cache_obs = Path(data_dir)/'gracefo'/str(year)/ds/f"GPS1B_{ds}_{gf}_04.pkl"
    cache_gnv = Path(data_dir)/'gracefo'/str(year)/ds/f"GNV1B_{ds}_{gf}_04.txt"
    cache_gnv.parent.mkdir(parents=True, exist_ok=True)

    if cache_obs.exists():
        try:
            obs = pickle.load(open(cache_obs,'rb'))
            gnv = parse_gnv1b(str(cache_gnv))
            n = sum(len(v) for v in obs.values())
            print(f"  [缓存] {ds}: {len(obs)} 历元 ({n} obs), GNV1B {len(gnv)} 历元")
            return obs, gnv
        except: pass

    fname = f"gracefo_1B_{ds}_RL04.ascii.noLRI.tgz"
    url = ISDC_TGZ.format(year=year) + fname
    print(f"  [下载] {fname}...", flush=True)
    ctx = ssl.create_default_context()
    ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE
    try:
        req = urllib.request.Request(url, headers={'User-Agent':'curl/7.88','Accept-Encoding':'gzip'})
        with urllib.request.urlopen(req, timeout=180, context=ctx) as r:
            data = r.read()
        print(f"  完成: {len(data)//1024//1024} MB")
    except Exception as e:
        print(f"  [错误] {e}"); return None, None

    try:
        tar = tarfile.open(fileobj=io.BytesIO(data), mode='r:*')
        gps_m = next((m for m in tar.getmembers() if 'GPS1B' in m.name and gf in m.name and m.name.endswith('.txt')), None)
        gnv_m = next((m for m in tar.getmembers() if 'GNV1B' in m.name and gf in m.name and m.name.endswith('.txt')), None)
        gps_text = tar.extractfile(gps_m).read().decode('ascii', errors='replace') if gps_m else ''
        gnv_text = tar.extractfile(gnv_m).read().decode('ascii', errors='replace') if gnv_m else ''
        tar.close()
    except Exception as e:
        print(f"  [解压错误] {e}"); return None, None

    if gps_text:
        obs = load_gps1b_text(gps_text, gf)
        pickle.dump(obs, open(cache_obs,'wb'))
        n = sum(len(v) for v in obs.values())
        print(f"  GPS1B: {len(obs)} 历元 ({n} obs), 缓存已存")
    else:
        print(f"  [错误] 无 GPS1B 数据"); return None, None

    if gnv_text:
        cache_gnv.write_text(gnv_text)
        print(f"  GNV1B: 已存")

    gnv = parse_gnv1b(str(cache_gnv)) if cache_gnv.exists() else {}
    return obs, gnv

def compute_mw(gps_obs, min_sats=6, max_ep=2000):
    """Melbourne-Wübbena 宽巷模糊度计算
    GPS1B: L1/L2 = 载波相位 (cycles), P1/P2 = 伪距 (m)
    MW = (L1 - L2)*λ_W - (P1 - P2)  [units must match]
    """
    LAM1 = C / F1   # ~0.1903 m
    LAM2 = C / F2   # ~0.2443 m
    LAM_W = C / (F1 - F2)  # ~0.8624 m (wide-lane wavelength)

    sv_vals = defaultdict(list)
    for sod in sorted(gps_obs.keys())[:max_ep]:
        if len(gps_obs[sod]) < min_sats: continue
        for sv, rec in gps_obs[sod].items():
            L1_cyc = rec.get('L1'); L2_cyc = rec.get('L2')
            P1_m = rec.get('P1'); P2_m = rec.get('P2')
            if None in (L1_cyc, L2_cyc, P1_m, P2_m): continue
            # 载波相位 (cycles) → 米
            L1_m = L1_cyc * LAM1
            L2_m = L2_cyc * LAM2
            # MW = (L1 - L2)*λ_W - (P1 - P2)  [cycles - meters/meters]
            phi_mw_cycles = (L1_cyc - L2_cyc)          # cycles
            p_mw_cycles   = (P1_m - P2_m) / LAM_W      # cycles
            b_mw = phi_mw_cycles - p_mw_cycles
            if -1.0 < b_mw < 1.0:
                sv_vals[sv].append(b_mw)

    results = {}
    for sv, vals in sv_vals.items():
        if len(vals) < 10: continue
        arr = np.array(vals)
        N_mean = float(np.nanmean(arr)); N_int = int(round(N_mean)); N_std = float(np.nanstd(arr))
        results[sv] = {'N_mean': N_mean, 'N_int': N_int, 'std': N_std, 'n': len(vals)}
        print(f"    MW {sv}: N={N_int} ({N_mean:+.4f}±{N_std:.4f} cycles, n={len(vals)})")
    return results

def apply_mw(gps_obs, mw):
    """MW校正：修正IF相位中的非整数模糊度
    
    GPS1B数据中：
    - L_if 是电离层-free组合，单位为载波相位"cycles"（L1/L2加权）
    - P_if 是电离层-free组合，单位为"m"（P1/P2加权）
    - 两者单位不一致！L_if是cycles，P_if是meters
    
    MW校正将L_if中的小数部分模糊度修正为整数：
    delta = (N_mean - N_int) * lambda_W  (lambda_W = 0.862 m)
    L_if_corr = L_if + delta  (cycles, 整数化)
    """
    corrected = {}
    for sod, sv_obs in gps_obs.items():
        corrected[sod] = {}
        for sv, rec in sv_obs.items():
            rec = dict(rec)
            if sv in mw:
                delta_cycles = mw[sv]['N_mean'] - mw[sv]['N_int']
                rec['L_if_corr'] = rec['L_if'] + delta_cycles  # 整数化L_if
                rec['mw_ok'] = True
            else:
                rec['L_if_corr'] = rec['L_if']
                rec['mw_ok'] = False
            corrected[sod][sv] = rec
    return corrected

class KFPPP:
    """Kalman Filter PPP: 状态 [X, Y, Z, clock, trop_wet]"""
    def __init__(self, pos0):
        self.x = np.array(list(pos0) + [0.0, 0.05], dtype=float)
        self.P = np.diag([1e4**2]*3 + [1e4**2] + [0.1**2])
        self.Q_pos = 1e-4**2  # m²/s
        self.Q_clk = 1e4**2   # m²/s
        self.Q_trop = 1e-8**2 # m²/s

    def predict(self, dt):
        F = np.eye(5)
        Q = np.diag([self.Q_pos*dt]*3 + [self.Q_clk*dt] + [self.Q_trop*dt])
        self.P = F @ self.P @ F.T + Q

    def update(self, obs_list):
        H_l, y_l, w_l = [], [], []
        for obs in obs_list:
            sv, sp, sv_v, Lif, Pif, el_d = obs
            el = math.radians(el_d)
            if el < 0.087: continue
            rho = math.sqrt(sum((self.x[:3][i]-sp[i])**2 for i in range(3)))
            if rho < 1e6 or rho > 1e8: continue
            e = np.array([(sp[i]-self.x[:3][i])/rho for i in range(3)])
            mf = 1.0 / max(math.sin(el), 0.05)
            # 载波相位 (MW校正后)
            mod_l = Lif - (rho + self.x[3] + self.x[4]*mf)
            w = 1.0 / (0.003**2 / (math.sin(el)**2 + 0.01))
            h = np.array([-e[0],-e[1],-e[2], 1.0, mf])
            H_l.append(h); y_l.append(mod_l); w_l.append(w)
            # 伪距
            mod_p = Pif - (rho + self.x[3] + self.x[4]*mf)
            wp = 1.0 / (0.300**2 / (math.sin(el)**2 + 0.01))
            hp = np.array([-e[0],-e[1],-e[2], 1.0, mf])
            H_l.append(hp); y_l.append(mod_p); w_l.append(wp)
        if not H_l: return False
        H = np.array(H_l); y = np.array(y_l); W = np.diag(w_l)
        try:
            PHT = self.P @ H.T
            S = H @ PHT + np.linalg.inv(W)
            K = PHT @ np.linalg.inv(S)
            self.x = self.x + K @ y
            self.P = (np.eye(5) - K @ H) @ self.P
            return True
        except: return False

def process_day(year, month, day, data_dir='./data',
                nhours=4.0, interval=30.0, gf='C', sp3_prod='FIN'):
    ds = f"{year:04d}-{month:02d}-{day:02d}"
    print(f"\n[{ds}] 处理中...")
    gps_obs, ref_orbit = download_gracefo(year, month, day, data_dir, gf)
    if not gps_obs or not ref_orbit: return [], {}

    # MW 模糊度
    mw = compute_mw(gps_obs)
    if not mw: print(f"  [警告] 无MW结果"); return [], {}
    gps_corr = apply_mw(gps_obs, mw)
    n_mw = sum(1 for sod in gps_corr for sv in gps_corr[sod] if gps_corr[sod][sv].get('mw_ok'))
    print(f"  MW应用到 {n_mw} 条")

    # IGS SP3
    doy = (datetime(year,month,day) - datetime(year,1,1)).days + 1
    print(f"  IGS {sp3_prod} SP3 DOY {doy}...")
    sp3 = load_igs_sp3(year, doy, data_dir, product=sp3_prod)
    if sp3: print(f"  SP3: {sp3.get('source','?')[:60]}")
    else: print(f"  [警告] SP3 加载失败")

    # 构建PPP记录
    t_start = datetime(year, month, day, 0, 0, 0)
    t_end = t_start + timedelta(hours=nhours)
    orbit_ts = sorted(ref_orbit.keys())
    records = []

    # λ_if: effective wavelength of ionospheric-free phase combination
    # λ_if = c / f_if where f_if = (f1² - f2²)/(f1 + f2) ≈ 2799.02 MHz
    # GPS L1=1575.42, L2=1227.60 → f_if ≈ 2799.02e6
    F1, F2 = 1575.42e6, 1227.60e6
    F_IF = (F1*F1 - F2*F2) / (F1 + F2)
    LAM_IF = C / F_IF   # ~0.1070 m

    for sod, sv_obs in sorted(gps_corr.items()):
        utc = gps_sod_to_utc(sod)
        if not (t_start <= utc <= t_end): continue
        dt_s = (utc - t_start).total_seconds()
        if abs(dt_s - round(dt_s/interval)*interval) > 2.0: continue

        # GRACE-FO位置
        t0 = t1 = None
        for j, ti in enumerate(orbit_ts):
            if ti >= utc: t1 = ti; t0 = orbit_ts[j-1] if j > 0 else None; break
            t0 = ti
        if t1 is None: t0 = t1 = orbit_ts[-1]
        if t0 is None: t0 = orbit_ts[0]
        dt0 = (utc - t0).total_seconds(); dtt = (t1 - t0).total_seconds()
        gr = ref_orbit[t0]*(1-dt0/dtt) + ref_orbit[t1]*(dt0/dtt) if dtt else ref_orbit[t0]

        lat, lon, _ = ecef_to_blh(gr)
        R = ecef_to_enu_matrix(lat, lon)

        for sv, rec in sv_obs.items():
            if sp3:
                sp, clk, sv_v = get_gps_pos_from_sp3(sp3, sv, utc)
            else: continue
            if sp is None: continue
            d = gr - sp; rng = float(np.linalg.norm(d))
            if not (2e7 < rng < 5e7): continue
            try:
                e_enu = R @ (d / rng)
                el = float(np.arcsin(np.clip(e_enu[2], -1.0, 1.0)))
                az = float(np.arctan2(e_enu[0], e_enu[1]))
                if az < 0: az += 2*math.pi
            except: continue
            if el < 0.087: continue
            records.append({
                'time': utc, 'sv': sv,
                'sp': sp.astype(float),
                'L_if': float(rec['L_if_corr']) * LAM_IF,  # cycles → meters
                'P_if': float(rec['P_if']),                # meters
                'el': float(np.degrees(el)),
            })

    print(f"  PPP records: {len(records)}")
    if len(records) < 30: return [], {}

    # KF PPP
    kf = KFPPP(list(ref_orbit[orbit_ts[0]]))
    results = []; prev_t = None
    records.sort(key=lambda r: r['time'].timestamp())
    groups = [(t, list(g)) for t, g in groupby(records,
               key=lambda r: round(r['time'].timestamp()/interval)*interval)]

    for t_key, grp in groups:
        dt = datetime.fromtimestamp(t_key)
        if not (t_start <= dt <= t_end): continue
        if prev_t:
            kf.predict(max((dt-prev_t).total_seconds(), 0))
        obs_list = [(r['sv'], r['sp'], None, r['L_if'], r['P_if'], r['el']) for r in grp]
        kf.update(obs_list)
        pos = kf.x[:3]

        # 参考位置
        rp0 = rp1 = None
        for j, ti in enumerate(orbit_ts):
            if ti >= dt: rp1 = ti; rp0 = orbit_ts[j-1] if j > 0 else None; break
            rp0 = ti
        if rp1 is None: rp0 = rp1 = orbit_ts[-1]
        if rp0 is None: rp0 = orbit_ts[0]
        dt0 = (dt-rp0).total_seconds(); dtt = (rp1-rp0).total_seconds()
        ref = ref_orbit[rp0]*(1-dt0/dtt) + ref_orbit[rp1]*(dt0/dtt) if dtt else ref_orbit[rp0]
        err = pos - ref
        enu = R @ err
        d3 = float(np.linalg.norm(err))
        results.append({'time': dt, 'dE': enu[0], 'dN': enu[1], 'dU': enu[2],
                        'd3': d3, 'n': len(grp)})
        prev_t = dt

    print(f"  KF: {len(results)} 历元收敛")
    return results, {'mw': mw, 'n_obs': len(records), 'n_mw': n_mw}

def rms(v):
    a = np.array(v, dtype=float)
    return float(np.sqrt(np.nanmean(a**2)))

def compute_stats(results):
    if not results: return None
    dE = np.array([r['dE']*100 for r in results])
    dN = np.array([r['dN']*100 for r in results])
    dU = np.array([r['dU']*100 for r in results])
    d3 = np.array([r['d3']*100 for r in results])
    good = d3 < 100
    def rg(v): return rms([v[i] for i in range(len(v)) if good[i]])
    times = [(r['time']-results[0]['time']).total_seconds()/3600 for r in results]
    return {
        'RMS_E': rg(dE), 'RMS_N': rg(dN), 'RMS_U': rg(dU), 'RMS_3D': rg(d3),
        'RMS_E_all': rms(dE), 'RMS_N_all': rms(dN), 'RMS_U_all': rms(dU), 'RMS_3D_all': rms(d3),
        'MAX_3D': float(np.nanmax(d3)),
        'n_total': len(results), 'n_good': int(good.sum()),
        'rate': float(good.sum())/len(results)*100,
        '_dE': dE.tolist(), '_dN': dN.tolist(), '_dU': dU.tolist(),
        '_d3': d3.tolist(), '_times': times,
    }

def make_svg(st, ds, out_path):
    W,H = 960,380; ML,MR,MT,MB = 72,20,40,48
    PW = W-ML-MR; ph = (H-MT-MB-24)//3; gap = 8
    dE,dN,dU,d3 = st['_dE'],st['_dN'],st['_dU'],st['_d3']
    times = st['_times']
    t0,t1 = times[0],times[-1] if times else 1
    ymx = max(0.5, np.nanpercentile(np.abs(np.array(d3)),95)) * 1.4

    def px(tv): return ML + (tv-t0)/(t1-t0+1e-9)*(PW-1)
    def py(val, top):
        v = float(val)
        if np.isnan(v): v = 0.0
        return top+ph-1-int((v+ymx)/(2*ymx)*(ph-1))
    def p0(top): return py(0, top)
    rows = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}">',
            f'<title>GRACE-FO PPP {ds}</title>',
            f'<rect width="{W}" height="{H}" fill="white"/>']

    for i,(data,col,label,rv) in enumerate([(dE,'#1565C0','E (东西)',st['RMS_E']),
                                             (dN,'#2E7D32','N (南北)',st['RMS_N']),
                                             (dU,'#C62828','U (垂直)',st['RMS_U'])]):
        top = MT+i*(ph+gap)
        rows += [f'<rect x="{ML}" y="{top}" width="{PW}" height="{ph}" fill="#f7f7f7"/>',
                 f'<line x1="{ML}" y1="{p0(top)}" x2="{ML+PW}" y2="{p0(top)}" stroke="#bbb"/>']
        pts = ' '.join(f'{px(times[j])},{py(data[j],top)}' for j in range(len(times)))
        rows.append(f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="1.5"/>')
        rows.append(f'<text x="{ML+4}" y="{top+14}" font-size="12" fill="#333" font-weight="bold">{label}</text>')
        rows.append(f'<text x="{ML+PW+4}" y="{top+14}" font-size="11" fill="{col}">RMS={rv:.2f} cm</text>')

    top = MT+3*(ph+gap)
    rows += [f'<rect x="{ML}" y="{top}" width="{PW}" height="{ph}" fill="#f7f7f7"/>',
             f'<line x1="{ML}" y1="{p0(top)}" x2="{ML+PW}" y2="{p0(top)}" stroke="#bbb"/>']
    for data,col,ls in [(dE,'#1565C0',''),(dN,'#2E7D32',''),(dU,'#C62828',''),(d3,'#333','4,4')]:
        pts = ' '.join(f'{px(times[j])},{py(data[j],top)}' for j in range(len(times)))
        rows.append(f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="2" stroke-dasharray="{ls}"/>')
    rows.append(f'<text x="{ML+4}" y="{top+14}" font-size="12" fill="#333" font-weight="bold">3D 位置误差</text>')
    rows.append(f'<text x="{ML+PW+4}" y="{top+14}" font-size="11" fill="#333">3D RMS={st["RMS_3D"]:.2f} cm</text>')
    rows.append(f'<text x="{W//2}" y="24" text-anchor="middle" font-size="13" font-weight="bold">GRACE-FO PPP v1.2.0 — {ds} (GPS1B+IGS Final+MW+KF)</text>')
    rows.append(f'<text x="{W//2}" y="38" text-anchor="middle" font-size="11" fill="#666">收敛 {st["n_good"]}/{st["n_total"]} ({st["rate"]:.1f}%) | E:{st["RMS_E"]:.2f} N:{st["RMS_N"]:.2f} U:{st["RMS_U"]:.2f} 3D:{st["RMS_3D"]:.2f} cm</text>')
    rows.append(f'<text x="{W//2}" y="{H-6}" text-anchor="middle" font-size="10" fill="#888">时间 (小时)</text>')
    rows.append('</svg>')
    Path(out_path).write_text('\n'.join(rows))

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description='GRACE-FO PPP v1.2.0')
    p.add_argument('--start', default='2024-05-01')
    p.add_argument('--end',   default='2024-06-30')
    p.add_argument('--data-dir', default='./data')
    p.add_argument('--out-dir',  default='./output_v12')
    p.add_argument('--nhours', type=float, default=4.0)
    p.add_argument('--interval', type=float, default=30.0)
    p.add_argument('--sp3', default='FIN')
    p.add_argument('--grace', default='C')
    p.add_argument('--redownload', type=int, default=0,
                   help='1=强制重新下载（忽略缓存）, 0=读缓存')
    args = p.parse_args()

    out_dir = Path(args.out_dir); out_dir.mkdir(exist_ok=True)
    t_start = datetime.fromisoformat(args.start)
    t_end   = datetime.fromisoformat(args.end)
    total_days = (t_end - t_start).days + 1

    print(f"{'='*56}")
    print(f"GRACE-FO PPP v1.2.0 批量处理")
    print(f"范围: {args.start} → {args.end} ({total_days} 天)")
    print(f"数据: GPS1B (ISDC) + IGS {args.sp3} SP3")
    print(f"算法: MW模糊度固定 + Kalman Filter PPP")
    print(f"{'='*56}")

    all_stats = {}
    cur = t_start
    while cur <= t_end:
        ds = cur.strftime('%Y-%m-%d')
        y,m,d = cur.year, cur.month, cur.day
        print(f"\n[{ds} ({(cur-t_start).days+1}/{total_days})]", flush=True)

        csv_f = out_dir/f"ppp_{ds}_{args.nhours:.0f}h.csv"
        svg_f = out_dir/f"ppp_{ds}_{args.nhours:.0f}h.svg"
        json_f = out_dir/f"ppp_{ds}_{args.nhours:.0f}h.json"

        # 跳过已完成
        if json_f.exists() and csv_f.exists() and not args.redownload:
            try:
                s = json.load(open(json_f))
                all_stats[ds] = s
                print(f"  [跳过] 已有: 3D={s.get('RMS_3D',0):.2f} cm ({s.get('rate',0):.1f}%)")
                cur += timedelta(days=1); continue
            except: pass

        res, meta = process_day(y,m,d, args.data_dir, args.nhours, args.interval, args.grace, args.sp3)
        if not res:
            print(f"  [无结果]"); cur += timedelta(days=1); continue

        st = compute_stats(res)
        # CSV
        with open(csv_f,'w') as f:
            f.write("time,dE_m,dN_m,dU_m,d3_m,n_sat\n")
            for r in res:
                f.write(f"{r['time'].isoformat()},{r['dE']:.4f},{r['dN']:.4f},{r['dU']:.4f},{r['d3']:.4f},{r['n']}\n")
        # SVG
        make_svg(st, ds, str(svg_f))
        # JSON
        st_ser = {k: (v.tolist() if hasattr(v,'tolist') else v)
                   for k,v in st.items() if not k.startswith('_')}
        st_ser['n_mw'] = meta.get('n_mw',0)
        json.dump(st_ser, open(json_f,'w'), indent=2, default=str)
        all_stats[ds] = st_ser

        print(f"  ✓ {ds}: E={st['RMS_E']:.2f} N={st['RMS_N']:.2f} U={st['RMS_U']:.2f} 3D={st['RMS_3D']:.2f} cm ({st['n_good']}/{st['n_total']} {st['rate']:.1f}%)")
        cur += timedelta(days=1)

    # 汇总
    if not all_stats:
        print("无结果"); sys.exit(0)
    all_e = [all_stats[d]['RMS_E'] for d in sorted(all_stats)]
    all_n = [all_stats[d]['RMS_N'] for d in sorted(all_stats)]
    all_u = [all_stats[d]['RMS_U'] for d in sorted(all_stats)]
    all_3 = [all_stats[d]['RMS_3D'] for d in sorted(all_stats)]
    print(f"\n{'='*56}")
    print(f"全时期统计 ({len(all_stats)}/{total_days} 天成功)")
    print(f"  E RMS: {np.mean(all_e):.2f} ± {np.std(all_e):.2f} cm  [{min(all_e):.2f}~{max(all_e):.2f}]")
    print(f"  N RMS: {np.mean(all_n):.2f} ± {np.std(all_n):.2f} cm  [{min(all_n):.2f}~{max(all_n):.2f}]")
    print(f"  U RMS: {np.mean(all_u):.2f} ± {np.std(all_u):.2f} cm  [{min(all_u):.2f}~{max(all_u):.2f}]")
    print(f"  3D RMS: {np.mean(all_3):.2f} ± {np.std(all_3):.2f} cm  [{min(all_3):.2f}~{max(all_3):.2f}]")

    # 汇总CSV
    sum_csv = out_dir/"summary_all_days.csv"
    with open(sum_csv,'w') as f:
        f.write("date,RMS_E_cm,RMS_N_cm,RMS_U_cm,RMS_3D_cm,n_good,n_total,rate_pct\n")
        for ds in sorted(all_stats):
            r = all_stats[ds]
            f.write(f"{ds},{r['RMS_E']:.3f},{r['RMS_N']:.3f},{r['RMS_U']:.3f},{r['RMS_3D']:.3f},{r['n_good']},{r['n_total']},{r['rate']:.1f}\n")
    json.dump(all_stats, open(out_dir/'summary.json','w'), indent=2, default=str)
    print(f"\n结果: {out_dir}")
    print(f"汇总: {sum_csv}")
