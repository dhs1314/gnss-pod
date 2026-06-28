#!/usr/bin/env python3
"""
GRACE-FO PPP 全月验证 — WHU IGS ultra-rapid SP3 精密星历
用法: python3 ppp_may2024_sp3.py
"""
import sys, os, csv, io, ssl, urllib.request, tarfile, datetime as dt, pickle
from pathlib import Path
from itertools import groupby
import numpy as np

# 复用 run_ppp.py 的所有函数（已含 --ephem sp3 支持）
import importlib.util
spec = importlib.util.spec_from_file_location('rp', '/workspace/gnss_pod/run_ppp.py')
rp = importlib.util.module_from_spec(spec)
sys.modules['rp_module'] = rp
spec.loader.exec_module(rp)

C = 299792458.0; MU_E = 3.986004418e14; OMEGA_E = 7.2921151467e-5
ISDC = "https://isdc-data.gfz.de/grace-fo/Level-1B/JPL/INSTRUMENT/RL04/{year}/"

def download_gnv1b(y, m, d):
    ds = f"{y:04d}-{m:02d}-{d:02d}"
    od = Path("data")/"gracefo"/str(y)/ds; od.mkdir(parents=True, exist_ok=True)
    gnv = od/f"GNV1B_{ds}_C_04.txt"
    if gnv.exists() and gnv.stat().st_size > 1000:
        print(f"  [缓存] {ds}"); return str(gnv)
    fname = f"gracefo_1B_{ds}_RL04.ascii.noLRI.tgz"
    url = ISDC.format(year=y) + fname
    ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'curl/7.88', 'Accept-Encoding': 'gzip'})
        with urllib.request.urlopen(req, timeout=180, context=ctx) as r: data = r.read()
        print(f"  GNV1B: {len(data)//1024}KB")
    except Exception as e:
        print(f"  错误: {e}"); return None
    try:
        tar = tarfile.open(fileobj=io.BytesIO(data), mode='r:*')
        for m in tar.getmembers():
            if 'GNV1B' in m.name and m.name.endswith('.txt'):
                fo = od/Path(m.name).name; f = tar.extractfile(m)
                if f: fo.write_bytes(f.read())
        tar.close()
    except: return None
    return str(gnv) if gnv.exists() else None

def parse_gnv1b(path):
    orbit = {}; gps0 = dt.datetime(2000, 1, 1, 12, 0, 0)
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            p = line.split()
            if len(p) < 6: continue
            try:
                tg = float(p[0]); flag = p[2]
                if flag not in ('C', 'E'): continue
                X, Y, Z = float(p[3]), float(p[4]), float(p[5])
                if abs(X) < 1e3: continue
                orbit[gps0 + dt.timedelta(seconds=tg)] = np.array([X, Y, Z])
            except: continue
    return orbit

def run_day(y, m, d, nhours=4.0, interval=30.0):
    t_start = dt.datetime(y, m, d, 0, 0, 0)
    t_end = t_start + dt.timedelta(hours=nhours)
    ds = f"{y:04d}-{m:02d}-{d:02d}"
    print(f"\n[{dt.datetime.now().strftime('%H:%M:%S')}] === {ds} [SP3] ===")
    gnv = download_gnv1b(y, m, d)
    if not gnv: return None
    ref_orbit = parse_gnv1b(gnv)
    print(f"  轨道: {len(ref_orbit)} epochs")
    doy = (t_start - dt.datetime(t_start.year, 1, 1)).days + 1
    print(f"  加载 WHU SP3 (DOY {doy})...")
    sp3_data = rp.load_whu_sp3_ultra_rapid(y, doy, "data")
    if sp3_data is None:
        print("  [错误] SP3 加载失败"); return None
    records = rp.generate_obs_from_orbit(ref_orbit, t_start, nhours, interval,
                                         ephem='sp3', sp3_data=sp3_data)
    print(f"  观测: {len(records)} obs")
    if not records: return None
    records.sort(key=lambda r: r['time'].timestamp())
    groups = [(t, list(g)) for t, g in groupby(records, key=lambda r: r['time'].timestamp())]
    t0 = min(ref_orbit.keys())
    x0 = np.concatenate([ref_orbit[t0], [0.0, 0.2]])
    results = []
    for i, (t, grp) in enumerate(groups):
        dt_obj = dt.datetime.fromtimestamp(t)
        if not (t_start <= dt_obj <= t_end): continue
        obs_list = [
            (r['sv'], r['sat_pos'], r['sat_vel'],
             r['L1'], r['L2'], r['P1'], r['P2'],
             np.radians(r['el']), np.radians(r['az']))
            for r in grp
        ]
        if len(obs_list) < 4: continue
        x = rp.ppp_single_epoch(obs_list, x0)
        pos_est = x[:3]
        ts = sorted(ref_orbit.keys()); ref_p0 = ref_p1 = None
        for j, ti in enumerate(ts):
            if ti >= dt_obj: ref_p1 = ti; ref_p0 = ts[j-1] if j > 0 else None; break
            ref_p0 = ti
        if ref_p1 is None: ref_p0 = ts[-1]; ref_p1 = ts[-1]
        if ref_p0 is None: ref_p0 = ts[0]
        ref_pos = ref_orbit[ref_p0]
        if ref_p0 != ref_p1:
            dt0 = (dt_obj - ref_p0).total_seconds()
            dt_tot = (ref_p1 - ref_p0).total_seconds()
            if dt_tot > 0:
                a = dt0 / dt_tot
                ref_pos = ref_orbit[ref_p0]*(1-a) + ref_orbit[ref_p1]*a
        err = pos_est - ref_pos
        lat, lon, _ = rp.ecef_to_blh(ref_pos)
        R = rp.ecef_to_enu_matrix(lat, lon)
        enu = R @ err
        results.append({
            'time': dt_obj,
            'dE': float(enu[0]), 'dN': float(enu[1]), 'dU': float(enu[2]),
            'n_sat': len(obs_list)
        })
        x0 = x.copy()
    print(f"  PPP: {len(results)}/{len(groups)} converged")
    if not results: return None
    dE = np.array([float(r['dE']) for r in results])
    dN = np.array([float(r['dN']) for r in results])
    dU = np.array([float(r['dU']) for r in results])
    d3 = np.sqrt(dE**2 + dN**2 + dU**2)
    ok = d3 < 1.0
    dE_ok, dN_ok, dU_ok = dE[ok], dN[ok], dU[ok]
    if len(dE_ok) == 0: return None
    re = float(np.sqrt(np.nanmean(dE_ok**2)) * 100)
    rn = float(np.sqrt(np.nanmean(dN_ok**2)) * 100)
    ru = float(np.sqrt(np.nanmean(dU_ok**2)) * 100)
    r3 = float(np.sqrt(np.nanmean(dE_ok**2 + dN_ok**2 + dU_ok**2)) * 100)
    print(f"  E={re:.2f}cm  N={rn:.2f}cm  U={ru:.2f}cm  3D={r3:.2f}cm  ({ok.sum()}/{len(results)} epochs)")
    csv_path = f"output/ppp_may_sp3_{ds}_{nhours:.0f}h.csv"
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['time','dE_m','dN_m','dU_m','d3_m','n_sat'])
        w.writeheader()
        for r in results:
            d3v = np.sqrt(float(r['dE'])**2+float(r['dN'])**2+float(r['dU'])**2)
            w.writerow({'time': r['time'].isoformat(), 'dE_m': f"{r['dE']:.4f}",
                        'dN_m': f"{r['dN']:.4f}", 'dU_m': f"{r['dU']:.4f}",
                        'd3_m': f"{d3v:.4f}", 'n_sat': r['n_sat']})
    return {'date': ds, 'csv': csv_path, 're': re, 'rn': rn, 'ru': ru, 'r3': r3,
            'n_ok': int(ok.sum()), 'n_total': len(results),
            'dE': list(dE_ok), 'dN': list(dN_ok), 'dU': list(dU_ok)}

def make_svg(dE, dN, dU, n, re, rn, ru, r3, label, W=900, H=520):
    times = [i * 30.0 / 3600.0 for i in range(n)]
    tmin, tmax = 0, max(times) if times else 1
    ymx = max(max(abs(v) for v in dE), max(abs(v) for v in dN), max(abs(v) for v in dU)) * 1.25
    ymx = max(ymx, 0.5)
    ML, MR, MT, MB, gl = 60, 20, 20, 36, 8
    pw = W - ML - MR
    ph = (H - MT - MB - 2*gl) // 3
    def px(tv): return ML + (tv - tmin) / (tmax - tmin + 1e-9) * (pw - 1)
    def py(tv, top): return top + ph - 1 - int((tv + ymx) / (2*ymx) * (ph - 1))
    def p0(top): return py(0.0, top)
    rows = [(dE, '#1565C0', 'E (东西)'), (dN, '#2E7D32', 'N (南北)'), (dU, '#C62828', 'U (垂直)')]
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" font-family="DejaVu Sans,sans-serif" font-size="11"><rect width="{W}" height="{H}" fill="white"/>']
    for i, (data, col, lbl) in enumerate(rows):
        ry = MT + gl + i * (ph + gl)
        svg.append(f'<rect x="{ML}" y="{ry}" width="{pw}" height="{ph}" fill="#f8f9fa"/>')
        for j in range(0, n, max(1, n//30)):
            xi = int(px(times[j]))
            if ML <= xi < ML + pw:
                svg.append(f'<line x1="{xi}" y1="{ry}" x2="{xi}" y2="{ry+ph}" stroke="#e0e0e0" stroke-width="0.5"/>')
        yz = p0(ry)
        svg.append(f'<line x1="{ML}" y1="{yz}" x2="{ML+pw}" y2="{yz}" stroke="#999" stroke-width="0.8"/>')
        pts = ' '.join(f'{px(times[j]):.1f},{py(data[j], ry):.1f}' for j in range(n) if -ymx <= data[j] <= ymx)
        if pts: svg.append(f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="1.2"/>')
        rv = {'E (东西)': re, 'N (南北)': rn, 'U (垂直)': ru}[lbl]
        svg.append(f'<text x="{ML+4}" y="{ry+14}" fill="{col}" font-weight="bold">{lbl}</text>')
        svg.append(f'<text x="{ML+pw-4}" y="{ry+14}" text-anchor="end" fill="#555">RMS={rv:.2f} cm</text>')
        svg.append(f'<rect x="{ML}" y="{ry}" width="{pw}" height="{ph}" fill="none" stroke="#ccc" stroke-width="1"/>')
    svg.append(f'<text x="{W//2}" y="{H-8}" text-anchor="middle" font-size="12" fill="#333">{label} [SP3] | 3D RMS={r3:.2f} cm ({n} epochs)</text>')
    svg.append('</svg>')
    return '\n'.join(svg)

CSS = """*{box-sizing:border-box;margin:0;padding:0}
body{font-family:DejaVu Sans,Arial,sans-serif;background:#f4f6fa;padding:16px}
.hdr{text-align:center;padding:24px;background:linear-gradient(135deg,#1a237e,#6A1B9A);color:white;border-radius:12px;margin-bottom:20px}
.hdr h1{font-size:22px;margin-bottom:6px}
.hdr p{font-size:13px;opacity:0.85}
.day-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(480px,1fr));gap:12px}
.card{background:white;border-radius:12px;padding:20px;box-shadow:0 2px 8px rgba(0,0,0,0.06)}
.card-hdr{font-size:15px;font-weight:bold;color:#6A1B9A;margin-bottom:10px;display:flex;justify-content:space-between;align-items:center}
.stats{display:grid;grid-template-columns:repeat(4,1fr);gap:8px;margin:10px 0}
.s{background:#f8f9fa;border-radius:6px;padding:10px;text-align:center}
.s .v{font-size:18px;font-weight:bold}
.s .l{font-size:10px;color:#888;margin-top:2px}
.e .v{color:#1565C0}.n .v{color:#2E7D32}.u .v{color:#C62828}.a .v{color:#6A1B9A}
.note{font-size:11px;color:#888;text-align:center;margin-top:6px}
.sum{background:white;border-radius:12px;padding:20px;margin-top:20px;box-shadow:0 2px 8px rgba(0,0,0,0.06)}
.sum h2{font-size:16px;color:#6A1B9A;margin-bottom:16px}
.sum-grid{display:grid;grid-template-columns:repeat(7,1fr);gap:8px}
.su{background:#f8f9fa;border-radius:6px;padding:8px;text-align:center;font-size:12px}
.su .v{font-size:16px;font-weight:bold}
.su .l{font-size:10px;color:#888}"""

def main():
    days = [(2024, 5, d) for d in range(1, 32)]
    all_results = []
    failed = []
    for y, m, d in days:
        r = run_day(y, m, d, nhours=4.0, interval=30.0)
        if r: all_results.append(r)
        else: failed.append(f"{y:04d}-{m:02d}-{d:02d}")

    print("\n" + "="*60)
    print("MAY 2024 — WHU ultra-rapid SP3 精密星历验证结果")
    print("="*60)
    print(f"{'日期':<12} {'E(cm)':>8} {'N(cm)':>8} {'U(cm)':>8} {'3D(cm)':>8} {'有效历元':>10}")
    print("-"*60)
    for r in all_results:
        print(f"{r['date']:<12} {r['re']:>8.2f} {r['rn']:>8.2f} {r['ru']:>8.2f} {r['r3']:>8.2f} {r['n_ok']:>6}/{r['n_total']:>3}")
    if failed: print(f"\n失败: {', '.join(failed)}")

    re_all = [r['re'] for r in all_results]
    rn_all = [r['rn'] for r in all_results]
    ru_all = [r['ru'] for r in all_results]
    r3_all = [r['r3'] for r in all_results]
    n_ok_all = sum(r['n_ok'] for r in all_results)
    n_tot_all = sum(r['n_total'] for r in all_results)
    print(f"\n月均值: E={np.mean(re_all):.2f}±{np.std(re_all):.2f}cm  N={np.mean(rn_all):.2f}±{np.std(rn_all):.2f}cm  U={np.mean(ru_all):.2f}±{np.std(ru_all):.2f}cm  3D={np.mean(r3_all):.2f}±{np.std(r3_all):.2f}cm")
    print(f"3D范围: {min(r3_all):.2f} ~ {max(r3_all):.2f} cm")
    print(f"成功率: {n_ok_all}/{n_tot_all} ({n_ok_all/n_tot_all*100:.1f}%)")

    html = ['<!DOCTYPE html><html lang="zh"><head>',
            '<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">',
            '<title>GRACE-FO PPP 2024年5月 — WHU SP3 精密星历验证报告</title>',
            '<style>' + CSS + '</style></head>',
            '<body><div class="hdr">',
            '<h1>GRACE-FO PPP — WHU ultra-rapid SP3 精密星历</h1>',
            '<p>GPS星历: WHU IGS ultra-rapid (~5 cm精度) | PPP vs GNV1B | 采样 30s | 每段 4 小时 | 过滤 3D 误差 &lt;1m</p>',
            '</div><div class="day-grid">']
    for r in all_results:
        svg = make_svg(r['dE'], r['dN'], r['dU'], r['n_ok'], r['re'], r['rn'], r['ru'], r['r3'], r['date'])
        html.append(f'<div class="card">'
                    f'<div class="card-hdr"><span>{r["date"]} [SP3]</span><span style="font-size:11px;color:#888">{r["n_ok"]}/{r["n_total"]} epochs</span></div>'
                    f'<div class="stats">'
                    f'<div class="s e"><div class="v">{r["re"]:.2f}</div><div class="l">E cm</div></div>'
                    f'<div class="s n"><div class="v">{r["rn"]:.2f}</div><div class="l">N cm</div></div>'
                    f'<div class="s u"><div class="v">{r["ru"]:.2f}</div><div class="l">U cm</div></div>'
                    f'<div class="s a"><div class="v">{r["r3"]:.2f}</div><div class="l">3D cm</div></div>'
                    f'</div>' + svg +
                    f'<div class="note">X轴=时间（h） Y轴=误差（cm） 蓝=E 绿=N 红=U | 参考: ISDC/GFZ GNV1B</div>'
                    f'</div>')
    html.append('</div>')
    html.append(f'<div class="sum"><h2>月统计 ({len(all_results)}天)</h2><div class="sum-grid">'
               f'<div class="su"><div class="v" style="color:#1565C0">{np.mean(re_all):.2f}±{np.std(re_all):.2f}</div><div>E (cm)</div></div>'
               f'<div class="su"><div class="v" style="color:#2E7D32">{np.mean(rn_all):.2f}±{np.std(rn_all):.2f}</div><div>N (cm)</div></div>'
               f'<div class="su"><div class="v" style="color:#C62828">{np.mean(ru_all):.2f}±{np.std(ru_all):.2f}</div><div>U (cm)</div></div>'
               f'<div class="su"><div class="v" style="color:#6A1B9A">{np.mean(r3_all):.2f}±{np.std(r3_all):.2f}</div><div>3D (cm)</div></div>'
               f'<div class="su"><div class="v">{min(r3_all):.1f}~{max(r3_all):.1f}</div><div>3D范围</div></div>'
               f'<div class="su"><div class="v">{n_ok_all/n_tot_all*100:.0f}%</div><div>成功率</div></div>'
               f'<div class="su"><div class="v">{len(failed)}</div><div>失败天</div></div>'
               f'</div></div>')
    html.append('</body></html>')
    out = 'output/report_may2024_sp3.html'
    Path(out).write_text(''.join(html))
    print(f"\n报告: {out} ({os.path.getsize(out)//1024}KB, {len(all_results)}天成功)")

if __name__ == '__main__':
    main()