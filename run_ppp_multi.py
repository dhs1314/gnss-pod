#!/usr/bin/env python3
"""
run_ppp_multi.py — 调用 run_ppp.py 验证过的函数，处理 Apr 29 - May 2 四天数据
用法: cd /workspace/gnss_pod && python3 run_ppp_multi.py
"""
import sys, os, ssl, tarfile, io
from pathlib import Path
from datetime import datetime as dt, timedelta
from itertools import groupby
import numpy as np

# 动态加载 run_ppp.py（避免 its __main__ 解析 CLI 参数）
import importlib.util
spec = importlib.util.spec_from_file_location('rp', '/workspace/gnss_pod/run_ppp.py')
rp = importlib.util.module_from_spec(spec)
sys.modules['rp_module'] = rp
spec.loader.exec_module(rp)

ISDC = "https://isdc-data.gfz.de/grace-fo/Level-1B/JPL/INSTRUMENT/RL04/{year}/"

def download_gnv1b(y, m, d):
    ds = f"{y:04d}-{m:02d}-{d:02d}"
    od = Path("data")/"gracefo"/str(y)/ds; od.mkdir(parents=True, exist_ok=True)
    gnv = od/f"GNV1B_{ds}_C_04.txt"
    if gnv.exists() and gnv.stat().st_size > 1000:
        print(f"  [缓存] {ds}"); return str(gnv)
    fname = f"gracefo_1B_{ds}_RL04.ascii.noLRI.tgz"
    url = ISDC.format(year=y) + fname
    print(f"  下载 {fname} ...")
    ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'curl/7.88', 'Accept-Encoding': 'gzip'})
        with urllib.request.urlopen(req, timeout=120, context=ctx) as r:
            data = r.read()
        print(f"  完成: {len(data)//1024}KB")
    except Exception as e:
        print(f"  错误: {e}"); return None
    try:
        tar = tarfile.open(fileobj=io.BytesIO(data), mode='r:*')
        for m in tar.getmembers():
            if 'GNV1B' in m.name and m.name.endswith('.txt'):
                fo = od/Path(m.name).name; f = tar.extractfile(m)
                if f: fo.write_bytes(f.read()); print(f"  解压: {fo.name}")
        tar.close()
    except Exception as e:
        print(f"  解压错误: {e}"); return None
    return str(gnv) if gnv.exists() else None

def parse_gnv1b(path):
    orbit = {}; gps0 = dt(2000, 1, 1, 12, 0, 0)
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            p = line.split()
            if len(p) < 6: continue
            try:
                tg = float(p[0]); flag = p[2]
                if flag not in ('C', 'E'): continue
                X, Y, Z = float(p[3]), float(p[4]), float(p[5])
                if abs(X) < 1e3: continue
                orbit[gps0 + timedelta(seconds=tg)] = np.array([X, Y, Z])
            except: continue
    return orbit

def run_day(y, m, d, nhours=4.0, interval=30.0):
    t_start = dt(y, m, d, 0, 0, 0)
    t_end = t_start + timedelta(hours=nhours)
    ds = f"{y:04d}-{m:02d}-{d:02d}"
    print(f"\n[{dt.now().strftime('%H:%M:%S')}] === {ds} ({nhours:.0f}h, {interval}s) ===")
    gnv = download_gnv1b(y, m, d)
    if not gnv: return None
    ref_orbit = parse_gnv1b(gnv)
    print(f"  轨道: {len(ref_orbit)} 历元")
    # 生成观测（调用 run_ppp.py 中验证过的函数）
    records = rp.generate_obs_from_orbit(ref_orbit, t_start, nhours, interval)
    print(f"  观测: {len(records)} 条")
    if not records:
        print("  [错误] 无观测"); return None
    # 按时间分组
    records.sort(key=lambda r: r['time'].timestamp())
    groups = [(t, list(g)) for t, g in groupby(records, key=lambda r: r['time'].timestamp())]
    print(f"  历元: {len(groups)}")
    # 初始位置
    t0 = min(ref_orbit.keys())
    x0 = np.concatenate([ref_orbit[t0], [0.0, 0.2]])  # [X,Y,Z,clock,trop]
    results = []
    for i, (t, grp) in enumerate(groups):
        dt_obj = dt.fromtimestamp(t)
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
        # 参考位置
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
    print(f"  PPP: {len(results)}/{len(groups)} 收敛")
    if not results: return None
    dE = np.array([r['dE'] for r in results])
    dN = np.array([r['dN'] for r in results])
    dU = np.array([r['dU'] for r in results])
    rms_e = float(np.sqrt(np.nanmean(dE**2)) * 100)
    rms_n = float(np.sqrt(np.nanmean(dN**2)) * 100)
    rms_u = float(np.sqrt(np.nanmean(dU**2)) * 100)
    rms_3 = float(np.sqrt(np.nanmean(dE**2 + dN**2 + dU**2)) * 100)
    print(f"  E={rms_e:.2f}cm  N={rms_n:.2f}cm  U={rms_u:.2f}cm  3D={rms_3:.2f}cm")
    # CSV
    csv_path = f"output/ppp_multi_{ds}_{nhours:.0f}h.csv"
    with open(csv_path, 'w') as f:
        f.write("time,dE_m,dN_m,dU_m,d3_m,n_sat\n")
        for r in results:
            f.write(f"{r['time'].isoformat()},{r['dE']:.4f},{r['dN']:.4f},{r['dU']:.4f},"
                    f"{np.sqrt(r['dE']**2+r['dN']**2+r['dU']**2):.4f},{r['n_sat']}\n")
    print(f"  CSV: {csv_path}")
    return {'date': ds, 'csv': csv_path, 're': rms_e, 'rn': rms_n,
            'ru': rms_u, 'r3': rms_3, 'n': len(results), 'results': results}

def make_svg(results, label, W=900, H=520):
    n = len(results)
    times = [(r['time'] - results[0]['time']).total_seconds() / 3600.0 for r in results]
    dE = [float(r['dE']) * 100 for r in results]
    dN = [float(r['dN']) * 100 for r in results]
    dU = [float(r['dU']) * 100 for r in results]
    re = float(np.sqrt(np.nanmean(np.array(dE)**2)))
    rn = float(np.sqrt(np.nanmean(np.array(dN)**2)))
    ru = float(np.sqrt(np.nanmean(np.array(dU)**2)))
    r3 = float(np.sqrt(np.nanmean(np.array(dE)**2 + np.array(dN)**2 + np.array(dU)**2)))
    tmin, tmax = min(times), max(times)
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
    svg.append(f'<text x="{W//2}" y="{H-8}" text-anchor="middle" font-size="12" fill="#333">{label} | 3D RMS={r3:.2f} cm ({n} epochs)</text>')
    svg.append('</svg>')
    return '\n'.join(svg)

CSS = """*{box-sizing:border-box;margin:0;padding:0}
body{font-family:DejaVu Sans,Arial,sans-serif;background:#f4f6fa;padding:16px}
.hdr{text-align:center;padding:24px;background:linear-gradient(135deg,#1a237e,#1565C0);color:white;border-radius:12px;margin-bottom:20px}
.hdr h1{font-size:22px;margin-bottom:6px}
.hdr p{font-size:13px;opacity:0.85}
.card{background:white;border-radius:12px;padding:20px;margin:0 auto 16px;max-width:960px;box-shadow:0 2px 8px rgba(0,0,0,0.06)}
.day-hdr{font-size:16px;font-weight:bold;color:#1a237e;margin-bottom:12px;display:flex;justify-content:space-between;align-items:center}
.stats{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin:12px 0}
.s{background:#f8f9fa;border-radius:8px;padding:12px;text-align:center}
.s .v{font-size:20px;font-weight:bold}
.s .l{font-size:11px;color:#888;margin-top:2px}
.e .v{color:#1565C0}.n .v{color:#2E7D32}.u .v{color:#C62828}.a .v{color:#6A1B9A}
.note{font-size:12px;color:#888;text-align:center;margin-top:8px}"""

def main():
    days = [(2024, 4, 29), (2024, 4, 30), (2024, 5, 1), (2024, 5, 2)]
    all_results = []
    for y, m, d in days:
        r = run_day(y, m, d, nhours=4.0, interval=30.0)
        if r: all_results.append(r)
    print("\n" + "="*50)
    print("SUMMARY — 四天验证结果")
    print("="*50)
    for r in all_results:
        print(f"  {r['date']}: E={r['re']:.2f}cm N={r['rn']:.2f}cm U={r['ru']:.2f}cm 3D={r['r3']:.2f}cm n={r['n']}")
    # HTML
    html = ['<!DOCTYPE html><html lang="zh"><head>',
            '<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">',
            '<title>GRACE-FO PPP 多日验证报告</title>',
            '<style>' + CSS + '</style></head>',
            '<body><div class="hdr">',
            '<h1>GRACE-FO PPP 精密定轨 — 多日验证报告</h1>',
            '<p>PPP vs GNV1B 参考轨道 | 光行时 + Sagnac + 相对论 + 天线 PCV 改正 | 采样 30s | 每段 4 小时</p>',
            '</div>']
    for r in all_results:
        svg = make_svg(r['results'], r['date'])
        html.append(f'<div class="card">'
                    f'<div class="day-hdr"><span>{r["date"]}</span><span style="font-size:12px;color:#888">{r["n"]} 历元</span></div>'
                    f'<div class="stats">'
                    f'<div class="s e"><div class="v">{r["re"]:.2f} cm</div><div class="l">E（东西）RMS</div></div>'
                    f'<div class="s n"><div class="v">{r["rn"]:.2f} cm</div><div class="l">N（南北）RMS</div></div>'
                    f'<div class="s u"><div class="v">{r["ru"]:.2f} cm</div><div class="l">U（垂直）RMS</div></div>'
                    f'<div class="s a"><div class="v">{r["r3"]:.2f} cm</div><div class="l">3D RMS</div></div>'
                    f'</div>' + svg +
                    f'<div class="note">X轴=时间（h） Y轴=误差（cm） 蓝=E 绿=N 红=U | 参考: ISDC/GFZ GNV1B</div>'
                    f'</div>')
    html.append('</body></html>')
    out = 'output/report.html'
    Path(out).write_text(''.join(html))
    print(f"\n报告: {out} ({os.path.getsize(out)//1024}KB)")

if __name__ == '__main__':
    main()
