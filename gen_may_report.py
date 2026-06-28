#!/usr/bin/env python3
"""生成5月全月报告"""
import csv, numpy as np, os
from pathlib import Path

od = Path('/workspace/gnss_pod/output')
files = sorted(od.glob('ppp_may_2024-05-*_4h.csv'))
print(f'找到 {len(files)} 个 CSV')

all_results = []
for f in files:
    rows = list(csv.DictReader(open(f)))
    valid = [r for r in rows if r['d3_m'] and not np.isnan(float(r['d3_m'])) and abs(float(r['d3_m'])) < 1.0]
    if not valid: continue
    dE = np.array([float(r['dE_m']) * 100 for r in valid])
    dN = np.array([float(r['dN_m']) * 100 for r in valid])
    dU = np.array([float(r['dU_m']) * 100 for r in valid])
    re = float(np.sqrt(np.nanmean(dE**2)))
    rn = float(np.sqrt(np.nanmean(dN**2)))
    ru = float(np.sqrt(np.nanmean(dU**2)))
    r3 = float(np.sqrt(np.nanmean(dE**2 + dN**2 + dU**2)))
    parts = f.stem.split('_')
    # 文件名格式: ppp_may_2024-05-01_4h.csv
    date = parts[2]  # '2024-05-01'
    all_results.append({'date': date, 're': re, 'rn': rn, 'ru': ru, 'r3': r3,
                        'n_ok': len(valid), 'n_total': len(rows),
                        'dE': list(dE), 'dN': list(dN), 'dU': list(dU)})

print(f'成功: {len(all_results)} 天')

re_all = [r['re'] for r in all_results]
rn_all = [r['rn'] for r in all_results]
ru_all = [r['ru'] for r in all_results]
r3_all = [r['r3'] for r in all_results]
n_ok_all = sum(r['n_ok'] for r in all_results)
n_tot_all = sum(r['n_total'] for r in all_results)
print(f'\n月统计:')
print(f'  E:  {np.mean(re_all):.2f} ± {np.std(re_all):.2f} cm  (range {min(re_all):.2f}~{max(re_all):.2f})')
print(f'  N:  {np.mean(rn_all):.2f} ± {np.std(rn_all):.2f} cm  (range {min(rn_all):.2f}~{max(rn_all):.2f})')
print(f'  U:  {np.mean(ru_all):.2f} ± {np.std(ru_all):.2f} cm  (range {min(ru_all):.2f}~{max(ru_all):.2f})')
print(f'  3D: {np.mean(r3_all):.2f} ± {np.std(r3_all):.2f} cm  (range {min(r3_all):.2f}~{max(r3_all):.2f})')
print(f'  成功率: {n_ok_all}/{n_tot_all} ({n_ok_all/n_tot_all*100:.1f}%)')

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
    rows_data = [(dE, '#1565C0', 'E (东西)'), (dN, '#2E7D32', 'N (南北)'), (dU, '#C62828', 'U (垂直)')]
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" font-family="DejaVu Sans,sans-serif" font-size="11"><rect width="{W}" height="{H}" fill="white"/>']
    for i, (data, col, lbl) in enumerate(rows_data):
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
.day-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(480px,1fr));gap:12px}
.card{background:white;border-radius:12px;padding:20px;box-shadow:0 2px 8px rgba(0,0,0,0.06)}
.card-hdr{font-size:15px;font-weight:bold;color:#1a237e;margin-bottom:10px;display:flex;justify-content:space-between;align-items:center}
.stats{display:grid;grid-template-columns:repeat(4,1fr);gap:8px;margin:10px 0}
.s{background:#f8f9fa;border-radius:6px;padding:10px;text-align:center}
.s .v{font-size:18px;font-weight:bold}
.s .l{font-size:10px;color:#888;margin-top:2px}
.e .v{color:#1565C0}.n .v{color:#2E7D32}.u .v{color:#C62828}.a .v{color:#6A1B9A}
.note{font-size:11px;color:#888;text-align:center;margin-top:6px}
.sum{background:white;border-radius:12px;padding:20px;margin-top:20px;box-shadow:0 2px 8px rgba(0,0,0,0.06)}
.sum h2{font-size:16px;color:#1a237e;margin-bottom:16px}
.sum-grid{display:grid;grid-template-columns:repeat(6,1fr);gap:8px}
.su{background:#f8f9fa;border-radius:6px;padding:8px;text-align:center;font-size:12px}
.su .v{font-size:16px;font-weight:bold}
.su .l{font-size:10px;color:#888}"""

html = ['<!DOCTYPE html><html lang="zh"><head>',
        '<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">',
        '<title>GRACE-FO PPP 2024年5月全月验证报告</title>',
        '<style>' + CSS + '</style></head>',
        '<body><div class="hdr">',
        '<h1>GRACE-FO PPP 精密定轨 — 2024年5月全月验证报告</h1>',
        '<p>PPP vs GNV1B 参考轨道 | 光行时 + Sagnac + 相对论 + 天线 PCV 改正 | 采样 30s | 每段 4 小时 | 过滤 3D 误差 &lt;1m</p>',
        '</div><div class="day-grid">']
for r in all_results:
    svg = make_svg(r['dE'], r['dN'], r['dU'], r['n_ok'], r['re'], r['rn'], r['ru'], r['r3'], r['date'])
    html.append(f'<div class="card">'
                f'<div class="card-hdr"><span>{r["date"]}</span><span style="font-size:11px;color:#888">{r["n_ok"]}/{r["n_total"]} epochs</span></div>'
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
           f'</div></div>')
html.append('</body></html>')
out = str(od / 'report_may2024.html')
Path(out).write_text(''.join(html))
print(f'\n报告: {out} ({os.path.getsize(out)//1024}KB)')