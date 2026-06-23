#!/usr/bin/env python3
"""生成过滤后的最终报告（<1m 阈值），只依赖已有 CSV"""
import csv, numpy as np, os
from pathlib import Path

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
.note{font-size:12px;color:#888;text-align:center;margin-top:8px}
.filtered{font-size:11px;color:#888;margin-top:4px}"""

def make_svg(dE, dN, dU, n, re, rn, ru, r3, label, W=900, H=520):
    times = [i * 30.0 / 3600.0 for i in range(n)]
    tmin, tmax = 0, max(times)
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
    svg.append(f'<text x="{W//2}" y="{H-8}" text-anchor="middle" font-size="12" fill="#333">{label} | 3D RMS={r3:.2f} cm ({n} epochs, filtered &lt;1m)</text>')
    svg.append('</svg>')
    return '\n'.join(svg)

def process_csv(path):
    rows = list(csv.DictReader(open(path)))
    # 过滤: 3D 误差 < 1m 且非 NaN
    valid = [r for r in rows if r['d3_m'] and not np.isnan(float(r['d3_m'])) and abs(float(r['d3_m'])) < 1.0]
    dE = [float(r['dE_m']) * 100 for r in valid]
    dN = [float(r['dN_m']) * 100 for r in valid]
    dU = [float(r['dU_m']) * 100 for r in valid]
    re = float(np.sqrt(np.nanmean(np.array(dE)**2)))
    rn = float(np.sqrt(np.nanmean(np.array(dN)**2)))
    ru = float(np.sqrt(np.nanmean(np.array(dU)**2)))
    r3 = float(np.sqrt(np.nanmean(np.array(dE)**2 + np.array(dN)**2 + np.array(dU)**2)))
    return dE, dN, dU, re, rn, ru, r3, len(valid), len(rows)

files = {
    '2024-04-29': 'output/ppp_multi_2024-04-29_4h.csv',
    '2024-04-30': 'output/ppp_multi_2024-04-30_4h.csv',
    '2024-05-01': 'output/ppp_multi_2024-05-01_4h.csv',
    '2024-05-02': 'output/ppp_multi_2024-05-02_4h.csv',
}

results = []
for date, path in files.items():
    if not os.path.exists(path):
        print(f"跳过 (无数据): {date}"); continue
    dE, dN, dU, re, rn, ru, r3, n_valid, n_total = process_csv(path)
    print(f"  {date}: E={re:.2f}cm N={rn:.2f}cm U={ru:.2f}cm 3D={r3:.2f}cm ({n_valid}/{n_total} epochs)")
    svg = make_svg(dE, dN, dU, n_valid, re, rn, ru, r3, date)
    results.append({'date': date, 're': re, 'rn': rn, 'ru': ru, 'r3': r3, 'n_valid': n_valid, 'n_total': n_total, 'svg': svg})

html = ['<!DOCTYPE html><html lang="zh"><head>',
        '<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">',
        '<title>GRACE-FO PPP 多日验证报告</title>',
        '<style>' + CSS + '</style></head>',
        '<body><div class="hdr">',
        '<h1>GRACE-FO PPP 精密定轨 — 多日验证报告</h1>',
        '<p>PPP vs GNV1B 参考轨道 | 光行时 + Sagnac + 相对论 + 天线 PCV 改正 | 采样 30s | 每段 4 小时 | 过滤 3D 误差 &lt;1m</p>',
        '</div>']
for r in results:
    html.append(f'<div class="card">'
                f'<div class="day-hdr"><span>{r["date"]}</span>'
                f'<span style="font-size:12px;color:#888">{r["n_valid"]}/{r["n_total"]} 历元（过滤&lt;1m）</span></div>'
                f'<div class="stats">'
                f'<div class="s e"><div class="v">{r["re"]:.2f} cm</div><div class="l">E（东西）RMS</div></div>'
                f'<div class="s n"><div class="v">{r["rn"]:.2f} cm</div><div class="l">N（南北）RMS</div></div>'
                f'<div class="s u"><div class="v">{r["ru"]:.2f} cm</div><div class="l">U（垂直）RMS</div></div>'
                f'<div class="s a"><div class="v">{r["r3"]:.2f} cm</div><div class="l">3D RMS</div></div>'
                f'</div>' + r['svg'] +
                f'<div class="note">X轴=时间（h） Y轴=误差（cm） 蓝=E 绿=N 红=U | 参考: ISDC/GFZ GNV1B</div>'
                f'</div>')
html.append('</body></html>')
Path('output/report.html').write_text(''.join(html))
print(f"\n报告: output/report.html ({os.path.getsize('output/report.html')//1024}KB)")
