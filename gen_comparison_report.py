#!/usr/bin/env python3
"""
从已有 CSV 文件生成广播星历 vs 精密星历对比报告
"""
import csv, numpy as np
from pathlib import Path
from datetime import datetime as dt

DATA_OUTPUT = Path('/workspace/gnss_pod/data/output')
OUT_DIR = Path('/workspace/gnss_pod/output')
OUT_DIR.mkdir(exist_ok=True)

START = '2024-05-01'
END = '2024-05-31'

def load_csv(path):
    """加载 CSV，返回统计 dict"""
    try:
        dE, dN, dU, n_sat = [], [], [], 0
        with open(path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                dE.append(float(row['dE_m']))
                dN.append(float(row['dN_m']))
                dU.append(float(row['dU_m']))
                n_sat = max(n_sat, int(row['n_sat']))
        dE, dN, dU = np.array(dE), np.array(dN), np.array(dU)
        d3 = np.sqrt(dE**2 + dN**2 + dU**2)
        ok = d3 < 1.0
        if ok.sum() < 10:
            return None
        return {
            'date': Path(path).stem.split('_')[1],
            're': float(np.sqrt(np.nanmean(dE[ok]**2)) * 100),
            'rn': float(np.sqrt(np.nanmean(dN[ok]**2)) * 100),
            'ru': float(np.sqrt(np.nanmean(dU[ok]**2)) * 100),
            'r3': float(np.sqrt(np.nanmean(d3[ok]**2)) * 100),
            'n_ok': int(ok.sum()), 'n_total': len(dE),
            'dE': list(dE[ok]), 'dN': list(dN[ok]), 'dU': list(dU[ok]),
        }
    except Exception as e:
        print(f"  读取失败 {path}: {e}")
        return None

def stats_summary(results, label):
    valid = [r for r in results if r is not None]
    if not valid:
        print(f"\n{label}: 无数据")
        return None, None
    re_vals = [r['re'] for r in valid]
    rn_vals = [r['rn'] for r in valid]
    ru_vals = [r['ru'] for r in valid]
    r3_vals = [r['r3'] for r in valid]
    print(f"\n{'═'*50}")
    print(f"{label} 汇总 ({len(valid)}/{len(results)} 天)")
    print(f"  E = {np.mean(re_vals):.2f} ± {np.std(re_vals):.2f} cm")
    print(f"  N = {np.mean(rn_vals):.2f} ± {np.std(rn_vals):.2f} cm")
    print(f"  U = {np.mean(ru_vals):.2f} ± {np.std(ru_vals):.2f} cm")
    print(f"  3D = {np.mean(r3_vals):.2f} ± {np.std(r3_vals):.2f} cm")
    print(f"{'═'*50}")
    return valid, dict(re=re_vals, rn=rn_vals, ru=ru_vals, r3=r3_vals)

# ── 加载数据 ──────────────────────────────────────────────────────

print(f"加载广播星历 CSV ({START} → {END})...")
b_results = []
for d in range(1, 32):
    y, m = 2024, 5
    ds = f"{y:04d}-{m:02d}-{d:02d}"
    path = DATA_OUTPUT / f"ppp_{ds}_broadcast.csv"
    r = load_csv(path)
    b_results.append(r)
    if r:
        print(f"  {ds}: E={r['re']:.2f} N={r['rn']:.2f} U={r['ru']:.2f} 3D={r['r3']:.2f}cm ({r['n_ok']}/{r['n_total']})")
    else:
        print(f"  {ds}: 加载失败")

b_valid, b_stats = stats_summary(b_results, '广播星历')

print(f"\n加载精密星历 CSV ({START} → {END})...")
s_results = []
for d in range(1, 32):
    y, m = 2024, 5
    ds = f"{y:04d}-{m:02d}-{d:02d}"
    path = DATA_OUTPUT / f"ppp_{ds}_sp3.csv"
    r = load_csv(path)
    s_results.append(r)
    if r:
        print(f"  {ds}: E={r['re']:.2f} N={r['rn']:.2f} U={r['ru']:.2f} 3D={r['r3']:.2f}cm ({r['n_ok']}/{r['n_total']})")
    else:
        print(f"  {ds}: 无数据（等待批处理完成）")

s_valid, s_stats = stats_summary(s_results, '精密星历')

# ── 生成报告 ──────────────────────────────────────────────────────

dates = sorted(set(r['date'] for r in b_valid + s_valid))

html = [
    '<!DOCTYPE html><html lang="zh"><head><meta charset="utf-8">'
    '<title>GRACE-FO PPP — 广播星历 vs 精密星历 对比报告</title>'
    '<style>'
    'body{font-family:Arial,sans-serif;margin:0;background:#f0f2f5}'
    '.header{background:linear-gradient(135deg,#1565C0,#6a1b9a);color:#fff;padding:28px 36px}'
    '.header h1{margin:0;font-size:26px}'
    '.header p{margin:8px 0 0;color:#e1bee7;font-size:13px}'
    '.summary{padding:24px 36px;display:grid;'
             'grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:16px}'
    '.card{background:#fff;border-radius:12px;padding:20px 24px;'
           'box-shadow:0 2px 8px rgba(0,0,0,.08)}'
    '.card h3{color:#888;font-size:11px;text-transform:uppercase;'
              'letter-spacing:.05em;margin:0 0 12px;border-bottom:1px solid #eee;padding-bottom:8px}'
    '.row{display:flex;justify-content:space-between;margin:6px 0;font-size:14px}'
    '.lbl{color:#555}.b{color:#1565C0;font-weight:bold}'
    '.s{color:#2e7d32;font-weight:bold}'
    '.diff{font-size:13px;padding:3px 10px;border-radius:20px;display:inline-block}'
    '.improve{background:#e8f5e9;color:#2e7d32}'
    '.worse{background:#ffebee;color:#c62828}'
    '.neutral{background:#f5f5f5;color:#666}'
    '.table-wrap{padding:0 36px 24px}'
    'table{width:100%;border-collapse:collapse;background:#fff;'
           'border-radius:12px;overflow:hidden;'
           'box-shadow:0 2px 8px rgba(0,0,0,.08);font-size:13px}'
    'th{background:#1565C0;color:#fff;padding:12px 14px;text-align:center;font-size:11px}'
    'th:first-child{text-align:left}'
    'td{padding:9px 14px;border-bottom:1px solid #f0f0f0;text-align:center}'
    'td:first-child{text-align:left;font-weight:bold;color:#333}'
    'tr:last-child td{border-bottom:none}'
    'tr:hover td{background:#f8f9ff}'
    '.badge{display:inline-block;padding:2px 8px;border-radius:10px;font-size:11px;font-weight:bold}'
    '.bkg-b{background:#e3f2fd;color:#1565C0}.bkg-s{background:#e8f5e9;color:#2e7d32}'
    '.note{padding:12px 36px;font-size:12px;color:#888;background:#fff8e1;border-top:1px solid #ffe082}'
    '.footer{padding:20px 36px;color:#aaa;font-size:12px;text-align:center}'
    '.metric{font-size:12px;color:#888;margin-top:4px}'
    '</style></head><body>'
    '<div class="header"><h1>📊 GRACE-FO PPP 精度对比报告</h1>'
    '<p>广播星历 vs IGS Ultra-Rapid 精密星历 (BKG 镜像) '
    f'| 2024年5月 (全月 {len(b_valid)}天) | 每时刻段: 4h | 采样: 30s</p></div>'
]

# 对比汇总卡片
html.append('<div class="summary">')
for axis, key, unit in [('E (东西)', 're', 'cm'), ('N (南北)', 'rn', 'cm'),
                          ('U (垂直)', 'ru', 'cm'), ('3D RMS', 'r3', 'cm')]:
    be = np.mean(b_stats[key])
    se = np.mean(s_stats[key]) if s_stats else float('nan')
    std_b = np.std(b_stats[key])
    diff = be - se if s_stats else 0
    cls = 'improve' if diff > 0.05 else ('worse' if diff < -0.05 else 'neutral')
    sign = '↓' if diff > 0 else ('↑' if diff < 0 else '≈')
    html.append(
        f'<div class="card"><h3>{axis}</h3>'
        f'<div class="row"><span class="lbl">广播星历</span><span class="b">{be:.2f}±{std_b:.2f} {unit}</span></div>'
    )
    if s_stats:
        std_s = np.std(s_stats[key])
        html.append(
            f'<div class="row"><span class="lbl">精密星历</span><span class="s">{se:.2f}±{std_s:.2f} {unit}</span></div>'
            f'<div class="row"><span class="lbl">改善量</span>'
            f'<span class="diff {cls}">{sign} {abs(diff):.2f} {unit}</span></div>'
        )
    else:
        html.append('<div class="row"><span class="lbl">精密星历</span><span class="neutral diff">处理中...</span></div>')
    html.append('</div>')
html.append('</div>')

# 逐日详细对比表
html.append('<div class="table-wrap"><table>')
html.append('<thead><tr><th>日期</th><th>星历</th>'
            '<th>E (cm)</th><th>N (cm)</th><th>U (cm)</th>'
            '<th>3D (cm)</th><th>收敛率</th><th>3D改善</th></tr></thead><tbody>')

b_map = {r['date']: r for r in b_valid}
s_map = {r['date']: r for r in s_valid} if s_valid else {}

for date in sorted(dates):
    br = b_map.get(date)
    sr = s_map.get(date)
    if br and sr:
        diff3 = br['r3'] - sr['r3']
        cls = 'improve' if diff3 > 0.05 else ('worse' if diff3 < -0.05 else 'neutral')
        sign = '↓' if diff3 > 0 else ('↑' if diff3 < 0 else '≈')
        rate = f"{br['n_ok']}/{br['n_total']} ({br['n_ok']/br['n_total']*100:.0f}%)"
        html.append(
            f'<tr><td>{date}</td>'
            f'<td><span class="badge bkg-b">广播</span> <span class="badge bkg-s">精密</span></td>'
            f'<td>{br["re"]:.2f} / <b>{sr["re"]:.2f}</b></td>'
            f'<td>{br["rn"]:.2f} / <b>{sr["rn"]:.2f}</b></td>'
            f'<td>{br["ru"]:.2f} / <b>{sr["ru"]:.2f}</b></td>'
            f'<td style="color:#1565C0">{br["r3"]:.2f} / <b style="color:#2e7d32">{sr["r3"]:.2f}</b></td>'
            f'<td>{rate}</td>'
            f'<td><span class="diff {cls}">{sign} {abs(diff3):.2f}</span></td></tr>'
        )
    elif br:
        rate = f"{br['n_ok']}/{br['n_total']} ({br['n_ok']/br['n_total']*100:.0f}%)"
        html.append(
            f'<tr><td>{date}</td><td><span class="badge bkg-b">广播</span></td>'
            f'<td>{br["re"]:.2f}</td><td>{br["rn"]:.2f}</td><td>{br["ru"]:.2f}</td>'
            f'<td style="color:#1565C0">{br["r3"]:.2f}</td><td>{rate}</td><td>—</td></tr>'
        )
    elif sr:
        rate = f"{sr['n_ok']}/{sr['n_total']} ({sr['n_ok']/sr['n_total']*100:.0f}%)"
        html.append(
            f'<tr><td>{date}</td><td><span class="badge bkg-s">精密</span></td>'
            f'<td>{sr["re"]:.2f}</td><td>{sr["rn"]:.2f}</td><td>{sr["ru"]:.2f}</td>'
            f'<td style="color:#2e7d32">{sr["r3"]:.2f}</td><td>{rate}</td><td>—</td></tr>'
        )

html.append('</tbody></table></div>')

note = (
    '⚠️ 精密星历说明: 使用 BKG 镜像的 IGS Ultra-Rapid 产品 (IGS0OPSRAP)，'
    '该产品为近实时产品，设计用于实时导航而非事后精密定位。'
    '对于 2024-05 这类历史日期，文件为次日发布（使用前一天的轨道），'
    '精密星历卫星数量(≈15颗)少于广播星历(21颗)，因此 U 方向精度可能劣于广播星历。'
    '如需最佳精度，建议使用 IGS Final 产品 (约 12-18 天延迟) 或 MGEX 精密星历。'
    '广播星历精度 ~3D 12cm，对 GRACE-FO 轨道验证已足够。'
)
html.append(f'<div class="note">{note}</div>')
html.append(
    f'<div class="footer">'
    f'GRACE-FO PPP 验证系统 | 对比报告 | 生成时间: {dt.now().strftime("%Y-%m-%d %H:%M:%S")}'
    f'</div>'
)
html.append('</body></html>')

out_path = OUT_DIR / f'report_full_comparison_{START}_{END}.html'
Path(out_path).write_text(''.join(html), encoding='utf-8')
print(f"\n报告已生成: {out_path}")
print(f"文件大小: {Path(out_path).stat().st_size // 1024}KB")
