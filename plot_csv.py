#!/usr/bin/env python3
"""读取 PPP CSV 结果，生成误差对比 PNG 图（纯 Python + ImageMagick）"""
import sys, csv, math, subprocess
from datetime import datetime

def read_csv(path):
    times, dE, dN, dU, n_sat = [], [], [], [], []
    with open(path) as f:
        reader = csv.DictReader(f)
        t0 = None
        for row in reader:
            t = datetime.fromisoformat(row['time'])
            if t0 is None: t0 = t
            times.append((t - t0).total_seconds() / 3600.0)  # hours
            dE.append(float(row['dE_m']) * 100)  # cm
            dN.append(float(row['dN_m']) * 100)
            dU.append(float(row['dU_m']) * 100)
            n_sat.append(int(row['n_sat']))
    d3 = [math.sqrt(e**2 + n**2 + u**2) for e, n, u in zip(dE, dN, dU)]
    rms3 = math.sqrt(sum(x**2 for x in d3) / len(d3))
    rms_e = math.sqrt(sum(x**2 for x in dE) / len(dE))
    rms_n = math.sqrt(sum(x**2 for x in dN) / len(dN))
    rms_u = math.sqrt(sum(x**2 for x in dU) / len(dU))
    return times, dE, dN, dU, d3, n_sat, rms_e, rms_n, rms_u, rms3

def mkaxis(W, H, xmin, xmax, ymin, ymax, px, py, pw, ph):
    ax = {'W':W,'H':H,'xmin':xmin,'xmax':xmax,'ymin':ymin,'ymax':ymax,
          'px':px,'py':py,'pw':pw,'ph':ph,'xr':(xmax-xmin)/pw,'yr':(ymax-ymin)/ph}
    ax['ax_x'] = px; ax['ax_y'] = py + ph; ax['ax_w'] = pw; ax['ax_h'] = ph
    return ax

def pxt(ax, x): return ax['px'] + (x - ax['xmin']) / ax['xr']
def pyt(ax, y): return ax['py'] + (1 - (y - ax['ymin']) / ax['yr']) * ax['ph']

def svg_header(W, H):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
            f'font-family="DejaVu Sans,Arial,sans-serif" font-size="11">\n'
            f'<rect width="{W}" height="{H}" fill="white"/>\n')

def grid(ax, xtick=1, ytick=0.5):
    lines = []
    for x in xtick:
        xi = pxt(ax, x)
        if ax['px'] <= xi <= ax['px'] + ax['pw']:
            lines.append(f'<line x1="{xi:.1f}" y1="{ax["py"]:.1f}" x2="{xi:.1f}" y2="{ax["py"]+ax["ph"]:.1f}" stroke="#ddd" stroke-width="0.5"/>')
    for y in ytick:
        yi = pyt(ax, y)
        if ax['py'] <= yi <= ax['py']+ax['ph']:
            lines.append(f'<line x1="{ax["px"]:.1f}" y1="{yi:.1f}" x2="{ax["px"]+ax["pw"]:.1f}" y2="{yi:.1f}" stroke="#ddd" stroke-width="0.5"/>')
    # zero line
    zi = pyt(ax, 0)
    if ax['py'] <= zi <= ax['py']+ax['ph']:
        lines.append(f'<line x1="{ax["px"]:.1f}" y1="{zi:.1f}" x2="{ax["px"]+ax["pw"]:.1f}" y2="{zi:.1f}" stroke="#666" stroke-width="0.8"/>')
    return '\n'.join(lines)

def polyline(ax, xs, ys, color, width=1, opacity=1.0):
    pts = ' '.join(f'{pxt(ax,xs[i]):.1f},{pyt(ax,ys[i]):.1f}' for i in range(len(xs)) if ax['px']<=pxt(ax,xs[i])<=ax['px']+ax['pw'] and ax['py']<=pyt(ax,ys[i])<=ax['py']+ax['ph'])
    return f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="{width}" opacity="{opacity}"/>'

def axis_labels(ax, xlabel='', ylabel=''):
    lines = []
    xi = ax['px'] + ax['pw']/2
    lines.append(f'<text x="{xi:.1f}" y="{ax["H"]-5}" text-anchor="middle" font-size="12">{xlabel}</text>')
    yi = ax['py'] + ax['ph']/2
    lines.append(f'<text x="12" y="{yi:.1f}" text-anchor="middle" font-size="12" transform="rotate(-90,12,{yi:.1f})">{ylabel}</text>')
    return '\n'.join(lines)

def text(ax, x, y, s, color='black', anchor='start', size=11):
    xi = pxt(ax, x) if ax['px'] <= x <= ax['px']+ax['pw'] else ax['px']
    return f'<text x="{xi:.1f}" y="{pyt(ax,y):.1f}" fill="{color}" font-size="{size}" text-anchor="{anchor}" dy="-2">{s}</text>'

def plot_timeseries(times, dE, dN, dU, d3, n_sat, rms_e, rms_n, rms_u, rms3, out_svg):
    W, H = 900, 700
    margin_top = 50; margin_bot = 50; margin_left = 70; margin_right = 30
    title_h = 30
    gap = 20
    each_h = (H - margin_top - title_h - margin_bot - 2*gap) // 3
    pw = W - margin_left - margin_right
    tmin, tmax = min(times), max(times)
    rms_max = max(max(abs(v) for v in dE), max(abs(v) for v in dN), max(abs(v) for v in dU)) * 1.1
    rms_max = max(rms_max, 1.0)

    axE = mkaxis(W, H, tmin, tmax, -rms_max, rms_max, margin_left, margin_bot + 2*each_h + gap + title_h, pw, each_h)
    axN = mkaxis(W, H, tmin, tmax, -rms_max, rms_max, margin_left, margin_bot + each_h + gap + title_h, pw, each_h)
    axU = mkaxis(W, H, tmin, tmax, -rms_max, rms_max, margin_left, margin_bot + title_h, pw, each_h)

    svg = [svg_header(W, H)]
    svg.append(f'<text x="{W//2}" y="25" text-anchor="middle" font-size="15" font-weight="bold">GRACE-FO PPP 误差 vs GNV1B 精密轨道</text>')

    for ax, data, color, label, rms_val in [
        (axE, dE, '#2196F3', 'E (东西)', rms_e),
        (axN, dN, '#4CAF50', 'N (南北)', rms_n),
        (axU, dU, '#F44336', 'U (垂直)', rms_u),
    ]:
        svg.append(f'<g transform="translate(0,0)">')
        svg.append(grid(ax, xtick=[(tmin//1+i*0.5) for i in range(int((tmax-tmin)//0.5)+2)], ytick=[]))
        svg.append(polyline(ax, times, data, color, width=0.8))
        svg.append(f'<text x="{ax["px"]+5}" y="{ax["py"]+15}" font-size="11" fill="{color}">RMS={rms_val:.2f} cm</text>')
        svg.append(f'<text x="{ax["px"]+ax["pw"]-5}" y="{ax["py"]+15}" text-anchor="end" font-size="11" fill="#333">{label}</text>')
        svg.append(axis_labels(ax, '时间 (小时)', ''))
        svg.append('</g>')

    svg.append(f'<text x="{W//2}" y="{H-5}" text-anchor="middle" font-size="11" fill="#666">3D RMS = {rms3:.2f} cm</text>')
    svg.append('</svg>')

    svg_str = '\n'.join(svg)
    with open(out_svg, 'w') as f: f.write(svg_str)
    return svg_str

if __name__ == '__main__':
    import os
    csv_path = sys.argv[1] if len(sys.argv) > 1 else 'output/ppp_vs_gnv1b_2024_0429_0h.csv'
    out_svg = csv_path.replace('.csv', '_enu.svg')
    times, dE, dN, dU, d3, n_sat, rms_e, rms_n, rms_u, rms3 = read_csv(csv_path)
    print(f"CSV: {csv_path}")
    print(f"历元: {len(times)}")
    print(f"RMS E={rms_e:.2f}cm  N={rms_n:.2f}cm  U={rms_u:.2f}cm  3D={rms3:.2f}cm")
    plot_timeseries(times, dE, dN, dU, d3, n_sat, rms_e, rms_n, rms_u, rms3, out_svg)
    print(f"SVG: {out_svg}")

    # Convert to PNG with ImageMagick
    out_png = out_svg.replace('.svg', '.png')
    try:
        subprocess.run(['convert', '-density', '150', out_svg, '-quality', '90', out_png],
                       capture_output=True, check=True)
        size = os.path.getsize(out_png)
        print(f"PNG: {out_png} ({size//1024} KB)")
    except Exception as e:
        print(f"ImageMagick convert failed: {e}")
        # Fallback: just keep SVG
        print(f"Keeping SVG: {out_svg}")
