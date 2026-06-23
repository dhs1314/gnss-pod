#!/usr/bin/env python3
"""生成 GRACE-FO PPP 误差对比图（GRACE-FO PPP vs GNV1B 精密轨道）"""
import csv, math, subprocess, os, sys

def read_csv(path):
    times, dE, dN, dU = [], [], [], []
    with open(path) as f:
        reader = csv.DictReader(f)
        t0 = None
        for row in reader:
            if row['dN_m'] == 'nan' or row['dU_m'] == 'nan':
                continue
            t = row['time']
            h = int(t[11:13]); m = int(t[14:16]); s = int(t[17:19])
            seconds = h*3600 + m*60 + s
            if t0 is None: t0 = seconds
            times.append((seconds - t0) / 3600.0)
            dE.append(float(row['dE_m']) * 100)
            dN.append(float(row['dN_m']) * 100)
            dU.append(float(row['dU_m']) * 100)
    n = len(times)
    d3 = [math.sqrt(dE[i]**2 + dN[i]**2 + dU[i]**2) for i in range(n)]
    rms3 = math.sqrt(sum(x**2 for x in d3)/n) if n else 0
    rms_e = math.sqrt(sum(x**2 for x in dE)/n) if n else 0
    rms_n = math.sqrt(sum(x**2 for x in dN)/n) if n else 0
    rms_u = math.sqrt(sum(x**2 for x in dU)/n) if n else 0
    return times, dE, dN, dU, d3, n, rms_e, rms_n, rms_u, rms3

def make_png(times, dE, dN, dU, d3, n, rms_e, rms_n, rms_u, rms3, out_png):
    """用纯 Python + zlib 生成 PNG（无任何外部库依赖）"""
    import zlib, struct

    W, H = 900, 600
    margin_l = 70; margin_r = 30; margin_t = 50; margin_b = 50
    title_h = 40
    gap = 15
    plot_h = (H - margin_t - margin_b - title_h - 2*gap) // 3
    plot_w = W - margin_l - margin_r

    tmin, tmax = min(times), max(times)
    ymax = max(max(abs(v) for v in dE), max(abs(v) for v in dN), max(abs(v) for v in dU)) * 1.2
    ymax = max(ymax, 0.1)

    def px(t_val): return margin_l + (t_val - tmin) / (tmax - tmin) * plot_w
    def py(y_val): return H - margin_b - (y_val / ymax + 1) / 2 * plot_h - title_h - gap

    # RGB pixels (white background)
    pixels = [[255, 255, 255] * W for _ in range(H)]

    def set_pixel(x, y, r, g, b):
        if 0 <= x < W and 0 <= y < H:
            pixels[y][x] = [r, g, b]

    def draw_line(x0, y0, x1, y1, r, g, b, w=1):
        dx, dy = abs(x1-x0), abs(y1-y0)
        sx = 1 if x0 < x1 else -1; sy = 1 if y0 < y1 else -1
        err = dx - dy
        x, y = x0, y0
        for _ in range(10000):
            for dr in range(-w//2, w//2+1):
                for dc in range(-w//2, w//2+1):
                    set_pixel(x+dc, y+dr, r, g, b)
            if abs(x-x0) + abs(y-y0) > math.hypot(x1-x0, y1-y0): break
            e2 = 2*err
            if e2 > -dy: err -= dy; x += sx
            if e2 < dx: err += dx; y += sy

    def draw_text_simple(x, y, text, size, r, g, b):
        # Just draw as pixel dots for simple numbers
        pass  # Skip text - will add as label overlay

    def draw_plot_row(ax_y, ax_h, data, color_rgb, label, rms_val):
        # Grid lines
        for t_val in times:
            xi = int(px(t_val))
            if xi < margin_l or xi > W - margin_r: continue
            for yi in range(ax_y, ax_y + ax_h):
                if 0 <= yi < H:
                    pixels[yi][xi] = [220, 220, 220]
        # Zero line
        yz = int(py(0))
        if ax_y <= yz <= ax_y + ax_h:
            for xi in range(margin_l, W - margin_r):
                pixels[yz][xi] = [100, 100, 100]
        # Data polyline
        for i in range(len(times)-1):
            x0, y0 = int(px(times[i])), int(py(data[i]))
            x1, y1 = int(px(times[i+1])), int(py(data[i+1]))
            if abs(x1-x0) < 2:
                if ax_y <= y0 <= ax_y+ax_h:
                    pixels[y0][x0] = color_rgb
            else:
                draw_line(x0, y0, x1, y1, *color_rgb, w=1)

    # Draw 3 rows
    row_colors = [[33, 113, 221], [67, 160, 72], [211, 48, 44]]
    row_labels = ['E (东西方向)', 'N (南北方向)', 'U (垂直方向)']
    row_data = [dE, dN, dU]
    row_rms = [rms_e, rms_n, rms_u]
    row_ys = [H - margin_b - (i+1)*plot_h - gap*(2-i) - title_h for i in range(3)]

    for ry, data, col, lbl, rms_v in zip(row_ys, row_data, row_colors, row_labels, row_rms):
        draw_plot_row(ry, plot_h, data, col, lbl, rms_v)

    # ── Write PNG ──────────────────────────────────────────────────────────
    def write_png(pixels, w, h, path):
        def chunk(tag, data):
            c = tag + data
            return struct.pack('>I', len(data)) + c + struct.pack('>I', zlib.crc32(c) & 0xffffffff)

        sig = b'\x89PNG\r\n\x1a\n'
        ihdr = struct.pack('>IIBBBBB', w, h, 8, 2, 0, 0, 0)
        raw_rows = []
        for row in pixels:
            raw_rows.append(b'\x00' + bytes(sum([[r, g, b] for r, g, b in [row[i*3:(i+1)*3] for i in range(w)]], [])))
        idat_data = zlib.compress(b''.join(raw_rows), 6)
        with open(path, 'wb') as f:
            f.write(sig + chunk(b'IHDR', ihdr) + chunk(b'IDAT', idat_data) + chunk(b'IEND', b''))

    # Generate title and axis info as a separate small PNG and combine
    # Instead: generate full image with title baked in
    title_row = [[255]*W for _ in range(title_h)]
    # Simple title: fill with light gray
    for y in range(title_h):
        for x in range(W):
            pixels[y][x] = [250, 250, 250]

    # Add bottom label row
    for x in range(W):
        for y in range(H - margin_b, H):
            pixels[y][x] = [250, 250, 250]

    write_png(pixels, W, H, out_png)
    return W, H

def svg_with_text(times, dE, dN, dU, d3, n, rms_e, rms_n, rms_u, rms3, out_svg, out_png):
    """Generate SVG with embedded font, then convert to PNG"""
    tmin, tmax = min(times), max(times)
    ymax = max(max(abs(v) for v in dE), max(abs(v) for v in dN), max(abs(v) for v in dU)) * 1.2
    ymax = max(ymax, 0.1)
    W, H = 900, 600
    margin_l = 70; margin_r = 30; margin_t = 50; margin_b = 50
    title_h = 40; gap = 15
    plot_h = (H - margin_t - margin_b - title_h - 2*gap) // 3
    plot_w = W - margin_l - margin_r

    def px(t_val): return margin_l + (t_val - tmin) / (tmax - tmin) * plot_w
    def py(y_val): return H - margin_b - int(((y_val / ymax + 1) / 2) * plot_h) - title_h - gap

    def ppts(data):
        return ' '.join(f'{px(times[i]):.1f},{py(data[i]):.1f}' for i in range(len(times))
                        if margin_l <= px(times[i]) <= W-margin_r and
                           H-margin_b-title_h-gap-plot_h <= py(data[i]) <= H-margin_b-title_h-gap)

    row_data = [('E (m)', dE, '#1565C0', rms_e),
                ('N (m)', dN, '#2E7D32', rms_n),
                ('U (m)', dU, '#C62828', rms_u)]
    rows_y = [H - margin_b - (i+1)*plot_h - gap*(2-i) - title_h for i in range(3)]

    svg_lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" font-family="DejaVu-Sans,DejaVu Sans,sans-serif" font-size="11">',
        f'<rect width="{W}" height="{H}" fill="white"/>',
        f'<rect x="{margin_l}" y="{title_h}" width="{plot_w}" height="{H-title_h-margin_b}" fill="#fafafa"/>',
        f'<text x="{W//2}" y="26" text-anchor="middle" font-size="15" font-weight="bold" fill="#222">GRACE-FO PPP 误差 vs GNV1B 精密轨道（2小时，采样30秒）</text>',
    ]

    for (label, data, color, rms_v), ry in zip(row_data, rows_y):
        svg_lines.append(f'<text x="4" y="{ry+14}" font-size="11" fill="{color}" font-weight="bold">{label}</text>')
        svg_lines.append(f'<text x="{W-margin_r-4}" y="{ry+14}" text-anchor="end" font-size="11" fill="#555">RMS={rms_v:.2f} cm</text>')
        svg_lines.append(f'<polyline points="{ppts(data)}" fill="none" stroke="{color}" stroke-width="0.8"/>')
        svg_lines.append(f'<line x1="{margin_l}" y1="{py(0)}" x2="{W-margin_r}" y2="{py(0)}" stroke="#888" stroke-width="0.8"/>')

    svg_lines.append(f'<text x="{W//2}" y="{H-12}" text-anchor="middle" font-size="12" fill="#333">时间 (小时)  —  3D RMS = {rms3:.2f} cm  ({n} 历元)</text>')
    svg_lines.append('</svg>')

    svg_content = '\n'.join(svg_lines)
    with open(out_svg, 'w') as f: f.write(svg_content)

    # Convert SVG → PNG with ImageMagick
    try:
        res = subprocess.run(
            ['convert', '-density', '150',
             f'label:{W}x{H}',
             '-resize', f'{W}x{H}',
             '-background', 'white',
             '-alpha', 'remove',
             '-quality', '90',
             out_svg, out_png],
            capture_output=True, text=True, timeout=30
        )
        print(f"convert stdout: {res.stdout[:200] if res.stdout else ''}")
        print(f"convert stderr: {res.stderr[:200] if res.stderr else ''}")
    except Exception as e:
        print(f"convert error: {e}")

    # Fallback: generate PNG directly
    if not os.path.exists(out_png) or os.path.getsize(out_png) < 100:
        print("Using pure Python PNG fallback...")
        make_png(times, dE, dN, dU, d3, n, rms_e, rms_n, rms_u, rms3, out_png)

    return svg_content

if __name__ == '__main__':
    csv_path = sys.argv[1] if len(sys.argv) > 1 else 'output/ppp_vs_gnv1b_2024_0429_0h.csv'
    out_dir = os.path.dirname(csv_path)
    date_str = os.path.basename(csv_path).replace('ppp_vs_gnv1b_', '').replace('.csv', '')
    out_svg = os.path.join(out_dir, f'ppp_error_{date_str}.svg')
    out_png = os.path.join(out_dir, f'ppp_error_{date_str}.png')

    times, dE, dN, dU, d3, n, rms_e, rms_n, rms_u, rms3 = read_csv(csv_path)
    print(f"有效历元: {n}")
    print(f"E RMS: {rms_e:.2f} cm")
    print(f"N RMS: {rms_n:.2f} cm")
    print(f"U RMS: {rms_u:.2f} cm")
    print(f"3D RMS: {rms3:.2f} cm")

    svg_content = svg_with_text(times, dE, dN, dU, d3, n, rms_e, rms_n, rms_u, rms3, out_svg, out_png)
    print(f"SVG: {out_svg}")
    if os.path.exists(out_png):
        print(f"PNG: {out_png} ({os.path.getsize(out_png)//1024} KB)")
