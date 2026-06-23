#!/usr/bin/env python3
"""
Pure Python PNG generator for PPP results (no third-party libs needed).
Uses only built-in struct + zlib to write valid PNG.
"""
import struct, zlib, sys, os
from pathlib import Path

def write_png(width, height, pixels, out_path):
    """Write a valid PNG file. pixels = flat RGB bytes (width*height*3)."""
    def chunk(name, data):
        name_bytes = name if isinstance(name, bytes) else name.encode('ascii')
        c = name_bytes + data
        return struct.pack('>I', len(data)) + c + struct.pack('>I', zlib.crc32(c) & 0xffffffff)
    sig = b'\x89PNG\r\n\x1a\n'
    ihdr = struct.pack('>IIBBBBB', width, height, 8, 2, 0, 0, 0)
    png = sig + chunk(b'IHDR', ihdr) + chunk(b'IDAT', zlib.compress(pixels, 6)) + chunk(b'IEND', b'')
    with open(out_path, 'wb') as f:
        f.write(png)


def hline(pixels, width, height, x, y, w, r, g, b):
    """Horizontal line at (x,y) for w pixels."""
    for dx in range(max(0,x), min(width, x+w)):
        if 0 <= y < height:
            i = (height-1-y)*width*3 + dx*3
            if 0 <= i < len(pixels)-2:
                pixels[i], pixels[i+1], pixels[i+2] = r, g, b


def vline(pixels, width, height, x, y, h, r, g, b):
    """Vertical line at (x,y) for h pixels."""
    for dy in range(max(0,y), min(height, y+h)):
        if 0 <= x < width:
            i = (height-1-dy)*width*3 + x*3
            if 0 <= i < len(pixels)-2:
                pixels[i], pixels[i+1], pixels[i+2] = r, g, b


def line(pixels, width, height, x0, y0, x1, y1, r, g, b, thick=1):
    """Bresenham line with thickness."""
    dx = abs(x1-x0); dy = -abs(y1-y0)
    sx = 1 if x0 < x1 else -1; sy = 1 if y0 < y1 else -1
    e = dx + dy
    while True:
        for tx in range(-thick//2, thick//2+1):
            for ty in range(-thick//2, thick//2+1):
                px, py = x0+tx, y0+ty
                if 0 <= px < width and 0 <= py < height:
                    i = (height-1-py)*width*3 + px*3
                    if 0 <= i < len(pixels)-2:
                        pixels[i], pixels[i+1], pixels[i+2] = r, g, b
        if x0 == x1 and y0 == y1: break
        e2 = 2*e
        if e2 >= dy: e += dy; x0 += sx
        if e2 <= dx: e += dx; y0 += sy


def rect(pixels, width, height, x, y, w, h, r, g, b, fill=False):
    """Rectangle at (x,y) size w*h."""
    if fill:
        for dy in range(max(0,y), min(height,y+h)):
            for dx in range(max(0,x), min(width,x+w)):
                i = (height-1-dy)*width*3 + dx*3
                if 0 <= i < len(pixels)-2:
                    pixels[i], pixels[i+1], pixels[i+2] = r, g, b
    else:
        for dx in range(max(0,x), min(width,x+w)):
            for ty in range(-1,2):
                py = y+ty
                if 0 <= py < height and 0 <= dx < width:
                    i = (height-1-py)*width*3 + dx*3
                    if 0 <= i < len(pixels)-2:
                        pixels[i], pixels[i+1], pixels[i+2] = r, g, b
        for dy in range(max(0,y), min(height,y+h)):
            for tx in range(-1,2):
                px = x+tx
                if 0 <= px < width and 0 <= dy < height:
                    i = (height-1-dy)*width*3 + px*3
                    if 0 <= i < len(pixels)-2:
                        pixels[i], pixels[i+1], pixels[i+2] = r, g, b


def text(pixels, width, height, text, x, y, r, g, b, scale=1):
    """Tiny 5x3 pixel font (only digits, '-', '.' supported)."""
    FONT = {
        '0': [[1,1,1],[1,0,1],[1,0,1],[1,0,1],[1,1,1]],
        '1': [[0,1,0],[1,1,0],[0,1,0],[0,1,0],[1,1,1]],
        '2': [[1,1,1],[0,0,1],[1,1,1],[1,0,0],[1,1,1]],
        '3': [[1,1,1],[0,0,1],[1,1,1],[0,0,1],[1,1,1]],
        '4': [[1,0,1],[1,0,1],[1,1,1],[0,0,1],[0,0,1]],
        '5': [[1,1,1],[1,0,0],[1,1,1],[0,0,1],[1,1,1]],
        '6': [[1,1,1],[1,0,0],[1,1,1],[1,0,1],[1,1,1]],
        '7': [[1,1,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1]],
        '8': [[1,1,1],[1,0,1],[1,1,1],[1,0,1],[1,1,1]],
        '9': [[1,1,1],[1,0,1],[1,1,1],[0,0,1],[1,1,1]],
        '-': [[0,0,0],[0,0,0],[1,1,1],[0,0,0],[0,0,0]],
        '.': [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[1,1,1]],
        'e': [[1,1,1],[1,0,0],[1,1,1],[1,0,0],[1,1,1]],
        'm': [[0,0,0],[1,1,1],[1,1,1],[1,0,1],[1,0,1]],
        'n': [[0,0,0],[1,1,1],[1,0,1],[1,0,1],[1,0,1]],
        'k': [[0,0,0],[1,0,1],[1,1,0],[1,0,1],[1,0,1]],
        's': [[1,1,1],[1,0,0],[1,1,1],[0,0,1],[1,1,1]],
        't': [[0,1,0],[1,1,1],[1,1,1],[0,1,0],[0,1,0]],
        'r': [[0,0,0],[1,1,0],[1,0,0],[1,0,0],[1,0,0]],
        ' ': [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]],
        'N': [[1,0,1],[1,1,1],[1,0,1],[1,0,1],[1,0,1]],
        'E': [[1,1,1],[1,0,0],[1,1,1],[1,0,0],[1,1,1]],
        'U': [[1,0,1],[1,0,1],[1,0,1],[1,0,1],[1,1,1]],
        'S': [[1,1,1],[1,0,0],[1,1,1],[0,0,1],[1,1,1]],
        'R': [[1,1,0],[1,0,1],[1,1,0],[1,0,1],[1,0,1]],
    }
    px = x
    for ch in text:
        glyph = FONT.get(ch, FONT.get(ch.lower(), FONT[' ']))
        for gy, row in enumerate(glyph):
            for gx, on in enumerate(row):
                if on:
                    for sx in range(scale):
                        for sy in range(scale):
                            wx = px + gx*scale + sx
                            wy = y - gy*scale + sy
                            if 0 <= wx < width and 0 <= wy < height:
                                i = (height-1-wy)*width*3 + wx*3
                                if 0 <= i < len(pixels)-2:
                                    pixels[i], pixels[i+1], pixels[i+2] = r, g, b
        px += 6 * scale


def timeseries_plot(vals_dict, labels, colors, out_png, title):
    """vals_dict: {field_name: [values]}, labels/colors for each."""
    W, H = 900, 300
    bg = (255, 255, 255)
    pixels = bytearray(W * H * 3)
    for i in range(W * H):
        pixels[i*3:i*3+3] = bytes(bg)

    # Count total points
    total = max(len(v) for v in vals_dict.values()) if vals_dict else 0
    if total == 0:
        print(f"  No data in {out_png}"); return

    # Determine x range
    MARGIN_LEFT = 70
    MARGIN_RIGHT = 20
    MARGIN_TOP = 20
    MARGIN_BOTTOM = 40
    plot_w = W - MARGIN_LEFT - MARGIN_RIGHT
    plot_h = H - MARGIN_TOP - MARGIN_BOTTOM

    # Border
    rect(pixels, W, H, MARGIN_LEFT, MARGIN_TOP, plot_w, plot_h, 0, 0, 0)
    # Grid lines
    for i in range(5):
        y = MARGIN_TOP + int(plot_h * i / 4)
        hline(pixels, W, H, MARGIN_LEFT, y, plot_w, 220, 220, 220)
        # Y label
        for digit in str(i*25):
            pass

    # Plot each series
    for (field, vals), color in zip(vals_dict.items(), colors):
        vals = list(vals)
        if not vals: continue
        valid = [(i, v) for i, v in enumerate(vals) if v == v]  # no NaN
        if not valid: continue
        mn = min(v for _, v in valid)
        mx = max(v for _, v in valid)
        if mx == mn:
            mx = mn + 0.001
        for j in range(len(valid)-1):
            i0, v0 = valid[j]
            i1, v1 = valid[j+1]
            x0 = int(MARGIN_LEFT + i0/(total-1)*plot_w)
            x1 = int(MARGIN_LEFT + i1/(total-1)*plot_w)
            y0 = int(MARGIN_TOP + plot_h - (v0-mn)/(mx-mn)*plot_h)
            y1 = int(MARGIN_TOP + plot_h - (v1-mn)/(mx-mn)*plot_h)
            line(pixels, W, H, x0, y0, x1, y1, *color, 1)

    # Legend
    for i, (lbl, color) in enumerate(zip(labels, colors)):
        lx = MARGIN_LEFT + 10 + i * 80
        ly = MARGIN_TOP + plot_h + 10
        rect(pixels, W, H, lx, ly, 20, 10, *color, True)
        text(pixels, W, H, lbl, lx+25, ly+8, *color, 1)

    # Title
    text(pixels, W, H, title, W//2 - 60, 12, 40, 40, 40, 1)

    write_png(W, H, bytes(pixels), out_png)
    sz = os.path.getsize(out_png)
    print(f"  PNG {W}x{H} ({sz//1024}KB): {Path(out_png).name}")


def bar_chart(stats, out_png, title):
    """Draw a simple bar chart PNG."""
    W, H = 500, 280
    bg = (255, 255, 255)
    pixels = bytearray(W * H * 3)
    for i in range(W * H):
        pixels[i*3:i*3+3] = bytes(bg)

    keys = list(stats.keys())
    vals = [v * 100 for v in stats.values()]
    n = len(keys)
    BAR_W = min(80, (W - 60) // max(n, 1) - 20)
    GAP = ((W - 60) - n * BAR_W) // (n + 1)

    rect(pixels, W, H, 50, 20, W-55, H-45, 0, 0, 0)
    bar_colors = [(220,60,60),(60,180,60),(60,60,220),(200,100,30)]

    for i, (k, v, col) in enumerate(zip(keys, vals, bar_colors)):
        bx = 50 + GAP + i * (BAR_W + GAP)
        v_clean = v if (v == v) else 0.0  # filter NaN
        max_v = max(vv for vv in vals if vv == vv) if any(vv == vv for vv in vals) else 0.001
        bh = min(int(v_clean / max(max_v, 0.001) * (H - 55)), H - 55)
        bh = max(bh, 2)
        rect(pixels, W, H, bx, H-30-bh, BAR_W, bh, *col, True)
        # Value label
        vtxt = f"{v:.1f}"
        text(pixels, W, H, vtxt, bx + BAR_W//2 - len(vtxt)*3, H-30-bh-12, *col, 1)
        # Key label
        text(pixels, W, H, str(k), bx + BAR_W//2 - len(str(k))*3, H-20, 60,60,60, 1)

    text(pixels, W, H, title, 55, 12, 30, 30, 30, 1)
    write_png(W, H, bytes(pixels), out_png)
    sz = os.path.getsize(out_png)
    print(f"  PNG {W}x{H} ({sz//1024}KB): {Path(out_png).name}")


def satcount_plot(sat_counts, out_png, title):
    W, H = 500, 180
    bg = (255, 255, 255)
    pixels = bytearray(W * H * 3)
    for i in range(W * H):
        pixels[i*3:i*3+3] = bytes(bg)

    MARGIN_LEFT = 60; MARGIN_TOP = 15; MARGIN_BOTTOM = 35
    plot_w = W - MARGIN_LEFT - 20
    plot_h = H - MARGIN_TOP - MARGIN_BOTTOM
    n = len(sat_counts)
    if n == 0:
        print(f"  No satcount data"); return

    rect(pixels, W, H, MARGIN_LEFT, MARGIN_TOP, plot_w, plot_h, 0, 0, 0)
    mn, mx = min(sat_counts), max(sat_counts)
    if mx == mn: mx = mn + 1

    for i in range(n-1):
        x0 = MARGIN_LEFT + int(i/(n-1)*plot_w)
        x1 = MARGIN_LEFT + int((i+1)/(n-1)*plot_w)
        y0 = MARGIN_TOP + int((mx - sat_counts[i])/(mx-mn)*plot_h)
        y1 = MARGIN_TOP + int((mx - sat_counts[i+1])/(mx-mn)*plot_h)
        line(pixels, W, H, x0, y0, x1, y1, 0, 100, 200, 1)

    avg_txt = f"avg={sum(sat_counts)/n:.0f}"
    text(pixels, W, H, title, MARGIN_LEFT + 5, 12, 30, 30, 30, 1)
    text(pixels, W, H, avg_txt, MARGIN_LEFT + 5, H-12, 0, 80, 160, 1)
    write_png(W, H, bytes(pixels), out_png)
    print(f"  PNG {W}x{H}: {Path(out_png).name}")


def orbit_plot(x_vals, y_vals, z_vals, out_png, title):
    """Plot ECEF orbit projection."""
    W, H = 480, 320
    bg = (255, 255, 255)
    pixels = bytearray(W * H * 3)
    for i in range(W * H):
        pixels[i*3:i*3+3] = bytes(bg)

    x, y, z = list(x_vals), list(y_vals), list(z_vals)
    cx = sum(x)/len(x); cy = sum(y)/len(y); cz = sum(z)/len(z)
    rx = (max(x)-min(x))/2 if max(x)!=min(x) else 1
    ry = (max(y)-min(y))/2 if max(y)!=min(y) else 1
    rz = (max(z)-min(z))/2 if max(z)!=min(z) else 1
    r = max(rx, ry, rz, 1)

    def proj(xx, yy, zz):
        return (int((xx-cx)/r*110+W//2), int((yy-cy)/r*110+H//2))

    # X axis
    px, py = proj(0, 0, 0)
    line(pixels, W, H, px, py, px+60, py, 220, 50, 50, 1)
    text(pixels, W, H, 'X+', px+62, py+4, 220, 50, 50, 1)
    line(pixels, W, H, px, py, px-40, py, 180, 80, 80, 1)
    # Y axis
    line(pixels, W, H, px, py, px, py+40, 50, 180, 50, 1)
    text(pixels, W, H, 'Y+', px+2, py+42, 50, 180, 50, 1)

    # Orbit trace
    for i in range(len(x)-1):
        x0, y0 = proj(x[i], y[i], z[i])
        x1, y1 = proj(x[i+1], y[i+1], z[i+1])
        line(pixels, W, H, x0, y0, x1, y1, 0, 60, 180, 1)

    # Start dot
    if x:
        sx, sy = proj(x[0], y[0], z[0])
        rect(pixels, W, H, sx-3, sy-3, 7, 7, 220, 80, 0, True)

    text(pixels, W, H, title, 5, 14, 30, 30, 30, 1)
    write_png(W, H, bytes(pixels), out_png)
    print(f"  PNG {W}x{H}: {Path(out_png).name}")


def convergence_plot(x_vals, y_vals, z_vals, out_png, title):
    """Plot position convergence."""
    W, H = 480, 200
    bg = (255, 255, 255)
    pixels = bytearray(W * H * 3)
    for i in range(W * H):
        pixels[i*3:i*3+3] = bytes(bg)

    x, y, z = list(x_vals), list(y_vals), list(z_vals)
    n = len(x)
    if n < 2:
        print(f"  Not enough data for convergence"); return

    dr = [((x[i]-x[0])**2+(y[i]-y[0])**2+(z[i]-z[0])**2)**0.5 for i in range(n)]
    MARGIN_LEFT = 60; MARGIN_TOP = 15; MARGIN_BOTTOM = 35
    plot_w = W - MARGIN_LEFT - 20
    plot_h = H - MARGIN_TOP - MARGIN_BOTTOM
    mx, mn = max(dr), min(dr)
    if mx == mn: mx = mn + 0.001

    rect(pixels, W, H, MARGIN_LEFT, MARGIN_TOP, plot_w, plot_h, 0, 0, 0)
    for i in range(n-1):
        x0 = MARGIN_LEFT + int(i/(n-1)*plot_w)
        x1 = MARGIN_LEFT + int((i+1)/(n-1)*plot_w)
        y0 = MARGIN_TOP + int((mx-dr[i])/(mx-mn)*plot_h)
        y1 = MARGIN_TOP + int((mx-dr[i+1])/(mx-mn)*plot_h)
        line(pixels, W, H, x0, y0, x1, y1, 200, 80, 0, 1)

    text(pixels, W, H, title, 5, 12, 30, 30, 30, 1)
    write_png(W, H, bytes(pixels), out_png)
    print(f"  PNG {W}x{H}: {Path(out_png).name}")


def main():
    import csv as csvmod

    out_dir = Path('/workspace/gnss_pod/output/plots')
    out_dir.mkdir(exist_ok=True)

    # Find latest GRACE-FO CSV
    csv_files = sorted(Path('/workspace/gnss_pod/output').glob('ppp_gracefo1*.csv'))
    if not csv_files:
        print("No GRACE-FO CSV found"); return
    csv_path = str(csv_files[-1])
    print(f"Reading: {csv_path}")

    times, dE, dN, dU, dX, dY, dZ, X, Y, Z, n_sats = [], [], [], [], [], [], [], [], [], [], []
    with open(csv_path) as f:
        for row in csvmod.DictReader(f):
            try:
                times.append(row.get('time', ''))
                dE.append(float(row.get('dE', 0)))
                dN.append(float(row.get('dN', 0)))
                dU.append(float(row.get('dU', 0)))
                dX.append(float(row.get('dX', 0)))
                dY.append(float(row.get('dY', 0)))
                dZ.append(float(row.get('dZ', 0)))
                X.append(float(row.get('X', 0)))
                Y.append(float(row.get('Y', 0)))
                Z.append(float(row.get('Z', 0)))
                n_sats.append(int(float(row.get('n_sats', 0))))
            except Exception:
                pass

    if not times:
        print("No valid rows"); return

    import numpy as np
    dE, dN, dU = np.array(dE), np.array(dN), np.array(dU)
    X, Y, Z = np.array(X), np.array(Y), np.array(Z)
    n_sats = np.array(n_sats)
    d3d = np.sqrt(dE**2 + dN**2 + dU**2)

    # 1. ENU time series
    timeseries_plot(
        {'dE (cm)': dE, 'dN (cm)': dN, 'dU (cm)': dU},
        ['E', 'N', 'U'],
        [(220,50,50),(50,180,50),(50,50,220)],
        str(out_dir / '01_enu_errors.png'),
        'ENU Position Errors — GRACE-FO PPP'
    )

    # 2. XYZ ECEF errors
    timeseries_plot(
        {'dX (m)': dX, 'dY (m)': dY, 'dZ (m)': dZ},
        ['X', 'Y', 'Z'],
        [(180,0,180),(0,150,150),(200,100,0)],
        str(out_dir / '02_xyz_errors.png'),
        'ECEF Coordinate Errors — GRACE-FO PPP'
    )

    # 3. Stats bar chart
    bar_chart(
        {'3D RMS': float(d3d.mean()),
         'E RMS': float(np.sqrt(np.nanmean(dE**2))),
         'N RMS': float(np.sqrt(np.nanmean(dN**2))),
         'U RMS': float(np.sqrt(np.nanmean(dU**2)))},
        str(out_dir / '04_stats.png'),
        'PPP RMS Error (cm)'
    )

    # 4. Satellite count
    satcount_plot(list(n_sats), str(out_dir / '05_satcount.png'), 'Visible GPS Sats')

    # 5. Orbit 3D projection
    orbit_plot(X, Y, Z, str(out_dir / '06_orbit.png'), 'GRACE-FO Orbit (ECEF)')

    # 6. Convergence
    convergence_plot(X, Y, Z, str(out_dir / '07_convergence.png'), 'Position Convergence (m)')

    print(f"\nAll PNG plots saved to {out_dir}/")
    for p in sorted(out_dir.glob('*.png')):
        print(f"  {p.name} ({p.stat().st_size//1024}KB)")


if __name__ == '__main__':
    main()
