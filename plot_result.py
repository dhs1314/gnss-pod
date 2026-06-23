#!/usr/bin/env python3
"""GRACE-FO PPP 误差 PNG — 纯 Python，无依赖"""
import csv, math, os, sys, zlib, struct

def rms(v):
    n = len(v)
    return math.sqrt(sum(x*x for x in v) / n) if n else 0.0

def read_csv(path):
    times, dE, dN, dU = [], [], [], []
    t0 = None
    with open(path) as f:
        for row in csv.DictReader(f):
            try:
                de, dn, du = float(row['dE_m']), float(row['dN_m']), float(row['dU_m'])
                if math.isnan(de + dn + du):
                    continue
                d3 = math.sqrt(de**2 + dn**2 + du**2)
                if d3 > 1.0:
                    continue
                h = int(row['time'][11:13])
                m = int(row['time'][14:16])
                s = int(row['time'][17:19])
                sec = h * 3600 + m * 60 + s
                if t0 is None:
                    t0 = sec
                times.append((sec - t0) / 3600.0)
                dE.append(de * 100)
                dN.append(dn * 100)
                dU.append(du * 100)
            except Exception:
                continue
    return times, dE, dN, dU

def write_png(W, H, pixels, path):
    def chunk(tag, data):
        c = tag + data
        return struct.pack('>I', len(data)) + c + struct.pack('>I', zlib.crc32(c) & 0xffffffff)
    sig = b'\x89PNG\r\n\x1a\n'
    ihdr = struct.pack('>IIBBBBB', W, H, 8, 2, 0, 0, 0)
    raw_rows = []
    for row in pixels:
        flat = b''
        for px in row:
            if isinstance(px, list):
                flat += bytes(px)
            else:
                flat += bytes([px])
        raw_rows.append(b'\x00' + flat)
    idat = zlib.compress(b''.join(raw_rows), 6)
    with open(path, 'wb') as f:
        f.write(sig + chunk(b'IHDR', ihdr) + chunk(b'IDAT', idat) + chunk(b'IEND', b''))

def plot(times, dE, dN, dU, out_png):
    W, H = 960, 640
    ML, MR, MT, MB = 72, 20, 56, 44
    title_h = 36
    gap = 12
    pw = W - ML - MR
    ph = (H - MT - MB - 2 * gap - title_h) // 3
    tmin, tmax = min(times), max(times)
    ymax = max(max(abs(v) for v in dE), max(abs(v) for v in dN), max(abs(v) for v in dU)) * 1.25
    ymax = max(ymax, 0.5)

    def px(tv):
        return ML + (tv - tmin) / (tmax - tmin + 1e-9) * (pw - 1)

    def py(tv, top):
        return top + ph - 1 - int((tv + ymax) / (2 * ymax) * (ph - 1))

    def py0(top):
        return py(0.0, top)

    def sp(x, y, c, img):
        if 0 <= x < W and 0 <= y < H:
            img[y][x] = c

    def hl(x0, x1, y, c, img, w=1):
        for xi in range(int(x0), int(x1) + 1):
            for dw in range(-w // 2, w // 2 + 1):
                sp(xi, y + dw, c, img)

    def drow(top, data, hexcol):
        r = int(hexcol[1:3], 16)
        g = int(hexcol[3:5], 16)
        b = int(hexcol[5:7], 16)
        for y in range(top, top + ph):
            for x in range(ML, ML + pw):
                img[y][x] = [242, 242, 242]
        for tv in times:
            xi = int(px(tv))
            if ML <= xi < ML + pw:
                for yi in range(top, top + ph):
                    img[yi][xi] = [225, 225, 225]
        yz = py0(top)
        hl(ML, ML + pw - 1, yz, [130, 130, 130], img, 1)
        px_ = None
        for i in range(len(times)):
            xi = int(px(times[i]))
            yi = py(data[i], top)
            yi = max(top, min(top + ph - 1, yi))
            if px_ is not None and abs(xi - px_) <= 3:
                hl(px_, xi, yi, [r, g, b], img, 1)
            else:
                sp(xi, yi, [r, g, b], img)
            px_ = xi

    img = [[255] * W for _ in range(H)]
    for y in range(MT - title_h, MT):
        for x in range(W):
            img[y][x] = [250, 250, 250]
    for y in range(H - MB, H):
        for x in range(W):
            img[y][x] = [250, 250, 250]

    rows = [
        (dE, '#1565C0'),
        (dN, '#2E7D32'),
        (dU, '#C62828'),
    ]
    for i, (data, col) in enumerate(rows):
        ry = MT + title_h + gap + i * (ph + gap)
        drow(ry, data, col)
        vl(ML, ry, ry + ph, [100, 100, 100], img, 2)
        vl(ML + pw, ry, ry + ph, [100, 100, 100], img, 2)
        hl(ML, ML + pw, ry, [100, 100, 100], img, 2)
        hl(ML, ML + pw, ry + ph, [100, 100, 100], img, 2)

    write_png(W, H, img, out_png)
    return W, H

def vl(x, y0, y1, c, img, w=1):
    for yi in range(int(y0), int(y1) + 1):
        for dw in range(-w // 2, w // 2 + 1):
            if 0 <= x + dw < 960 and 0 <= yi < 640:
                img[yi][x + dw] = c

if __name__ == '__main__':
    p = sys.argv[1] if len(sys.argv) > 1 else '/workspace/gnss_pod/output/ppp_vs_gnv1b_2024_0429_2h.csv'
    out = p.replace('.csv', '_enu_clean.png')
    times, dE, dN, dU = read_csv(p)
    n = len(times)
    print(f'有效历元: {n}')
    r3 = rms([math.sqrt(dE[i]**2 + dN[i]**2 + dU[i]**2) for i in range(n)])
    print(f'E RMS:  {rms(dE):.2f} cm')
    print(f'N RMS:  {rms(dN):.2f} cm')
    print(f'U RMS:  {rms(dU):.2f} cm')
    print(f'3D RMS: {r3:.2f} cm')
    W, H = plot(times, dE, dN, dU, out)
    sz = os.path.getsize(out)
    print(f'PNG: {out} ({sz // 1024} KB, {W}x{H})')
