#!/usr/bin/env python3
"""
自动探测 WHU SP3 正确文件名并下载
用法: python3 src/sp3_probe.py [year month day]
"""
import subprocess, gzip, pickle, os, re, sys, time
from pathlib import Path
from datetime import datetime, timedelta

HOST = 'igs.gnsswhu.cn'
BASE = f'ftp://{HOST}'

def gps_wd(year, m, d):
    t = datetime(year, m, d)
    g0 = datetime(1980, 1, 6)
    tot = (t - g0).days
    return tot // 7, tot % 7, (t - datetime(t.year,1,1)).days + 1

def probe_and_download(gw, gps_d, doy, timeout_per=15, max_workers=6):
    """
    探测正确文件名，然后下载。
    """
    # 生成候选文件名（WHU ultra-rapid + IGS 标准格式）
    cands = []
    for hh in [0, 6, 12, 18]:
        # WHU ultra-rapid (最可能): whug{wwww}{doy}{hh}.sp3.gz
        cands.append((f'whug{gw:04d}{doy:03d}{hh:02d}.sp3.gz',
                      f'/pub/gps/products/{gw:04d}/whug{gw:04d}{doy:03d}{hh:02d}.sp3.gz'))
        # WHU 00h daily
        cands.append((f'whug{gw:04d}{doy:03d}00.sp3.gz',
                      f'/pub/gps/products/{gw:04d}/whug{gw:04d}{doy:03d}00.sp3.gz'))
        # Without .gz
        cands.append((f'whug{gw:04d}{doy:03d}{hh:02d}.sp3',
                      f'/pub/gps/products/{gw:04d}/whug{gw:04d}{doy:03d}{hh:02d}.sp3'))
    # IGS ultra-rapid: wuh0igs{s}.sp3.gz (Bernese format)
    for s in [f'{doy}0', f'{doy}6', f'{doy}12', f'{doy}18', f'{doy}']:
        cands.append((f'wuh0igs{s}.sp3.gz', f'/pub/gps/products/{gw:04d}/wuh0igs{s}.sp3.gz'))
        cands.append((f'WUH0IGS{s}.SP3.gz', f'/pub/gps/products/{gw:04d}/WUH0IGS{s}.SP3.gz'))
    # WHU products dir alternative
    for hh in [0, 6, 12, 18]:
        cands.append((f'whug{gw}{doy:03d}{hh:02d}.sp3.gz',
                      f'/pub/gps/products/{gw:04d}/whug{gw}{doy:03d}{hh:02d}.sp3.gz'))
    # Try without /products/ prefix
    cands.append((f'whug{gw:04d}{doy:03d}00.sp3.gz', f'/pub/gps/{gw:04d}/whug{gw:04d}{doy:03d}00.sp3.gz'))

    # 去重
    seen = set(); unique = []
    for name, path in cands:
        if path not in seen:
            seen.add(path); unique.append((name, path))

    print(f"  探测 {len(unique)} 个候选文件名...")
    for name, path in unique:
        tmp = f'/tmp/whu_probe_{name}'
        r = subprocess.run(
            ['curl', '-s', '--insecure', '--ftp-ssl', '--epsv', '-P', '-',
             '-o', tmp, '-w', '%{http_code}', '-m', str(timeout_per),
             f'{BASE}{path}'],
            capture_output=True, text=True, timeout=timeout_per + 5
        )
        code = r.stdout.strip()
        sz = os.path.getsize(tmp) if os.path.exists(tmp) else 0
        if code == '226' and sz > 5000:
            print(f"  ✅ 成功: {name} ({sz}B)")
            # 验证内容
            try:
                with open(tmp, 'rb') as f: raw = f.read()
                if raw[:2] == b'\x1f\x8b':
                    text = gzip.decompress(raw).decode('ascii', errors='replace')
                else:
                    text = raw.decode('ascii', errors='replace')
                if 'SP3' in text or '+' in text:
                    print(f"  ✅ 有效 SP3! 前3行:")
                    for l in text.split('\n')[:3]: print(f"     {l[:80]}")
                    return tmp, name, path, text
            except Exception as e:
                print(f"  ❌ 解析失败: {e}")
        else:
            print(f"  ❌ {code} {sz}B: {name}")
    return None, None, None, None

def parse_sp3_text(text):
    epochs = {}; cur_sats = []; cur_ep = None
    for line in text.split('\n'):
        if line.startswith('+'):
            cur_sats = [s.strip() for s in line[9:].split() if s.strip()]; continue
        if line.startswith('*') and len(line) > 1:
            p = line[1:].split()
            if len(p) >= 6:
                try: cur_ep = datetime(int(p[0]),int(p[1]),int(p[2]),int(p[3]),int(p[4]),int(float(p[5])))
                except: cur_ep = None
                if cur_ep and cur_ep not in epochs: epochs[cur_ep] = {}
            continue
        if line.startswith('EOF') or not line.strip() or not cur_sats or cur_ep is None: continue
        parts = line.split()
        if len(parts) < 4: continue
        sv = parts[0]
        if sv not in cur_sats: continue
        try:
            x = float(parts[1])*1000; y = float(parts[2])*1000; z = float(parts[3])*1000
            clk = float(parts[4])*1e-6*299792458.0 if len(parts)>4 and parts[4] not in ('0','*') else 0.0
            epochs[cur_ep][sv] = [x,y,z,clk]
        except: continue
    return sorted(epochs.keys()), epochs

def load_whu_sp3_ftp(year, month, day, data_dir='./data'):
    gw, gps_d, doy = gps_wd(year, month, day)
    cache = Path(data_dir)/str(year)/f'{doy:03d}'/'whu_sp3.pkl'
    if cache.exists():
        try:
            with open(cache,'rb') as f: d=pickle.load(f)
            print(f"  [缓存] WHU SP3 ({len(d['ts'])} epochs)")
            return d
        except: pass
    print(f"  WHU FTP探测: week={gw} DOY={doy}")
    tmp, name, path, text = probe_and_download(gw, gps_d, doy)
    if not tmp:
        print(f"  [错误] 所有文件名探测失败"); return None
    ts, epochs = parse_sp3_text(text)
    if not ts:
        print(f"  [错误] SP3 解析失败"); return None
    result = {'epochs': epochs, 'ts': ts}
    cache.parent.mkdir(parents=True, exist_ok=True)
    with open(cache, 'wb') as f: pickle.dump(result, f)
    try: os.remove(tmp)
    except: pass
    n_gps = sum(1 for s in (list(epochs.values())[0] if epochs else {}) if str(s).startswith('G'))
    print(f"  ✅ SP3 OK: {len(ts)} epochs, {n_gps} GPS sats, {ts[0]} -> {ts[-1]}")
    return result

if __name__ == '__main__':
    y, m, d = 2024, 5, 1
    if len(sys.argv) >= 2: y = int(sys.argv[1])
    if len(sys.argv) >= 3: m = int(sys.argv[2])
    if len(sys.argv) >= 4: d = int(sys.argv[3])
    r = load_whu_sp3_ftp(y, m, d)
    print('结果:', 'OK' if r else 'FAILED')
