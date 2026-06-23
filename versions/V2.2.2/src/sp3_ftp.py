"""
WHU IGS ultra-rapid SP3 下载器
数据源: ftp://igs.gnsswhu.cn/pub/gps/products/{gps_week}/
工具: wget (目录索引) + curl EPSV (文件下载)
"""
import subprocess, gzip, pickle, os, re
from pathlib import Path
from datetime import datetime, timedelta

HOST = 'igs.gnsswhu.cn'
BASE = f'ftp://{HOST}'

def gps_wd(year, month, day):
    t = datetime(year, month, day)
    g0 = datetime(1980, 1, 6)
    tot = (t - g0).days
    return tot // 7, tot % 7, (t - datetime(t.year, 1, 1)).days + 1

def _wget(path, local, timeout=90):
    r = subprocess.run(
        ['wget', '-q', f'--timeout={timeout}',
         '--ftp-user=anonymous', '--ftp-password=ftp@x.com',
         '-O', local, f'{BASE}{path}'],
        capture_output=True, timeout=timeout + 10
    )
    return r.returncode == 0

def _curl_epsv(remote, local, timeout=120):
    r = subprocess.run(
        ['curl', '-s', '--insecure', '--ftp-ssl', '--epsv', '-P', '-',
         '-o', local, '-w', '%{http_code}',
         f'{BASE}{remote}'],
        capture_output=True, timeout=timeout
    )
    ok = r.stdout.strip() == '226'
    sz = os.path.getsize(local) if os.path.exists(local) else 0
    return ok, sz

def _get_dir_html(gps_week):
    html_file = f'/tmp/whu_week{gps_week}.html'
    ok = _wget(f'/pub/gps/products/{gps_week:04d}/', html_file, timeout=90)
    if not ok or not os.path.exists(html_file) or os.path.getsize(html_file) < 1000:
        return ''
    return open(html_file, errors='replace').read()

def _parse_sp3_files(html):
    links = re.findall(r'href="([^"]+\.SP3\.gz)"', html, re.I)
    files = sorted(set(n.split('/')[-1] for n in links if '.SP3.gz' in n.upper()))
    return files

def _select_best(files, gw, doy):
    whu_ultra = [f for f in files if 'WUH0' in f.upper() and 'ULT' in f.upper()]
    whu_mgex  = [f for f in files if 'WUH0' in f.upper() and 'FIN' in f.upper()]
    wuh_any   = [f for f in files if 'WUH0' in f.upper()]
    any_sp3   = [f for f in files if 'SP3' in f.upper()]

    def doy_score(f):
        m = re.search(r'_(\d{11})_', f)
        if m:
            file_doy = int(m.group(1)[4:7])
            return abs(file_doy - doy)
        return 999

    for lst in [whu_ultra, whu_mgex, wuh_any]:
        if lst:
            return min(lst, key=doy_score)
    igr = [f for f in any_sp3 if 'IGR' in f.upper() or 'IGU' in f.upper()]
    if igr: return min(igr, key=doy_score)
    mgex = [f for f in any_sp3 if 'MGX' in f.upper()]
    if mgex: return min(mgex, key=doy_score)
    return any_sp3[0] if any_sp3 else None

def _parse_sp3(text):
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

    print(f"  WHU FTP: week={gw} DOY={doy}")
    html = _get_dir_html(gw)
    if not html:
        print(f"  [错误] 目录索引获取失败"); return None

    files = _parse_sp3_files(html)
    print(f"  目录 SP3: {len(files)} 个")
    for f in files[:3]: print(f"    {f}")
    if len(files) > 3: print(f"    ... (+{len(files)-3})")

    best = _select_best(files, gw, doy)
    if not best:
        print(f"  [错误] 无可用 SP3"); return None
    print(f"  选择: {best}")

    remote = f'/pub/gps/products/{gw:04d}/{best}'
    tmp = f'/tmp/whu_sp3_{best}'
    print(f"  下载...")
    ok, sz = _curl_epsv(remote, tmp)
    if not ok or sz < 1000:
        print(f"  [错误] 下载失败: ok={ok} sz={sz}"); return None
    print(f"  完成: {sz} bytes")

    try:
        with open(tmp,'rb') as f: raw = f.read()
        if raw[:2] == b'\x1f\x8b':
            text = gzip.decompress(raw).decode('ascii', errors='replace')
        elif raw[:2] == b'PK':
            import zipfile
            with zipfile.ZipFile(tmp) as zf:
                text = zf.read(zf.namelist()[0]).decode('ascii', errors='replace')
        else:
            text = raw.decode('ascii', errors='replace')
    except Exception as e:
        print(f"  [错误] 解压: {e}"); return None

    ts, epochs = _parse_sp3(text)
    if not ts:
        print(f"  [错误] SP3 解析失败"); return None
    result = {'epochs': epochs, 'ts': ts}
    cache.parent.mkdir(parents=True, exist_ok=True)
    with open(cache,'wb') as f: pickle.dump(result, f)
    try: os.remove(tmp)
    except: pass
    n_gps = sum(1 for s in (list(epochs.values())[0] if epochs else {}) if str(s).startswith('G'))
    print(f"  ✅ SP3: {len(ts)} epochs, {n_gps} GPS sats, {ts[0]} -> {ts[-1]}")
    return result

if __name__ == '__main__':
    import sys
    y, m, d = 2024, 5, 1
    if len(sys.argv) >= 2: y = int(sys.argv[1])
    if len(sys.argv) >= 3: m = int(sys.argv[2])
    if len(sys.argv) >= 4: d = int(sys.argv[3])
    r = load_whu_sp3_ftp(y, m, d)
    print('结果:', 'OK' if r else 'FAILED')
