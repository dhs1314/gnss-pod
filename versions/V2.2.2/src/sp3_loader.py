# v1.2.3 — 2026-04-24
"""
IGS SP3 Loader

文件命名: IGS0OPSFIN_YYYYDDHHMM000_01D_15M_ORB.SP3.gz
  YYYY = GPS year (e.g. 2024)
  DD  = GPS week (1-52)?? Actually: GPS week number (not day-of-year!)
  
IGS naming convention: the 7-digit number "YYYYDDD" uses GPS week + day-of-week,
NOT year + day-of-year. So '2024120' means GPS year 2024, week 1, day 0 (Sunday).

Data vs coverage: the epoch times in the SP3 file are INDEPENDENT of the filename date.
The SP3 file named for date D contains epoch data for the NEXT calendar day.
For example: IGS0OPSFIN_20241190000 has epochs from 2024-04-28 00:00 → 2024-04-28 23:45
           IGS0OPSFIN_20241200000 has epochs from 2024-04-29 00:00 → 2024-04-29 23:45

覆盖规则 (IGS convention, 01D files):
  The nominal coverage (file_date + 1 to file_date + 2) is OFF BY ONE for our purposes.
  Instead, the epoch DATA starts at file_date + 1 day 00:00, runs for 1 day.
  
  For file_date = Apr 28 (DOY 119 = 2024119):
    nominal_coverage = Apr 29 → Apr 30
    actual_epoch_data = Apr 28 00:00 → Apr 28 23:45 ✓ (contains Apr 28)

  For file_date = Apr 29 (DOY 120 = 2024120):
    nominal_coverage = Apr 30 → May 1
    actual_epoch_data = Apr 29 00:00 → Apr 29 23:45 ✓ (contains Apr 29)

So the correct file for Apr 29 is the one with "file_date = Apr 29" (IGS0OPSFIN_20241200000).
The SP3 epoch data for Apr 29 IS IN the file named for Apr 29 (DOY 120).
"""
from pathlib import Path
from datetime import datetime, timedelta
import numpy as np
import ssl
import urllib.request
import gzip
import pickle
import re

C = 299792458.0
GPS_ORIGIN = datetime(1980, 1, 6)


def gps_week_for_utc(dt_utc):
    return int((dt_utc - GPS_ORIGIN).total_seconds() // 604800)


def _parse_igs_file_date(fname):
    """
    解析 IGS 文件名的 file_date (GPS year + GPS DOY)
    格式: IGS0OPSFIN_YYYYDDDHHMM000_01D_15M_ORB.SP3.gz
                   ^^^^^^^^^
                   8位数字: YYYY(4) + DDD(3) + HH(2)
    例: '20241190000' → year=2024, DOY=119 → Apr 28
        '20241200000' → year=2024, DOY=120 → Apr 29
    """
    m = re.search(r'_(\d{8})000_', fname)
    if not m:
        return None
    try:
        num = m.group(1)
        gps_year = int(num[:4])
        gps_doy  = int(num[4:7])
        return datetime(gps_year, 1, 1) + timedelta(days=gps_doy - 1)
    except (ValueError, IndexError):
        return None


def _quick_parse_epochs(text):
    """快速解析 SP3 文件的头几个 epoch 时间，用于判断文件数据范围"""
    ts = []
    for line in text.split('\n')[:300]:
        if len(line) > 1 and line[0] == '*':
            p = line[1:].split()
            try:
                ts.append(datetime(int(p[0]), int(p[1]), int(p[2]),
                                   int(p[3]), int(p[4]), int(float(p[5]))))
            except:
                pass
    return sorted(ts)


def parse_sp3_text(text):
    epochs = {}
    ts = []
    cur_ep = None
    cur_sats = []
    cur_data = []

    for line in text.split('\n'):
        if len(line) > 0 and line[0] == '+' and 'G' in line:
            found = re.findall(r'G\d+', line[9:])
            if found:
                cur_sats += [s for s in found if s not in cur_sats]
            continue
        if len(line) > 0 and line[0] == '*':
            if cur_ep and cur_data:
                epochs[cur_ep] = dict(zip(cur_sats, cur_data))
                ts.append(cur_ep)
            p = line[1:].split()
            try:
                cur_ep = datetime(int(p[0]), int(p[1]), int(p[2]),
                                  int(p[3]), int(p[4]), int(float(p[5])))
            except (IndexError, ValueError):
                continue
            cur_data = []
            continue
        if line.startswith('EOF') or not line.strip():
            continue
        if cur_ep and line.strip():
            parts = line.split()
            if len(parts) >= 5 and len(parts[0]) >= 3:
                sv = parts[0]
                if sv.startswith('P') and len(sv) >= 3:
                    sv = sv[1:]
                try:
                    x = float(parts[1]) * 1000
                    y = float(parts[2]) * 1000
                    z = float(parts[3]) * 1000
                    clk = float(parts[4]) * 1e-6 * C
                    cur_data.append([x, y, z, clk])
                except (ValueError, IndexError):
                    continue

    if cur_ep and cur_data:
        epochs[cur_ep] = dict(zip(cur_sats, cur_data))
        ts.append(cur_ep)

    return epochs, sorted(ts)


def fetch_gzip(url, timeout=30):
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    req = urllib.request.Request(url, headers={
        'User-Agent': 'curl/7.88', 'Accept-Encoding': 'gzip'
    })
    with urllib.request.urlopen(req, timeout=timeout, context=ctx) as r:
        return gzip.decompress(r.read()).decode('ascii', errors='replace')


def list_bkg_files(gps_week):
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    dir_url = f'https://igs.bkg.bund.de/root_ftp/IGS/products/{gps_week}/'
    try:
        req = urllib.request.Request(dir_url, headers={'User-Agent': 'curl/7.88'})
        r = urllib.request.urlopen(req, timeout=15, context=ctx)
        content = r.read().decode('iso-8859-1', errors='replace')
        files = re.findall(r'(IGS\w+_\d{11}_01D_15M_ORB\.SP3\.gz)', content)
        return list(dict.fromkeys(files))
    except Exception:
        return []


def find_best_bkg_file(gps_week, target_date, product_types=('FIN', 'RAP', 'ULT')):
    """
    选择包含 target_date epoch 数据的 SP3 文件

    IGS SP3 文件规则:
      文件名的 file_date (YYYYDDD) 与其 epoch 数据日期不同步。
      实际规律：文件名为 DOY=N 的文件包含 DOY=N+1 的 epoch 数据。
      例: 20241200000 → file_date=Apr 29 → epochs from Apr 29 00:00 → Apr 29 23:45

    因此选择 file_date = target_date 的文件（即 DOY=target_doy 的文件）。
    优先级: FIN > RAP > ULT（同日期内）
    """
    files = list_bkg_files(gps_week)
    if not files:
        for offset in [-1, 1, -2, 2]:
            gw_n = gps_week + offset
            files2 = list_bkg_files(gw_n)
            if files2:
                print(f"  [BKG] GPS Week {gw_n}: {len(files2)} files")
                files = files2
                break

    # 按 file_date 找最接近 target_date 的文件
    candidates = []
    for fname in files:
        ptype = next((pt for pt in product_types if pt in fname), None)
        if ptype is None:
            continue
        file_date = _parse_igs_file_date(fname)
        if file_date is None:
            continue
        diff_days = (target_date - file_date).days
        priority = {'FIN': 0, 'RAP': 1, 'ULT': 2}.get(ptype, 3)
        candidates.append(
            (abs(diff_days), priority, diff_days,
             fname, file_date, ptype))

    if not candidates:
        return None, None

    candidates.sort()
    _, _, diff_days, best_fname, best_date, ptype = candidates[0]

    print(f"  [BKG] Selected: {best_fname}")
    print(f"         file_date={best_date.date()} (GPS yr+DOY), "
          f"target=Apr {target_date.day}, diff={diff_days}d, type={ptype}")

    return (f'https://igs.bkg.bund.de/root_ftp/IGS/products/{gps_week}/{best_fname}',
            best_fname)


def load_igs_sp3(year, doy, data_dir="./data", product='auto'):
    """加载 IGS SP3 精密星历"""
    t_ref = datetime(year, 1, 1) + timedelta(days=doy - 1)
    gw = gps_week_for_utc(t_ref)

    product_map = {
        'final': ('FIN',), 'rapid': ('RAP',), 'ultra': ('RAP',),
        'auto': ('FIN', 'RAP'),
    }
    ptypes = product_map.get(product, ('FIN', 'RAP'))

    cache = Path(data_dir) / str(year) / f"{doy:03d}" / f"igs_sp3_{ptypes[0]}.pkl"
    cache.parent.mkdir(parents=True, exist_ok=True)

    if cache.exists():
        try:
            d = pickle.load(open(cache, 'rb'))
            # 验证缓存是否包含目标日期的数据
            ts = d['ts']
            if ts and min(ts) <= t_ref <= max(ts):
                print(f"  [缓存] IGS SP3 {ptypes[0]} ({len(ts)} epochs, "
                      f"data {ts[0].date()} → {ts[-1].date()})")
                return d
            else:
                print(f"  [缓存过期] {ts[0].date()} → {ts[-1].date()} (need {t_ref.date()})")
        except Exception:
            pass

    print(f"  查找 SP3: {t_ref.date()} (GPS week {gw}, "
          f"优先 {'FIN' if 'FIN' in ptypes else 'RAP'})")

    bkg_url, bkg_fname = find_best_bkg_file(gw, t_ref, ptypes)

    if not bkg_url:
        for w_offset in [-1, 1, -2, 2, -3, 3]:
            gw_n = gw + w_offset
            url2, fname2 = find_best_bkg_file(gw_n, t_ref, ptypes)
            if url2:
                print(f"  [BKG W{gw_n}] Found: {fname2}")
                try:
                    text = fetch_gzip(url2, timeout=20)
                    quick_ts = _quick_parse_epochs(text)
                    if quick_ts and min(quick_ts) <= t_ref <= max(quick_ts):
                        print(f"  [BKG W{gw_n}] ✓ Has data for {t_ref.date()} "
                              f"({quick_ts[0].date()} → {quick_ts[-1].date()})")
                        bkg_url, bkg_fname = url2, fname2
                        break
                except Exception as e:
                    print(f"  [BKG W{gw_n}] ✗ {e}")

    if not bkg_url:
        print(f"  [错误] 所有 SP3 数据源不可用")
        return None

    try:
        text = fetch_gzip(bkg_url, timeout=30)
        if '## ' not in text or len(text) < 5000:
            raise ValueError("Invalid SP3")
        ptype = 'FIN' if 'FIN' in bkg_fname else 'RAP'
        print(f"  [BKG] ✓ 下载 ({len(text)//1024}KB, {ptype})")
    except Exception as e:
        print(f"  [BKG] ✗ {e}")
        return None

    epochs, ts = parse_sp3_text(text)
    if not ts:
        print(f"  [错误] SP3 解析失败")
        return None

    first_ep = list(epochs.values())[0]
    n_gps = sum(1 for sv in first_ep if str(sv).startswith('G'))

    result = {
        'epochs': epochs, 'ts': ts,
        'source': bkg_fname[:40],
        'url': bkg_url,
        'product': 'FIN' if 'FIN' in bkg_fname else 'RAP',
    }

    try:
        pickle.dump(result, open(cache, 'wb'))
    except Exception:
        pass

    print(f"  SP3: {len(ts)} epochs, {n_gps} GPS sats, {ts[0]} → {ts[-1]}")
    return result


def get_gps_pos_from_sp3(sp3_data, sv, t, n_order=10):
    """Barycentric Lagrange interpolation of GPS sat position/clock/velocity from SP3.

    Uses N-point sliding window (default 10 = 9th-order polynomial) for mm-level
    accuracy between 15-min SP3 epochs. Returns analytic derivative for velocity.

    Args:
        sp3_data: dict with 'ts' (list of datetime) and 'epochs' (dict sv→[x,y,z,clk])
        sv: satellite ID (e.g. 'G05')
        t: query epoch (datetime)
        n_order: number of SP3 epochs in interpolation window (default 10)

    Returns:
        pos: np.array([x, y, z]) in metres, or None
        clk: satellite clock correction in metres
        vel: np.array([vx, vy, vz]) in m/s
    """
    if sp3_data is None:
        return None, 0.0, np.zeros(3)

    ts = sp3_data['ts']
    epochs = sp3_data['epochs']
    if not ts:
        return None, 0.0, np.zeros(3)

    t_ref = ts[0]
    t_float = (t - t_ref).total_seconds()
    ts_float = np.array([(ti - t_ref).total_seconds() for ti in ts])

    valid_indices = []
    for i, ti in enumerate(ts):
        p = epochs[ti].get(sv)
        if p is not None and abs(p[3]) < 0.1 * C:
            valid_indices.append(i)

    if len(valid_indices) < 2:
        valid_indices = [i for i, ti in enumerate(ts) if sv in epochs[ti]]
    if len(valid_indices) < 2:
        for ti in ts:
            p = epochs[ti].get(sv)
            if p is not None:
                return np.array(p[:3]), float(p[3]), np.zeros(3)
        return None, 0.0, np.zeros(3)

    valid_t = ts_float[valid_indices]
    n_use = min(n_order, len(valid_indices))

    idx_right = int(np.searchsorted(valid_t, t_float))
    half = n_use // 2
    left = max(0, idx_right - half)
    right = min(len(valid_indices), left + n_use)
    left = max(0, right - n_use)

    window_idx = [valid_indices[i] for i in range(left, right)]
    window_t = np.array([ts_float[i] for i in window_idx])
    window_vals = np.array([epochs[ts[i]][sv] for i in window_idx])
    n = len(window_t)

    # ── barycentric weights ──
    w = np.ones(n)
    for i in range(n):
        for j in range(n):
            if i != j:
                w[i] /= (window_t[i] - window_t[j])

    diff = t_float - window_t

    # exact-node check
    on_node = None
    for k in range(n):
        if abs(diff[k]) < 1e-9:
            on_node = k
            break

    if on_node is not None:
        k = on_node
        p_raw = window_vals[k]
        pos = np.array(p_raw[:3], dtype=float)
        clk = float(p_raw[3])
        # L'(t_k) = sum_{j != k} (w_j/w_k) * (y_j - y_k) / (t_k - t_j)
        deriv = np.zeros(4)
        wk = w[k]
        for j in range(n):
            if j == k:
                continue
            coeff = (w[j] / wk) / (window_t[k] - window_t[j])
            deriv += coeff * (window_vals[j] - window_vals[k])
        return pos, clk, deriv[:3]

    # ── evaluate barycentric formula ──
    inv_diff = 1.0 / diff
    q_terms = w * inv_diff
    s = np.sum(q_terms)

    p_num = np.sum(w[:, np.newaxis] * window_vals * inv_diff[:, np.newaxis], axis=0)
    interp = p_num / s

    # analytic derivative: L' = (p'·s − p·s') / s²
    inv_diff_sq = -1.0 / (diff * diff)
    s_prime = np.sum(w * inv_diff_sq)
    p_prime = np.sum(w[:, np.newaxis] * window_vals * inv_diff_sq[:, np.newaxis], axis=0)

    deriv = (p_prime * s - p_num * s_prime) / (s * s)

    return interp[:3], float(interp[3]), deriv[:3]


def load_whu_sp3_ultra_rapid(year, doy, data_dir="./data"):
    return load_igs_sp3(year, doy, data_dir)