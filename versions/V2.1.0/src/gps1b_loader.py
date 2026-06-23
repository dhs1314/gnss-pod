"""
GRACE-FO GPS Level-1B 真实观测数据解析器 — v1.2.0
支持 ISDC 下载和本地 tgz 归档两种方式

数据来源: ISDC @ GFZ
URL: https://isdc-data.gfz.de/grace-fo/Level-1B/JPL/INSTRUMENT/RL04/{year}/
文件: gracefo_1B_{date}_RL04.ascii.noLRI.tgz → GPS1B_{date}_C_04.txt

格式: ASCII NetCDF (含 YAML 头 + 固定格式数据行)
采样率: 0.1 Hz (10秒间隔)

关键说明:
  - prod_flag: 16位二进制串，LSB (最右位) = bit 0
  - 数据字段顺序与 YAML 变量定义顺序一致，与 flag_meanings 的 bit 编号不同
  - prod_flag=4095 (bits 0-11) 表示前12个字段存在（12个数据值）
"""
from pathlib import Path
from datetime import datetime, timedelta
import numpy as np
import ssl
import urllib.request
import tarfile
import io
import re
import pickle

GPS_ORIGIN = datetime(2000, 1, 1, 12, 0, 0)
C = 299792458.0
F1 = 1575.42e6
F2 = 1227.60e6
F1_SQ = F1 * F1
F2_SQ = F2 * F2

ISDC_TGZ_URL = (
    "https://isdc-data.gfz.de/grace-fo/"
    "Level-1B/JPL/INSTRUMENT/RL04/{year}/"
    "gracefo_1B_{date}_RL04.ascii.noLRI.tgz"
)

# GPS1B 文件中固定 header 字段数量（prod_flag 之前的列）
N_HEADER_FIELDS = 7  # rcvtime_intg, rcvtime_frac, GRACEFO_id, prn_id, ant_id, prod_flag, qualflg

# 按数据出现顺序排列的产品字段名（与 YAML variables 列表顺序一致）
# 每个字段对应 prod_flag 的一个 bit（bit 0 = 第一个字段，bit 1 = 第二个字段，...）
YAML_PROD_FIELDS_ORDER = [
    'CA_range',   # bit 0
    'L1_range',   # bit 1
    'L2_range',   # bit 2
    'CA_phase',   # bit 3
    'L1_phase',   # bit 4
    'L2_phase',   # bit 5
    'CA_SNR',     # bit 6
    'L1_SNR',     # bit 7
    'L2_SNR',     # bit 8
    'CA_chan',    # bit 9
    'L1_chan',    # bit 10
    'L2_chan',    # bit 11
    'L2_raw',     # bit 12
    'Ka_phase',   # bit 13
    'K_SNR',      # bit 14
    'Ka_SNR',     # bit 15
]

# ─────────────────────────────────────────────────────────────────
# GPS 时间转换
# ─────────────────────────────────────────────────────────────────
_LEAP_TABLE = [
    (datetime(2000, 1, 1),  32),
    (datetime(2006, 1, 1),  33),
    (datetime(2009, 1, 1),  34),
    (datetime(2012, 7, 1),  35),
    (datetime(2015, 7, 1),  36),
    (datetime(2017, 1, 1),  37),
    (datetime(2024, 1, 1),  18),
]


def gps_sod_to_utc(gps_sod):
    """GPS 秒（自 2000-01-01 12:00 GPS）→ UTC datetime"""
    gps_dt = GPS_ORIGIN + timedelta(seconds=gps_sod)
    for d, leap in reversed(_LEAP_TABLE):
        if gps_dt >= d:
            return gps_dt - timedelta(seconds=leap)
    return gps_dt


def parse_prod_flag(s):
    """解析 prod_flag 二进制字符串（如 '0000111111111111'）→ int"""
    return int(s.strip(), 2)


def parse_qual_flag(s):
    """解析 qualflg 二进制字符串 → int"""
    return int(s.strip(), 2)


# ─────────────────────────────────────────────────────────────────
# GPS1B 单行解析
# ─────────────────────────────────────────────────────────────────
def parse_gps1b_record(parts, grace_filter='C'):
    """解析单条 GPS1B 数据行

    Args:
        parts: 空格分割的数据字段（不含 # 注释）
        grace_filter: 卫星 ID 过滤（'C' 或 'D'）
    Returns:
        rec dict 或 None
    """
    if len(parts) < N_HEADER_FIELDS:
        return None

    try:
        gps_sod = int(parts[0]) + int(parts[1]) * 1e-6
        grace_id = parts[2].strip()
        prn = int(parts[3])
        prod_flag_int = parse_prod_flag(parts[5])
        qual_flag_int = parse_qual_flag(parts[6])
    except (ValueError, IndexError):
        return None

    if prn < 1 or prn > 32:
        return None
    if grace_filter and grace_id != grace_filter:
        return None

    # 提取有效产品字段（按数据顺序）
    prod_values = parts[N_HEADER_FIELDS:]
    rec = {
        'sv': f"G{prn:02d}",
        'grace_id': grace_id,
        'gps_sod': gps_sod,
        'prod_flag_int': prod_flag_int,
        'qual_flag_int': qual_flag_int,
    }

    # 按 YAML 定义顺序 + prod_flag 动态提取
    for bit_pos, field_name in enumerate(YAML_PROD_FIELDS_ORDER):
        if (prod_flag_int >> bit_pos) & 1:
            if bit_pos < len(prod_values):
                try:
                    rec[field_name] = float(prod_values[bit_pos])
                except ValueError:
                    rec[field_name] = None
            else:
                rec[field_name] = None
        # else: 字段不存在，不设置

    # ── 派生量 ───────────────────────────────────────────────────
    # L1/L2 载波相位（来自 iono-smoothed range 或原始相位）
    L1_raw = rec.get('L1_range') or rec.get('L1_phase')
    L2_raw = rec.get('L2_range') or rec.get('L2_phase')
    P1_raw = rec.get('CA_range') or rec.get('L1_range')
    P2_raw = rec.get('L2_range')

    if L1_raw is None or L2_raw is None:
        return None

    rec['L1'] = L1_raw
    rec['L2'] = L2_raw
    rec['P1'] = P1_raw
    rec['P2'] = P2_raw

    # SNR V/V → dB-Hz
    snr1 = rec.get('L1_SNR')
    rec['SNR1_dB'] = (20 * np.log10(max(snr1, 1e-6))
                      if snr1 is not None and snr1 > 0 else None)

    # 质量过滤：L1 SNR < 5 V/V（约 14 dB-Hz）
    if snr1 is not None and snr1 < 5:
        return None

    # 电离层自由组合（IF）
    rec['L_if'] = (F1_SQ * L1_raw - F2_SQ * L2_raw) / (F1_SQ - F2_SQ)
    if P1_raw is not None and P2_raw is not None:
        rec['P_if'] = (F1_SQ * P1_raw - F2_SQ * P2_raw) / (F1_SQ - F2_SQ)
    else:
        rec['P_if'] = P1_raw

    return rec


# ─────────────────────────────────────────────────────────────────
# 解析 GPS1B 文本文件
# ─────────────────────────────────────────────────────────────────
def parse_gps1b_text(text, grace_filter='C'):
    """解析 GPS1B ASCII 文件

    Returns:
        {gps_sod: {sv: rec}}  按时钟面（SOW）索引
    """
    lines = text.split('\n')

    # 找 YAML 头结束位置
    header_end = 0
    for i, l in enumerate(lines):
        if l.strip() == '# End of YAML header':
            header_end = i + 1
            break

    print(f"  GPS1B YAML: {header_end} header lines, "
          f"{len(YAML_PROD_FIELDS_ORDER)} product fields, filter={grace_filter}")

    gps_obs: dict = {}
    n_total = n_valid = 0
    sat_counts: dict = {}

    for line in lines[header_end:]:
        if not line.strip() or line.startswith('#'):
            continue
        n_total += 1
        parts = line.split()
        rec = parse_gps1b_record(parts, grace_filter)
        if rec is None:
            continue
        gps_sod = rec['gps_sod']
        sv = rec['sv']
        if gps_sod not in gps_obs:
            gps_obs[gps_sod] = {}
        gps_obs[gps_sod][sv] = rec
        sat_counts[sv] = sat_counts.get(sv, 0) + 1
        n_valid += 1

    print(f"  GPS1B: {n_total} rows -> {n_valid} valid, "
          f"{len(gps_obs)} epochs, {len(sat_counts)} sats")
    for sv in sorted(sat_counts.keys()):
        print(f"    {sv}: {sat_counts[sv]} epochs")

    return gps_obs


# ─────────────────────────────────────────────────────────────────
# 下载
# ─────────────────────────────────────────────────────────────────
def download_gps1b(year, month, day, data_dir="./data", grace_filter='C'):
    """下载并解析 GRACE-FO GPS1B 数据（自动缓存 .pkl）

    Returns:
        (gps_obs, utc_t_start) 或 (None, None)
    """
    date_str = f"{year:04d}-{month:02d}-{day:02d}"
    cache = Path(data_dir) / "gracefo" / str(year) / date_str / f"gps1b_{grace_filter}.pkl"
    cache.parent.mkdir(parents=True, exist_ok=True)

    # 读缓存
    if cache.exists():
        try:
            gps_obs = pickle.load(open(cache, 'rb'))
            n_obs = sum(len(v) for v in gps_obs.values())
            print(f"  [GPS1B cache] {date_str} ({grace_filter}): "
                  f"{len(gps_obs)} epochs, {n_obs} obs")
            return gps_obs, datetime(year, month, day)
        except Exception:
            pass

    # 找本地 tgz（可能在 data/ 根目录）
    tgz_path = Path(data_dir) / f"gracefo_1B_{date_str}_RL04.ascii.noLRI.tgz"
    if not tgz_path.exists():
        tgz_path = Path(data_dir) / "gracefo" / str(year) / f"gracefo_1B_{date_str}_RL04.ascii.noLRI.tgz"

    if tgz_path.exists():
        print(f"  [tgz 缓存 ✓] {tgz_path.name}")
    else:
        tgz_url = ISDC_TGZ_URL.format(year=year, date=date_str)
        print(f"  [GPS1B 下载] {tgz_url.split('/')[-1]}")
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        try:
            req = urllib.request.Request(tgz_url, headers={
                'User-Agent': 'curl/7.88', 'Accept-Encoding': 'gzip'
            })
            with urllib.request.urlopen(req, timeout=300, context=ctx) as r:
                data = r.read()
            tgz_path.parent.mkdir(parents=True, exist_ok=True)
            tgz_path.write_bytes(data)
            print(f"  [完成] {len(data) // 1024 // 1024} MB")
        except Exception as e:
            print(f"  [下载失败] {e}")
            return None, None

    # 从 tgz 提取 GPS1B 文件
    try:
        tar = tarfile.open(tgz_path)
        names = tar.getnames()
        gps_member = [n for n in names if 'GPS1B' in n and n.endswith('.txt')]
        if not gps_member:
            print(f"  [错误] tgz 中无 GPS1B .txt: {names}")
            return None, None
        gps_member = gps_member[0]
        print(f"  [GPS1B 解包] {gps_member}")
        f = tar.extractfile(gps_member)
        text = f.read().decode('ascii', errors='replace')
        f.close()
        tar.close()
    except Exception as e:
        print(f"  [解包失败] {e}")
        return None, None

    gps_obs = parse_gps1b_text(text, grace_filter=grace_filter)
    if gps_obs:
        pickle.dump(gps_obs, open(cache, 'wb'))
        print(f"  [缓存已保存] {cache}")
    return gps_obs, datetime(year, month, day)


# ─────────────────────────────────────────────────────────────────
# 构建 PPP 观测记录
# ─────────────────────────────────────────────────────────────────
def build_ppp_records(gps_obs, ref_orbit, t_start, nhours, interval,
                     ephem, sp3_data=None, satpos_func=None,
                     gps_sv_plan=None):
    """用真实 GPS1B 观测数据驱动 PPP

    Args:
        gps_obs:      {gps_sod: {sv: rec}} 来自 parse_gps1b_text
        ref_orbit:   {utc_datetime: [X,Y,Z]} GNV1B 参考轨道
        t_start:     处理开始 UTC datetime
        nhours:      处理时长（小时）
        interval:    采样间隔（秒）
        ephem:       'broadcast' 或 'sp3'
        sp3_data:    SP3 精密星历
        satpos_func: (sv, utc_dt) → (pos, vel)  自定义卫星位置函数
                     若为 None，使用 gps_sv_plan（广播星历）
        gps_sv_plan: GPS_SV_PLAN 常量（satpos_from_sv 使用）
    Returns:
        [{time, sv, sat_pos, sat_vel, L_if, P_if, el, az}, ...]
    """
    from run_ppp import satpos_from_sv, ecef_to_blh, ecef_to_enu_matrix
    from src.sp3_loader import get_gps_pos_from_sp3

    t_end = t_start + timedelta(hours=nhours)
    records: list = []
    orbit_ts = sorted(ref_orbit.keys())
    n_skip_rng = n_skip_el = n_no_pos = 0

    for gps_sod, sv_obs in sorted(gps_obs.items()):
        utc_dt = gps_sod_to_utc(gps_sod)
        if not (t_start <= utc_dt <= t_end):
            continue

        # 降采样：接受每个 interval 秒段的第一个可用历元
        # GPS1B 原始数据 10s 间隔；使用 round-to-nearest 匹配目标时刻
        dt_s = (utc_dt - t_start).total_seconds()
        if interval > 5:
            nearest = int((dt_s + interval / 2) / interval) * interval
            if abs(dt_s - nearest) > 2.0:
                continue

        # GRACE-FO 位置（GNV1B 参考轨道线性插值）
        t0_orb = t1_orb = None
        for j, ti in enumerate(orbit_ts):
            if ti >= utc_dt:
                t1_orb = ti
                t0_orb = orbit_ts[j - 1] if j > 0 else None
                break
            t0_orb = ti
        if t1_orb is None:
            t0_orb = t1_orb = orbit_ts[-1]
        if t0_orb is None:
            t0_orb = orbit_ts[0]

        dt_frac = (utc_dt - t0_orb).total_seconds()
        dt_tot = (t1_orb - t0_orb).total_seconds()
        if dt_tot == 0:
            grace_pos = np.array(ref_orbit[t0_orb], dtype=float)
        else:
            a = dt_frac / dt_tot
            grace_pos = (np.array(ref_orbit[t0_orb], dtype=float) * (1 - a) +
                         np.array(ref_orbit[t1_orb], dtype=float) * a)

        # 合理性检查
        grace_r = float(np.linalg.norm(grace_pos))
        if not (6e6 < grace_r < 8e6):
            continue

        # ECEF → BLH → ENU 旋转矩阵（用于天顶角/方位角）
        try:
            lat, lon, _ = ecef_to_blh(grace_pos)
            R_enu = ecef_to_enu_matrix(lat, lon)
        except Exception:
            continue

        for sv, rec in sv_obs.items():
            # GPS 卫星位置
            sat_pos = sat_vel = None
            if satpos_func is not None:
                result = satpos_func(sv, utc_dt)
                if result is None:
                    n_no_pos += 1; continue
                sat_pos, sat_vel = result[0], result[1]
            elif ephem == 'sp3' and sp3_data is not None:
                sat_pos, _clk, sat_vel = get_gps_pos_from_sp3(sp3_data, sv, utc_dt)
                if sat_pos is None:
                    n_no_pos += 1; continue
            else:
                # 广播星历：搜索匹配的 SV 记录
                if gps_sv_plan is not None:
                    for sv_rec in gps_sv_plan:
                        if f"G{sv_rec[0]:02d}" == sv:
                            sat_pos, sat_vel = satpos_from_sv(sv_rec, utc_dt)
                            break
                if sat_pos is None:
                    n_no_pos += 1; continue

            # 卫地距
            delta = grace_pos - sat_pos
            rng = float(np.linalg.norm(delta))
            if not (2e7 < rng < 5e7):
                n_skip_rng += 1; continue

            # 高度角 / 方位角
            try:
                e_enu = R_enu @ (delta / rng)
                el = float(np.arcsin(np.clip(e_enu[2], -1.0, 1.0)))
                az = float(np.arctan2(e_enu[0], e_enu[1]))
                if az < 0:
                    az += 2 * np.pi
            except Exception:
                continue

            if el < 0.087:   # < 5°
                n_skip_el += 1; continue

            records.append({
                'time': utc_dt,
                'sv': sv,
                'sat_pos': sat_pos,
                'sat_vel': sat_vel,
                'L_if': float(rec['L_if']),
                'P_if': float(rec['P_if']),
                'L1': float(rec['L1']),   # 原始 L1 载波相位 (m)
                'L2': float(rec['L2']),   # 原始 L2 载波相位 (m)
                'P1': float(rec['P1']) if rec.get('P1') is not None else float(rec['L1']),
                'P2': float(rec['P2']) if rec.get('P2') is not None else float(rec['L2']),
                'el': float(np.degrees(el)),
                'az': float(np.degrees(az)),
                'rng': rng,
                'SNR1': rec.get('SNR1_dB'),
            })

    print(f"  PPP records: {len(records)} "
          f"(rng={n_skip_rng}, el={n_skip_el}, no_pos={n_no_pos})")
    return records
