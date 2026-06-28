#!/usr/bin/env python3
"""
GRACE-FO PPP 批量验证脚本 - v1.1.0
支持广播星历 / 精密星历 / 两者对比

用法示例:
  python3 run_batch.py --ephem both --start 2024-05-01 --end 2024-05-31
  python3 run_batch.py --ephem sp3 --sp3-product final --start 2024-05-01 --end 2024-05-31
  python3 run_batch.py --ephem broadcast --obs-source gps1b --start 2024-04-29 --end 2024-04-29
"""
from __future__ import annotations
import sys, os, argparse, ssl, urllib.request, tarfile, io, csv
import importlib.util
from pathlib import Path
from datetime import datetime as dt, timedelta
from itertools import groupby
from typing import Optional, List, Dict, Any

import numpy as np

WORKDIR = Path(__file__).parent
sys.path.insert(0, str(WORKDIR))

# Import PPP engine
spec = importlib.util.spec_from_file_location('_rp', str(WORKDIR / 'run_ppp.py'))
_rp = importlib.util.module_from_spec(spec)
sys.modules['_rp'] = _rp
spec.loader.exec_module(_rp)

ISDC = "https://isdc-data.gfz.de/grace-fo/Level-1B/JPL/INSTRUMENT/RL04/{year}/"

# ─────────────────────────────────────────────────────────────────────────────
# 数据下载
# ─────────────────────────────────────────────────────────────────────────────
def download_gnv1b(year: int, month: int, day: int, data_dir: str = "./data") -> Optional[str]:
    """下载 GRACE-FO GNV1B 精密轨道文件（缓存）"""
    ds = f"{year:04d}-{month:02d}-{day:02d}"
    out_dir = Path(data_dir) / "gracefo" / str(year) / ds
    out_dir.mkdir(parents=True, exist_ok=True)
    gnv_path = out_dir / f"GNV1B_{ds}_C_04.txt"
    if gnv_path.exists() and gnv_path.stat().st_size > 1000:
        print(f"  [缓存] {ds}")
        return str(gnv_path)
    fname = f"gracefo_1B_{ds}_RL04.ascii.noLRI.tgz"
    url = ISDC.format(year=year) + fname
    print(f"  下载 {fname} ...", end="", flush=True)
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    try:
        req = urllib.request.Request(url, headers={
            "User-Agent": "curl/7.88",
            "Accept-Encoding": "gzip"
        })
        with urllib.request.urlopen(req, timeout=180, context=ctx) as r:
            data = r.read()
        print(f" {len(data) // 1024} KB")
    except Exception as e:
        print(f"  error: {e}")
        return None
    try:
        tar = tarfile.open(fileobj=io.BytesIO(data), mode="r:*")
        for member in tar.getmembers():
            if "GNV1B" in member.name and member.name.endswith(".txt"):
                fo = out_dir / Path(member.name).name
                f = tar.extractfile(member)
                if f:
                    fo.write_bytes(f.read())
        tar.close()
    except Exception as e:
        print(f"  extract error: {e}")
        return None
    return str(gnv_path) if gnv_path.exists() else None


def parse_gnv1b(path: str) -> Dict[dt, np.ndarray]:
    """解析 GNV1B 精密轨道文件，返回 {datetime: [X,Y,Z]} dict"""
    orbit: Dict[dt, np.ndarray] = {}
    gps_origin = dt(2000, 1, 1, 12, 0, 0)
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                tg = float(parts[0])
                flag = parts[2]
                if flag not in ("C", "E"):
                    continue
                X, Y, Z = float(parts[3]), float(parts[4]), float(parts[5])
                if abs(X) < 1e3:
                    continue
                orbit[gps_origin + timedelta(seconds=tg)] = np.array([X, Y, Z])
            except Exception:
                continue
    return orbit


# ─────────────────────────────────────────────────────────────────────────────
# SP3 精密星历
# ─────────────────────────────────────────────────────────────────────────────
def load_sp3_for_date(year: int, month: int, day: int,
                       data_dir: str = "./data", product: str = "auto") -> Any:
    """加载指定日期的 SP3 精密星历"""
    from src.sp3_loader import load_igs_sp3
    doy = int((dt(year, month, day) - dt(year, 1, 1)).total_seconds() // 86400) + 1
    return load_igs_sp3(year, int(doy), data_dir, product=product)


# ─────────────────────────────────────────────────────────────────────────────
# PPP 单日处理
# ─────────────────────────────────────────────────────────────────────────────
def run_ppp_day(year: int, month: int, day: int, ephem: str,
                product: str = "auto",
                nhours: float = 4.0,
                interval: float = 30.0,
                data_dir: str = "./data",
                grace_satellite: str = "gracefo1",
                obs_source: str = "simulated") -> Optional[Dict[str, Any]]:
    """
    对单日数据运行 PPP

    obs_source:
      simulated  - 用 GNV1B 轨道正向模拟 GPS 双频观测（默认）
      gps1b      - 用 GRACE-FO 真实 GPS1B L1B 观测数据

    Returns: dict with keys re, rn, ru, r3, n_ok, n_total or None
    """
    t_start = dt(year, month, day, 0, 0, 0)
    t_end = t_start + timedelta(hours=nhours)
    ds = f"{year:04d}-{month:02d}-{day:02d}"
    print(f"\n[{dt.now().strftime('%H:%M:%S')}] {ds} [{ephem}] obs={obs_source}")

    # 1. 参考轨道
    gnv = download_gnv1b(year, month, day, data_dir)
    if not gnv:
        return None
    ref_orbit = parse_gnv1b(gnv)
    print(f"  轨道: {len(ref_orbit)} epochs  ({min(ref_orbit)} -> {max(ref_orbit)})")

    # 2. 精密星历（可选）
    sp3_data = None
    if ephem == "sp3":
        sp3_data = load_sp3_for_date(year, month, day, data_dir, product=product)
        if sp3_data is None:
            print(f"  [警告] SP3 不可用，降级为广播星历")
            ephem = "broadcast"

    # 3. 生成或加载观测
    if obs_source == "gps1b":
        from src.gps1b_loader import download_gps1b, build_ppp_records
        grace_char = "C" if grace_satellite in ("gracefo1", "gracefo") else "D"
        gps1b_obs, _ = download_gps1b(year, month, day, data_dir, grace_filter=grace_char)
        if gps1b_obs is None:
            print(f"  [错误] GPS1B 数据不可用")
            return None
        records = build_ppp_records(
            gps1b_obs, ref_orbit, t_start, nhours, interval,
            ephem=ephem, sp3_data=sp3_data,
            satpos_func=None,
            gps_sv_plan=_rp.GPS_SV_PLAN)
        if not records:
            print(f"  [错误] GPS1B 观测为空")
            return None
        print(f"  GPS1B 观测: {len(records)} 条")
    else:
        records = _rp.generate_obs_from_orbit(
            ref_orbit, t_start, nhours, interval,
            ephem=ephem, sp3_data=sp3_data)

    if not records:
        print(f"  [错误] 无有效观测")
        return None

    # 4. 按时间排序、分组
    records.sort(key=lambda r: r["time"].timestamp())
    groups = [(t, list(g)) for t, g in groupby(records, key=lambda r: r["time"].timestamp())]

    # 5. PPP 求解
    ts = sorted(ref_orbit.keys())
    t0_orb = ts[0]
    x0 = np.concatenate([ref_orbit[t0_orb], [0.0, 0.2]])  # [x, y, z, clock, sigma]

    results: List[Dict[str, Any]] = []
    for i, (t, grp) in enumerate(groups):
        dt_obj = dt.fromtimestamp(t)
        if not (t_start <= dt_obj <= t_end):
            continue
        obs_list = [
            (r["sv"], r["sat_pos"], r["sat_vel"],
             r["L1"], r["L2"], r["P1"], r["P2"],
             np.radians(r["el"]), np.radians(r["az"]))
            for r in grp
        ]
        if len(obs_list) < 4:
            continue
        x = _rp.ppp_single_epoch(obs_list, x0)
        pos_est = x[:3]

        # 参考轨道插值
        t0_orb = t1_orb = None
        for j, ti in enumerate(ts):
            if ti >= dt_obj:
                t1_orb = ti
                t0_orb = ts[j - 1] if j > 0 else None
                break
            t0_orb = ti
        if t1_orb is None:
            t0_orb = t1_orb = ts[-1]
        if t0_orb is None:
            t0_orb = ts[0]
        ref_pos = ref_orbit[t0_orb]
        if t0_orb != t1_orb:
            dt0 = (dt_obj - t0_orb).total_seconds()
            dt_tot = (t1_orb - t0_orb).total_seconds()
            if dt_tot > 0:
                alpha = dt0 / dt_tot
                ref_pos = ref_orbit[t0_orb] * (1 - alpha) + ref_orbit[t1_orb] * alpha

        err = pos_est - ref_pos
        lat, lon, _ = _rp.ecef_to_blh(ref_pos)
        R = _rp.ecef_to_enu_matrix(lat, lon)
        enu = R @ err
        results.append({
            "time": dt_obj,
            "dE": float(enu[0]),
            "dN": float(enu[1]),
            "dU": float(enu[2]),
            "n_sat": len(obs_list),
        })
        x0 = x.copy()

    print(f"  PPP: {len(results)}/{len(groups)} epochs 收敛")

    if not results:
        return None

    # 6. 统计（过滤异常值 < 1m）
    dE = np.array([float(r["dE"]) for r in results])
    dN = np.array([float(r["dN"]) for r in results])
    dU = np.array([float(r["dU"]) for r in results])
    d3 = np.sqrt(dE**2 + dN**2 + dU**2)
    ok = d3 < 1.0
    if ok.sum() < 10:
        return None

    re = float(np.sqrt(np.nanmean(dE[ok] ** 2)) * 100)
    rn = float(np.sqrt(np.nanmean(dN[ok] ** 2)) * 100)
    ru = float(np.sqrt(np.nanmean(dU[ok] ** 2)) * 100)
    r3 = float(np.sqrt(np.nanmean(d3[ok] ** 2)) * 100)
    print(f"  E={re:.2f}cm N={rn:.2f}cm U={ru:.2f}cm 3D={r3:.2f}cm  ({ok.sum()}/{len(results)} epochs)")

    # 7. 保存 CSV
    out_dir = Path(data_dir) / "output"
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"ppp_{ds}_{obs_source}.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["time", "dE_m", "dN_m", "dU_m", "d3_m", "n_sat"])
        writer.writeheader()
        for r in results:
            d3r = np.sqrt(float(r["dE"])**2 + float(r["dN"])**2 + float(r["dU"])**2)
            writer.writerow({
                "time": r["time"].isoformat(),
                "dE_m": f"{r['dE']:.4f}",
                "dN_m": f"{r['dN']:.4f}",
                "dU_m": f"{r['dU']:.4f}",
                "d3_m": f"{d3r:.4f}",
                "n_sat": r["n_sat"],
            })

    return {
        "date": ds,
        "csv": str(csv_path),
        "re": re, "rn": rn, "ru": ru, "r3": r3,
        "n_ok": int(ok.sum()),
        "n_total": len(results),
    }


# ─────────────────────────────────────────────────────────────────────────────
# 汇总统计
# ─────────────────────────────────────────────────────────────────────────────
def _summarize(label: str, results: List[Optional[Dict[str, Any]]]) -> None:
    valid = [r for r in results if r is not None]
    if not valid:
        print(f"\n{label}: 无数据")
        return
    re = np.array([r["re"] for r in valid])
    rn = np.array([r["rn"] for r in valid])
    ru = np.array([r["ru"] for r in valid])
    r3 = np.array([r["r3"] for r in valid])
    print(f"\n{'='*50}\n{label} 汇总 ({len(valid)}/{len(results)} 天)")
    for lbl, arr in [("E", re), ("N", rn), ("U", ru), ("3D", r3)]:
        print(f"  {lbl} = {np.mean(arr):.2f} +/- {np.std(arr):.2f} cm")
    print("=" * 50)


# ─────────────────────────────────────────────────────────────────────────────
# 日期工具
# ─────────────────────────────────────────────────────────────────────────────
def parse_date(s: str) -> dt:
    return dt.strptime(s, "%Y-%m-%d")


def date_range(start: dt, end: dt):
    cur = start
    while cur <= end:
        yield cur
        cur += timedelta(days=1)


# ─────────────────────────────────────────────────────────────────────────────
# 主程序
# ─────────────────────────────────────────────────────────────────────────────
def main() -> None:
    parser = argparse.ArgumentParser(
        description="GRACE-FO PPP 批量验证 - 支持广播/精密星历对比")
    parser.add_argument("--ephem", required=True,
        choices=["broadcast", "sp3", "both"],
        help="broadcast=广播星历, sp3=精密星历, both=两者对比")
    parser.add_argument("--sp3-product", default="auto",
        choices=["final", "rapid", "ultra", "auto"],
        help="final=事后精密(推荐), rapid=快速精密, ultra=Ultra-Rapid, auto=自动")
    parser.add_argument("--obs-source", default="simulated",
        choices=["simulated", "gps1b"],
        help="simulated=模拟观测(默认), gps1b=真实GRACE-FO GPS1B观测数据")
    parser.add_argument("--start", required=True,
        help="开始日期 (YYYY-MM-DD)")
    parser.add_argument("--end", required=True,
        help="结束日期 (YYYY-MM-DD)")
    parser.add_argument("--hours", type=float, default=4.0,
        help="每日处理时长（小时，默认 4）")
    parser.add_argument("--interval", type=float, default=30.0,
        help="观测采样间隔（秒，默认 30）")
    parser.add_argument("--data-dir", default="./data",
        help="数据目录（默认 ./data）")
    parser.add_argument("--output-dir", default="./output",
        help="输出目录（默认 ./output）")
    args = parser.parse_args()

    start_dt = parse_date(args.start)
    end_dt = parse_date(args.end)
    dates = list(date_range(start_dt, end_dt))
    n_days = len(dates)

    print(f"""
{'='*60}
GRACE-FO PPP 批量验证 v1.1.0
日期范围: {args.start} -> {args.end}  ({n_days} 天)
星历模式: {args.ephem}
观测数据: {args.obs_source}
SP3 产品: {args.sp3_product}
每日处理: {args.hours}h @ {args.interval}s 间隔
{'='*60}
""")

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    def run_for_ephem(ephem_type: str) -> List[Optional[Dict[str, Any]]]:
        results: List[Optional[Dict[str, Any]]] = []
        for i, cur_date in enumerate(dates):
            y, m, d = cur_date.year, cur_date.month, cur_date.day
            r = run_ppp_day(
                y, m, d, ephem_type,
                product=args.sp3_product,
                nhours=args.hours,
                interval=args.interval,
                data_dir=args.data_dir,
                obs_source=args.obs_source)
            results.append(r)
            print(f"  进度: {i+1}/{n_days}  ({(i+1)*100//n_days}%)")
        return results

    if args.ephem == "both":
        print(">>> 运行广播星历模式")
        b_results = run_for_ephem("broadcast")
        _summarize("广播星历", b_results)
        print(">>> 运行精密星历模式")
        s_results = run_for_ephem("sp3")
        _summarize("SP3 精密", s_results)
        print(f"\n对比完成! 输出目录: {args.output_dir}")
    elif args.ephem == "broadcast":
        results = run_for_ephem("broadcast")
        _summarize("广播星历", results)
    else:
        results = run_for_ephem("sp3")
        _summarize("SP3 精密", results)

    print(f"\n全部完成! 日期范围: {args.start} -> {args.end}")


if __name__ == "__main__":
    main()
