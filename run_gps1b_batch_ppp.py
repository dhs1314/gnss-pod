#!/usr/bin/env python3
"""
GRACE-FO Batch Least Squares PPP with Basic Dynamics (Phase 2)

Estimates initial orbit state [r0, v0, Cd] plus ambiguity parameters
over a short arc using batch least squares with orbit integration.

This is the key methodological improvement from kinematic PPP (~1.2m)
toward reduced-dynamic POD (~0.1-0.2m).

Usage:
  py -3.12 run_gps1b_batch_ppp.py --date 2024-04-29 --hours 0.5 --interval 30
"""
import sys, os, pickle, csv, math, argparse
from pathlib import Path
from datetime import datetime, timedelta
from collections import defaultdict
import numpy as np

sys.path.insert(0, str(Path(__file__).parent / 'src'))
from gps1b_rnx_loader import load_gps1b_rnx
from sp3_loader import get_gps_pos_from_sp3 as _sp3_get
from orbit_integrator import integrate_orbit, integrate_orbit_with_stm
from orbit_dynamics import total_acc
from batch_lsq import solve_batch_lsq

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
F1_SQ, F2_SQ = F1 * F1, F2 * F2
ALPHA = F1_SQ / (F1_SQ - F2_SQ)
BETA = -F2_SQ / (F1_SQ - F2_SQ)
OMEGA_E = 7.2921151467e-5
J2000 = datetime(2000, 1, 1, 12, 0, 0)


# -- Coordinate transforms --
def ecef_to_blh(pos):
    X, Y, Z = float(pos[0]), float(pos[1]), float(pos[2])
    p = np.sqrt(X**2 + Y**2)
    if p < 1e-6:
        return np.array([np.pi/2 if Z >= 0 else -np.pi/2, 0.0, 0.0])
    lat = np.arctan2(Z, p)
    for _ in range(10):
        sinL = np.sin(lat)
        N = 6378137.0 / np.sqrt(1 - 0.00669437999014 * sinL**2)
        lat_new = np.arctan2(Z + 0.00669437999014 * N * sinL, p)
        if abs(lat_new - lat) < 1e-15: lat = lat_new; break
        lat = lat_new
    lon = np.arctan2(Y, X)
    return np.array([lat, lon, 0.0])


def ecef_to_enu_matrix(lat, lon):
    sl, cl = np.sin(lat), np.cos(lat)
    sn, cn = np.sin(lon), np.cos(lon)
    return np.array([[-sn, cn, 0],
                     [-sl*cn, -sl*sn, cl],
                     [cl*cn, cl*sn, sl]])


# -- SP3 geometry --
def get_sat_geometry(sp3, sv, utc_dt, rcv_pos):
    pos, clk, vel = _sp3_get(sp3, sv, utc_dt)
    if pos is None or abs(clk) > 0.1 * C:
        return None, None, None
    rho = float(np.linalg.norm(pos - rcv_pos))
    for _ in range(5):
        travel_time = rho / C
        tx_dt = utc_dt - timedelta(seconds=travel_time)
        pos_tx, clk_tx, vel_tx = _sp3_get(sp3, sv, tx_dt)
        if pos_tx is None: break
        rho_new = float(np.linalg.norm(pos_tx - rcv_pos))
        if abs(rho_new - rho) < 1e-8:
            pos, clk = pos_tx, clk_tx; rho = rho_new; break
        pos, clk = pos_tx, clk_tx; rho = rho_new
    if not (1.8e7 < rho < 2.8e7):
        return None, None, None
    sag = (OMEGA_E / C) * (pos[0] * rcv_pos[1] - pos[1] * rcv_pos[0])
    return pos, clk, rho + sag


# -- GNV1B reference orbit --
def load_gnv1b(filepath):
    orbit = {}
    with open(filepath, encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            parts = line.split()
            if len(parts) < 6: continue
            try:
                t = float(parts[0]); flag = parts[2]
                if flag not in ('C', 'E'): continue
                X, Y, Z = float(parts[3]), float(parts[4]), float(parts[5])
                if abs(X) < 1e3: continue
                orbit[t] = np.array([X, Y, Z])
            except (ValueError, IndexError): continue
    return orbit


def interpolate_ref(ref_orbit, gps_sod):
    ts = sorted(ref_orbit.keys())
    t0 = t1 = None
    for i, ti in enumerate(ts):
        if ti >= gps_sod: t1 = ti; t0 = ts[i-1] if i > 0 else None; break
        t0 = ti
    if t1 is None: t0 = t1 = ts[-1]
    if t0 is None: t0 = ts[0]
    if t0 == t1: return ref_orbit[t0]
    a = (gps_sod - t0) / (t1 - t0)
    return ref_orbit[t0] * (1-a) + ref_orbit[t1] * a


def compute_initial_velocity(ref_orbit, gps_sod):
    """Compute approximate ECEF velocity from consecutive GNV1B positions."""
    ts = sorted(ref_orbit.keys())
    # Find closest epoch
    i = 0
    for j, t in enumerate(ts):
        if t >= gps_sod:
            i = j; break
        i = j
    # Use positions at i and i+1 or i-1 and i
    if i + 1 < len(ts):
        dt = ts[i+1] - ts[i]
        if dt > 0.01:
            return (ref_orbit[ts[i+1]] - ref_orbit[ts[i]]) / dt
    if i > 0:
        dt = ts[i] - ts[i-1]
        if dt > 0.01:
            return (ref_orbit[ts[i]] - ref_orbit[ts[i-1]]) / dt
    return np.zeros(3)


def process_day_batch(date_str, nhours=0.5, data_dir='./data', grace_id='C', interval=30.0):
    y, m, d = [int(x) for x in date_str.split('-')]
    doy = (datetime(y, m, d) - datetime(y, 1, 1)).days + 1
    dp = Path(data_dir)

    # Load GNV1B reference
    gnv_path = dp / 'gracefo' / str(y) / date_str / f'GNV1B_{date_str}_{grace_id}_04.txt'
    ref_orbit = load_gnv1b(str(gnv_path))
    print(f"[GNV1B] {len(ref_orbit)} epochs")

    # Load GPS1B RINEX
    rnx_path = dp / f'GPS1B_{date_str}_{grace_id}_04.rnx'
    if not rnx_path.exists():
        print(f"[FATAL] RINEX not found: {rnx_path}")
        return None
    gps1b_raw = load_gps1b_rnx(str(rnx_path))
    if gps1b_raw is None:
        return None
    print(f"[RINEX] {len(gps1b_raw)} epochs")

    # Load SP3
    sp3_pkl = dp / str(y) / f'{doy:03d}' / 'igs_sp3_FIN.pkl'
    if not sp3_pkl.exists():
        sp3_pkl = dp / str(y) / f'{doy:03d}' / 'igs_sp3.pkl'
    sp3 = pickle.load(open(str(sp3_pkl), 'rb'))
    print(f"[SP3] {len(sp3['ts'])} epochs")

    # Select epochs
    gps_sod_start = min(gps1b_raw.keys())
    gps_sod_end = gps_sod_start + nhours * 3600
    epochs = []
    for gps_sod in sorted(gps1b_raw.keys()):
        if not (gps_sod_start <= gps_sod <= gps_sod_end): continue
        dt = gps_sod - gps_sod_start
        nearest = round(dt / interval) * interval
        if abs(dt - nearest) > max(2.0, interval * 0.1): continue
        epochs.append(gps_sod)
    print(f"[EPOCH] {len(epochs)} selected (dt={interval}s)")

    # -- Pre-compute geometry --
    N_BIAS = min(60, len(epochs))
    print(f"\n-- Computing geometry ({len(epochs)} epochs) --")
    epoch_data_list = []
    sv_p_residuals = defaultdict(list)

    for gps_sod in epochs:
        utc_dt = J2000 + timedelta(seconds=gps_sod)
        ref_pos = interpolate_ref(ref_orbit, gps_sod)
        ep_data = []

        for sv_id, rec in gps1b_raw[gps_sod].items():
            if 'L_if' not in rec: continue
            sat_pos, sat_clk, rho_corr = get_sat_geometry(sp3, sv_id, utc_dt, ref_pos)
            if sat_pos is None: continue
            el = math.asin(abs(sat_pos[2] - ref_pos[2]) / rho_corr)
            if el < 0.087: continue

            L_if_raw = float(rec['L_if'])
            P_if_raw = float(rec['P_if'])
            L_r = L_if_raw + sat_clk - rho_corr
            P_r = P_if_raw + sat_clk - rho_corr

            ep_data.append({
                'sv': sv_id, 'sat_pos': sat_pos, 'sat_clk': sat_clk,
                'rho_corr': rho_corr, 'el': el,
                'L_r': L_r, 'P_r': P_r,
                'L_if_raw': L_if_raw, 'P_if_raw': P_if_raw,
            })

            if len(epoch_data_list) < N_BIAS:
                sv_p_residuals[sv_id].append(P_r)

        epoch_data_list.append({
            'gps_sod': gps_sod,
            'utc_dt': utc_dt,
            'ep_data': ep_data,
        })

    # Estimate per-SV code biases
    sv_bias = {}
    for sv, p_vals in sv_p_residuals.items():
        if len(p_vals) >= 10:
            sv_bias[sv] = float(np.median(p_vals))
    print(f"  Per-SV code biases: {len(sv_bias)} SVs, "
          f"range [{min(sv_bias.values()):.1f}, {max(sv_bias.values()):.1f}]m")

    # -- Initial state for batch LSQ --
    r0_init = interpolate_ref(ref_orbit, epochs[0])
    v0_init = compute_initial_velocity(ref_orbit, epochs[0])
    print(f"\n-- Initial orbit: r0=[{r0_init[0]/1000:.1f},{r0_init[1]/1000:.1f},{r0_init[2]/1000:.1f}]km "
          f"v0=[{v0_init[0]:.1f},{v0_init[1]:.1f},{v0_init[2]:.1f}]m/s")

    # -- Batch LSQ --
    print(f"\n-- Batch LSQ ({len(epochs)} epochs, {len(epochs)*10:.0f} obs approx) --")
    result = solve_batch_lsq(
        epoch_data_list=epoch_data_list,
        ref_orbit=ref_orbit,
        sp3=sp3,
        sv_bias=sv_bias,
        gps_sods=epochs,
        r0_init=r0_init,
        v0_init=v0_init,
        Cd_init=2.2,
        max_iter=5,
        sigma_phase=0.01,
        sigma_code=0.30,
    )

    if result is None:
        print("[FATAL] Batch LSQ failed")
        return None

    # -- Compute position errors at each epoch --
    print(f"\n-- Computing post-fit errors --")
    t0 = epochs[0]
    integ = integrate_orbit(result['r0'], result['v0'], (0, epochs[-1] - t0),
                            Cd=result['Cd'], dt=10.0)

    results = []
    for gps_sod in epochs:
        t_rel = gps_sod - t0
        i_step = int(np.round(t_rel / 10.0))
        if i_step < 0 or i_step >= len(integ['t']):
            continue

        r_est = integ['r'][i_step]
        ref_pos = interpolate_ref(ref_orbit, gps_sod)
        err = r_est - ref_pos
        d3 = float(np.linalg.norm(err))

        lat, lon, _ = ecef_to_blh(ref_pos)
        R = ecef_to_enu_matrix(lat, lon)
        enu = R @ err

        results.append({
            'time': J2000 + timedelta(seconds=gps_sod),
            'gps_sod': gps_sod,
            'dE': float(enu[0]), 'dN': float(enu[1]), 'dU': float(enu[2]),
            'd3': d3,
            'n_sat': len([d for ed in epoch_data_list if ed['gps_sod'] == gps_sod
                         for d in ed['ep_data']]),
        })

    if len(results) < 10:
        print("[FATAL] Too few valid epochs")
        return None

    # -- Statistics --
    dE = np.array([r['dE'] for r in results])
    dN = np.array([r['dN'] for r in results])
    dU = np.array([r['dU'] for r in results])
    d3 = np.array([r['d3'] for r in results])

    def rms(a): return float(np.sqrt(np.nanmean(a**2)))

    stats = {}
    for label, thresh in [('all', 1e9), ('<5m', 5.0), ('<2m', 2.0), ('<1m', 1.0)]:
        mask = d3 < thresh
        n_v = int(mask.sum())
        if n_v > 2:
            stats[label] = {
                'n': n_v,
                'rms_e': rms(dE[mask])*100, 'rms_n': rms(dN[mask])*100,
                'rms_u': rms(dU[mask])*100, 'rms_3d': rms(d3[mask])*100,
                'mean_3d': float(np.nanmean(d3[mask]))*100,
                'max_3d': float(np.nanmax(d3[mask]))*100 if n_v > 0 else 0,
            }

    print("\n" + "="*65)
    print(f"  GRACE-FO Batch LSQ PPP -- Two-Body+J2+Drag Dynamics")
    print(f"  Date: {date_str}  Hours: {nhours}h  dt: {interval}s")
    print("="*65)
    hdr = f"  {'Filter':<10s} {'N':>6s} {'E(cm)':>8s} {'N(cm)':>8s} {'U(cm)':>8s} {'3D(cm)':>8s} {'Mean':>8s}"
    print(hdr); print("  " + "-"*55)
    for label, _ in [('all', 1e9), ('<5m', 5.0), ('<2m', 2.0), ('<1m', 1.0)]:
        if label in stats:
            s = stats[label]
            print(f"  {label:<10s} {s['n']:>6d} {s['rms_e']:>8.2f} {s['rms_n']:>8.2f} "
                  f"{s['rms_u']:>8.2f} {s['rms_3d']:>8.2f} {s['mean_3d']:>8.2f}")
    print("="*65)

    return {
        'stats': stats, 'results': results, 'date': date_str,
        'grace_id': grace_id,
        'orbit': result,
    }


# -- Main --
def main():
    parser = argparse.ArgumentParser(description='GRACE-FO Batch LSQ PPP with Dynamics')
    parser.add_argument('--date', default='2024-04-29')
    parser.add_argument('--hours', type=float, default=0.5)
    parser.add_argument('--data-dir', default='./data')
    parser.add_argument('--output-dir', default='./output')
    parser.add_argument('--grace-id', default='C')
    parser.add_argument('--interval', type=float, default=30.0)
    args = parser.parse_args()

    print(f"\n{'='*65}")
    print(f"  GRACE-FO Batch LSQ PPP -- Two-Body+J2+Drag Dynamics")
    print(f"  Date: {args.date}  Hours: {args.hours}h  dt: {args.interval}s")
    print(f"{'='*65}")

    output = process_day_batch(date_str=args.date, nhours=args.hours,
                               data_dir=args.data_dir, grace_id=args.grace_id,
                               interval=args.interval)
    if output is None:
        print("[FATAL] Processing failed")
        sys.exit(1)

    # Save CSV
    out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)
    tag = output['date'].replace('-', '') + '_' + output.get('grace_id', 'C')
    csv_path = out / f'ppp_batch_{tag}.csv'
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['time', 'gps_sod', 'dE_m', 'dN_m', 'dU_m', 'd3_m', 'n_sat'])
        for r in output['results']:
            w.writerow([r['time'].isoformat(), f"{r['gps_sod']:.6f}",
                        f"{r['dE']:.4f}", f"{r['dN']:.4f}", f"{r['dU']:.4f}",
                        f"{r['d3']:.4f}", r['n_sat']])
    print(f"[CSV] {csv_path}")
    print(f"[DONE] {args.date}")


if __name__ == '__main__':
    main()
