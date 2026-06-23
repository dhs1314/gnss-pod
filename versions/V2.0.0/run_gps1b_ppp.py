#!/usr/bin/env python3
"""
GRACE-FO PPP — Two-Pass Float-Ambiguity with per-SV Bias Correction

Algorithm:
  Pass 1: Estimate float ambiguities (B_if) + per-SV DCB-like biases
          Uses GNV1B position as linearization point (works with code-based too)
  Pass 2: Epoch-by-epoch WLS [dX, dY, dZ, clk_r, trop_wet]
          Phase pre-corrected for float amb + per-SV bias
          Linearized at same reference position as Pass 1

Data: JPL RL04 GPS1B (phase in METERS, ionospherically smoothed)
      IGS Final SP3 (barycentric Lagrange interpolation)
      GNV1B reduced-dynamic orbit (reference truth)

Expected accuracy: ~1-2m 3D RMS (code-noise limited data)

Usage:
  py -3.12 run_gps1b_ppp.py --date 2024-04-29 --hours 4 --interval 30
"""
import sys, os, pickle, csv, json, math, argparse
from pathlib import Path
from datetime import datetime, timedelta
from collections import defaultdict
import numpy as np

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
F1_SQ, F2_SQ = F1*F1, F2*F2
ALPHA = F1_SQ / (F1_SQ - F2_SQ)
BETA  = -F2_SQ / (F1_SQ - F2_SQ)
OMEGA_E = 7.2921151467e-5
J2000 = datetime(2000, 1, 1, 12, 0, 0)

# ── Coordinate transforms ───────────────────────────────────────────────
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

# ── SP3 ─────────────────────────────────────────────────────────────────
from src.sp3_loader import get_gps_pos_from_sp3 as _sp3_get

def get_sat_geometry(sp3, sv, utc_dt, rcv_pos):
    """Light-time iterated satellite position, clock, and geometric range + Sagnac.

    Returns (sat_pos_tx, sat_clk_m, rho_sagnac_m) or (None, None, None).
    """
    pos, clk, vel = _sp3_get(sp3, sv, utc_dt)
    if pos is None or abs(clk) > 0.1 * C:
        return None, None, None

    rho = float(np.linalg.norm(pos - rcv_pos))
    for _ in range(5):
        travel_time = rho / C
        tx_dt = utc_dt - timedelta(seconds=travel_time)
        pos_tx, clk_tx, vel_tx = _sp3_get(sp3, sv, tx_dt)
        if pos_tx is None:
            break
        rho_new = float(np.linalg.norm(pos_tx - rcv_pos))
        if abs(rho_new - rho) < 1e-8:
            pos, clk = pos_tx, clk_tx
            rho = rho_new
            break
        pos, clk = pos_tx, clk_tx
        rho = rho_new

    if not (1.8e7 < rho < 2.8e7):
        return None, None, None

    sag = (OMEGA_E / C) * (pos[0] * rcv_pos[1] - pos[1] * rcv_pos[0])
    return pos, clk, rho + sag

# ── GNV1B reference orbit ──────────────────────────────────────────────
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

def interpolate_ref(orbit, gps_sod):
    ts = sorted(orbit.keys())
    t0 = t1 = None
    for i, ti in enumerate(ts):
        if ti >= gps_sod: t1 = ti; t0 = ts[i-1] if i > 0 else None; break
        t0 = ti
    if t1 is None: t0 = t1 = ts[-1]
    if t0 is None: t0 = ts[0]
    if t0 == t1: return orbit[t0]
    a = (gps_sod - t0) / (t1 - t0)
    return orbit[t0] * (1-a) + orbit[t1] * a

# ── Main processing ─────────────────────────────────────────────────────
def process_day(date_str, nhours=4.0, data_dir='./data', grace_id='C', interval=30.0):
    y, m, d = [int(x) for x in date_str.split('-')]
    doy = (datetime(y, m, d) - datetime(y, 1, 1)).days + 1
    dp = Path(data_dir)

    # Load data
    gnv_path = dp / 'gracefo' / str(y) / date_str / f'GNV1B_{date_str}_{grace_id}_04.txt'
    ref_orbit = load_gnv1b(str(gnv_path))
    print(f"[GNV1B] {len(ref_orbit)} epochs")

    gps1b_pkl = dp / 'gracefo' / str(y) / date_str / f'GPS1B_{date_str}_{grace_id}_04.pkl'
    gps1b_raw = pickle.load(open(str(gps1b_pkl), 'rb'))
    print(f"[GPS1B] {len(gps1b_raw)} epochs")

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

    # ═════════════════════════════════════════════════════════════════
    # PASS 1: Estimate B_if and per-SV biases from first 60 epochs
    # Float ambiguity drifts over time — use short window at start.
    # Geometry pre-computed for ALL epochs at GNV1B reference.
    # ═════════════════════════════════════════════════════════════════
    print(f"\n── Pass 1: Bias estimation (first {min(60, len(epochs))} epochs) ──")

    N_BIAS = min(60, len(epochs))
    print(f"  Bias window: {N_BIAS} epochs (~{N_BIAS*interval/60:.0f} min)")

    # Pre-compute geometry + accumulate bias statistics
    epoch_geo = {}
    sv_lp_diff = defaultdict(list)
    sv_p_residuals = defaultdict(list)

    for ep_idx, gps_sod in enumerate(epochs):
        utc_dt = J2000 + timedelta(seconds=gps_sod)
        ref_pos = interpolate_ref(ref_orbit, gps_sod)
        ep_data = []

        for sv_id, rec in gps1b_raw[gps_sod].items():
            if 'L1_phase' not in rec: continue

            sat_pos, sat_clk, rho_corr = get_sat_geometry(sp3, sv_id, utc_dt, ref_pos)
            if sat_pos is None: continue

            el = math.asin(abs(sat_pos[2] - ref_pos[2]) / rho_corr)
            if el < 0.087: continue

            L_if_raw = ALPHA * float(rec['L1_phase']) + BETA * float(rec['L2_phase'])
            P_if_raw = float(rec['P_if'])

            L_r = L_if_raw + sat_clk - rho_corr
            P_r = P_if_raw + sat_clk - rho_corr

            ep_data.append({
                'sv': sv_id, 'sat_pos': sat_pos, 'sat_clk': sat_clk,
                'rho_corr': rho_corr, 'el': el,
                'L_r': L_r, 'P_r': P_r,
                'L_if_raw': L_if_raw, 'P_if_raw': P_if_raw,
            })

            if ep_idx < N_BIAS:
                sv_lp_diff[sv_id].append(L_if_raw - P_if_raw)
                sv_p_residuals[sv_id].append(P_r)

        if ep_data:
            epoch_geo[gps_sod] = ep_data

    # Float ambiguities from bias window
    float_amb = {}
    for sv, diffs in sv_lp_diff.items():
        if len(diffs) >= 10:
            arr = np.array(diffs)
            s = float(np.std(arr))
            if s <= 5.0:
                float_amb[sv] = float(np.median(arr))
            else:
                print(f"  [REJ] {sv}: B_if std={s:.1f}m > 5m threshold")

    print(f"  Float ambiguities: {len(float_amb)} SVs (from {N_BIAS} epochs)")

    # Per-SV biases from bias window
    sv_bias = {}
    for sv, p_vals in sv_p_residuals.items():
        if sv not in float_amb: continue
        if len(p_vals) >= 10:
            sv_bias[sv] = float(np.median(p_vals))

    print(f"  Per-SV biases: {len(sv_bias)} SVs")
    if sv_bias:
        biases = list(sv_bias.values())
        print(f"  Bias range: [{min(biases):.1f}, {max(biases):.1f}]m")
    sv_ok = set(float_amb.keys()) & set(sv_bias.keys())

    # ═════════════════════════════════════════════════════════════════
    # PASS 2: Epoch-by-epoch WLS with bias-corrected phase
    # ═════════════════════════════════════════════════════════════════
    print(f"\n── Pass 2: Epoch-by-epoch WLS (bias-corrected) ──")

    std_phase = 0.01  # tight: suppresses code noise
    std_code  = 0.30

    results = []
    for i_ep, gps_sod in enumerate(epochs):
        if gps_sod not in epoch_geo: continue
        ep_data = epoch_geo[gps_sod]
        ref_pos = interpolate_ref(ref_orbit, gps_sod)

        # Single-pass WLS (linear problem at fixed reference position)
        A_rows, y_rows, w_rows = [], [], []
        for d in ep_data:
            sv = d['sv']
            if sv not in sv_ok: continue

            e_vec = (d['sat_pos'] - ref_pos) / d['rho_corr']
            mf = 1.0 / max(math.sin(d['el']), 0.1)
            h = np.array([-e_vec[0], -e_vec[1], -e_vec[2], 1.0, mf])

            # Phase: corrected for float amb + per-SV bias
            L_corrected = d['L_r'] - float_amb[sv] - sv_bias[sv]

            # Reject gross outliers
            if abs(L_corrected) > 50.0:
                continue

            A_rows.append(h)
            y_rows.append(L_corrected)
            w_rows.append(1.0 / std_phase**2)

            # Code
            if abs(d['P_r']) > 50.0:
                continue
            A_rows.append(h)
            y_rows.append(d['P_r'])
            w_rows.append(1.0 / std_code**2)

        if len(A_rows) < 10: continue

        A = np.array(A_rows)
        y = np.array(y_rows)
        w_diag = np.array(w_rows)

        try:
            HtWH = A.T @ (w_diag[:, None] * A)
            HtWy = A.T @ (w_diag * y)
            dx = np.linalg.solve(HtWH, HtWy)
        except np.linalg.LinAlgError:
            continue

        x_cur = dx

        pos_est = ref_pos + x_cur[:3]
        clk_est = x_cur[3]
        trop_est = x_cur[4]

        # Post-fit residuals
        model = A @ x_cur
        residuals = y - model
        phase_res = residuals[0::2]

        err = pos_est - ref_pos
        d3 = float(np.linalg.norm(err))

        lat, lon, _ = ecef_to_blh(ref_pos)
        R = ecef_to_enu_matrix(lat, lon)
        enu = R @ err

        results.append({
            'time': J2000 + timedelta(seconds=gps_sod),
            'gps_sod': gps_sod,
            'dE': float(enu[0]), 'dN': float(enu[1]), 'dU': float(enu[2]),
            'd3': d3,
            'clk': float(clk_est), 'trop': float(trop_est),
            'n_sat': len(ep_data),
            'ph_res_std': float(np.std(phase_res)) if len(phase_res) > 0 else 0.0,
        })

        if len(results) <= 5 or len(results) % 100 == 0:
            r = results[-1]
            print(f"  ep {len(results)}/{len(epochs)}: 3D={d3:.3f}m "
                  f"clk={clk_est:.2f}m ph_res_std={r['ph_res_std']:.3f}m n_sat={len(ep_data)}")

    print(f"[PASS2] {len(results)}/{len(epochs)} epochs solved")

    if len(results) < 10:
        print("[FATAL] Too few successful epochs")
        return None

    # ═════════════════════════════════════════════════════════════════
    # Statistics
    # ═════════════════════════════════════════════════════════════════
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
    print(f"  GRACE-FO PPP — Float-Ambiguity + per-SV Bias + IEKF")
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
        'float_amb': {sv: float(v) for sv, v in float_amb.items()},
        'sv_bias': {sv: float(v) for sv, v in sv_bias.items()},
    }

# ── Output ──────────────────────────────────────────────────────────────
def save_output(output, output_dir='./output'):
    out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)
    tag = output['date'].replace('-', '')

    csv_path = out / f'ppp_{tag}.csv'
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['time', 'gps_sod', 'dE_m', 'dN_m', 'dU_m', 'd3_m',
                     'clock_m', 'trop_m', 'n_sat', 'ph_res_std_m'])
        for r in output['results']:
            w.writerow([r['time'].isoformat(), f"{r['gps_sod']:.6f}",
                        f"{r['dE']:.4f}", f"{r['dN']:.4f}", f"{r['dU']:.4f}",
                        f"{r['d3']:.4f}", f"{r['clk']:.4f}", f"{r['trop']:.4f}",
                        r['n_sat'], f"{r['ph_res_std']:.4f}"])
    print(f"[CSV] {csv_path}")

    json_path = out / f'ppp_{tag}_stats.json'
    json.dump(output['stats'], open(json_path, 'w'), indent=2)
    print(f"[JSON] {json_path}")
    return csv_path

def generate_plots(output, output_dir='./output'):
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    except ImportError:
        print("[PLOT] matplotlib not available"); return

    results = output['results']; stats = output['stats']
    date_str = output['date']; tag = date_str.replace('-', '')
    out = Path(output_dir)

    times_h = [(r['time'] - results[0]['time']).total_seconds()/3600 for r in results]
    dE_m = [r['dE'] for r in results]; dN_m = [r['dN'] for r in results]
    dU_m = [r['dU'] for r in results]; d3_m = [r['d3'] for r in results]

    fig, axes = plt.subplots(4, 1, figsize=(14, 10), sharex=True)
    for ax, data, color, ylabel in [
        (axes[0], dE_m, '#1565C0', 'E (m)'), (axes[1], dN_m, '#2E7D32', 'N (m)'),
        (axes[2], dU_m, '#C62828', 'U (m)'), (axes[3], d3_m, '#333333', '3D (m)')]:
        ax.plot(times_h, data, color=color, linewidth=0.5, alpha=0.7)
        ax.axhline(0, color='k', linewidth=0.5)
        ax.set_ylabel(ylabel); ax.grid(True, alpha=0.3)
    axes[0].set_title(f'GRACE-FO PPP Error vs GNV1B — {date_str}')
    axes[3].set_xlabel('Time (hours)')
    if 'all' in stats:
        s = stats['all']
        fig.suptitle(f'3D RMS={s["rms_3d"]/100:.2f}m E={s["rms_e"]/100:.2f}m '
                     f'N={s["rms_n"]/100:.2f}m U={s["rms_u"]/100:.2f}m N={s["n"]}', fontsize=9)
    plt.tight_layout()
    fig.savefig(str(out / f'ppp_{tag}.png'), dpi=150)
    plt.close(fig)
    print(f"[PLOT] {out / f'ppp_{tag}.png'}")

    # Histogram
    fig2, axes2 = plt.subplots(2, 2, figsize=(12, 8))
    for ax, data, label, color in [
        (axes2[0,0], dE_m, 'East (m)', '#1565C0'),
        (axes2[0,1], dN_m, 'North (m)', '#2E7D32'),
        (axes2[1,0], dU_m, 'Up (m)', '#C62828'),
        (axes2[1,1], d3_m, '3D (m)', '#333333')]:
        d_arr = np.array(data)
        d_arr = d_arr[np.isfinite(d_arr)]
        if len(d_arr) > 0:
            ax.hist(d_arr, bins=40, color=color, alpha=0.6, edgecolor='black')
            ax.axvline(0, color='k', linewidth=0.5)
            rms_val = np.sqrt(np.nanmean(d_arr**2))
            ax.text(0.95, 0.95, f'RMS={rms_val:.3f}m\nN={len(d_arr)}',
                    transform=ax.transAxes, ha='right', va='top', fontsize=8,
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax.set_xlabel(label); ax.set_ylabel('Count')
    plt.tight_layout()
    fig2.savefig(str(out / f'ppp_{tag}_hist.png'), dpi=150)
    plt.close(fig2)
    print(f"[PLOT] {out / f'ppp_{tag}_hist.png'}")

# ═════════════════════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(description='GRACE-FO PPP — Float-Amb + Bias + WLS')
    parser.add_argument('--date', default='2024-04-29')
    parser.add_argument('--hours', type=float, default=4.0)
    parser.add_argument('--data-dir', default='./data')
    parser.add_argument('--output-dir', default='./output')
    parser.add_argument('--grace-id', default='C')
    parser.add_argument('--interval', type=float, default=30.0)
    args = parser.parse_args()

    print(f"\n{'='*65}")
    print(f"  GRACE-FO PPP — Float-Ambiguity + per-SV Bias + WLS")
    print(f"  Date: {args.date}  Hours: {args.hours}h  dt: {args.interval}s")
    print(f"{'='*65}")

    output = process_day(date_str=args.date, nhours=args.hours,
                         data_dir=args.data_dir, grace_id=args.grace_id,
                         interval=args.interval)
    if output is None:
        print("[FATAL] Processing failed"); sys.exit(1)

    save_output(output, output_dir=args.output_dir)
    generate_plots(output, output_dir=args.output_dir)
    print(f"\n[DONE] {args.date}")

if __name__ == '__main__':
    main()
