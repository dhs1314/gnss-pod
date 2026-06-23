"""Multi-date / multi-satellite POD validation (Phase 13.0).

Decoupled from GRACE-FO specifics: auto-detects satellite config,
loads standard RINEX data, and runs EKF POD with consistent metrics.

Usage:
    py -3.12 run_multiday.py --mission GRACE-FO --sat-id C \
      --dates 2024-04-29 2024-04-30 --hours 0.17,0.5 --interval 30

    py -3.12 run_multiday.py --list-sats    (list all supported satellites)
"""
import sys, os, pickle, math, argparse
from pathlib import Path
from datetime import datetime, timedelta
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from satellite_config import get_config, get_ekf_params, list_satellites
from pod_io import (find_rinex_obs, find_reference_orbit,
                     load_precision_products, compute_metrics)
from gravity_model import read_icgem_gfc

# ── Constants ──
C_LIGHT = 299792458.0
SEC_PER_DAY = 86400.0
OMEGA_E = 7.2921151467e-5
J2000 = datetime(2000, 1, 1, 12, 0, 0)
MJD_J2000 = 51544.5
GPS_UTC_OFFSET = 18.0


def parse_args():
    p = argparse.ArgumentParser(description='Multi-date multi-satellite POD')
    p.add_argument('--mission', default='GRACE-FO', help='Mission name')
    p.add_argument('--sat-id', default='C', help='Satellite ID')
    p.add_argument('--dates', nargs='+', default=['2024-04-29'],
                   help='Dates to process (YYYY-MM-DD)')
    p.add_argument('--hours', default='0.17,0.5',
                   help='Arc lengths [h], comma-separated')
    p.add_argument('--interval', type=float, default=30.0)
    p.add_argument('--data-dir', default='d:/prj/gnss_pod/data')
    p.add_argument('--gravity-nmax', type=int, default=90)
    p.add_argument('--list-sats', action='store_true', default=False)
    p.add_argument('--output-dir', default='results/multiday')
    return p.parse_args()


def run_ekf_arc(gps1b_raw, ref_orbit, ref_vel, sp3, clk_data,
                ekf_config, sv_bias_global, dcb_if,
                epochs, mjd_start, r0_eci, v0_eci,
                doy, out_dir=None):
    """Run EKF for one arc segment."""
    from sequential_filter import SequentialEKF
    from coordinates import ecef_to_eci, eci_to_ecef

    ekf = SequentialEKF(ekf_config)
    state = ekf.initialize(r0_eci, v0_eci, mjd_start, epochs[0])

    epoch_dr = []
    r_ecef_list = []
    for i_ep, gps_sod in enumerate(epochs):
        mjd_utc = MJD_J2000 + gps_sod / 86400.0
        mjd_tt = mjd_utc + 69.184 / 86400.0
        if i_ep > 0:
            mjd_prev = MJD_J2000 + epochs[i_ep - 1] / 86400.0
            state = ekf.predict(state, gps_sod, mjd_prev,
                                 mjd_prev + 69.184 / 86400.0)

        # Simplified geometry computation
        rcv_ecef, _ = eci_to_ecef(state.r_eci, state.v_eci, mjd_utc)
        r_ecef_list.append(rcv_ecef)

        # Use the main pipeline's compute_epoch_geometry
        from run_sequential_pod import compute_epoch_geometry, interpolate_ref
        ep_data = compute_epoch_geometry(gps_sod, gps1b_raw, sp3, rcv_ecef, clk_data)
        if not ep_data:
            continue

        sv_bias = {}; sv_bias_ref = 0.0
        state, stats = ekf.process_epoch(state, ep_data, sp3, sv_bias,
                                          sv_bias_ref, mjd_utc, mjd_tt, doy)
        r_ecef, _ = eci_to_ecef(state.r_eci, state.v_eci, mjd_utc)
        r_gnv = interpolate_ref(ref_orbit, gps_sod)
        if r_gnv is not None:
            epoch_dr.append(np.linalg.norm(r_ecef - r_gnv))

    metrics = compute_metrics([], [], [])
    if epoch_dr:
        r_ecef_for_metrics = r_ecef_list[:len(epoch_dr)]
        ref_for_metrics = [interpolate_ref(ref_orbit, ep)
                           for ep in epochs[:len(epoch_dr)]]
        metrics = compute_metrics(r_ecef_for_metrics, ref_for_metrics, epochs)
    return metrics


def main():
    args = parse_args()
    if args.list_sats:
        list_satellites()
        return

    # ── Satellite config ──
    sat_cfg = get_config(args.mission, args.sat_id)
    ekf_params = get_ekf_params(args.mission, args.sat_id)
    arc_hours_list = [float(h) for h in args.hours.split(',')]

    print(f"=" * 80)
    print(f"Multi-Date POD Validation — {args.mission}-{args.sat_id}")
    print(f"  {sat_cfg['description']}")
    print(f"  Altitude: {sat_cfg['altitude_km']}km  "
          f"Mass: {sat_cfg['mass_kg']}kg  "
          f"GNSS: {', '.join(sat_cfg['gnss_systems'])}")
    print(f"  Dates: {len(args.dates)}  Arcs: {args.hours}h  "
          f"Step: {args.interval}s")
    print(f"=" * 80)

    dp = Path(args.data_dir)

    # ── Load shared products ──
    y, m, d_ = [int(x) for x in args.dates[0].split('-')]
    doy_first = (datetime(y, m, d_) - datetime(y, 1, 1)).days + 1
    products = load_precision_products(args.data_dir, y, doy_first)
    sp3 = products['sp3']
    clk_data = products['clk_data']
    print(f"  SP3: {'OK' if sp3 else 'MISSING'}  "
          f"CLK: {'OK' if clk_data else 'MISSING'}")

    # ── Load gravity (once) ──
    gfc_path = dp / 'gravity' / 'GGM05C.gfc'
    Cnm, Snm, Nmax, GM_grav, R_grav = read_icgem_gfc(str(gfc_path))
    Nmax = min(Nmax, args.gravity_nmax)
    print(f"  Gravity: GGM05C Nmax={Nmax}")

    # ── Build EKF config (satellite-specific) ──
    ekf_config = {
        **ekf_params,
        'bodies': ['Sun', 'Moon'],
        'Cnm': Cnm, 'Snm': Snm, 'GM_grav': GM_grav, 'R_grav': R_grav,
        'gravity_nmax': Nmax,
        'sigma_acc_process': 1e-3, 'tau_emp': 600.0, 'sigma_emp_ss': 1e-8,
        'sigma_zwd_rw': 1e-9, 'sigma_phase': 0.20, 'sigma_code': 0.30,
        'chi2_threshold': 100.0, 'el_min': 0.087,
        'use_phase_windup': True, 'use_relativity': True, 'use_cycle_slip': False,
        'ar_min_epochs': 6,
        'antex_data': products.get('antex_data'),
        'dcb_data': products.get('dcb_data'),
        'elev_exp_phase': 1.0, 'elev_exp_code': 0.70,
        'mw_max_epochs': 200,
    }

    if products.get('antex_data'):
        print(f"  ANTEX: OK")
    if products.get('dcb_data'):
        print(f"  DCB: OK")

    # ── Process all dates ──
    all_results = defaultdict(dict)  # date → arc_h → metrics

    for date_str in args.dates:
        print(f"\n--- {date_str} ---")
        y, m, d_ = [int(x) for x in date_str.split('-')]
        doy = (datetime(y, m, d_) - datetime(y, 1, 1)).days + 1

        # Find RINEX
        rinex_path = find_rinex_obs(args.data_dir, args.mission, args.sat_id, date_str)
        if not rinex_path:
            print(f"  No RINEX data found for {date_str}")
            continue
        print(f"  RINEX: {rinex_path}")

        # Load GPS1B
        if str(rinex_path).endswith('.pkl'):
            gps1b_raw = pickle.load(open(str(rinex_path), 'rb'))
        else:
            from gps1b_rnx_loader import load_gps1b_rnx
            gps1b_raw = load_gps1b_rnx(str(rinex_path))

        # Find reference orbit
        ref_path, ref_type = find_reference_orbit(args.data_dir, args.mission, args.sat_id, date_str)
        if not ref_path:
            print(f"  No reference orbit found for {date_str}")
            continue
        print(f"  Ref orbit ({ref_type}): {ref_path}")

        # Load reference
        from run_sequential_pod import load_gnv1b, interpolate_ref
        ref_orbit, ref_vel = load_gnv1b(str(ref_path))

        for arc_hours in arc_hours_list:
            gps_sod_start = min(gps1b_raw.keys())
            gps_sod_end = gps_sod_start + arc_hours * 3600

            # Epoch selection
            epochs = sorted(set(
                gps_sod for gps_sod in gps1b_raw.keys()
                if gps_sod_start <= gps_sod <= gps_sod_end
                and abs((gps_sod - gps_sod_start) % args.interval) <= 2.0
            ))
            if len(epochs) < 6:
                print(f"  {arc_hours}h: only {len(epochs)} epochs, skip")
                continue

            # Initial state from ref orbit
            r0_ecef = interpolate_ref(ref_orbit, epochs[0])
            v0_ecef = interpolate_ref(ref_vel, epochs[0])
            if r0_ecef is None or v0_ecef is None:
                print(f"  {arc_hours}h: no ref at epoch 0, skip")
                continue

            from coordinates import ecef_to_eci
            mjd_start = MJD_J2000 + epochs[0] / 86400.0
            r0_eci, v0_eci = ecef_to_eci(r0_ecef, v0_ecef, mjd_start)

            # DCB computation
            dcb_if = {}
            if products.get('dcb_data'):
                from precision_products import compute_dcb_if_correction
                for prn in products['dcb_data']:
                    dcb_if[prn] = compute_dcb_if_correction(products['dcb_data'], prn)

            # Run EKF
            metrics = run_ekf_arc(gps1b_raw, ref_orbit, ref_vel, sp3, clk_data,
                                   ekf_config, {}, dcb_if,
                                   epochs, mjd_start, r0_eci, v0_eci, doy)

            rms = metrics.get('rms_3d', 0)
            all_results[date_str][arc_hours] = metrics
            print(f"  {arc_hours}h: RMS={rms:.3f}m  "
                  f"(mean={metrics.get('mean_3d',0):.3f}  "
                  f"max={metrics.get('max_3d',0):.3f}  "
                  f"n={metrics.get('n_epochs',0)})")

    # ── Summary ──
    print(f"\n{'='*80}")
    print(f"Summary — {args.mission}-{args.sat_id}")
    print(f"{'='*80}")

    rms_by_arc = defaultdict(list)
    for date_str, arcs in all_results.items():
        for arc_h, metrics in arcs.items():
            rms_by_arc[arc_h].append(metrics['rms_3d'])

    for arc_h, vals in sorted(rms_by_arc.items()):
        v = np.array(vals)
        print(f"  {arc_h}h: n={len(v)}  mean={np.mean(v):.3f}m  "
              f"median={np.median(v):.3f}m  std={np.std(v):.3f}m  "
              f"best={np.min(v):.3f}m  worst={np.max(v):.3f}m")

    # ── Save ──
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"multiday_{args.mission}_{args.sat_id}.pkl"
    pickle.dump({
        'results': dict(all_results),
        'args': vars(args),
        'sat_config': sat_cfg,
    }, open(str(out_path), 'wb'))
    print(f"\n  Saved: {out_path}")

    # ── Plot ──
    if len(all_results) >= 2:
        fig, axes = plt.subplots(len(arc_hours_list), 1, figsize=(12, 4*len(arc_hours_list)))
        if len(arc_hours_list) == 1:
            axes = [axes]

        for idx, arc_h in enumerate(sorted(arc_hours_list)):
            ax = axes[idx]
            dates_plt = sorted(all_results.keys())
            rms_plt = [all_results[d].get(arc_h, {}).get('rms_3d', 0) for d in dates_plt]
            ax.plot(range(len(dates_plt)), rms_plt, 'o-', color='#2196F3', markersize=6)
            ax.axhline(y=np.mean(rms_plt), color='red', linestyle='--', alpha=0.5,
                       label=f'Mean {np.mean(rms_plt):.3f}m')
            ax.set_xticks(range(len(dates_plt)))
            ax.set_xticklabels([d[5:] for d in dates_plt], rotation=45, fontsize=8)
            ax.set_ylabel(f'{arc_h}h 3D RMS [m]')
            ax.set_title(f'{args.mission}-{args.sat_id} {arc_h}h Arc POD Accuracy')
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)

        plt.tight_layout()
        png_path = out_dir / f"multiday_{args.mission}_{args.sat_id}.png"
        plt.savefig(str(png_path), dpi=150, bbox_inches='tight')
        print(f"  Saved plot: {png_path}")


if __name__ == '__main__':
    main()
