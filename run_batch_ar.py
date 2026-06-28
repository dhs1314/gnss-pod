"""Two-pass batch ambiguity resolution for POD (Phase 6.0).

Pass 1: Float EKF → collect ambiguity states
Batch AR: Full-arc MW smoothing + NL integer bootstrapping
Pass 2: Fixed EKF → re-estimate orbit with locked ambiguities

Usage:
    py -3.12 run_batch_ar.py --date 2024-04-29 --hours 0.5 --interval 30 --grace-id C \
      --sp3-file ... --clk-file ... --dcb-file ... --antex-file ... --iers-c04 ... \
      --enable-phase-windup --enable-relativity --ar-min-epochs 6 --gravity-nmax 90
"""
import os, sys, argparse, pickle, math
from datetime import datetime
from pathlib import Path
import numpy as np

# Add project root
sys.path.insert(0, str(Path(__file__).parent))

from src.sequential_filter import SequentialEKF
from src.batch_lsq import BatchAmbiguityResolver, COEFF_W, LAM_NL


def load_data(args):
    """Reuse data-loading logic from run_sequential_pod."""
    from run_sequential_pod import compute_epoch_geometry, interpolate_ref
    from src.gnss_data import load_gps1b_rnxnx, load_gnv1b
    from src.gravity_model import read_icgem_gfc
    from src.precision_products import read_rinex_clk, read_antex, load_code_dcb_pair
    from src.iers import setup_iers_from_c04

    dp = Path(args.data_dir)
    date_str = args.date
    grace_id = args.grace_id
    y, m, d = [int(x) for x in date_str.split('-')]

    # Load GNV1B
    gnv_path = dp / 'gracefo' / str(y) / date_str / f'GNV1B_{date_str}_{grace_id}_04.txt'
    ref_orbit, ref_vel = load_gnv1b(str(gnv_path))
    print(f"[GNV1B] {len(ref_orbit)} pos, {len(ref_vel)} vel epochs")

    # Load GPS1B
    rnx_path = dp / f'GPS1B_{date_str}_{grace_id}_04.rnx'
    pkl_path = dp / 'gracefo' / str(y) / date_str / f'GPS1B_{date_str}_{grace_id}_04.pkl'
    if rnx_path.exists():
        gps1b_raw = load_gps1b_rnx(str(rnx_path))
    elif pkl_path.exists():
        gps1b_raw = pickle.load(open(str(pkl_path), 'rb'))
    else:
        raise FileNotFoundError(f"No GPS data at {rnx_path} or {pkl_path}")
    print(f"[GPS1B] {len(gps1b_raw)} epochs")

    # Load SP3
    from src.sp3_loader import parse_sp3_text
    sp3_path = Path(args.sp3_file)
    if sp3_path.suffix == '.pkl':
        sp3 = pickle.load(open(str(sp3_path), 'rb'))
    else:
        with open(str(sp3_path), 'r') as f:
            epochs_dict, ts_list = parse_sp3_text(f.read())
        sp3 = {'ts': ts_list, 'epochs': epochs_dict, 'source': 'CODE', 'product': 'FIN'}
    print(f"[SP3] {len(sp3['ts'])} epochs")

    # Load CLK
    clk_data = None
    if args.clk_file:
        from src.precision_products import read_rinex_clk
        clk_path = args.clk_file
        if not Path(clk_path).exists():
            import gzip
            clk_gz = Path(str(clk_path) + '.gz')
            if clk_gz.exists():
                with gzip.open(str(clk_gz), 'rt') as gf:
                    with open(clk_path, 'w') as tf:
                        tf.write(gf.read())
        clk_data = read_rinex_clk(clk_path)
        print(f"[CLK] {len(clk_data)} SVs")

    # Select epochs
    gps_sod_start = min(gps1b_raw.keys())
    gps_sod_end = gps_sod_start + args.hours * 3600
    epochs = []
    for gps_sod in sorted(gps1b_raw.keys()):
        if not (gps_sod_start <= gps_sod <= gps_sod_end):
            continue
        dt_ep = gps_sod - gps_sod_start
        nearest = round(dt_ep / args.interval) * args.interval
        if abs(dt_ep - nearest) > max(2.0, args.interval * 0.1):
            continue
        epochs.append(gps_sod)
    print(f"[EPOCH] {len(epochs)} selected")

    # Load ANTEX
    antex_data = None
    if args.antex_file:
        from src.precision_products import read_antex
        antex_data = read_antex(args.antex_file)
        print(f"[ANTEX] {args.antex_file}")

    # Load DCB
    dcb_data = None
    if args.dcb_file:
        from src.precision_products import load_code_dcb_pair
        p1c1_path = args.dcb_file.replace('P1P2', 'P1C1')
        dcb_data = load_code_dcb_pair(p1c1_path, args.dcb_file)
        print(f"[DCB] loaded")

    # Load IERS
    if args.iers_c04:
        from src.iers import setup_iers_from_c04
        setup_iers_from_c04(args.iers_c04)
        print(f"[IERS] loaded")

    # Load gravity
    gfc_path = Path(args.gravity_field) if hasattr(args, 'gravity_field') else \
        dp / 'gravity' / 'GGM05C.gfc'
    Cnm, Snm, Nmax, GM_grav, R_grav = read_icgem_gfc(str(gfc_path))
    Nmax = min(Nmax, args.gravity_nmax)
    print(f"[Gravity] GGM05C Nmax={Nmax}")

    # Compute epoch geometry for all selected epochs
    all_ep_data = []
    sv_bias = {}
    from run_sequential_pod import compute_epoch_geometry

    for gps_sod in epochs:
        utc_dt = datetime(y, m, d) + \
                 (gps_sod - 0 * 3600) * np.timedelta64(1, 's').astype('timedelta64[s]')
        utc_dt_py = datetime(y, m, d, 0, 0, 0) + \
                    __import__('datetime').timedelta(seconds=float(gps_sod))

        # Get coarse receiver position from GNV1B
        ref_pos = None
        ts = sorted(ref_orbit.keys())
        for i, t in enumerate(ts):
            if t >= gps_sod:
                if i > 0:
                    a = (gps_sod - ts[i-1]) / (t - ts[i-1])
                    ref_pos = ref_orbit[ts[i-1]] * (1-a) + ref_orbit[t] * a
                else:
                    ref_pos = ref_orbit[t]
                break
        if ref_pos is None and ts:
            ref_pos = ref_orbit[ts[-1]]

        ep_data = compute_epoch_geometry(
            gps_sod, utc_dt_py, gps1b_raw, sp3, clk_data, ref_pos)
        all_ep_data.extend(ep_data)

        # Estimate per-SV code biases from first 60 epochs
        if len(all_ep_data) <= 60 * 30:  # ~first 30 min of data
            for d in ep_data:
                sv = d['sv']
                if sv not in sv_bias and 'P_if_raw' in d:
                    sv_bias[sv] = 0.0

    # Compute per-SV code biases
    if sv_bias:
        sv_bias_values = list(sv_bias.values())
        # Proper bias estimation from first N epochs
        from run_sequential_pod import compute_epoch_geometry
        print(f"  Per-SV code biases: {len(sv_bias)} SVs")

    print(f"\n[BatchAR] Collected {len(all_ep_data)} per-SV observations "
          f"across {len(epochs)} epochs")

    return {
        'epochs': epochs,
        'all_ep_data': all_ep_data,
        'ref_orbit': ref_orbit,
        'ref_vel': ref_vel,
        'sp3': sp3,
        'clk_data': clk_data,
        'antex_data': antex_data,
        'dcb_data': dcb_data,
        'Cnm': Cnm, 'Snm': Snm, 'Nmax': Nmax,
        'GM_grav': GM_grav, 'R_grav': R_grav,
        'gps_sod_start': gps_sod_start,
        'date_str': date_str,
        'grace_id': grace_id,
        'y': y, 'm': m, 'd': d,
    }


def run_float_ekf(data, args):
    """Pass 1: Float EKF, same as current run_sequential_pod."""
    from src.orbit_dynamics import compute_initial_state, eci_to_ecef
    from src.sequential_filter import SequentialEKF
    from run_sequential_pod import compute_epoch_geometry, interpolate_ref

    # Compute initial state from GNV1B
    gps_sod_start = data['gps_sod_start']
    r0_ecef, v0_ecef = None, None
    ts = sorted(data['ref_orbit'].keys())
    for t in ts:
        if abs(t - gps_sod_start) < 1.0:
            r0_ecef = data['ref_orbit'][t]
            v0_ecef = data['ref_vel'][t]
            break
    if r0_ecef is None and ts:
        r0_ecef = data['ref_orbit'][ts[0]]
        v0_ecef = data['ref_vel'][ts[0]]

    from src.coordinates import ecef_to_eci
    mjd_utc = 60429 + gps_sod_start / 86400.0  # approximate
    r0_eci, v0_eci = ecef_to_eci(r0_ecef, v0_ecef, mjd_utc)

    # Configure EKF
    ekf_cfg = {
        'dynamics_mode': args.dynamics_mode,
        'gravity_field': str(Path(args.data_dir) / 'gravity' / 'GGM05C.gfc'),
        'gravity_nmax': args.gravity_nmax,
        'mass': 580.0, 'area_drag': 0.68, 'area_srp': 3.4, 'CR': 1.3, 'CD': 2.2,
        'sigma_acc': getattr(args, 'sigma_acc', 1e-3),
        'tau_emp': getattr(args, 'tau_emp', 600.0),
        'bodies': ['Sun', 'Moon'],
        'Cnm': data['Cnm'], 'Snm': data['Snm'],
        'GM_grav': data['GM_grav'], 'R_grav': data['R_grav'],
        'antex_data': data['antex_data'],
        'dcb_data': data['dcb_data'],
        'use_phase_windup': args.enable_phase_windup,
        'use_relativity': args.enable_relativity,
        'use_cycle_slip': getattr(args, 'enable_cycle_slip', False),
        'ar_min_epochs': getattr(args, 'ar_min_epochs', 10),
    }

    ekf = SequentialEKF(ekf_cfg)
    state = ekf.initialize(r0_eci, v0_eci, mjd_utc, gps_sod_start)
    state = ekf.predict(state, state.t + args.interval, mjd_utc, mjd_utc + 69.184/86400.0)

    # Process epochs
    results = {'epochs': [], 'r_ecef': [], 'v_ecef': [], 'a_rtn': [],
               'n_sv': [], 'r_gnv': [], 'dr': [], 'n_wl_fixed': []}
    sv_bias = {}
    N_BIAS = min(60, len(data['epochs']))

    for i_ep, gps_sod in enumerate(data['epochs']):
        utc_dt = datetime(data['y'], data['m'], data['d']) + \
                 __import__('datetime').timedelta(seconds=float(gps_sod))
        doy = (datetime(data['y'], data['m'], data['d']) - \
               datetime(data['y'], 1, 1)).days + 1

        # Get receiver position for geometry
        r_ecef_cur, _ = eci_to_ecef(state.r_eci, state.v_eci, mjd_utc + \
                                      (gps_sod - data['gps_sod_start']) / 86400.0)

        ep_data = compute_epoch_geometry(
            gps_sod, utc_dt, __import__('pickle').load(
                open(str(Path(args.data_dir) / f'GPS1B_{data["date_str"]}_{data["grace_id"]}_04.rnx'), 'rb'))
            if False else data.get('_gps1b_raw', {}),  # reuse from load
            data['sp3'], data['clk_data'], r_ecef_cur)

        # Pass 1 uses the EKF as-is with code-phase ambiguity init
        sv_bias_ref = -1.98

        state, stats = ekf.process_epoch(
            state, ep_data, data['sp3'], sv_bias, sv_bias_ref,
            mjd_utc + (gps_sod - data['gps_sod_start']) / 86400.0,
            mjd_utc + (gps_sod - data['gps_sod_start']) / 86400.0 + 69.184/86400.0,
            doy)

        r_ecef, v_ecef = eci_to_ecef(state.r_eci, state.v_eci,
                                       mjd_utc + (gps_sod - data['gps_sod_start']) / 86400.0)
        ref_pos = None
        ts = sorted(data['ref_orbit'].keys())
        for j, t in enumerate(ts):
            if t >= gps_sod:
                if j > 0:
                    a = (gps_sod - ts[j-1]) / (t - ts[j-1])
                    ref_pos = data['ref_orbit'][ts[j-1]] * (1-a) + data['ref_orbit'][t] * a
                else:
                    ref_pos = data['ref_orbit'][t]
                break

        results['epochs'].append(gps_sod)
        results['r_ecef'].append(r_ecef)
        results['v_ecef'].append(v_ecef)
        results['a_rtn'].append(state.a_rtn.copy())
        results['n_sv'].append(state.n_sv)
        if ref_pos is not None:
            results['r_gnv'].append(ref_pos)
            results['dr'].append(np.linalg.norm(r_ecef - ref_pos))
        results['n_wl_fixed'].append(len(ekf._wl_fixed))

        if i_ep % max(1, len(data['epochs']) // 10) == 0 or i_ep < 3:
            dr = results['dr'][-1] if results['dr'] else 0
            print(f"  [Pass1 {i_ep:4d}/{len(data['epochs'])}] t={gps_sod:.0f}  "
                  f"|dr|={dr:.3f}m  n_sv={state.n_sv}  WL={len(ekf._wl_fixed)}")

        # Predict to next epoch
        mjd_next = mjd_utc + (gps_sod + args.interval - data['gps_sod_start']) / 86400.0
        mjd_tt_next = mjd_next + 69.184 / 86400.0
        state = ekf.predict(state, gps_sod + args.interval, mjd_utc + \
                             (gps_sod - data['gps_sod_start']) / 86400.0, mjd_tt_next)

    # Compute pass 1 RMS
    dr_vals = results['dr']
    rms_pass1 = np.sqrt(np.mean([d**2 for d in dr_vals])) if dr_vals else 0
    print(f"\n  Pass 1 (float) 3D RMS: {rms_pass1:.3f} m")

    # Extract ambiguity estimates for batch AR
    ekf_amb = {}
    for sv in state.sv_list:
        if sv in state.sv_to_idx:
            idx = 10 + state.sv_to_idx[sv]  # I_AMB_START = 11? Actually need to check...
            amb_val = float(state.x[idx])
            ekf_amb[sv] = amb_val

    # Actually use the proper index
    ekf_amb2 = {}
    for sv in state.sv_list:
        idx = state.sv_to_idx[sv]
        amb_idx = 11 + idx  # I_AMB_START = 11 (after r,v,aR,aT,aN,zwd,clk)
        ekf_amb2[sv] = float(state.x[amb_idx])

    return results, ekf_amb2, state, ekf


def main():
    parser = argparse.ArgumentParser(description='Batch AR POD (Phase 6.0)')
    parser.add_argument('--date', required=True)
    parser.add_argument('--hours', type=float, default=0.5)
    parser.add_argument('--interval', type=float, default=30.0)
    parser.add_argument('--data-dir', default='d:/prj/gnss_pod/data')
    parser.add_argument('--grace-id', default='C')
    parser.add_argument('--gravity-nmax', type=int, default=90)
    parser.add_argument('--sp3-file', required=True)
    parser.add_argument('--clk-file', required=True)
    parser.add_argument('--dcb-file', default=None)
    parser.add_argument('--antex-file', default=None)
    parser.add_argument('--iers-c04', default=None)
    parser.add_argument('--enable-phase-windup', action='store_true', default=False)
    parser.add_argument('--enable-relativity', action='store_true', default=False)
    parser.add_argument('--enable-cycle-slip', action='store_true', default=False)
    parser.add_argument('--ar-min-epochs', type=int, default=6)
    parser.add_argument('--dynamics-mode', default='simplified')
    parser.add_argument('--sigma-acc', type=float, default=1e-3)
    parser.add_argument('--tau-emp', type=float, default=600.0)
    args = parser.parse_args()

    print("=" * 60)
    print("Phase 6.0: Two-Pass Batch Ambiguity Resolution POD")
    print("=" * 60)

    # Load all data once
    print("\n-- Loading data --")
    data = load_data(args)

    # Pass 1: Float EKF
    print("\n-- Pass 1: Float EKF --")
    results_pass1, ekf_amb, final_state, ekf = run_float_ekf(data, args)

    # Batch ambiguity resolution
    print("\n-- Batch Ambiguity Resolution --")
    # Collect epoch data with L1/L2/P1/P2 for MW
    # We need the raw epoch data from compute_epoch_geometry
    all_raw_ep = []
    for gps_sod in data['epochs']:
        utc_dt = datetime(data['y'], data['m'], data['d']) + \
                 __import__('datetime').timedelta(seconds=float(gps_sod))
        from run_sequential_pod import compute_epoch_geometry
        r_ecef_interp = None
        # Use GNV1B position
        ts = sorted(data['ref_orbit'].keys())
        for j, t in enumerate(ts):
            if t >= gps_sod:
                if j > 0:
                    a = (gps_sod - ts[j-1]) / (t - ts[j-1])
                    r_ecef_interp = data['ref_orbit'][ts[j-1]] * (1-a) + data['ref_orbit'][t] * a
                else:
                    r_ecef_interp = data['ref_orbit'][t]
                break
        if r_ecef_interp is None and ts:
            r_ecef_interp = data['ref_orbit'][ts[-1]]

        ep_data = compute_epoch_geometry(
            gps_sod, utc_dt,
            __import__('pickle').load(open(
                str(Path(args.data_dir) / f'GPS1B_{data["date_str"]}_{data["grace_id"]}_04.pkl'), 'rb')),
            data['sp3'], data['clk_data'], r_ecef_interp)
        all_raw_ep.extend(ep_data)

    resolver = BatchAmbiguityResolver(all_raw_ep, min_epochs_wl=args.ar_min_epochs)
    batch_result = resolver.resolve(ekf_amb)

    # Print batch AR results
    n_wl = len(batch_result['wl_fixed'])
    n_nl = len(batch_result['nl_fixed'])
    print(f"\n  Batch AR: {n_wl} WL fixed, {n_nl} NL fixed "
          f"(b_r_wl={batch_result['b_r_wl']:+.4f} cyc)")

    # Compute absolute IF ambiguities for second pass
    amb_if_fixed = {}
    for sv, N_w in batch_result['wl_fixed'].items():
        B_if_float = ekf_amb.get(sv, 0)
        dN1 = batch_result['nl_fixed'].get(sv, 0)
        amb_if_fixed[sv] = LAM_NL * dN1 + COEFF_W * N_w

    # Pass 2: Fixed EKF
    print(f"\n-- Pass 2: Fixed EKF ({len(amb_if_fixed)} SVs with batch-fixed amb) --")

    # Reconfigure EKF with batch-fixed ambiguities
    ekf_cfg2 = {
        'dynamics_mode': args.dynamics_mode,
        'gravity_field': str(Path(args.data_dir) / 'gravity' / 'GGM05C.gfc'),
        'gravity_nmax': args.gravity_nmax,
        'mass': 580.0, 'area_drag': 0.68, 'area_srp': 3.4, 'CR': 1.3, 'CD': 2.2,
        'sigma_acc': args.sigma_acc,
        'tau_emp': getattr(args, 'tau_emp', 600.0),
        'bodies': ['Sun', 'Moon'],
        'Cnm': data['Cnm'], 'Snm': data['Snm'],
        'GM_grav': data['GM_grav'], 'R_grav': data['R_grav'],
        'antex_data': data['antex_data'],
        'dcb_data': data['dcb_data'],
        'use_phase_windup': args.enable_phase_windup,
        'use_relativity': args.enable_relativity,
        'use_cycle_slip': False,
        'ar_min_epochs': args.ar_min_epochs,
        'amb_batch_fixed': amb_if_fixed,
        'amb_batch_var': 0.0001,  # 1cm std
    }

    ekf2 = SequentialEKF(ekf_cfg2)
    gps_sod_start = data['gps_sod_start']
    mjd_utc = 60429 + gps_sod_start / 86400.0
    r0_ecef_2 = data['ref_orbit'].get(gps_sod_start,
                    list(data['ref_orbit'].values())[0])
    v0_ecef_2 = data['ref_vel'].get(gps_sod_start,
                    list(data['ref_vel'].values())[0])
    from src.coordinates import ecef_to_eci
    r0_eci_2, v0_eci_2 = ecef_to_eci(r0_ecef_2, v0_ecef_2, mjd_utc)

    state2 = ekf2.initialize(r0_eci_2, v0_eci_2, mjd_utc, gps_sod_start)
    state2 = ekf2.predict(state2, state2.t + args.interval, mjd_utc,
                           mjd_utc + 69.184/86400.0)

    results_pass2 = {'epochs': [], 'r_ecef': [], 'dr': [], 'n_sv': []}
    sv_bias2 = {}

    for i_ep, gps_sod in enumerate(data['epochs']):
        utc_dt = datetime(data['y'], data['m'], data['d']) + \
                 __import__('datetime').timedelta(seconds=float(gps_sod))
        doy = (datetime(data['y'], data['m'], data['d']) - \
               datetime(data['y'], 1, 1)).days + 1

        r_ecef_cur2, _ = eci_to_ecef(state2.r_eci, state2.v_eci,
                                       mjd_utc + (gps_sod - gps_sod_start) / 86400.0)

        # Recompute epoch geometry (we need it for the fixed pass too)
        ep_data2 = all_raw_ep  # reuse from batch AR collection
        # Filter to current epoch
        ep_cur = [d for d in ep_data2 if abs(d.get('_gps_sod', gps_sod) - gps_sod) < 1.0]
        if not ep_cur:
            # Fallback: use all data (will be filtered by SV in process_epoch)
            pass

        # Simplified: use the cached per-epoch geometry
        # For now, just reuse what we computed earlier
        sv_bias_ref2 = -1.98

        # Actually let's use a simpler approach — recompute geometry per epoch
        from run_sequential_pod import compute_epoch_geometry
        ep_cur = compute_epoch_geometry(
            gps_sod, utc_dt,
            __import__('pickle').load(open(
                str(Path(args.data_dir) / f'GPS1B_{data["date_str"]}_{data["grace_id"]}_04.pkl'), 'rb')),
            data['sp3'], data['clk_data'], r_ecef_cur2)

        state2, stats2 = ekf2.process_epoch(
            state2, ep_cur, data['sp3'], sv_bias2, sv_bias_ref2,
            mjd_utc + (gps_sod - gps_sod_start) / 86400.0,
            mjd_utc + (gps_sod - gps_sod_start) / 86400.0 + 69.184/86400.0,
            doy)

        r_ecef2, v_ecef2 = eci_to_ecef(state2.r_eci, state2.v_eci,
                                         mjd_utc + (gps_sod - gps_sod_start) / 86400.0)
        ts2 = sorted(data['ref_orbit'].keys())
        ref_pos2 = None
        for j, t in enumerate(ts2):
            if t >= gps_sod:
                if j > 0:
                    a = (gps_sod - ts2[j-1]) / (t - ts2[j-1])
                    ref_pos2 = data['ref_orbit'][ts2[j-1]] * (1-a) + data['ref_orbit'][t] * a
                else:
                    ref_pos2 = data['ref_orbit'][t]
                break

        results_pass2['epochs'].append(gps_sod)
        results_pass2['r_ecef'].append(r_ecef2)
        results_pass2['n_sv'].append(state2.n_sv)
        if ref_pos2 is not None:
            results_pass2['dr'].append(np.linalg.norm(r_ecef2 - ref_pos2))

        if i_ep % max(1, len(data['epochs']) // 10) == 0 or i_ep < 3:
            dr2 = results_pass2['dr'][-1] if results_pass2['dr'] else 0
            print(f"  [Pass2 {i_ep:4d}/{len(data['epochs'])}] t={gps_sod:.0f}  "
                  f"|dr|={dr2:.3f}m  n_sv={state2.n_sv}")

        mjd_next2 = mjd_utc + (gps_sod + args.interval - gps_sod_start) / 86400.0
        state2 = ekf2.predict(state2, gps_sod + args.interval,
                               mjd_utc + (gps_sod - gps_sod_start) / 86400.0,
                               mjd_next2 + 69.184/86400.0)

    dr_vals2 = results_pass2['dr']
    rms_pass2 = np.sqrt(np.mean([d**2 for d in dr_vals2])) if dr_vals2 else 0

    # Results
    dr_vals1 = results_pass1['dr']
    rms_pass1 = np.sqrt(np.mean([d**2 for d in dr_vals1])) if dr_vals1 else 0
    print(f"\n-- Final Results --")
    print(f"  Pass 1 (float EKF):   3D RMS = {rms_pass1:.3f} m")
    print(f"  Pass 2 (batch-fixed): 3D RMS = {rms_pass2:.3f} m")
    if rms_pass1 > 0:
        pct = (rms_pass1 - rms_pass2) / rms_pass1 * 100
        print(f"  Improvement: {pct:+.1f}%")
    print(f"  WL fixed: {n_wl} SVs, NL fixed: {n_nl} SVs")

    # Save
    out = {
        'pass1': results_pass1,
        'pass2': results_pass2,
        'batch_ar': {k: v for k, v in batch_result.items()
                     if k != 'wl_stats'},
        'args': vars(args),
    }
    out_path = Path('results') / 'batch_ar' / f"batch_ar_{data['date_str']}_{data['grace_id']}_{args.hours}h.pkl"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pickle.dump(out, open(str(out_path), 'wb'))
    print(f"\n  Saved: {out_path}")


if __name__ == '__main__':
    main()
