#!/usr/bin/env python3
"""
Batch Least Squares POD for GRACE-FO.

Estimates initial state [r0, v0] and force parameters (Cd, CR) from
GNSS observations over an arc using ECI orbit integration with the
full force model (GGM05C gravity + third-body + drag + SRP).

Usage:
  py -3.12 run_batch_pod.py --date 2024-04-29 --hours 2 --interval 30
"""
import sys, os, pickle, csv, math, argparse
from pathlib import Path
from datetime import datetime, timedelta
from collections import defaultdict
import numpy as np

sys.path.insert(0, str(Path(__file__).parent / 'src'))
from gps1b_rnx_loader import load_gps1b_rnx
from sp3_loader import get_gps_pos_from_sp3 as _sp3_get
from orbit_dynamics import GM, OMEGA_E
from coordinates import ecef_to_eci, eci_to_ecef
from gravity_model import read_icgem_gfc
from batch_estimator import run_batch_lsq

C_LIGHT = 299792458.0
SEC_PER_DAY = 86400.0
J2000 = datetime(2000, 1, 1, 12, 0, 0)
MJD_J2000 = 51544.5
GPS_UTC_OFFSET = 18.0  # seconds, valid since 2017-01-01


def load_gnv1b(filepath):
    """Load GNV1B reference orbit (position + velocity).

    Returns (pos_dict, vel_dict) where each maps gps_sod -> np.array([x,y,z]) or [vx,vy,vz].
    """
    pos_orbit, vel_orbit = {}, {}
    with open(filepath, encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 12:
                continue
            try:
                t = float(parts[0])
                flag = parts[2]
                if flag not in ('C', 'E'):
                    continue
                X, Y, Z = float(parts[3]), float(parts[4]), float(parts[5])
                VX, VY, VZ = float(parts[9]), float(parts[10]), float(parts[11])
                if abs(X) < 1e3:
                    continue
                pos_orbit[t] = np.array([X, Y, Z])
                vel_orbit[t] = np.array([VX, VY, VZ])
            except (ValueError, IndexError):
                continue
    return pos_orbit, vel_orbit


def interpolate_ref(ref_orbit, gps_sod):
    ts = sorted(ref_orbit.keys())
    t0 = t1 = None
    for i, ti in enumerate(ts):
        if ti >= gps_sod:
            t1 = ti
            t0 = ts[i - 1] if i > 0 else None
            break
        t0 = ti
    if t1 is None:
        t0 = t1 = ts[-1]
    if t0 is None:
        t0 = ts[0]
    if t0 == t1:
        return ref_orbit[t0]
    a = (gps_sod - t0) / (t1 - t0)
    return ref_orbit[t0] * (1 - a) + ref_orbit[t1] * a


def compute_initial_velocity(ref_orbit, gps_sod):
    ts = sorted(ref_orbit.keys())
    i = 0
    for j, t in enumerate(ts):
        if t >= gps_sod:
            i = j
            break
        i = j
    if i + 1 < len(ts):
        dt = ts[i + 1] - ts[i]
        if dt > 0.01:
            return (ref_orbit[ts[i + 1]] - ref_orbit[ts[i]]) / dt
    if i > 0:
        dt = ts[i] - ts[i - 1]
        if dt > 0.01:
            return (ref_orbit[ts[i]] - ref_orbit[ts[i - 1]]) / dt
    return np.zeros(3)


def get_sat_geometry(sp3, sv, utc_dt, rcv_pos):
    pos, clk, vel = _sp3_get(sp3, sv, utc_dt)
    if pos is None or abs(clk) > 0.1 * C_LIGHT:
        return None, None, None
    rho = float(np.linalg.norm(pos - rcv_pos))
    for _ in range(5):
        travel_time = rho / C_LIGHT
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
    sag = (OMEGA_E / C_LIGHT) * (pos[0] * rcv_pos[1] - pos[1] * rcv_pos[0])
    return pos, clk, rho + sag


def main():
    parser = argparse.ArgumentParser(description='Batch LSQ POD')
    parser.add_argument('--date', required=True, help='YYYY-MM-DD')
    parser.add_argument('--hours', type=float, default=2.0)
    parser.add_argument('--interval', type=float, default=30.0)
    parser.add_argument('--data-dir', default='./data')
    parser.add_argument('--grace-id', default='C')
    parser.add_argument('--max-iter', type=int, default=6)
    parser.add_argument('--gravity-nmax', type=int, default=90)
    parser.add_argument('--seg-duration', type=float, default=900.0,
                        help='Empirical RTN segment duration [s] (default: 900 = 15 min)')
    args = parser.parse_args()

    date_str = args.date
    y, m, d = [int(x) for x in date_str.split('-')]
    doy = (datetime(y, m, d) - datetime(y, 1, 1)).days + 1
    dp = Path(args.data_dir)
    grace_id = args.grace_id

    # -- Load GNV1B reference orbit --
    gnv_path = dp / 'gracefo' / str(y) / date_str / f'GNV1B_{date_str}_{grace_id}_04.txt'
    ref_orbit, ref_vel = load_gnv1b(str(gnv_path))
    print(f"[GNV1B] {len(ref_orbit)} pos, {len(ref_vel)} vel epochs")

    # -- Load GPS1B observations --
    rnx_path = dp / f'GPS1B_{date_str}_{grace_id}_04.rnx'
    pkl_path = dp / 'gracefo' / str(y) / date_str / f'GPS1B_{date_str}_{grace_id}_04.pkl'
    if rnx_path.exists():
        gps1b_raw = load_gps1b_rnx(str(rnx_path))
        print(f"[RINEX] {len(gps1b_raw)} epochs")
    elif pkl_path.exists():
        gps1b_raw = pickle.load(open(str(pkl_path), 'rb'))
        print(f"[GPS1B PKL] {len(gps1b_raw)} epochs")
    else:
        print(f"[FATAL] No GPS data found")
        return

    # -- Load SP3 --
    sp3_pkl = dp / str(y) / f'{doy:03d}' / 'igs_sp3_FIN.pkl'
    if not sp3_pkl.exists():
        sp3_pkl = dp / str(y) / f'{doy:03d}' / 'igs_sp3.pkl'
    sp3 = pickle.load(open(str(sp3_pkl), 'rb'))
    print(f"[SP3] {len(sp3['ts'])} epochs")

    # -- Select epochs --
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
    print(f"[EPOCH] {len(epochs)} selected (dt={args.interval}s)")

    # -- Pre-compute geometry at GNV1B reference --
    N_BIAS = min(60, len(epochs))
    print(f"\n-- Computing geometry ({len(epochs)} epochs) --")
    epoch_geo = {}
    sv_p_residuals = defaultdict(list)

    for gps_sod in epochs:
        utc_dt = J2000 + timedelta(seconds=gps_sod)
        ref_pos = interpolate_ref(ref_orbit, gps_sod)
        ep_data = []

        recs = gps1b_raw.get(gps_sod, {})
        for sv_id, rec in recs.items():
            if 'L_if' not in rec:
                continue
            sat_pos, sat_clk, rho_corr = get_sat_geometry(sp3, sv_id, utc_dt, ref_pos)
            if sat_pos is None:
                continue
            el = math.asin(abs(sat_pos[2] - ref_pos[2]) / rho_corr)
            if el < 0.087:
                continue

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

            if len(epoch_geo) < N_BIAS:
                sv_p_residuals[sv_id].append(P_r)

        if ep_data:
            epoch_geo[gps_sod] = ep_data

    # Estimate per-SV code biases
    sv_bias = {}
    for sv, p_vals in sv_p_residuals.items():
        if len(p_vals) >= 10:
            sv_bias[sv] = float(np.median(p_vals))
    print(f"  Per-SV code biases: {len(sv_bias)} SVs")

    # -- Load gravity model --
    GRAVITY_NMAX = args.gravity_nmax
    gravity_path = dp / 'gravity' / 'GGM05C.gfc'
    if gravity_path.exists():
        Cnm, Snm, _, GM_grav, R_grav = read_icgem_gfc(str(gravity_path))
        print(f"  Gravity: GGM05C Nmax={GRAVITY_NMAX}")
    else:
        print(f"  [FATAL] Gravity model not found: {gravity_path}")
        return

    # -- Initial state from GNV1B --
    r0_ecef = interpolate_ref(ref_orbit, epochs[0])
    v0_ecef = interpolate_ref(ref_vel, epochs[0])
    v0_fd = compute_initial_velocity(ref_orbit, epochs[0])
    mjd_start = MJD_J2000 + (epochs[0] - GPS_UTC_OFFSET) / SEC_PER_DAY
    print(f"  Initial: r0=[{r0_ecef[0]/1000:.1f},{r0_ecef[1]/1000:.1f},{r0_ecef[2]/1000:.1f}]km "
          f"|v0|={np.linalg.norm(v0_ecef)/1000:.3f}km/s")
    print(f"  File v0=[{v0_ecef[0]:.3f}, {v0_ecef[1]:.3f}, {v0_ecef[2]:.3f}] m/s")
    print(f"  FD   v0=[{v0_fd[0]:.3f}, {v0_fd[1]:.3f}, {v0_fd[2]:.3f}] m/s")
    print(f"  |dv| file-FD = {np.linalg.norm(v0_ecef - v0_fd):.4f} m/s")

    # -- Run batch LSQ --
    print(f"\n-- Batch LSQ ({args.hours}h arc, {args.interval}s interval, "
          f"seg={args.seg_duration:.0f}s) --")
    result = run_batch_lsq(
        epoch_geo, ref_orbit, sp3, sv_bias, epochs,
        r0_ecef, v0_ecef, mjd_start,
        Cnm, Snm, GRAVITY_NMAX, GM_grav, R_grav,
        Cd_init=2.2, CR_init=1.3,
        max_iter=args.max_iter,
        dt_integ=10.0,
        sigma_phase=0.01, sigma_code=0.30,
        param_names=['aR', 'aT', 'aN'],
        seg_duration=args.seg_duration,
    )

    if result is None:
        print("[FAIL] Batch LSQ failed")
        return

    # -- Compare with GNV1B --
    print(f"\n-- Results --")
    print(f"  Iterations: {result['iterations']}, Converged: {result['converged']}")
    print(f"  Post-fit RMS: {result['postfit_rms']:.4f}")
    print(f"  Cd: {result['Cd']:.4f}, CR: {result['CR']:.4f}")
    print(f"  aR={result['aR']:.3e} aT={result['aT']:.3e} aN={result['aN']:.3e} m/s^2 (mean over {result['n_segments']} segments)")
    if result['n_segments'] <= 4:
        print(f"    aR per seg: {[f'{v:.2e}' for v in result['aR_seg']]}")
        print(f"    aT per seg: {[f'{v:.2e}' for v in result['aT_seg']]}")
        print(f"    aN per seg: {[f'{v:.2e}' for v in result['aN_seg']]}")
    print(f"  r0 ECEF: [{result['r0'][0]:.3f}, {result['r0'][1]:.3f}, {result['r0'][2]:.3f}] m")
    print(f"  v0 ECEF: [{result['v0'][0]:.3f}, {result['v0'][1]:.3f}, {result['v0'][2]:.3f}] m/s")

    # Compute 3D RMS vs GNV1B
    gnv_r0 = interpolate_ref(ref_orbit, epochs[0])
    pos_diff = result['r0'] - gnv_r0
    print(f"  Position diff vs GNV1B at t0: |dr|={np.linalg.norm(pos_diff):.3f} m")

    # Integrate the estimated orbit and compare at each epoch
    print(f"\n-- Orbit comparison (estimated vs GNV1B) --")
    from src.orbit_dynamics import total_acc_eci
    from src.orbit_integrator import integrate_orbit_eci_with_stm

    def force_model(pos_eci, vel_eci, **kwargs):
        return total_acc_eci(pos_eci, vel_eci,
                             Cnm=Cnm, Snm=Snm, Nmax=GRAVITY_NMAX,
                             GM_gravity=GM_grav, R_gravity=R_grav,
                             **kwargs)

    # Segment-by-segment integration for comparison
    arc_dur = epochs[-1] - epochs[0]
    n_seg = result['n_segments']
    seg_dur = args.seg_duration

    r_cur, v_cur = result['r0_eci'].copy(), result['v0_eci'].copy()
    all_t, all_r = [], []

    for i_seg in range(n_seg):
        t_start = i_seg * seg_dur
        t_end = min((i_seg + 1) * seg_dur, arc_dur)
        if t_end <= t_start:
            break

        emp_seg = np.array([
            result['aR_seg'][i_seg] if i_seg < len(result['aR_seg']) else 0.0,
            result['aT_seg'][i_seg] if i_seg < len(result['aT_seg']) else 0.0,
            result['aN_seg'][i_seg] if i_seg < len(result['aN_seg']) else 0.0,
        ])

        mjd_seg = mjd_start + t_start / SEC_PER_DAY
        integ_seg = integrate_orbit_eci_with_stm(
            r_cur, v_cur, (0.0, t_end - t_start), force_model,
            Cd=result['Cd'], CR=result['CR'],
            area_drag=0.68, area_srp=3.4, mass=580.0,
            empirical_acc_rtn=emp_seg,
            param_names=['aR', 'aT', 'aN'],
            dt=10.0,
            mjd_tt=mjd_seg + 69.184 / SEC_PER_DAY,
            mjd_utc=mjd_seg,
            bodies=['Sun', 'Moon'],
        )

        start_idx = 0 if i_seg == 0 else 1
        for k in range(start_idx, len(integ_seg['t'])):
            all_t.append(t_start + integ_seg['t'][k])
            all_r.append(integ_seg['r'][k])

        r_cur = integ_seg['r'][-1].copy()
        v_cur = integ_seg['v'][-1].copy()

    pos_diffs = []
    for gps_sod in epochs:
        t_rel = gps_sod - epochs[0]
        # Find closest integration point
        i_step = int(round(t_rel / 10.0))
        if 0 <= i_step < len(all_t):
            r_eci_est = all_r[i_step]
            mjd_ep = mjd_start + t_rel / SEC_PER_DAY
            r_ecef_est, _ = eci_to_ecef(r_eci_est, np.zeros(3), mjd_ep)
            r_gnv = interpolate_ref(ref_orbit, gps_sod)
            pos_diffs.append(np.linalg.norm(r_ecef_est - r_gnv))

    if pos_diffs:
        pos_diffs = np.array(pos_diffs)
        rms_3d = np.sqrt(np.mean(pos_diffs**2))
        print(f"  3D RMS vs GNV1B: {rms_3d:.3f} m "
              f"(mean={np.mean(pos_diffs):.3f}, max={np.max(pos_diffs):.3f})")

    # Save results
    out_dir = Path('results') / 'batch_lsq'
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"batch_{date_str}_{grace_id}_{args.hours}h.pkl"
    pickle.dump({
        'result': result,
        'pos_diffs_vs_gnv': pos_diffs,
        'epochs': epochs,
        'config': vars(args),
    }, open(str(out_path), 'wb'))
    print(f"\n  Saved: {out_path}")


if __name__ == '__main__':
    main()
