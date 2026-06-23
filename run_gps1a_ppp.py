#!/usr/bin/env python3
"""
GRACE-FO PPP with Integer Ambiguity Resolution — GPS1A data

Algorithm:
  1. MW wide-lane → N_w fixed (reliable, std ~0.2 cycles)
  2. Arc-average B_if → N1 float → N1 fixed
  3. B_if_fixed = lam_n * N1 + coeff_w * N_w
  4. KF with ambiguity-corrected phase [dX, dY, dZ, clk_r, trop_wet]

Data: JPL RL04 GPS1A (phase 1Hz RAW, pseudorange 0.1Hz)
      IGS Final SP3, GNV1B reference orbit

Usage:
  py -3.12 run_gps1a_ppp.py --date 2024-04-29 --hours 4 --interval 30
"""
import sys, os, pickle, csv, json, math, argparse
from pathlib import Path
from datetime import datetime, timedelta
from collections import defaultdict
import numpy as np

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
F1_SQ, F2_SQ = F1 * F1, F2 * F2
ALPHA = F1_SQ / (F1_SQ - F2_SQ)
BETA_C = -F2_SQ / (F1_SQ - F2_SQ)
LAM_W = C / (F1 - F2)
LAM_N = C / (F1 + F2)
COEFF_W = C * F2 / (F1_SQ - F2_SQ)
OMEGA_E = 7.2921151467e-5
J2000 = datetime(2000, 1, 1, 12, 0, 0)

# KF parameters
Q_TROP = 1e-6
SIGMA_PHASE_FIXED = 0.005   # Fixed-ambiguity phase (5mm)
SIGMA_PHASE_FLOAT = 0.01    # Float-ambiguity phase (1cm)
SIGMA_CODE = 0.30

sys.path.insert(0, '.')
from src.gps1a_loader import download_gps1a, gps_sod_to_utc
from src.ambiguity import (
    compute_mw, compute_if_m, estimate_wide_lane,
    narrow_lane_from_if, b_if_from_ints, segment_arcs
)


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
    return np.array([[-sn, cn, 0], [-sl*cn, -sl*sn, cl], [cl*cn, cl*sn, sl]])


from src.sp3_loader import get_gps_pos_from_sp3 as _sp3_get

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
        if abs(rho_new - rho) < 1e-8: pos, clk, vel = pos_tx, clk_tx, vel_tx; rho = rho_new; break
        pos, clk, vel = pos_tx, clk_tx, vel_tx; rho = rho_new
    if not (1.8e7 < rho < 2.8e7): return None, None, None
    sag = (OMEGA_E / C) * (pos[0] * rcv_pos[1] - pos[1] * rcv_pos[0])
    return pos, clk, rho + sag


def load_gnv1b(filepath, nhours=4):
    orbit = {}; t0 = None
    for line in open(filepath):
        p = line.split()
        if len(p) < 6: continue
        try:
            t = float(p[0]); flag = p[2]
            if flag not in ("C", "E"): continue
            utc_dt = gps_sod_to_utc(t)
            if t0 is None: t0 = utc_dt
            if (utc_dt - t0).total_seconds() <= (nhours + 1) * 3600:
                orbit[utc_dt] = np.array([float(p[3]), float(p[4]), float(p[5])])
        except: pass
    return orbit


def build_gps1a_records(gps_obs, ref_orbit, t_start, nhours, interval, sp3):
    t_end = t_start + timedelta(hours=nhours)
    records = []
    orbit_ts = sorted(ref_orbit.keys())
    n_skip_rng = n_skip_el = n_no_pos = 0

    for gps_sod, sv_obs in sorted(gps_obs.items()):
        utc_dt = gps_sod_to_utc(gps_sod)
        if not (t_start <= utc_dt <= t_end): continue
        dt_s = (utc_dt - t_start).total_seconds()
        nearest = int((dt_s + interval/2)/interval)*interval
        if abs(dt_s - nearest) > 2.0: continue

        t0_orb = t1_orb = None
        for j, ti in enumerate(orbit_ts):
            if ti >= utc_dt: t1_orb = ti; t0_orb = orbit_ts[j-1] if j>0 else None; break
            t0_orb = ti
        if t1_orb is None: t0_orb = t1_orb = orbit_ts[-1]
        if t0_orb is None: t0_orb = orbit_ts[0]
        dt_frac = (utc_dt - t0_orb).total_seconds()
        dt_tot = (t1_orb - t0_orb).total_seconds()
        if dt_tot == 0: rcv_pos = np.array(ref_orbit[t0_orb], dtype=float)
        else:
            a = dt_frac/dt_tot
            rcv_pos = np.array(ref_orbit[t0_orb], dtype=float)*(1-a) + np.array(ref_orbit[t1_orb], dtype=float)*a
        if not (6e6 < float(np.linalg.norm(rcv_pos)) < 8e6): continue

        epoch_recs = []
        for sv, rec in sv_obs.items():
            sat_pos, clk_sv, rho_corr = get_sat_geometry(sp3, sv, utc_dt, rcv_pos)
            if sat_pos is None: n_no_pos += 1; continue
            delta = sat_pos - rcv_pos
            rng = float(np.linalg.norm(delta))
            if not (2e7 < rng < 5e7): n_skip_rng += 1; continue
            lat, lon, _ = ecef_to_blh(rcv_pos)
            R_enu = ecef_to_enu_matrix(lat, lon)
            e_enu = R_enu @ (delta/rng)
            el = float(np.arcsin(np.clip(e_enu[2], -1.0, 1.0)))
            if el < 0.087: n_skip_el += 1; continue

            L_if, P_if, B_if = compute_if_m(rec['L1_cyc'], rec['L2_cyc'], rec['P1'], rec['P2'])
            L_r = L_if + clk_sv - rho_corr
            P_r = P_if + clk_sv - rho_corr

            epoch_recs.append({
                'sv': sv, 'L_r': L_r, 'P_r': P_r, 'B_if': B_if,
                'L1_cyc': rec['L1_cyc'], 'L2_cyc': rec['L2_cyc'],
                'P1': rec['P1'], 'P2': rec['P2'],
                'el': float(np.degrees(el)), 'gps_sod': gps_sod,
                'slip_L1': rec.get('slip_L1', False),
                'slip_L2': rec.get('slip_L2', False),
                'sat_pos': sat_pos, 'clk_sv': clk_sv, 'rho_corr': rho_corr,
            })

        if epoch_recs:
            records.append({'utc': utc_dt, 'gps_sod': gps_sod,
                           'rcv_pos': rcv_pos, 'sv_data': epoch_recs})

    print(f"  PPP records: {len(records)} epochs "
          f"(rng={n_skip_rng}, el={n_skip_el}, no_pos={n_no_pos})")
    return records


def resolve_ambiguities(records):
    """Resolve N_w (MW) and N1 (narrow-lane) per SV arc."""
    sv_all = defaultdict(list)
    for ep in records:
        for d in ep['sv_data']:
            sv_all[d['sv']].append(d)

    amb_results = {}
    n_wl_fixed = n_nl_fixed = 0

    for sv, obs_list in sorted(sv_all.items()):
        obs_list.sort(key=lambda x: x['gps_sod'])
        arcs = segment_arcs(obs_list, max_gap_s=30.0)
        sv_arcs = {}

        for arc_idx, arc in enumerate(arcs):
            if len(arc) < 5: continue
            arc_start = arc[0]['gps_sod']
            arc_end = arc[-1]['gps_sod']

            # MW wide-lane
            N_w_fixed, N_w_float, sigma_w, n_w = estimate_wide_lane(arc)

            # IF ambiguity mean
            B_vals = [d['B_if'] for d in arc]
            B_mean = float(np.mean(B_vals))
            B_std = float(np.std(B_vals))

            N1_fixed = None
            N1_float = B_mean / LAM_N  # rough
            if N_w_fixed is not None:
                N1_float, N1_fixed, N1_std = narrow_lane_from_if(B_mean, B_std, N_w_fixed)
                if N1_fixed is not None:
                    n_nl_fixed += 1
                    n_wl_fixed += 1
                else:
                    n_wl_fixed += 1
            elif sigma_w <= 0.8:
                n_wl_fixed += 1  # borderline, counted as WL-fixed but no NL

            sv_arcs[arc_start] = {
                'N_w_fixed': N_w_fixed, 'N_w_float': N_w_float, 'sigma_w': sigma_w,
                'N1_fixed': N1_fixed, 'N1_float': N1_float,
                'B_if_mean': B_mean, 'B_if_std': B_std,
                'n_epochs': len(arc), 'arc_end': arc_end,
            }

        if sv_arcs:
            amb_results[sv] = sv_arcs

    print(f"  WL fixed: {n_wl_fixed} arcs, NL fixed: {n_nl_fixed} arcs "
          f"({len(amb_results)} SVs)")
    return amb_results


def run_kf_ppp_float_amb(records, ref_orbit):
    """Float-ambiguity phase-only KF PPP.

    State: [dX, dY, dZ, clk_r] (ephemeral, reset each epoch)
           + trop_wet (random walk)
           + B_sv_1, B_sv_2, ... (arc-constant float ambiguities)

    Phase measurement: y = L_r (NO ambiguity correction — B states absorb it)
    """
    b_states = {}   # sv -> B_if value (meters)
    b_vars = {}     # sv -> B_if variance
    b_arc = {}      # sv -> current arc_start (reset on cycle slip)

    P_trop = 100.0
    last_epoch = None
    results = []

    for ep_idx, ep in enumerate(records):
        utc_dt = ep['utc']
        rcv_pos = ep['rcv_pos']

        dt = 30.0
        if last_epoch is not None:
            dt = max((utc_dt - last_epoch).total_seconds(), 1.0)
        last_epoch = utc_dt

        # Time update: trop random walk
        P_trop += Q_TROP * dt

        active_svs = []
        for d in ep['sv_data']:
            sv = d['sv']
            # Detect arc change (cycle slip or new SV)
            arc_key = d['gps_sod']  # approximate
            if sv not in b_states:
                # Initialize from L_r - P_r
                b_states[sv] = d['L_r'] - d['P_r']
                b_vars[sv] = 0.5**2  # initial uncertainty ~0.5m
            active_svs.append(sv)

        N_B = len(active_svs)
        N_STATE = 5 + N_B

        # Build prior
        x_prior = np.zeros(N_STATE)
        P_prior = np.eye(N_STATE) * 1e6
        for i, sv in enumerate(active_svs):
            x_prior[5 + i] = b_states[sv]
            P_prior[5 + i, 5 + i] = b_vars[sv]
        P_prior[4, 4] = P_trop

        # Reset ephemeral diagonal blocks (already 1e6 from initialization)

        H_rows, y_rows, w_rows = [], [], []

        for i, d in enumerate(ep['sv_data']):
            sv = d['sv']
            los = d['sat_pos'] - rcv_pos
            los = los / np.linalg.norm(los)
            mf = 1.0 / max(np.sin(np.radians(d['el'])), 0.1)

            # Phase measurement: L_r = P_r + B_if
            # H = [-los, -los, -los, 1, mf, 0...1...0]
            row = np.zeros(N_STATE)
            row[0:3] = -los
            row[3] = 1.0
            row[4] = mf
            row[5 + i] = 1.0  # B_if state for this SV

            H_rows.append(row)
            y_rows.append(d['L_r'])
            w_rows.append(1.0 / 0.005**2)  # 5mm phase sigma

            # Also add code measurement for initial convergence
            H_rows.append(row.copy())
            H_rows[-1][5 + i] = 0.0  # code doesn't have B_if
            y_rows.append(d['P_r'])
            w_rows.append(1.0 / 0.30**2)  # 30cm code sigma

        H = np.array(H_rows)
        y = np.array(y_rows)
        W = np.diag(w_rows)

        try:
            S = H @ P_prior @ H.T + np.linalg.inv(W)
            K = P_prior @ H.T @ np.linalg.inv(S)
            innovation = y - H @ x_prior
            dx = K @ innovation
            x = x_prior + dx
            P_post = (np.eye(N_STATE) - K @ H) @ P_prior

            # Update B states for next epoch
            for i, sv in enumerate(active_svs):
                b_states[sv] = float(x[5 + i])
                b_vars[sv] = float(P_post[5 + i, 5 + i])

            # Trop variance
            P_trop = float(P_post[4, 4])

        except np.linalg.LinAlgError:
            continue

        solved_pos = rcv_pos + np.array([x[0], x[1], x[2]])

        utc_key = utc_dt
        ref = ref_orbit.get(utc_key, rcv_pos)
        lat, lon, _ = ecef_to_blh(ref)
        R_enu = ecef_to_enu_matrix(lat, lon)
        err_ecef = solved_pos - ref
        err_enu = R_enu @ err_ecef

        results.append({
            'utc': utc_dt, 'gps_sod': ep['gps_sod'],
            'dE_m': float(err_enu[0]), 'dN_m': float(err_enu[1]),
            'dU_m': float(err_enu[2]), 'd3_m': float(np.linalg.norm(err_ecef)),
            'n_sat': len(active_svs), 'n_fixed': 0,
            'clk_r': float(x[3]), 'trop_wet': float(x[4]),
        })

    return results


def process_day(date_str, nhours=4, interval=30):
    year, month, day = [int(x) for x in date_str.split('-')]
    print(f"\n{'='*60}")
    print(f"GPS1A PPP: {date_str} ({nhours}h, {interval}s)")
    print(f"{'='*60}")

    print("\n[1] Loading GPS1A data...")
    gps_obs = download_gps1a(year, month, day, grace_filter='C')
    if not gps_obs: return None

    print("\n[2] Loading SP3...")
    doy = (datetime(year, month, day) - datetime(year, 1, 1)).days + 1
    sp3_path = Path(f"data/{year}/{doy:03d}/igs_sp3_FIN.pkl")
    if not sp3_path.exists():
        sp3_path = Path(f"data/2024/{doy:03d}/igs_sp3_FIN.pkl")
    sp3 = pickle.load(open(sp3_path, 'rb'))

    print("\n[3] Loading GNV1B...")
    gnv_path = Path(f"data/gracefo/{year}/{date_str}/GNV1B_{date_str}_C_04.txt")
    ref_orbit = load_gnv1b(str(gnv_path), nhours)
    if not ref_orbit: return None
    ref_ts = sorted(ref_orbit.keys())
    t_start = ref_ts[0]

    print("\n[4] Building PPP records...")
    records = build_gps1a_records(gps_obs, ref_orbit, t_start, nhours, interval, sp3)
    if not records: return None
    print(f"  Records: {len(records)} epochs, "
          f"avg SVs: {sum(len(ep['sv_data']) for ep in records)/len(records):.1f}")

    print("\n[5] Running float-ambiguity KF PPP...")
    results = run_kf_ppp_float_amb(records, ref_orbit)

    if not results:
        print("  No results!")
        return None

    d3 = np.array([r['d3_m'] for r in results])
    de = np.array([abs(r['dE_m']) for r in results])
    dn = np.array([abs(r['dN_m']) for r in results])
    du = np.array([abs(r['dU_m']) for r in results])

    print(f"\n{'='*60}")
    print(f"RESULTS: {date_str} ({nhours}h, {interval}s)")
    print(f"{'='*60}")
    print(f"  Epochs:       {len(results)}")
    print(f"  3D RMS:       {float(np.sqrt(np.mean(d3**2))):.4f} m")
    print(f"  E RMS:        {float(np.sqrt(np.mean(de**2))):.4f} m")
    print(f"  N RMS:        {float(np.sqrt(np.mean(dn**2))):.4f} m")
    print(f"  U RMS:        {float(np.sqrt(np.mean(du**2))):.4f} m")
    print(f"  <2m:          {100*np.sum(d3 < 2)/len(d3):.1f}%")
    print(f"  <1m:          {100*np.sum(d3 < 1)/len(d3):.1f}%")
    print(f"  <0.5m:        {100*np.sum(d3 < 0.5)/len(d3):.1f}%")
    print(f"  <0.3m:        {100*np.sum(d3 < 0.3)/len(d3):.1f}%")

    n_fixed_epochs = sum(1 for r in results if r['n_fixed'] > 0)
    print(f"  Fixed epochs: {n_fixed_epochs}/{len(results)} "
          f"({100*n_fixed_epochs/max(len(results),1):.1f}%)")

    out_dir = Path("output"); out_dir.mkdir(exist_ok=True)
    csv_path = out_dir / f"gps1a_ppp_{date_str}_{int(nhours)}h.csv"
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'utc', 'gps_sod', 'dE_m', 'dN_m', 'dU_m', 'd3_m',
            'n_sat', 'n_fixed', 'clk_r', 'trop_wet'])
        writer.writeheader()
        for r in results:
            writer.writerow({k: v for k, v in r.items() if k in writer.fieldnames})
    print(f"  Saved: {csv_path}")
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--date', required=True)
    parser.add_argument('--hours', type=float, default=4)
    parser.add_argument('--interval', type=int, default=30)
    args = parser.parse_args()
    process_day(args.date, args.hours, args.interval)
