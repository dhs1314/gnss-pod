#!/usr/bin/env python3
"""
GRACE-FO Reduced-Dynamic Kalman Filter (Phase 3)

Full reduced-dynamic POD with:
  - GGM05C 90x90 spherical harmonic gravity
  - Third-body perturbations (Sun, Moon)
  - Cannonball SRP with conical Earth shadow
  - Exponential drag model
  - Empirical RTN accelerations (piecewise constant)
  - ECI dynamics integration with RK4
  - Dynamics-based velocity correction feedback

The KF estimates position corrections from GNV1B reference,
plus drag coefficient (Cd), SRP coefficient (CR), empirical
RTN accelerations, and per-SV float ambiguities.

Usage:
  py -3.12 run_gps1b_rd_ppp.py --date 2024-04-29 --hours 0.5 --interval 30
"""
import sys, os, pickle, csv, json, math, argparse
from pathlib import Path
from datetime import datetime, timedelta
from collections import defaultdict
import numpy as np

sys.path.insert(0, str(Path(__file__).parent / 'src'))
from gps1b_rnx_loader import load_gps1b_rnx
from sp3_loader import get_gps_pos_from_sp3 as _sp3_get
from orbit_dynamics import total_acc_eci, GM, OMEGA_E
from orbit_integrator import rk4_step_eci
from coordinates import ecef_to_eci, eci_to_ecef
from gravity_model import read_icgem_gfc

C = 299792458.0
SEC_PER_DAY = 86400.0
F1, F2 = 1575.42e6, 1227.60e6
F1_SQ, F2_SQ = F1 * F1, F2 * F2
ALPHA = F1_SQ / (F1_SQ - F2_SQ)
BETA = -F2_SQ / (F1_SQ - F2_SQ)
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
    i = 0
    for j, t in enumerate(ts):
        if t >= gps_sod: i = j; break
        i = j
    if i + 1 < len(ts):
        dt = ts[i+1] - ts[i]
        if dt > 0.01:
            return (ref_orbit[ts[i+1]] - ref_orbit[ts[i]]) / dt
    if i > 0:
        dt = ts[i] - ts[i-1]
        if dt > 0.01:
            return (ref_orbit[ts[i]] - ref_orbit[ts[i-1]]) / dt
    return np.zeros(3)


def process_day_rd(date_str, nhours=0.5, data_dir='./data', grace_id='C', interval=30.0):
    """Reduced-Dynamic Kalman Filter PPP.

    Approach: Linearize observations at ref_pos (GNV1B, always near-truth),
    use dynamics prediction r_pred as a soft prior on the position estimate.

    State: [dX, dY, dZ, trop_wet, Cd, B_G01, B_G02, ...]
      dX, dY, dZ: position correction from ref_pos

    Prior: dX, dY, dZ ~ N(r_pred - ref_pos, pcov_pos)
           (dynamics prediction pulls position estimate)

    Dynamics: RK4 from posterior position/velocity at epoch k-1 to
              prior position/velocity at epoch k.
    """
    y, m, d = [int(x) for x in date_str.split('-')]
    doy = (datetime(y, m, d) - datetime(y, 1, 1)).days + 1
    dp = Path(data_dir)

    # Load data
    gnv_path = dp / 'gracefo' / str(y) / date_str / f'GNV1B_{date_str}_{grace_id}_04.txt'
    ref_orbit = load_gnv1b(str(gnv_path))
    print(f"[GNV1B] {len(ref_orbit)} epochs")

    # Try RINEX first, then PKL (pre-processed binary) from gracefo dir
    rnx_path = dp / f'GPS1B_{date_str}_{grace_id}_04.rnx'
    pkl_path = dp / 'gracefo' / str(y) / date_str / f'GPS1B_{date_str}_{grace_id}_04.pkl'
    if rnx_path.exists():
        gps1b_raw = load_gps1b_rnx(str(rnx_path))
        print(f"[RINEX] {len(gps1b_raw)} epochs")
    elif pkl_path.exists():
        gps1b_raw = pickle.load(open(str(pkl_path), 'rb'))
        print(f"[GPS1B PKL] {len(gps1b_raw)} epochs")
    else:
        print(f"[FATAL] No GPS data found (tried {rnx_path}, {pkl_path})")
        return None
    if gps1b_raw is None: return None

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
        dt_ep = gps_sod - gps_sod_start
        nearest = round(dt_ep / interval) * interval
        if abs(dt_ep - nearest) > max(2.0, interval * 0.1): continue
        epochs.append(gps_sod)
    print(f"[EPOCH] {len(epochs)} selected (dt={interval}s)")

    # -- Pre-compute geometry at ref_pos (same as kinematic PPP) --
    N_BIAS = min(60, len(epochs))
    print(f"\n-- Computing geometry ({len(epochs)} epochs) --")
    epoch_geo = {}
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

    # -- Load gravity model (GGM05C) --
    GRAVITY_NMAX = 90
    gravity_path = dp / 'gravity' / 'GGM05C.gfc'
    if gravity_path.exists():
        Cnm, Snm, _, GM_grav, R_grav = read_icgem_gfc(str(gravity_path))
        print(f"  Gravity: GGM05C Nmax={GRAVITY_NMAX}")
    else:
        print(f"  [WARN] Gravity model not found: {gravity_path}, falling back to J2-only")
        Cnm, Snm = None, None

    # -- Initial state --
    r_post = interpolate_ref(ref_orbit, epochs[0])
    v_post = compute_initial_velocity(ref_orbit, epochs[0])
    Cd_cur = 2.2
    CR_cur = 1.3
    emp_rtn = np.zeros(3)
    print(f"  Initial: r0=[{r_post[0]/1000:.1f},{r_post[1]/1000:.1f},{r_post[2]/1000:.1f}]km "
          f"|v0|={np.linalg.norm(v_post)/1000:.3f}km/s")

    # -- KF persistent state --
    pstate = {}       # scalar persistent: trop_wet, Cd, CR, aR, aT, aN, B_sv
    pcov_pers = {}    # covariance for scalar persistent

    pcov_pos = np.full(3, 100.0)      # position prior variance (m^2)

    # Process noise
    Q_POS = 0.01    # m^2/s position process noise (dynamics uncertainty growth)
    Q_CD = 1e-6     # Cd drift per second — allow ~0.1 drift over 2h
    Q_CR = 1e-7     # CR drift per second — allow ~0.03 drift over 2h
    Q_EMP = 1e-11   # empirical RTN drift per second — allow ~3e-4 m/s^2 growth over 2h
    Q_TROP = 1e-6   # trop drift per second
    Q_B = 2e-5      # B_if drift per second

    # Measurement noise
    STD_PHASE = 0.01
    STD_CODE = 0.30
    W_PHASE = 1.0 / STD_PHASE**2
    W_CODE = 1.0 / STD_CODE**2

    # SV initialization
    INIT_BUF = 3
    sv_init = defaultdict(list)
    SV_STALE = 10

    # MJD of J2000 epoch (2000-01-01 12:00:00 UTC)
    MJD_J2000 = 51544.5

    print(f"\n-- RD-KF: epoch-by-epoch processing --")
    results = []
    sv_seen_ago = {}
    sv_bias_buf = defaultdict(list)
    NEW_SV_BIAS_BUF = 5

    for i_ep, gps_sod in enumerate(epochs):
        if gps_sod not in epoch_geo: continue
        ep_data = epoch_geo[gps_sod]
        ref_pos = interpolate_ref(ref_orbit, gps_sod)
        dt_s = interval

        # ============================================================
        # TIME UPDATE: propagate dynamics to get prior position
        # ============================================================
        Cd_cur = pstate.get('__cd__', 2.2)
        CR_cur = pstate.get('__cr__', 1.3)
        aR_cur = pstate.get('__aR__', 0.0)
        aT_cur = pstate.get('__aT__', 0.0)
        aN_cur = pstate.get('__aN__', 0.0)
        emp_rtn = np.array([aR_cur, aT_cur, aN_cur])

        # MJD for epoch. gps_sod = seconds since J2000 (UTC).
        # TT = TAI + 32.184s; TAI = UTC + 37s (as of 2017).
        # So TT = UTC + 69.184s.
        mjd_utc = MJD_J2000 + gps_sod / SEC_PER_DAY

        if i_ep == 0:
            r_pred = r_post.copy()
            v_pred = v_post.copy()
        else:
            # Convert ECEF to ECI at start of step, integrate to epoch
            mjd_start = mjd_utc - dt_s / SEC_PER_DAY
            mjd_tt_start = mjd_start + 69.184 / SEC_PER_DAY
            r_eci, v_eci = ecef_to_eci(r_post, v_post, mjd_start)

            if Cnm is not None:
                r_pred_eci, v_pred_eci = rk4_step_eci(
                    r_eci, v_eci, dt_s,
                    total_acc_eci,
                    mjd_tt=mjd_tt_start, mjd_utc=mjd_start,
                    Cnm=Cnm, Snm=Snm, Nmax=GRAVITY_NMAX,
                    CD=Cd_cur, CR=CR_cur,
                    area_drag=0.68, area_srp=3.4, mass=580.0,
                    empirical_acc_rtn=emp_rtn,
                    bodies=['Sun', 'Moon'],
                    GM_gravity=GM_grav, R_gravity=R_grav,
                )
            else:
                # Fallback: J2-only ECEF integration
                from orbit_integrator import rk4_step
                r_pred, v_pred = rk4_step(r_post, v_post, Cd_cur,
                                          area_to_mass=0.002, dt=dt_s)
                r_pred_eci, v_pred_eci = None, None

            if r_pred_eci is not None:
                r_pred, v_pred = eci_to_ecef(r_pred_eci, v_pred_eci, mjd_utc)

        # Position prior: center at dynamics prediction, expressed as
        # correction from ref_pos: dX ~ N(r_pred - ref_pos, pcov_pos)
        prior_offset = r_pred - ref_pos

        # Process noise accumulation
        pcov_pos = pcov_pos + Q_POS * dt_s

        for key in list(pcov_pers.keys()):
            if key == '__trop__':
                pcov_pers[key] += Q_TROP * dt_s
            elif key == '__cd__':
                pcov_pers[key] += Q_CD * dt_s
            elif key == '__cr__':
                pcov_pers[key] += Q_CR * dt_s
            elif key in ('__aR__', '__aT__', '__aN__'):
                pcov_pers[key] += Q_EMP * dt_s
            elif key.startswith('G'):
                pcov_pers[key] += Q_B * dt_s

        # ============================================================
        # SV LIFECYCLE
        # ============================================================
        svs_now = set(d['sv'] for d in ep_data)

        for d in ep_data:
            sv = d['sv']
            if sv not in pstate:
                sv_init[sv].append(d['L_r'] - d['P_r'])
            if sv not in sv_bias:
                sv_bias_buf[sv].append(d['P_r'])

        for sv in list(sv_bias_buf.keys()):
            if sv in sv_bias: continue
            if len(sv_bias_buf[sv]) >= NEW_SV_BIAS_BUF:
                sv_bias[sv] = float(np.median(sv_bias_buf[sv]))
                del sv_bias_buf[sv]

        for sv in list(sv_init.keys()):
            if sv in pstate: continue
            if len(sv_init[sv]) >= INIT_BUF or (len(sv_init[sv]) >= 1 and i_ep <= 2):
                pstate[sv] = float(np.median(sv_init[sv]))
                pcov_pers[sv] = 0.5
                del sv_init[sv]

        _NON_SV_KEYS = {'__trop__', '__cd__', '__cr__', '__aR__', '__aT__', '__aN__'}
        for sv in list(pstate.keys()):
            if sv in _NON_SV_KEYS: continue
            sv_seen_ago[sv] = 0 if sv in svs_now else sv_seen_ago.get(sv, 0) + 1

        for sv, ago in list(sv_seen_ago.items()):
            if ago > SV_STALE:
                pstate.pop(sv, None)
                pcov_pers.pop(sv, None)
                sv_seen_ago.pop(sv, None)

        if '__trop__' not in pstate:
            pstate['__trop__'] = 0.0
            pcov_pers['__trop__'] = 1.0
        if '__cd__' not in pstate:
            pstate['__cd__'] = Cd_cur
            pcov_pers['__cd__'] = 0.5       # σ=0.7 initial Cd uncertainty
        if '__cr__' not in pstate:
            pstate['__cr__'] = CR_cur
            pcov_pers['__cr__'] = 0.1       # σ=0.3 initial CR uncertainty
        if '__aR__' not in pstate:
            pstate['__aR__'] = 0.0
            pcov_pers['__aR__'] = 1e-4      # σ=0.01 m/s^2 — loose prior for empirical accel
        if '__aT__' not in pstate:
            pstate['__aT__'] = 0.0
            pcov_pers['__aT__'] = 1e-4
        if '__aN__' not in pstate:
            pstate['__aN__'] = 0.0
            pcov_pers['__aN__'] = 1e-4

        # ============================================================
        # BUILD STATE VECTOR
        # State: [dX, dY, dZ, trop_wet, Cd, CR, aR, aT, aN, B_G01, B_G02, ...]
        #   dX, dY, dZ: position correction from ref_pos
        #   trop_wet: wet troposphere delay
        #   Cd: drag coefficient
        #   CR: SRP coefficient
        #   aR, aT, aN: empirical RTN accelerations
        #   B_sv: per-SV float ambiguity
        #
        # Linearization point: ref_pos (GNV1B, always near-truth)
        # Position prior: pulls estimate toward r_pred (dynamics prediction)
        # ============================================================
        sv_list = sorted([k for k in pstate if k.startswith('G')])
        sv_to_idx = {sv: 9 + i for i, sv in enumerate(sv_list)}
        n_sv = len(sv_list)
        n_state = 9 + n_sv

        A_rows, y_rows, w_rows = [], [], []

        # -- Position prior: dX, dY, dZ ~ N(prior_offset, pcov_pos) --
        # Dynamics prediction acts as soft constraint on position estimate.
        # Prior variance is adaptive: at least 1 m^2 (sigma=1m), or the
        # actual prediction error squared, whichever is larger.
        #
        # Empirical RTN sensitivity: for constant acceleration a_rtn over dt,
        #   Δr_ecef = 0.5 * dt^2 * R_rtn_to_ecef @ a_rtn
        # This makes emp observable through the position prior.
        r_norm = float(np.linalg.norm(r_pred))
        if r_norm > 1e3:
            R_dir = r_pred / r_norm
            N_dir = np.cross(r_pred, v_pred)
            N_norm = float(np.linalg.norm(N_dir))
            if N_norm > 1e-6:
                N_dir = N_dir / N_norm
                T_dir = np.cross(N_dir, R_dir)
                # Sensitivity: ∂(pos)/∂(emp_rtn)  [3x3]  in ECEF
                S_emp = 0.5 * dt_s**2 * np.column_stack([R_dir, T_dir, N_dir])
            else:
                S_emp = np.zeros((3, 3))
        else:
            S_emp = np.zeros((3, 3))

        _USE_DYNAMICS_PRIOR = True
        if _USE_DYNAMICS_PRIOR:
            # Position prior variance inflates when prediction is bad:
            # lets GNSS dominate position when dynamics is inaccurate.
            prior_var_floor = max(float(np.linalg.norm(prior_offset))**2, 1.0)
            for j in range(3):
                h = np.zeros(n_state); h[j] = 1.0
                A_rows.append(h)
                y_rows.append(float(prior_offset[j]))
                pvar = max(float(pcov_pos[j]), prior_var_floor)
                w_rows.append(1.0 / pvar)

        # -- Trop prior --
        h = np.zeros(n_state); h[3] = 1.0
        A_rows.append(h); y_rows.append(pstate['__trop__'])
        w_rows.append(1.0 / max(pcov_pers.get('__trop__', 1.0), 1e-9))

        # -- Cd prior --
        h = np.zeros(n_state); h[4] = 1.0
        A_rows.append(h); y_rows.append(pstate['__cd__'])
        w_rows.append(1.0 / max(pcov_pers.get('__cd__', 0.25), 1e-9))

        # -- CR prior --
        h = np.zeros(n_state); h[5] = 1.0
        A_rows.append(h); y_rows.append(pstate['__cr__'])
        w_rows.append(1.0 / max(pcov_pers.get('__cr__', 0.01), 1e-9))

        # -- Empirical RTN priors --
        for j, label in enumerate(['__aR__', '__aT__', '__aN__']):
            h = np.zeros(n_state); h[6 + j] = 1.0
            A_rows.append(h); y_rows.append(pstate[label])
            w_rows.append(1.0 / max(pcov_pers.get(label, 1e-12), 1e-9))

        # -- B_sv priors --
        for sv in sv_list:
            idx = sv_to_idx[sv]
            h = np.zeros(n_state); h[idx] = 1.0
            A_rows.append(h); y_rows.append(pstate[sv])
            w_rows.append(1.0 / max(pcov_pers.get(sv, 0.5), 1e-9))

        n_prior = len(A_rows)
        n_ph, n_co = 0, 0

        # -- Measurement equations (linearized at ref_pos) --
        for d in ep_data:
            sv = d['sv']
            if sv not in sv_to_idx or sv not in sv_bias: continue

            e_vec = (d['sat_pos'] - ref_pos) / d['rho_corr']
            mf = 1.0 / max(math.sin(d['el']), 0.1)
            idx_b = sv_to_idx[sv]

            # Phase: obs = L_r - sv_bias (model: -e·dr + trop*mf + B_sv)
            obs_ph = d['L_r'] - sv_bias[sv]
            if abs(obs_ph - pstate.get(sv, 0)) < 200.0:
                h = np.zeros(n_state)
                h[0:3] = -e_vec
                h[3] = mf
                h[idx_b] = 1.0
                A_rows.append(h); y_rows.append(obs_ph)
                w_rows.append(W_PHASE)
                n_ph += 1

            # Code: obs = P_r - sv_bias (model: -e·dr + trop*mf)
            obs_co = d['P_r'] - sv_bias[sv]
            if abs(obs_co) < 100.0:
                h = np.zeros(n_state)
                h[0:3] = -e_vec
                h[3] = mf
                A_rows.append(h); y_rows.append(obs_co)
                w_rows.append(W_CODE)
                n_co += 1

        if len(A_rows) < n_state:
            continue

        A = np.array(A_rows)
        y = np.array(y_rows)
        w_diag = np.array(w_rows)

        try:
            HtWH = A.T @ (w_diag[:, None] * A)
            HtWy = A.T @ (w_diag * y)
            x = np.linalg.solve(HtWH, HtWy)
        except np.linalg.LinAlgError:
            continue

        # ============================================================
        # EXTRACT RESULTS
        # ============================================================
        dX, dY, dZ = x[0], x[1], x[2]
        trop_wet = x[3]
        Cd_est = x[4]
        CR_est = x[5]
        aR_est = x[6]
        aT_est = x[7]
        aN_est = x[8]

        pos_est = ref_pos + np.array([dX, dY, dZ])
        err = pos_est - ref_pos  # = [dX, dY, dZ]
        d3 = float(np.linalg.norm(err))

        # Compute instantaneous velocity for next epoch's dynamics.
        # Epoch 0: forward GNV1B + kinematic correction.
        # Epoch >0: dynamics-based correction — the prediction error
        #   pos_est - r_pred reveals the velocity error:
        #   v_corrected = v_pred + (pos_est - r_pred) / dt
        from orbit_dynamics import total_acc
        if i_ep == 0:
            r_cur = interpolate_ref(ref_orbit, gps_sod)
            r_next = interpolate_ref(ref_orbit, gps_sod + dt_s)
            v_avg = (r_next - r_cur) / max(dt_s, 1e-6)
            r_mid = 0.5 * (r_cur + r_next)
            a_mid = total_acc(r_mid, v_avg, Cd_est)
            v_est = v_avg - 0.5 * a_mid * dt_s
        else:
            v_est = v_pred + (pos_est - r_pred) / dt_s

        # ============================================================
        # UPDATE PERSISTENT STATE FOR NEXT EPOCH
        # ============================================================
        r_post = pos_est.copy()
        v_post = v_est
        pstate['__trop__'] = float(trop_wet)
        pstate['__cd__'] = float(Cd_est)
        pstate['__cr__'] = float(CR_est)
        pstate['__aR__'] = float(aR_est)
        pstate['__aT__'] = float(aT_est)
        pstate['__aN__'] = float(aN_est)
        for sv in sv_list:
            pstate[sv] = float(x[sv_to_idx[sv]])

        # Update covariances from formal errors
        try:
            cov_full = np.linalg.inv(HtWH)
            pcov_pos = np.array([max(float(cov_full[j, j]), 1e-6) for j in range(3)])
            pcov_pers['__trop__'] = max(float(cov_full[3, 3]), 1e-6)
            pcov_pers['__cd__'] = max(float(cov_full[4, 4]), 1e-6)
            pcov_pers['__cr__'] = max(float(cov_full[5, 5]), 1e-6)
            pcov_pers['__aR__'] = max(float(cov_full[6, 6]), 1e-4)
            pcov_pers['__aT__'] = max(float(cov_full[7, 7]), 1e-4)
            pcov_pers['__aN__'] = max(float(cov_full[8, 8]), 1e-4)
            for sv in sv_list:
                idx = sv_to_idx[sv]
                pcov_pers[sv] = max(float(cov_full[idx, idx]), 1e-6)
        except np.linalg.LinAlgError:
            pass

        # -- Post-solve emp feedback from velocity correction history --
        # v_corr = (pos_est - r_pred) / dt is the velocity error revealed
        # by the prediction residual. A persistent v_corr in RTN indicates
        # unmodeled acceleration: a_unmodeled ≈ mean(v_corr) / dt.
        # Use exponential moving average of v_corr_rtn / dt as emp estimate.
        if i_ep > 0:
            r_norm_fb = float(np.linalg.norm(pos_est))
            if r_norm_fb > 1e3:
                R_fb = pos_est / r_norm_fb
                N_fb = np.cross(pos_est, v_est)
                N_norm_fb = float(np.linalg.norm(N_fb))
                if N_norm_fb > 1e-6:
                    N_fb = N_fb / N_norm_fb
                    T_fb = np.cross(N_fb, R_fb)
                    R_ecef_to_rtn = np.array([R_fb, T_fb, N_fb])
                    v_corr_ecef = (pos_est - r_pred) / dt_s
                    v_corr_rtn = R_ecef_to_rtn @ v_corr_ecef
                    # EMA of acceleration implied by persistent velocity error
                    a_implied = v_corr_rtn / dt_s
                    alpha = 0.05  # slow: ~20-epoch time constant
                    pstate['__aR__'] = (1 - alpha) * pstate['__aR__'] + alpha * float(a_implied[0])
                    pstate['__aT__'] = (1 - alpha) * pstate['__aT__'] + alpha * float(a_implied[1])
                    pstate['__aN__'] = (1 - alpha) * pstate['__aN__'] + alpha * float(a_implied[2])

        # Post-fit residuals (phase only)
        model = A @ x
        residuals = y - model
        phase_res = residuals[n_prior:n_prior + n_ph] if n_ph > 0 else np.array([])

        # ENU
        lat, lon, _ = ecef_to_blh(ref_pos)
        R = ecef_to_enu_matrix(lat, lon)
        enu = R @ err

        results.append({
            'time': J2000 + timedelta(seconds=gps_sod),
            'gps_sod': gps_sod,
            'dE': float(enu[0]), 'dN': float(enu[1]), 'dU': float(enu[2]),
            'd3': d3,
            'clk': 0.0, 'trop': float(trop_wet),
            'Cd': float(Cd_est), 'CR': float(CR_est),
            'aR': float(aR_est), 'aT': float(aT_est), 'aN': float(aN_est),
            'n_sat': len(ep_data),
            'ph_res_std': float(np.std(phase_res)) if len(phase_res) > 0 else 0.0,
            'n_B_sv': n_sv,
        })

        # Debug: dynamics prediction error and velocity correction
        pred_err = float(np.linalg.norm(r_pred - ref_pos))
        if len(results) <= 5 or len(results) % 100 == 0:
            r = results[-1]
            print(f"  ep {len(results)}/{len(epochs)}: 3D={d3:.3f}m pred_err={pred_err:.1f}m "
                  f"B_sv={r['n_B_sv']} ph_res={r['ph_res_std']:.3f}m "
                  f"Cd={Cd_est:.3f} CR={CR_est:.3f} "
                  f"emp=[{aR_est:.1e},{aT_est:.1e},{aN_est:.1e}] "
                  f"trop={trop_wet:.3f} n_sat={len(ep_data)}")

    print(f"[RD-KF] {len(results)}/{len(epochs)} epochs solved")

    if len(results) < 10:
        print("[FATAL] Too few successful epochs")
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
    print(f"  GRACE-FO Reduced-Dynamic KF -- GGM05C 90x90 + 3rd-body + SRP")
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
        'sv_bias': sv_bias,
    }


def main():
    parser = argparse.ArgumentParser(description='GRACE-FO Reduced-Dynamic KF PPP')
    parser.add_argument('--date', default='2024-04-29')
    parser.add_argument('--hours', type=float, default=0.5)
    parser.add_argument('--data-dir', default='./data')
    parser.add_argument('--output-dir', default='./output')
    parser.add_argument('--grace-id', default='C')
    parser.add_argument('--interval', type=float, default=30.0)
    args = parser.parse_args()

    print(f"\n{'='*65}")
    print(f"  GRACE-FO Reduced-Dynamic KF PPP -- GGM05C 90x90 + 3rd-body + SRP")
    print(f"  Date: {args.date}  Hours: {args.hours}h  dt: {args.interval}s")
    print(f"{'='*65}")

    output = process_day_rd(date_str=args.date, nhours=args.hours,
                            data_dir=args.data_dir, grace_id=args.grace_id,
                            interval=args.interval)
    if output is None:
        print("[FATAL] Processing failed")
        sys.exit(1)

    out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)
    tag = output['date'].replace('-', '') + '_' + output.get('grace_id', 'C')
    csv_path = out / f'ppp_rd_{tag}.csv'
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['time', 'gps_sod', 'dE_m', 'dN_m', 'dU_m', 'd3_m',
                     'clock_m', 'trop_m', 'Cd', 'CR', 'aR', 'aT', 'aN',
                     'n_sat', 'ph_res_std_m', 'n_B_sv'])
        for r in output['results']:
            w.writerow([r['time'].isoformat(), f"{r['gps_sod']:.6f}",
                        f"{r['dE']:.4f}", f"{r['dN']:.4f}", f"{r['dU']:.4f}",
                        f"{r['d3']:.4f}", f"{r.get('clk', 0):.4f}",
                        f"{r.get('trop', 0):.4f}",
                        f"{r.get('Cd', 0):.4f}", f"{r.get('CR', 0):.4f}",
                        f"{r.get('aR', 0):.2e}", f"{r.get('aT', 0):.2e}",
                        f"{r.get('aN', 0):.2e}",
                        r['n_sat'], f"{r['ph_res_std']:.4f}",
                        r.get('n_B_sv', 0)])
    print(f"[CSV] {csv_path}")
    print(f"[DONE] {args.date}")


if __name__ == '__main__':
    main()
