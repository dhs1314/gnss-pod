"""
Batch Least Squares solver for reduced-dynamic POD.

Estimates initial conditions [r0, v0] and force parameters (Cd) from
GNSS observations over an arc, using orbit integration with force models.

State vector: [rx, ry, rz, vx, vy, vz, Cd, clk_r, trop_wet, B_Gxx...]
"""
import numpy as np
from datetime import datetime, timedelta

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
F1_SQ, F2_SQ = F1 * F1, F2 * F2
ALPHA = F1_SQ / (F1_SQ - F2_SQ)
BETA = -F2_SQ / (F1_SQ - F2_SQ)
OMEGA_E = 7.2921151467e-5

from src.orbit_integrator import integrate_orbit_with_stm
from src.orbit_dynamics import total_acc


def solve_batch_lsq(epoch_data_list, ref_orbit, sp3, sv_bias, gps_sods,
                    r0_init, v0_init, Cd_init=2.2, max_iter=5,
                    sigma_phase=0.01, sigma_code=0.30):
    """Batch least squares for a short arc.

    Args:
        epoch_data_list: list of dicts, each with keys:
            'gps_sod', 'utc_dt', 'ep_data' (list of per-SV dicts)
        ref_orbit: dict {gps_sod: ECEF_position} for linearization reference
        sp3: SP3 data dict
        sv_bias: dict {sv: bias_meters}
        gps_sods: sorted list of GPS seconds-of-day
        r0_init: initial ECEF position (3,)
        v0_init: initial ECEF velocity (3,)
        Cd_init: initial drag coefficient
        max_iter: maximum iterations
        sigma_phase, sigma_code: measurement noise (m)

    Returns:
        dict with keys: 'r0', 'v0', 'Cd', 'converged', 'iterations',
                        'postfit_rms_phase', 'postfit_rms_code',
                        'state', 'cov', 'results'
    """
    N_EPOCH = len(gps_sods)
    if N_EPOCH < 3:
        return None

    # Determine which SVs appear in the arc
    all_svs = set()
    for ed in epoch_data_list:
        for d in ed['ep_data']:
            if d['sv'] in sv_bias:
                all_svs.add(d['sv'])
    sv_list = sorted(all_svs)
    N_SV = len(sv_list)
    sv_to_idx = {sv: i for i, sv in enumerate(sv_list)}

    # State: r0(3), v0(3), Cd, clk_r, trop_wet, B_sv(N_SV)
    N_STATE = 9 + N_SV
    # Parameter indices
    I_R0, I_V0 = slice(0, 3), slice(3, 6)
    I_CD = 6
    I_CLK = 7
    I_TROP = 8
    I_BSV = slice(9, 9 + N_SV)

    r0 = np.array(r0_init, dtype=float).copy()
    v0 = np.array(v0_init, dtype=float).copy()
    Cd = float(Cd_init)
    clk_r = 0.0
    trop_wet = 0.0
    B_sv = np.zeros(N_SV)

    t0 = gps_sods[0]
    t_end = gps_sods[-1]

    W_PHASE = 1.0 / sigma_phase**2
    W_CODE = 1.0 / sigma_code**2

    results = None

    for iteration in range(max_iter):
        # -- Integrate reference orbit --
        integ = integrate_orbit_with_stm(r0, v0, (0, t_end - t0), Cd=Cd, dt=10.0)
        t_integ = integ['t']
        r_integ = integ['r']
        v_integ = integ['v']
        phi_integ = integ['phi']
        s_cd_integ = integ['s_cd']

        # -- Build normal equations --
        N = np.zeros((N_STATE, N_STATE))
        b = np.zeros(N_STATE)

        n_obs = 0
        n_phase = 0
        n_code = 0
        total_res_sq = 0.0

        for i_ep, gps_sod in enumerate(gps_sods):
            ed = epoch_data_list[i_ep]
            if gps_sod != ed['gps_sod']:
                continue

            # Interpolate reference position at this epoch
            # Use GNV1B as linearization point (can be replaced later)
            ref_pos = _interp_ref(ref_orbit, gps_sod)
            if ref_pos is None:
                continue

            # Find integration step for this epoch
            t_rel = gps_sod - t0
            i_step = int(np.round(t_rel / 10.0))
            if i_step < 0 or i_step >= len(t_integ):
                continue

            r_ref = r_integ[i_step]
            v_ref = v_integ[i_step]
            phi_step = phi_integ[i_step]
            s_cd_step = s_cd_integ[i_step]

            # STM sub-blocks
            phi_r_r0 = phi_step[0:3, 0:3]
            phi_r_v0 = phi_step[0:3, 3:6]
            phi_v_r0 = phi_step[3:6, 0:3]
            phi_v_v0 = phi_step[3:6, 3:6]
            s_r_cd = s_cd_step[0:3]
            s_v_cd = s_cd_step[3:6]

            ep_data = ed['ep_data']

            for d in ep_data:
                sv = d['sv']
                if sv not in sv_to_idx or sv not in sv_bias:
                    continue

                i_sv = sv_to_idx[sv]
                sat_pos = d['sat_pos']
                sat_clk = d['sat_clk']
                rho_corr = d['rho_corr']
                el = d['el']
                e_vec = (sat_pos - ref_pos) / rho_corr
                mf = 1.0 / max(np.sin(el), 0.1)

                # -- Phase observation --
                if abs(d['L_r'] - sv_bias.get(sv, 0) - B_sv[i_sv]) < 200.0:
                    # Compute modeled range
                    rho_model = float(np.linalg.norm(sat_pos - r_ref))
                    sag = (OMEGA_E / C) * (sat_pos[0] * r_ref[1] - sat_pos[1] * r_ref[0])
                    rho_corr_model = rho_model + sag

                    # Observation residual (O - C)
                    obs = d['L_if_raw'] + sat_clk - rho_corr_model - sv_bias[sv]
                    # Model: -e·dr + clk_r + trop·mf + B_sv
                    # dr = r_ref + δr - r0 is approximated through STM

                    # Partial w.r.t. r0: ∂r/∂r0 = Φ_r_r0, ∂r/∂v0 = Φ_r_v0
                    # ∂rho/∂r = -e_vec (unit vector from receiver to satellite)
                    # ∂rho/∂r0 = ∂rho/∂r · ∂r/∂r0 = -e_vec · Φ_r_r0
                    h = np.zeros(N_STATE)
                    h[I_R0] = -e_vec @ phi_r_r0
                    h[I_V0] = -e_vec @ phi_r_v0
                    h[I_CD] = -e_vec @ s_r_cd
                    h[I_CLK] = 1.0
                    h[I_TROP] = mf
                    h[I_BSV][i_sv] = 1.0

                    N += W_PHASE * np.outer(h, h)
                    b += W_PHASE * obs * h
                    n_phase += 1
                    total_res_sq += W_PHASE * obs**2

                # -- Code observation --
                if abs(d['P_r'] - sv_bias.get(sv, 0)) < 100.0:
                    rho_model = float(np.linalg.norm(sat_pos - r_ref))
                    sag = (OMEGA_E / C) * (sat_pos[0] * r_ref[1] - sat_pos[1] * r_ref[0])
                    rho_corr_model = rho_model + sag

                    obs = d['P_if_raw'] + sat_clk - rho_corr_model - sv_bias[sv]

                    h = np.zeros(N_STATE)
                    h[I_R0] = -e_vec @ phi_r_r0
                    h[I_V0] = -e_vec @ phi_r_v0
                    h[I_CD] = -e_vec @ s_r_cd
                    h[I_CLK] = 1.0
                    h[I_TROP] = mf
                    # No B_sv for code

                    N += W_CODE * np.outer(h, h)
                    b += W_CODE * obs * h
                    n_code += 1
                    total_res_sq += W_CODE * obs**2

        if n_phase + n_code < N_STATE:
            print(f"    iter {iteration}: insufficient observations ({n_phase + n_code} < {N_STATE})")
            return None

        # -- Constrain weak parameters --
        # Cd prior: 2.2 ± 1.0
        N[I_CD, I_CD] += 1.0 / 1.0
        b[I_CD] += (2.2 - Cd) / 1.0

        # trop prior: 0 ± 2.0
        N[I_TROP, I_TROP] += 1.0 / 2.0
        b[I_TROP] += (0.0 - trop_wet) / 2.0

        # B_sv prior: 0 ± 10.0
        for i_sv in range(N_SV):
            idx = 9 + i_sv
            N[idx, idx] += 1.0 / 100.0  # weak prior
            # No b prior (zero-mean)

        # -- Solve --
        try:
            dx = np.linalg.solve(N, b)
        except np.linalg.LinAlgError:
            print(f"    iter {iteration}: singular normal matrix")
            return None

        # -- Update state --
        r0 += dx[I_R0]
        v0 += dx[I_V0]
        Cd += dx[I_CD]
        clk_r += dx[I_CLK]
        trop_wet += dx[I_TROP]
        B_sv += dx[I_BSV]

        # -- Convergence check --
        pos_change = float(np.linalg.norm(dx[I_R0]))
        vel_change = float(np.linalg.norm(dx[I_V0]))
        rms_res = np.sqrt(total_res_sq / max(n_phase + n_code, 1))
        print(f"    iter {iteration}: dr={pos_change:.3f}m dv={vel_change:.4f}m/s "
              f"dCd={dx[I_CD]:.4f} obs={n_phase}+{n_code} rms_res={rms_res:.3f}")

        if pos_change < 0.001 and vel_change < 0.0001:
            print(f"    Converged after {iteration + 1} iterations")
            try:
                cov = np.linalg.inv(N)
            except np.linalg.LinAlgError:
                cov = np.zeros((N_STATE, N_STATE))
            results = {
                'r0': r0.copy(), 'v0': v0.copy(), 'Cd': Cd,
                'clk_r': clk_r, 'trop_wet': trop_wet,
                'B_sv': {sv: B_sv[i] for sv, i in sv_to_idx.items()},
                'converged': True, 'iterations': iteration + 1,
                'postfit_rms': rms_res,
                'n_phase': n_phase, 'n_code': n_code,
                'state': np.concatenate([r0, v0, [Cd, clk_r, trop_wet], B_sv]),
                'cov': cov,
                'sv_list': sv_list,
            }
            break

    if results is None:
        print(f"    Did not converge after {max_iter} iterations")
        try:
            cov = np.linalg.inv(N)
        except np.linalg.LinAlgError:
            cov = np.zeros((N_STATE, N_STATE))
        results = {
            'r0': r0.copy(), 'v0': v0.copy(), 'Cd': Cd,
            'clk_r': clk_r, 'trop_wet': trop_wet,
            'B_sv': {sv: B_sv[i] for sv, i in sv_to_idx.items()},
            'converged': False, 'iterations': max_iter,
            'postfit_rms': rms_res,
            'n_phase': n_phase, 'n_code': n_code,
            'state': np.concatenate([r0, v0, [Cd, clk_r, trop_wet], B_sv]),
            'cov': cov,
            'sv_list': sv_list,
        }

    return results


def _interp_ref(ref_orbit, gps_sod):
    """Interpolate reference position at given GPS second of day."""
    ts = sorted(ref_orbit.keys())
    if not ts:
        return None

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
