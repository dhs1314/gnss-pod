"""Batch least-squares estimator for reduced-dynamic POD.

Integrates a full arc in ECI with STM + parameter sensitivity propagation,
builds the design matrix by mapping observation partials back to the initial
state via the STM, accumulates normal equations with clock pre-elimination,
and solves for initial state + force parameters.

Reference: Montenbruck & Gill, "Satellite Orbits", Chapter 8.
"""
import numpy as np
from src.orbit_integrator import integrate_orbit_eci_with_stm
from src.coordinates import eci_to_ecef, ecef_to_eci

C_LIGHT = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
F1_SQ, F2_SQ = F1 * F1, F2 * F2
ALPHA = F1_SQ / (F1_SQ - F2_SQ)
BETA = -F2_SQ / (F1_SQ - F2_SQ)
OMEGA_E = 7.2921151467e-5
SEC_PER_DAY = 86400.0
MJD_J2000 = 51544.5


def run_batch_lsq(epoch_geo, ref_orbit, sp3, sv_bias, epochs,
                  r0_ecef, v0_ecef, mjd_start,
                  Cnm, Snm, Nmax, GM_grav=None, R_grav=None,
                  Cd_init=2.2, CR_init=1.3,
                  max_iter=6, dt_integ=5.0,
                  sigma_phase=0.01, sigma_code=0.30,
                  param_names=None,
                  seg_duration=900.0):
    """Run batch least squares with piecewise-constant empirical RTN.

    Args:
        epoch_geo: dict {gps_sod: [per_sv_dicts]} with pre-computed geometry
        ref_orbit: GNV1B reference orbit for linearization
        sp3: SP3 data dict
        sv_bias: per-SV code bias dict
        epochs: sorted list of GPS seconds-of-day
        r0_ecef, v0_ecef: initial ECEF state (from GNV1B)
        mjd_start: MJD(UTC) at arc start epoch
        Cnm, Snm: gravity field coefficients
        Nmax: maximum gravity degree
        GM_grav, R_grav: gravity model constants
        Cd_init, CR_init: initial force parameters
        max_iter: maximum iterations
        dt_integ: integration step size [s]
        sigma_phase, sigma_code: measurement noise [m]
        seg_duration: empirical RTN segment duration [s] (default 900s = 15min)

    Returns:
        dict with keys: 'r0', 'v0', 'Cd', 'CR', 'converged', 'iterations',
            'postfit_rms', 'residuals', 'cov', 'state', 'param_names'
    """
    N_EPOCH = len(epochs)
    if N_EPOCH < 3:
        return None

    # Identify SVs in the arc
    all_svs = set()
    for gps_sod in epochs:
        for d in epoch_geo.get(gps_sod, []):
            if d['sv'] in sv_bias:
                all_svs.add(d['sv'])
    sv_list = sorted(all_svs)
    N_SV = len(sv_list)
    sv_to_idx = {sv: i for i, sv in enumerate(sv_list)}

    t0 = epochs[0]
    t_end = epochs[-1]
    arc_duration = t_end - t0

    # Build per-segment empirical parameter names
    n_segments = max(1, int(np.ceil(arc_duration / seg_duration)))
    if param_names is None:
        param_names = ['Cd', 'CR', 'aR', 'aT', 'aN']
    # Global params: those without per-segment variants
    user_global_names = [p for p in param_names if p in ('Cd', 'CR')]
    user_emp_names = [p for p in param_names if p not in ('Cd', 'CR')]
    # Build parameter list: global + per-segment empirical
    seg_param_names = list(user_global_names)
    for i in range(n_segments):
        for base in user_emp_names:
            seg_param_names.append(f'{base}_{i}')

    Np = len(seg_param_names)
    N_ORB = 6 + Np + 1   # r0(3) + v0(3) + force(Np) + trop_wet
    I_R0 = slice(0, 3)
    I_V0 = slice(3, 6)
    I_TROP = 6 + Np
    param_idx = {name: 6 + k for k, name in enumerate(seg_param_names)}  # state index
    param_col = {name: k for k, name in enumerate(seg_param_names)}       # S-matrix column

    N_FULL = N_ORB + N_EPOCH + N_SV

    # Convert initial state to ECI
    r0_eci, v0_eci = ecef_to_eci(np.array(r0_ecef), np.array(v0_ecef), mjd_start)
    trop_wet = 0.0
    clk = np.zeros(N_EPOCH)
    B_sv = np.zeros(N_SV)
    # Force param values (tracked by name)
    force_val = {'Cd': float(Cd_init), 'CR': float(CR_init)}
    for pname in seg_param_names:
        if pname not in force_val:
            force_val[pname] = 0.0

    W_PHASE = 1.0 / sigma_phase**2
    W_CODE = 1.0 / sigma_code**2

    for iteration in range(max_iter):
        # -- Segment-by-segment orbit integration in ECI --
        mjd_now = mjd_start
        Cd_cur = force_val.get('Cd', Cd_init)
        CR_cur = force_val.get('CR', CR_init)

        force_model = _make_force_model(Cnm, Snm, Nmax, GM_grav, R_grav)

        # Integrate segment by segment with chained STM/S
        r_cur, v_cur = r0_eci.copy(), v0_eci.copy()
        Phi_cur = np.eye(6)   # STM from arc start to current segment start
        S_cur = np.zeros((6, Np))  # S from arc start to current segment start

        all_t, all_r, all_v, all_phi, all_S = [], [], [], [], []

        for i_seg in range(n_segments):
            t_start = i_seg * seg_duration
            t_end = min((i_seg + 1) * seg_duration, arc_duration)
            if t_end <= t_start:
                break

            seg_dur = t_end - t_start

            # Empirical acceleration for this segment
            emp_seg = np.array([
                force_val.get(f'aR_{i_seg}', 0.0) if f'aR_{i_seg}' in force_val else 0.0,
                force_val.get(f'aT_{i_seg}', 0.0) if f'aT_{i_seg}' in force_val else 0.0,
                force_val.get(f'aN_{i_seg}', 0.0) if f'aN_{i_seg}' in force_val else 0.0,
            ])

            # Active params for this segment: global + this segment's empiricals
            # Use BASE names for integrator (it knows 'Cd','CR','aR','aT','aN')
            seg_active_base = list(user_global_names) + list(user_emp_names)
            # Map from S-matrix column (base name) to full param name (with segment suffix)
            seg_col_map = list(user_global_names) + [f'{base}_{i_seg}' for base in user_emp_names]

            # Advance MJD for this segment's start time
            mjd_seg = mjd_start + t_start / SEC_PER_DAY

            # Integrate this segment
            integ_seg = integrate_orbit_eci_with_stm(
                r_cur, v_cur, (0.0, seg_dur), force_model,
                Cd=Cd_cur, CR=CR_cur,
                area_drag=0.68, area_srp=3.4, mass=580.0,
                empirical_acc_rtn=emp_seg,
                param_names=seg_active_base,
                dt=dt_integ,
                mjd_tt=mjd_seg + 69.184 / SEC_PER_DAY,
                mjd_utc=mjd_seg,
                bodies=['Sun', 'Moon'],
            )

            seg_t = integ_seg['t']
            seg_r = integ_seg['r']
            seg_v = integ_seg['v']
            seg_phi = integ_seg['phi']  # STM from segment start
            seg_S = integ_seg['S']      # S from segment start, (N,6,5)

            # Map segment results to full parameter space
            for k in range(len(seg_t)):
                # Full STM: chain from arc start
                Phi_full = seg_phi[k] @ Phi_cur

                # Full S: homogeneous part (past params carried forward)
                S_full = seg_phi[k] @ S_cur
                # Add current segment's active parameter contributions
                # seg_S[k] columns: [Cd, CR, aR_i, aT_i, aN_i]
                for col, pname in enumerate(seg_col_map):
                    if pname in param_col:
                        S_full[:, param_col[pname]] += seg_S[k, :, col]

                # Skip duplicate point at segment boundary
                if i_seg == 0 or k > 0:
                    all_t.append(t_start + seg_t[k])
                    all_r.append(seg_r[k])
                    all_v.append(seg_v[k])
                    all_phi.append(Phi_full)
                    all_S.append(S_full)

            # Update state for next segment
            r_cur = seg_r[-1].copy()
            v_cur = seg_v[-1].copy()
            Phi_cur = seg_phi[-1] @ Phi_cur
            # Carry forward S: homogeneous propagation + accumulated active params
            S_cur = seg_phi[-1] @ S_cur
            for col, pname in enumerate(seg_col_map):
                if pname in param_col:
                    S_cur[:, param_col[pname]] += seg_S[-1, :, col]

        # Convert to arrays
        t_integ = np.array(all_t)
        r_eci = np.array(all_r)
        v_eci = np.array(all_v)
        phi_grid = np.array(all_phi)
        S_grid = np.array(all_S)

        # -- Build normal equations --
        # We'll accumulate into orbital + ambiguity params first,
        # then pre-eliminate clocks
        N_orb_amb = N_ORB + N_SV
        N_full = np.zeros((N_FULL, N_FULL))
        b_full = np.zeros(N_FULL)

        n_phase = 0
        n_code = 0
        total_res_sq = 0.0
        all_residuals = []

        for i_ep, gps_sod in enumerate(epochs):
            ed_list = epoch_geo.get(gps_sod, [])
            if not ed_list:
                continue

            # Interpolate integrated orbit to this epoch
            t_rel = gps_sod - t0
            i_step = int(np.round(t_rel / dt_integ))
            if i_step < 0 or i_step >= len(t_integ):
                continue

            # Reference state in ECI and ECEF
            r_ref_eci = r_eci[i_step]
            v_ref_eci = v_eci[i_step]
            mjd_ep = mjd_start + t_rel / SEC_PER_DAY
            r_ref_ecef, v_ref_ecef = eci_to_ecef(r_ref_eci, v_ref_eci, mjd_ep)

            phi_ep = phi_grid[i_step]  # (6, 6)
            S_ep = S_grid[i_step]      # (6, Np)

            # STM sub-blocks for position
            phi_rr = phi_ep[0:3, 0:3]   # ∂r/∂r0
            phi_rv = phi_ep[0:3, 3:6]   # ∂r/∂v0

            # Sensitivity sub-blocks (position only, for range partials)
            # S_r_p = position partials w.r.t. parameter p
            S_r = S_ep[0:3, :]  # (3, Np)

            for d in ed_list:
                sv = d['sv']
                if sv not in sv_to_idx or sv not in sv_bias:
                    continue

                i_sv = sv_to_idx[sv]
                sat_pos = d['sat_pos']
                sat_clk = d['sat_clk']
                rho_corr = d['rho_corr']
                el = d['el']
                e_vec = (sat_pos - r_ref_ecef) / rho_corr  # ECEF LOS
                mf = 1.0 / max(np.sin(el), 0.1)

                # Compute ECI LOS for STM-based partials.
                # phi_rr, phi_rv, S_r are ∂r_ECI/∂(r0,v0,p) — must be
                # dotted with the LOS expressed in ECI, not ECEF.
                sat_pos_eci, _ = ecef_to_eci(sat_pos, np.zeros(3), mjd_ep)
                e_vec_eci = (sat_pos_eci - r_ref_eci)
                e_vec_eci = e_vec_eci / float(np.linalg.norm(e_vec_eci))

                # Compute modeled range + Sagnac correction (ECEF)
                rho_model = float(np.linalg.norm(sat_pos - r_ref_ecef))
                sag = (OMEGA_E / C_LIGHT) * (sat_pos[0] * r_ref_ecef[1] - sat_pos[1] * r_ref_ecef[0])
                rho_corr_model = rho_model + sag

                # Common partials for position, velocity, trop
                h_base = np.zeros(N_FULL)
                h_base[I_R0] = -e_vec_eci @ phi_rr
                h_base[I_V0] = -e_vec_eci @ phi_rv
                h_base[I_TROP] = mf
                # Force param partials (dynamic, from S matrix)
                for k, pname in enumerate(seg_param_names):
                    h_base[param_idx[pname]] = -e_vec_eci @ S_r[:, k]

                # -- Phase observation --
                phase_residual = (d['L_if_raw'] + sat_clk - rho_corr_model
                                  - sv_bias[sv] - clk[i_ep] - trop_wet * mf - B_sv[i_sv])
                if abs(phase_residual) < 200.0:
                    obs = phase_residual

                    h = h_base.copy()
                    h[N_ORB + i_ep] = 1.0          # Δclk
                    h[N_ORB + N_EPOCH + i_sv] = 1.0  # ΔB

                    N_full += W_PHASE * np.outer(h, h)
                    b_full += W_PHASE * obs * h
                    n_phase += 1
                    total_res_sq += W_PHASE * obs**2
                    all_residuals.append(obs)

                # -- Code observation --
                code_residual = (d['P_if_raw'] + sat_clk - rho_corr_model
                                 - sv_bias[sv] - clk[i_ep] - trop_wet * mf)
                if abs(code_residual) < 100.0:
                    obs = code_residual

                    h = h_base.copy()
                    h[N_ORB + i_ep] = 1.0
                    # No ambiguity for code

                    N_full += W_CODE * np.outer(h, h)
                    b_full += W_CODE * obs * h
                    n_code += 1
                    total_res_sq += W_CODE * obs**2
                    all_residuals.append(obs)

        if n_phase + n_code < N_ORB:
            print(f"  iter {iteration}: insufficient obs ({n_phase + n_code} < {N_ORB})")
            return None

        # -- Prior constraints --
        prior_cfg = {}
        if 'Cd' in seg_param_names:
            prior_cfg['Cd'] = (2.2, 1.0)
        if 'CR' in seg_param_names:
            prior_cfg['CR'] = (1.3, 0.25)
        # Per-segment empirical priors
        # Dynamics error ~20-30m over 30min → need ~5e-5 m/s^2 to compensate
        # Set prior sigma to allow this range
        EMP_PRIOR_SIGMA = 1e-5  # m/s^2 (~50m pos effect over 900s segment)
        for i_seg in range(n_segments):
            for base in user_emp_names:
                pname = f'{base}_{i_seg}'
                prior_cfg[pname] = (0.0, EMP_PRIOR_SIGMA**2)

        for pname in seg_param_names:
            if pname in param_idx and pname in prior_cfg:
                idx = param_idx[pname]
                ref_val, pvar = prior_cfg[pname]
                cur_val = force_val.get(pname, ref_val)
                N_full[idx, idx] += 1.0 / pvar
                b_full[idx] += (ref_val - cur_val) / pvar

        # Trop wet: 0 ± 0.5m
        N_full[I_TROP, I_TROP] += 1.0 / 0.25
        b_full[I_TROP] += (0.0 - trop_wet) / 0.25

        # Ambiguity: 0 ± 100m (weak)
        for i_sv in range(N_SV):
            idx = N_ORB + N_EPOCH + i_sv
            N_full[idx, idx] += 1.0 / 10000.0

        # -- Solve full system (skip clock pre-elimination for debugging) --
        # Add small regularization for numerical stability
        N_reg = N_full + np.eye(N_FULL) * 1e-10
        try:
            dx_full = np.linalg.solve(N_reg, b_full)
        except np.linalg.LinAlgError:
            print(f"  iter {iteration}: singular normal matrix")
            return None

        dx_reduced = np.zeros(N_ORB + N_SV)
        dx_reduced[:N_ORB] = dx_full[:N_ORB]
        dx_reduced[N_ORB:] = dx_full[N_ORB + N_EPOCH:N_ORB + N_EPOCH + N_SV]
        dx_clock = dx_full[N_ORB:N_ORB + N_EPOCH]

        # -- Update state --
        dr0_eci = dx_reduced[I_R0]
        dv0_eci = dx_reduced[I_V0]
        dtrop = dx_reduced[I_TROP]

        r0_eci += dr0_eci
        v0_eci += dv0_eci
        trop_wet += dtrop
        clk += dx_clock

        # Update force params
        d_force = {}
        for pname in seg_param_names:
            if pname in param_idx:
                idx = param_idx[pname]
                dval = float(dx_reduced[idx])
                force_val[pname] += dval
                d_force[pname] = dval

        # Update ambiguities
        for i_sv in range(N_SV):
            idx_red = N_ORB + i_sv
            B_sv[i_sv] += dx_reduced[idx_red]

        # -- Convergence check --
        pos_change = float(np.linalg.norm(dr0_eci))
        vel_change = float(np.linalg.norm(dv0_eci))
        rms_res = np.sqrt(total_res_sq / max(n_phase + n_code, 1))
        df_summary = ' '.join(f"d{pn}={d_force.get(pn, 0.0):.3e}"
                              for pn in seg_param_names[:8])  # first 8 params
        print(f"  iter {iteration}: dr={pos_change:.4f}m dv={vel_change:.5f}m/s "
              f"dtrop={dtrop:.4f} {df_summary} obs={n_phase}p+{n_code}c "
              f"rms_res={rms_res:.4f}")

        if pos_change < 0.001 and vel_change < 0.0001:
            print(f"  Converged after {iteration + 1} iterations")
            break

    # Convert final state back to ECEF
    r0_ecef, v0_ecef = eci_to_ecef(r0_eci, v0_eci, mjd_start)

    # Build aR/aT/aN arrays for backward compatibility (first segment or mean)
    aR_vals = [force_val.get(f'aR_{i}', 0.0) for i in range(n_segments)]
    aT_vals = [force_val.get(f'aT_{i}', 0.0) for i in range(n_segments)]
    aN_vals = [force_val.get(f'aN_{i}', 0.0) for i in range(n_segments)]

    return {
        'r0': r0_ecef,
        'v0': v0_ecef,
        'r0_eci': r0_eci,
        'v0_eci': v0_eci,
        'Cd': force_val.get('Cd', Cd_init),
        'CR': force_val.get('CR', CR_init),
        'aR': float(np.mean(aR_vals)) if aR_vals else 0.0,
        'aT': float(np.mean(aT_vals)) if aT_vals else 0.0,
        'aN': float(np.mean(aN_vals)) if aN_vals else 0.0,
        'aR_seg': aR_vals,
        'aT_seg': aT_vals,
        'aN_seg': aN_vals,
        'n_segments': n_segments,
        'trop_wet': trop_wet,
        'clk': clk,
        'B_sv': {sv: B_sv[i] for sv, i in sv_to_idx.items()},
        'converged': pos_change < 0.001 and vel_change < 0.0001,
        'iterations': min(iteration + 1, max_iter),
        'postfit_rms': rms_res,
        'n_phase': n_phase,
        'n_code': n_code,
        'residuals': all_residuals,
        'state': np.concatenate([r0_ecef, v0_ecef,
                                 [force_val.get(pn, 0.0) for pn in seg_param_names],
                                 [trop_wet], clk, B_sv]),
        'sv_list': sv_list,
        'param_names': seg_param_names,
    }


def _make_force_model(Cnm, Snm, Nmax, GM_grav, R_grav):
    """Create a closure for the ECI force model with fixed gravity coefficients."""
    from src.orbit_dynamics import total_acc_eci

    def force_model(pos_eci, vel_eci, **kwargs):
        return total_acc_eci(
            pos_eci, vel_eci,
            Cnm=Cnm, Snm=Snm, Nmax=Nmax,
            GM_gravity=GM_grav, R_gravity=R_grav,
            **kwargs,
        )

    return force_model


def _pre_eliminate_clocks(N_full, b_full, N_ORB, N_EPOCH, N_SV):
    """Pre-eliminate epoch-wise clock parameters via Schur complement.

    Each clock parameter appears only in its own epoch's observations.
    This makes the pre-elimination exact and efficient.

    Partition:
        N = [[N_aa, N_ab],    a = orbital + ambiguity params
             [N_ba, N_bb]]    b = clock params (diagonal-dominant)

    Reduced system:
        N_reduced = N_aa - N_ab @ N_bb^{-1} @ N_ba
        b_reduced = b_a - N_ab @ N_bb^{-1} @ b_b

    Returns:
        N_reduced: (N_ORB + N_SV, N_ORB + N_SV)
        b_reduced: (N_ORB + N_SV,)
        clock_idx_map: mapping from old clock index to new index
    """
    N_RED = N_ORB + N_SV
    N_reduced = np.zeros((N_RED, N_RED))
    b_reduced = np.zeros(N_RED)

    # Indices for "a" block (orbital + ambiguity)
    a_indices = list(range(N_ORB)) + list(range(N_ORB + N_EPOCH, N_ORB + N_EPOCH + N_SV))

    # Copy N_aa and b_a
    for i, ai in enumerate(a_indices):
        for j, aj in enumerate(a_indices):
            N_reduced[i, j] = N_full[ai, aj]
        b_reduced[i] = b_full[ai]

    # For each clock parameter, do the Schur complement
    for k in range(N_EPOCH):
        clk_idx = N_ORB + k
        n_bb = N_full[clk_idx, clk_idx]

        if n_bb < 1e-15:
            continue

        # N_ab column for this clock
        n_ab = np.zeros(N_RED)
        for i, ai in enumerate(a_indices):
            n_ab[i] = N_full[ai, clk_idx]

        b_b = b_full[clk_idx]

        # Schur complement: N_reduced -= n_ab * n_ab^T / n_bb
        N_reduced -= np.outer(n_ab, n_ab) / n_bb
        b_reduced -= n_ab * b_b / n_bb

    return N_reduced, b_reduced, None


def _recover_clocks(N_full, b_full, N_ORB, N_EPOCH, dx_reduced):
    """Recover clock parameter updates via back-substitution.

    dx_clk[k] = (b_clk - N_clk,a · dx_a) / N_clk,clk
    """
    dx_clock = np.zeros(N_EPOCH)
    a_indices = list(range(N_ORB)) + list(range(N_ORB + N_EPOCH, N_ORB + N_EPOCH + len(dx_reduced) - N_ORB))

    for k in range(N_EPOCH):
        clk_idx = N_ORB + k
        n_bb = N_full[clk_idx, clk_idx]

        if n_bb < 1e-15:
            dx_clock[k] = 0.0
            continue

        # Sum of N_clk,a_j · dx_a_j
        n_dot_dx = 0.0
        for j, aj in enumerate(a_indices):
            if j < len(dx_reduced):
                n_dot_dx += N_full[clk_idx, aj] * dx_reduced[j]

        dx_clock[k] = (b_full[clk_idx] - n_dot_dx) / n_bb

    return dx_clock
