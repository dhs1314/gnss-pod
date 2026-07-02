"""Code-only batch orbit determination (Phase 26.0).

Pseudorange-only Gauss-Newton orbit solver for self-consistent initial orbit.
No GNV1B dependency, no ambiguity parameters.

Pipeline:
  Step 0: Kinematic WLS single-epoch positioning (~5-10m)
  Step 1: Code-only GN over long arc (~0.15-0.20m with Orekit dynamics)
  Step 2: Outlier detection from code residuals

Usage:
    from src.code_orbit import kinematic_wls_single_epoch, CodeOnlyOrbitSolver

    # Step 0
    r0_kin, clk_kin = kinematic_wls_single_epoch(gps1b, gps_sod, sp3)

    # Step 1
    solver = CodeOnlyOrbitSolver(geometry, t_epochs, mjd_start, mjd_tt, orekit_prop)
    result = solver.solve(r0_kin, np.zeros(3))

    # Step 2
    bad_svs = solver.detect_outliers(result, threshold=5.0)
"""

import numpy as np
from datetime import datetime, timedelta

C_LIGHT = 299792458.0
OMEGA_E = 7.2921151467e-5
SEC_PER_DAY = 86400.0
MJD0 = 51544.5  # J2000.0


# ══════════════════════════════════════════════════════════════════════
# Step 0: Kinematic WLS single-epoch positioning
# ══════════════════════════════════════════════════════════════════════

def kinematic_wls_single_epoch(gps1b, gps_sod, sp3,
                                max_iter=8, tol=0.01):
    """Single-epoch weighted least squares from pseudorange only.

    Solves: P_if = |sat_pos - rcv| + clk + noise
    Unknowns: rcv(3) in ECEF, clk(1) in meters

    Args:
        gps1b: dict {gps_sec: {sv: {obs...}}}
        gps_sod: GPS seconds of day
        sp3: SP3 data dict
        max_iter: Newton iterations
        tol: convergence tolerance [m]

    Returns:
        (r_ecef, clk): np.array(3), float — ECEF position [m] and clock [m]
    """
    J2000 = datetime(2000, 1, 1, 12, 0, 0)
    utc_dt = J2000 + timedelta(seconds=gps_sod)
    recs = gps1b.get(int(gps_sod), gps1b.get(gps_sod, {}))

    # Collect visible SVs
    from run_sequential_pod import get_sat_geometry

    sv_obs = []
    for sv_id, rec in recs.items():
        if 'P_if' not in rec:
            continue
        sat_pos, sat_clk, _ = get_sat_geometry(sp3, sv_id, utc_dt,
                                                np.array([0, 0, 6371000.0]),
                                                clk_data=None)
        if sat_pos is None:
            continue
        sv_obs.append({
            'sv': sv_id,
            'sat_pos': np.asarray(sat_pos, dtype=float),
            'sat_clk': float(sat_clk),
            'P_if': float(rec['P_if']),
        })

    n_sv = len(sv_obs)
    if n_sv < 4:
        return np.array([0.0, 0.0, 6371000.0]), 0.0

    # Initial guess: use SVs' average range to estimate altitude
    # GRACE-FO altitude ~490km
    h_guess = 490000.0  # m above equatorial radius
    r_est = np.array([6378000.0 + h_guess, 0.0, 0.0])  # start near equator
    clk_est = 0.0

    best_r = r_est.copy()
    best_cost = 1e30

    for it in range(max_iter):
        A = np.zeros((n_sv, 4))
        b = np.zeros(n_sv)

        for i, obs in enumerate(sv_obs):
            sat_pos = obs['sat_pos']
            los = r_est - sat_pos
            rho = float(np.linalg.norm(los))

            # Sag correction
            sag = (OMEGA_E / C_LIGHT) * (sat_pos[0] * r_est[1]
                                          - sat_pos[1] * r_est[0])

            # Model: P_if = rho + sag + clk + sat_clk
            model = rho + sag + clk_est + obs['sat_clk']

            # Jacobian: dr/∂x = (x - x_sat)/rho
            A[i, 0:3] = los / rho
            A[i, 3] = 1.0
            b[i] = obs['P_if'] - model

        # Solve: dx = (A^T W A)^{-1} A^T W b
        W = np.eye(n_sv) / (0.30**2)
        try:
            dx = np.linalg.solve(A.T @ W @ A, A.T @ W @ b)
        except np.linalg.LinAlgError:
            dx = np.linalg.lstsq(A, b, rcond=1e-6)[0]

        r_est += dx[0:3]
        clk_est += dx[3]

        # Reject solutions inside Earth
        h = float(np.linalg.norm(r_est)) - 6378137.0
        if h < 100000:  # < 100km altitude = likely wrong
            r_est = best_r.copy()
            for i in range(3):
                r_est[i] += np.random.normal(0, 1e6)
            continue

        # Cost check
        resid = b - A @ dx
        cost = float(np.sum(resid**2))
        if cost < best_cost:
            best_cost = cost
            best_r = r_est.copy()

        if float(np.linalg.norm(dx[0:3])) < tol:
            break

    # Validate final result
    h_final = float(np.linalg.norm(best_r)) - 6378137.0
    if h_final < 100000 or h_final > 2000000:
        # Failed — fall back to approximate orbit
        best_r = np.array([6378137.0 + 490000.0, 0.0, 0.0])

    return best_r.copy(), float(clk_est)


# ══════════════════════════════════════════════════════════════════════
# Step 1: Code-only GN orbit solver
# ══════════════════════════════════════════════════════════════════════

class CodeOnlyOrbitSolver:
    """Pseudorange-only Gauss-Newton orbit determination.

    Uses ONLY P_if observations. No phase, no ambiguity parameters.
    Orekit dynamics provides orbit constraint.

    Parameters (GN loop):
      r0(3), v0(3), aR, aT, aN = 9 nonlinear parameters
      clk_i(N_epoch) + zwd_i(N_epoch) = solved analytically per GN iteration

    Usage:
        solver = CodeOnlyOrbitSolver(geometry, t_epochs, mjd_start, mjd_tt, orekit_prop)
        result = solver.solve(r0_kinematic, v0=np.zeros(3))
    """

    def __init__(self, geometry, t_epochs, mjd_utc_start, mjd_tt_start,
                 orekit_prop, sigma_code=0.30, max_iter=4, damping=0.5,
                 prior_r0=5.0, prior_v0=1.0, prior_emp=1e-7,
                 force_fn=None):
        """
        Args:
            geometry: list of lists, geometry[i_ep] = [ep_dict, ...]
                      each ep_dict must have: sv, sat_pos, sat_clk, el,
                      _obs_code, _geo_full (or equivalent)
            t_epochs: array of epoch times [s] from t0
            mjd_utc_start: MJD(UTC) at first epoch
            mjd_tt_start: MJD(TT) at first epoch
            orekit_prop: OrekitPropagator instance (or None for Python-only)
            sigma_code: pseudorange measurement noise [m]
            max_iter: max GN iterations
            damping: Levenberg-Marquardt damping factor
            prior_r0: prior sigma for r0 [m]
            prior_v0: prior sigma for v0 [m/s]
            prior_emp: prior sigma for empirical RTN [m/s²]
            force_fn: Python force function for dynamics (total_acc_eci closure)
        """
        self.geometry = geometry
        self.t_epochs = np.asarray(t_epochs, dtype=float)
        self.mjd_utc_start = mjd_utc_start
        self.mjd_tt_start = mjd_tt_start
        self.orekit_prop = orekit_prop
        self.force_fn = force_fn
        self.sigma_code = sigma_code
        self.max_iter = max_iter
        self.damping = damping
        self.prior_r0 = prior_r0
        self.prior_v0 = prior_v0
        self.prior_emp = prior_emp
        self.N_epoch = len(geometry)
        self._use_orekit = orekit_prop is not None
        self._index_svs()
        self._precompute_geo()

    def _index_svs(self):
        """Index unique SVs in geometry."""
        sv_set = set()
        for ep_list in self.geometry:
            for d in ep_list:
                if '_obs_code' in d:
                    sv_set.add(d['sv'])
        self.sv_list = sorted(sv_set)
        self.N_sv = len(self.sv_list)

    def _precompute_geo(self):
        """Pre-compute _geo_full for code observations if missing."""
        for ep_list in self.geometry:
            for d in ep_list:
                if '_geo_full' not in d and '_obs_code' in d:
                    # Approximate geo from satellite geometry — caller
                    # should have filled this in
                    d['_geo_full'] = 0.0

    # ─── GN outer loop ───

    def solve(self, r0_guess, v0_guess):
        """Gauss-Newton iteration for r0, v0, aRTN.

        Args:
            r0_guess: initial ECEF position [m] (from kinematic WLS)
            v0_guess: initial ECEF velocity [m/s] (zeros is fine)

        Returns:
            dict with keys:
                r_eci: (N_epoch, 3) ECI positions
                v_eci: (N_epoch, 3) ECI velocities
                clk: (N_epoch,) clock [m]
                zwd: (N_epoch,) ZWD [m]
                rms_code: float
                per_sv_rms: {sv: float} code RMS per SV
                converged: bool
                iterations: int
                r0: (3,) final initial position
                v0: (3,) final initial velocity
                a_emp: (3,) empirical RTN [m/s²]
                residuals: (N_obs,) code residuals [m]
        """
        # Convert ECEF guess → ECI
        from src.coordinates import ecef_to_eci, eci_to_ecef
        mjd_start = self.mjd_utc_start
        r0i, v0i = ecef_to_eci(np.asarray(r0_guess, dtype=float),
                                np.asarray(v0_guess, dtype=float),
                                mjd_start)

        # State vector: r0(3), v0(3), aR,aT,aN(3)
        x_nl = np.zeros(9)
        x_nl[0:3] = r0i
        x_nl[3:6] = v0i
        # aRTN starts at zero

        converged = False
        prev_cost = 1e30

        for it in range(self.max_iter):
            # 1. Propagation (Orekit preferred, Python fallback)
            if self._use_orekit:
                try:
                    r_eci, v_eci = self.orekit_prop.propagate_continuous_arc(
                        x_nl[0:3], x_nl[3:6],
                        self.t_epochs, x_nl[6:9],
                        self.mjd_utc_start,
                    )
                    # Build Phi/S from Python for Jacobian
                    _, _, phi_ep, S_ep = self._propagate_with_stm(x_nl)
                except Exception:
                    if it == 0:
                        print(f"  [CodeGN] Orekit propagation failed, using Python dynamics")
                    self._use_orekit = False
                    r_eci, v_eci, phi_ep, S_ep = self._propagate_with_stm(x_nl)
            else:
                r_eci, v_eci, phi_ep, S_ep = self._propagate_with_stm(x_nl)

            # 2. Rebuild geometry for current orbit
            geo_cur = self._rebuild_geometry(r_eci, v_eci)

            # 3. Solve linear sub-problem: clk + zwd
            clk, zwd, code_res, per_sv_res = self._solve_linear(geo_cur)

            code_res_arr = np.asarray(code_res, dtype=float)
            cost = float(0.5 * np.sum(code_res_arr**2))
            rms_code = float(np.sqrt(np.mean(code_res_arr**2)))

            if it == 0:
                print(f"  [CodeGN {it}] cost={cost:.1f} rms={rms_code:.3f}m "
                      f"a0=[{x_nl[6]:.2e},{x_nl[7]:.2e},{x_nl[8]:.2e}]")

            # 4. Check convergence
            dr_norm = float(np.linalg.norm(x_nl[0:3] - r0i))
            if it > 0:
                if dr_norm < 0.005 and abs(cost - prev_cost) / max(abs(cost), 1) < 1e-4:
                    converged = True

            if converged or it == self.max_iter - 1:
                # Final output
                r_ecef = np.zeros((self.N_epoch, 3))
                v_ecef = np.zeros((self.N_epoch, 3))
                for i_ep in range(self.N_epoch):
                    mjd_u = mjd_start + self.t_epochs[i_ep] / SEC_PER_DAY
                    r_ecef[i_ep], v_ecef[i_ep] = eci_to_ecef(
                        r_eci[i_ep], v_eci[i_ep], mjd_u)
                return {
                    'r_eci': r_eci, 'v_eci': v_eci,
                    'r_ecef': r_ecef, 'v_ecef': v_ecef,
                    'clk': clk, 'zwd': zwd,
                    'rms_code': rms_code,
                    'per_sv_rms': {sv: float(np.sqrt(np.mean(res**2)))
                                   for sv, res in per_sv_res.items()},
                    'converged': converged, 'iterations': it + 1,
                    'r0': x_nl[0:3].copy(), 'v0': x_nl[3:6].copy(),
                    'a_emp': x_nl[6:9].copy(),
                    'residuals': code_res,
                }

            prev_cost = cost

            # 5. Build Jacobian (chain rule through Phi STM)
            _, _, phi_ep, S_ep = self._propagate_with_stm(x_nl)

            J = self._build_jacobian(r_eci, phi_ep, S_ep, code_res, geo_cur)

            # 6. Normal equations
            H = J.T @ J
            g = J.T @ code_res

            # Priors
            H[0, 0] += 1.0 / self.prior_r0**2
            H[1, 1] += 1.0 / self.prior_r0**2
            H[2, 2] += 1.0 / self.prior_r0**2
            H[3, 3] += 1.0 / self.prior_v0**2
            H[4, 4] += 1.0 / self.prior_v0**2
            H[5, 5] += 1.0 / self.prior_v0**2
            for j in range(6, 9):
                H[j, j] += 1.0 / self.prior_emp**2

            # Damping
            H += self.damping * np.diag(np.diag(H) + 1e-6)

            try:
                dx = -np.linalg.solve(H, g)
            except np.linalg.LinAlgError:
                dx = -np.linalg.lstsq(H, g, rcond=1e-6)[0]

            # 7. Line search
            alpha = 1.0
            for _ in range(4):
                x_try = x_nl + alpha * dx
                if self._use_orekit:
                    try:
                        r_try, v_try = self.orekit_prop.propagate_continuous_arc(
                            x_try[0:3], x_try[3:6],
                            self.t_epochs, x_try[6:9],
                            self.mjd_utc_start,
                        )
                    except Exception:
                        r_try, v_try, _, _ = self._propagate_with_stm(x_try)
                else:
                    r_try, v_try, _, _ = self._propagate_with_stm(x_try)
                geo_try = self._rebuild_geometry(r_try, v_try)
                _, _, res_try, _ = self._solve_linear(geo_try)
                if 0.5 * np.sum(res_try**2) < cost:
                    break
                alpha *= 0.5

            x_nl += alpha * dx
            dr_n = float(np.linalg.norm(alpha * dx[0:3]))
            dv_n = float(np.linalg.norm(alpha * dx[3:6]))
            da_n = float(np.linalg.norm(alpha * dx[6:9]))
            ls_tag = f" a={alpha:.1f}" if alpha < 1.0 else ""
            print(f"  [CodeGN {it}] cost={cost:.1f} rms={rms_code:.3f}m "
                  f"dr={dr_n:.4f}m dv={dv_n:.6f}m/s da={da_n:.2e}{ls_tag}"
                  f"  a0=[{x_nl[6]:.2e},{x_nl[7]:.2e},{x_nl[8]:.2e}]"
                  + (" CONV" if converged else ""))

        # Should not reach here
        rms_code = float(r0i[0])  # dummy
        return {'rms_code': rms_code, 'converged': False,
                'iterations': self.max_iter}

    # ─── Linear sub-problem ───

    def _solve_linear(self, geo_cur):
        """Solve clk_i + zwd_i per epoch independently.

        Model: _obs_code - _geo_full = clk_i + zwd_i * mf + noise
        Each epoch is 2x2 independent (no amb linking epochs).

        Returns:
            clk: (N_epoch,), zwd: (N_epoch,),
            residuals: (N_obs,) all code residuals,
            per_sv_res: {sv: [residuals]}
        """
        from src.batch_solver import _trop_mf

        W = 1.0 / self.sigma_code**2
        clk = np.zeros(self.N_epoch)
        zwd = np.zeros(self.N_epoch)
        all_res = []
        per_sv_res = {}

        for i_ep, ep_list in enumerate(geo_cur):
            # Build 2x2 normal equations for this epoch
            ATA = np.zeros((2, 2))
            ATb = np.zeros(2)

            for d in ep_list:
                sv = d['sv']
                geo = float(d.get('_geo_full', 0.0))
                obs = float(d.get('_obs_code', 0.0))
                el = float(d.get('el', 0.5))
                mf = _trop_mf(el)

                ATA[0, 0] += W
                ATA[0, 1] += W * mf
                ATA[1, 0] += W * mf
                ATA[1, 1] += W * mf * mf

                res_pre = obs - geo
                ATb[0] += W * res_pre
                ATb[1] += W * res_pre * mf

            # Regularize (ZWD weakly)
            ATA[1, 1] += 1e-4

            try:
                sol = np.linalg.solve(ATA, ATb)
            except np.linalg.LinAlgError:
                # Fallback: clk only
                clk[i_ep] = ATb[0] / max(ATA[0, 0], 1e-6)
                zwd[i_ep] = 0.0
            else:
                clk[i_ep] = float(sol[0])
                zwd[i_ep] = float(sol[1])

            # Compute residuals
            for d in ep_list:
                sv = d['sv']
                geo = float(d.get('_geo_full', 0.0))
                obs = float(d.get('_obs_code', 0.0))
                el = float(d.get('el', 0.5))
                mf = _trop_mf(el)
                model = geo + clk[i_ep] + zwd[i_ep] * mf
                res = obs - model
                all_res.append(res)
                per_sv_res.setdefault(sv, []).append(res)

        return clk, zwd, np.array(all_res), per_sv_res

    # ─── Geometry rebuild ───

    def _rebuild_geometry(self, r_eci, v_eci):
        """Rebuild epoch_geometry with _geo_full for current orbit.

        _geo_full = rho + sag - sat_clk + zhd*mf
        where rho = |sat_eci - rcv_eci|
        """
        from src.coordinates import ecef_to_eci, eci_to_ecef

        new_geo = []
        for i_ep, ep_list in enumerate(self.geometry):
            mjd_u = self.mjd_utc_start + self.t_epochs[i_ep] / SEC_PER_DAY
            r_ecef, _ = eci_to_ecef(r_eci[i_ep], v_eci[i_ep], mjd_u)

            ep_new = []
            for d in ep_list:
                d_new = dict(d)
                sat_ecef = np.asarray(d['sat_pos'], dtype=float)
                sat_clk = float(d.get('sat_clk', 0.0))
                rho = float(np.linalg.norm(sat_ecef - r_ecef))
                sag = (OMEGA_E / C_LIGHT) * (sat_ecef[0] * r_ecef[1]
                                              - sat_ecef[1] * r_ecef[0])
                d_new['_geo_full'] = rho + sag - sat_clk
                ep_new.append(d_new)
            new_geo.append(ep_new)
        return new_geo

    # ─── Jacobian (chain rule through STM) ───

    def _propagate_with_stm(self, x_nl):
        """Propagate with Python integrator to get orbit, Phi, and S at each epoch."""
        from src.orbit_integrator import _rk4_step_eci_with_stm
        from src.orbit_dynamics import total_acc_eci

        force = self.force_fn if self.force_fn is not None else total_acc_eci

        r0 = x_nl[0:3].copy()
        v0 = x_nl[3:6].copy()
        a_rtn = x_nl[6:9].copy()

        N = len(self.t_epochs)
        r = np.zeros((N, 3))
        v = np.zeros((N, 3))
        phi_all = []
        S_all = []

        state_r = r0.copy()
        state_v = v0.copy()
        Phi = np.eye(6)
        S = np.zeros((6, 3))
        param_names = ['aR', 'aT', 'aN']
        t_cur = 0.0

        for i_ep in range(N):
            r[i_ep] = state_r.copy()
            v[i_ep] = state_v.copy()
            phi_all.append(Phi.copy())
            S_all.append(S.copy())

            if i_ep == N - 1:
                break

            dt = self.t_epochs[i_ep + 1] - t_cur
            n_sub = max(1, int(np.ceil(dt / 10.0)))
            h = dt / n_sub

            for _ in range(n_sub):
                mjd_utc = self.mjd_utc_start + t_cur / SEC_PER_DAY
                mjd_tt = self.mjd_tt_start + t_cur / SEC_PER_DAY
                state_r, state_v, Phi, S = _rk4_step_eci_with_stm(
                    state_r, state_v, Phi, S, force,
                    Cd=2.2, CR=1.3, area_drag=0.68, area_srp=3.4, mass=580.0,
                    empirical_acc_rtn=a_rtn, dt=h, param_names=param_names,
                    mjd_utc=mjd_utc, mjd_tt=mjd_tt)
                t_cur += h

        return r, v, phi_all, S_all

    def _build_jacobian(self, r_eci, phi_ep, S_ep, residuals, geo_cur):
        """Build Jacobian J = dres/d(r0,v0,aRTN).

        Chain rule:
          dres/dr0 = dres/dr_ep × dr_ep/dr0 = LOS × Phi_rr
          dres/dv0 = LOS × Phi_rv
          dres/daRTN = LOS × S(6×3)

        Code residual = obs - (geo + clk + zwd*mf)
        For code: dres/dpos = -LOS vector (unit from receiver to satellite)
        """
        from src.coordinates import eci_to_ecef

        n_obs = len(residuals)
        N_nl = 9
        J = np.zeros((n_obs, N_nl))
        obs_idx = 0

        for i_ep, ep_list in enumerate(geo_cur):
            Phi = phi_ep[i_ep]
            S_i = S_ep[i_ep]
            mjd_u = self.mjd_utc_start + self.t_epochs[i_ep] / SEC_PER_DAY
            r_ecef_ep, _ = eci_to_ecef(r_eci[i_ep], np.zeros(3), mjd_u)

            for d in ep_list:
                sat_ecef = np.asarray(d['sat_pos'], dtype=float)
                los_ecef = (r_ecef_ep - sat_ecef)
                rho = float(np.linalg.norm(los_ecef))
                if rho > 1e-6:
                    los_ecef = los_ecef / rho

                # Convert LOS from ECEF to ECI
                from src.coordinates import ecef_to_eci
                los_eci, _ = ecef_to_eci(los_ecef, np.zeros(3), mjd_u)

                # dr_ep/dr0 = Phi[0:3, 0:3]
                J[obs_idx, 0:3] = -Phi[0:3, 0:3].T @ los_eci
                # dr_ep/dv0 = Phi[0:3, 3:6]
                J[obs_idx, 3:6] = -Phi[0:3, 3:6].T @ los_eci
                # dr_ep/daRTN = S[0:3, :] (6x3 matrix, first 3 rows for position)
                J[obs_idx, 6:9] = -S_i[0:3, :].T @ los_eci

                obs_idx += 1

        return J

    # ─── Outlier detection ───

    def detect_outliers(self, result, threshold=5.0, max_bad_pct=0.30):
        """Detect outlier SVs from code residuals.

        Args:
            result: dict from solve()
            threshold: |residual| threshold in meters (default 5.0m ≈ 17σ)
            max_bad_pct: max fraction of bad epochs before rejecting SV

        Returns:
            dict {sv: {'n_bad': int, 'rms': float, 'reject': bool}}
        """
        per_sv_rms = result.get('per_sv_rms', {})
        residuals = result.get('residuals', np.array([]))

        # Reconstruct per-SV residual arrays
        per_sv_res = {}
        for ep_list in self.geometry:
            for d in ep_list:
                sv = d['sv']
                if '_obs_code' not in d:
                    continue
                geo = float(d.get('_geo_full', 0.0))
                obs = float(d['_obs_code'])
                resid = obs - geo  # approximate — should come from solve
                per_sv_res.setdefault(sv, []).append(resid)

        # Re-compute using the solved clk/zwd
        if result.get('clk') is not None:
            clk = result['clk']
            zwd = result.get('zwd', np.zeros_like(clk))
            from src.batch_solver import _trop_mf
            per_sv_res = {}
            for i_ep, ep_list in enumerate(self.geometry):
                for d in ep_list:
                    sv = d['sv']
                    geo = float(d.get('_geo_full', 0.0))
                    obs = float(d.get('_obs_code', 0.0))
                    el = float(d.get('el', 0.5))
                    mf = _trop_mf(el)
                    model = geo + clk[i_ep] + zwd[i_ep] * mf
                    per_sv_res.setdefault(sv, []).append(obs - model)

        n_total_ep = self.N_epoch
        decisions = {}

        for sv, res_list in sorted(per_sv_res.items()):
            if len(res_list) < 3:
                decisions[sv] = {'n_bad': 0, 'rms': 0.0, 'reject': False}
                continue
            res_arr = np.array(res_list)
            n_bad = int(np.sum(np.abs(res_arr) > threshold))
            rms = float(np.sqrt(np.mean(res_arr**2)))
            pct_bad = n_bad / max(len(res_list), 1)

            reject = (pct_bad > max_bad_pct) or (rms > 3.0)
            decisions[sv] = {
                'n_bad': n_bad,
                'n_total': len(res_list),
                'pct_bad': float(pct_bad),
                'rms': rms,
                'mean_resid': float(np.mean(res_arr)),
                'reject': reject,
            }

        n_reject = sum(1 for v in decisions.values() if v['reject'])
        if n_reject > 0:
            bad_list = [sv for sv, v in decisions.items() if v['reject']]
            print(f"  [CodeOutlier] {n_reject}/{len(decisions)} SVs flagged: {bad_list}")
            for sv in bad_list:
                d = decisions[sv]
                print(f"    {sv}: rms={d['rms']:.2f}m mean={d['mean_resid']:.2f}m "
                      f"bad={d['n_bad']}/{d['n_total']} ({d['pct_bad']*100:.0f}%)")

        return decisions
