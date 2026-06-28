"""Full-arc batch least-squares orbit determination (Phase 6.0).

Jointly estimates:
  Nonlinear: r0(3), v0(3), aR, aT, aN  (9 params)
  Linear:    clk(N), zwd(N), amb(M)     (via BatchLinearSolver)

Uses finite-difference Jacobian + Gauss-Newton with the existing
orbit integrator. Solves the clock-ambiguity cross-reference problem
by estimating everything in one consistent system.

Usage:
    solver = BatchOrbitLSQ(pass1_geometry, force_fn, t_epochs,
                            mjd_utc_start, mjd_tt_start)
    solution = solver.solve(r0_prior, v0_prior)
"""

import numpy as np
from src.orbit_integrator import integrate_orbit_eci_with_stm
from src.coordinates import eci_to_ecef
from src.batch_solver import BatchLinearSolver

C_LIGHT = 299792458.0
OMEGA_E = 7.2921151467e-5
SEC_PER_DAY = 86400.0


def _propagate_and_interp_segments(r0, v0, a_emp_segments, segment_dt,
                                    t_epochs, force_fn,
                                    mjd_utc_start, mjd_tt_start):
    """Propagate orbit with piecewise-constant RTN per segment.

    Integrates from segment 0 to segment K-1, switching empirical
    acceleration at each segment boundary. Interpolates to observation
    epochs.

    Args:
        a_emp_segments: (K, 3) RTN per segment [m/s²]
        segment_dt: segment duration [s]

    Returns:
        r_eci: (N_epoch, 3), v_eci: (N_epoch, 3), phi: (N_epoch, 6, 6)
    """
    K = len(a_emp_segments)
    N = len(t_epochs)

    r_arr = np.zeros((N, 3))
    v_arr = np.zeros((N, 3))
    phi_arr = np.zeros((N, 6, 6))

    r = np.asarray(r0, dtype=float).copy()
    v = np.asarray(v0, dtype=float).copy()
    t_cur = 0.0
    obs_idx = 0

    for k in range(K):
        a_rtn = np.asarray(a_emp_segments[k], dtype=float)
        t_seg_end = (k + 1) * segment_dt
        t_end_k = min(t_seg_end, max(t_epochs) + 30.0)

        # Integrate this segment
        integ = integrate_orbit_eci_with_stm(
            r, v, (0, t_end_k - t_cur), force_fn,
            Cd=2.2, CR=1.3,
            area_drag=0.68, area_srp=3.4, mass=580.0,
            empirical_acc_rtn=a_rtn,
            param_names=['Cd', 'CR', 'aR', 'aT', 'aN'],
            dt=10.0,
            mjd_utc=mjd_utc_start + t_cur / SEC_PER_DAY,
            mjd_tt=mjd_tt_start + t_cur / SEC_PER_DAY,
            bodies=['Sun', 'Moon'],
        )

        # Interpolate observation epochs within this segment
        t_integ = integ['t'] + t_cur
        while obs_idx < N and t_epochs[obs_idx] <= t_end_k + 10.0:
            t_obs = t_epochs[obs_idx]
            if t_obs >= t_cur - 1.0 and t_obs <= t_integ[-1] + 1.0:
                idx = np.argmin(np.abs(t_integ - t_obs))
                r_arr[obs_idx] = integ['r'][idx]
                v_arr[obs_idx] = integ['v'][idx]
                phi_arr[obs_idx] = integ['phi'][idx]
            obs_idx += 1

        # Update initial conditions for next segment
        r = integ['r'][-1].copy()
        v = integ['v'][-1].copy()
        t_cur = t_end_k

    return r_arr, v_arr, phi_arr


class BatchOrbitLSQ:
    """Full-arc batch least-squares orbit determination.

    Inner loop: BatchLinearSolver (clock, zwd, amb) — linear
    Outer loop: Gauss-Newton (r0, v0, a_emp per 5-min segment) — nonlinear via FD

    Empirical accelerations are piecewise-constant over K 5-minute segments.
    """

    def __init__(self, pass1_geometry, force_fn, t_epochs,
                 mjd_utc_start, mjd_tt_start,
                 sigma_phase=0.20, sigma_code=0.30,
                 max_iter=6, dx_tol_pos=0.01, dx_tol_vel=0.001,
                 damping=10.0, segment_minutes=5):
        self.geometry = pass1_geometry
        self.force_fn = force_fn
        self.t_epochs = np.asarray(t_epochs, dtype=float)
        self.mjd_utc_start = mjd_utc_start
        self.mjd_tt_start = mjd_tt_start
        self.sigma_phase = sigma_phase
        self.sigma_code = sigma_code
        self.max_iter = max_iter
        self.dx_tol_pos = dx_tol_pos
        self.dx_tol_vel = dx_tol_vel
        self.damping = damping

        self.N_epoch = len(pass1_geometry)
        t_max = max(t_epochs) - min(t_epochs) + 30.0
        self.segment_dt = segment_minutes * 60.0
        self.K_seg = max(1, int(np.ceil(t_max / self.segment_dt)))
        self.segment_dt = t_max / self.K_seg
        self.N_nl = 6 + 3 * self.K_seg  # r0(3), v0(3), a_rtn(3*K)
        print(f"  [BatchLSQ] {self.N_nl} NL params: r0+v0 + {self.K_seg}x3 emp "
              f"({self.segment_dt/60:.0f}min seg)")
        self._index_svs()

    def _index_svs(self):
        """Index SVs across all epochs."""
        sv_set = set()
        for ep_list in self.geometry:
            for d in ep_list:
                sv_set.add(d['sv'])
        self.sv_list = sorted(sv_set)
        self.sv_to_idx = {sv: i for i, sv in enumerate(self.sv_list)}

    def _forward_full(self, x_nl):
        """Full forward model: orbit propagation + linear solver.

        Args:
            x_nl: [r0(3), v0(3), aR1,aT1,aN1, ..., aRK,aTK,aNK]

        Returns:
            residuals_norm: (N_obs,) normalized observation residuals
            geo_epochs: list of lists with updated _geo_full for each dict
        """
        r0 = x_nl[0:3]
        v0 = x_nl[3:6]
        a_emp_segments = x_nl[6:].reshape(self.K_seg, 3)

        # Propagate orbit with piecewise-constant empirical RTN
        r_eci, v_eci, _ = _propagate_and_interp_segments(
            r0, v0, a_emp_segments, self.segment_dt, self.t_epochs,
            self.force_fn, self.mjd_utc_start, self.mjd_tt_start)

        # Convert to ECEF and compute geometric range per observation
        geo_epochs = []
        for i_ep, ep_list in enumerate(self.geometry):
            mjd_utc = self.mjd_utc_start + self.t_epochs[i_ep] / SEC_PER_DAY
            r_ecef, _ = eci_to_ecef(r_eci[i_ep], v_eci[i_ep], mjd_utc)
            new_ep = []
            for d in ep_list:
                dd = dict(d)  # shallow copy
                sat_ecef = np.asarray(d['sat_pos'], dtype=float)
                sat_clk = float(d.get('sat_clk', 0))
                rho_ecef = np.linalg.norm(sat_ecef - r_ecef)
                sag = (OMEGA_E / C_LIGHT) * (sat_ecef[0] * r_ecef[1]
                                              - sat_ecef[1] * r_ecef[0])
                dd['_geo_full'] = rho_ecef + sag - sat_clk
                new_ep.append(dd)
            geo_epochs.append(new_ep)

        # Run batch linear solver with updated geometry
        bls = BatchLinearSolver(geo_epochs,
                                 sigma_phase=self.sigma_phase,
                                 sigma_code=self.sigma_code)

        x_lin, res_norm, _ = bls._solve_internal()

        return np.array(res_norm), geo_epochs, bls, x_lin

    def solve(self, r0_prior, v0_prior):
        """Full Gauss-Newton batch LSQ.

        Args:
            r0_prior, v0_prior: initial ECI state guess

        Returns:
            dict with: r_eci, v_eci, clk, zwd, amb, amb_dict, rms, converged, iterations
        """
        x_nl = np.zeros(self.N_nl)
        x_nl[0:3] = np.asarray(r0_prior, dtype=float)
        x_nl[3:6] = np.asarray(v0_prior, dtype=float)
        # empirical accels start at 0

        cost_history = []
        bls = None
        x_lin = None

        # Perturb sizes for FD Jacobian
        perturbs = np.zeros(self.N_nl)
        perturbs[0:3] = 0.1       # pos: 0.1 m
        perturbs[3:6] = 0.001     # vel: 0.001 m/s
        perturbs[6:] = 1e-9       # emp: 1e-9 m/s²

        for iteration in range(self.max_iter):
            # --- Forward model ---
            res_norm, geo_epochs, bls, x_lin = self._forward_full(x_nl)
            cost = 0.5 * np.sum(res_norm**2)
            cost_history.append(cost)

            # --- FD Jacobian ---
            J = np.zeros((len(res_norm), self.N_nl))

            for i in range(self.N_nl):
                x_pert = x_nl.copy()
                x_pert[i] += perturbs[i]
                r_pert, _, _, _ = self._forward_full(x_pert)
                J[:, i] = (r_pert - res_norm) / perturbs[i]

            # --- Gauss-Newton ---
            H = J.T @ J
            g = J.T @ res_norm

            if self.damping > 0:
                H += self.damping * np.diag(np.diag(H) + 1e-6)

            try:
                dx = -np.linalg.solve(H, g)
            except np.linalg.LinAlgError:
                dx = -np.linalg.lstsq(H, g, rcond=1e-6)[0]

            x_nl += dx

            dr_norm = np.linalg.norm(dx[0:3])
            dv_norm = np.linalg.norm(dx[3:6])
            da_norm = np.linalg.norm(dx[6:9])

            converged = dr_norm < self.dx_tol_pos and dv_norm < self.dx_tol_vel
            tag = " ✓" if converged and iteration >= 2 else ""
            print(f"  [GN {iteration}] cost={cost:.1f} dr={dr_norm:.4f}m "
                  f"dv={dv_norm:.6f}m/s da={da_norm:.2e}{tag}")

            if converged and iteration >= 2:
                break

            self.damping = max(1e-4, self.damping * 0.4)

        # --- Final solution ---
        _, geo_final, bls, x_lin = self._forward_full(x_nl)
        clk = x_lin[0:self.N_epoch]
        zwd = x_lin[self.N_epoch:2*self.N_epoch]
        amb = x_lin[2*self.N_epoch:2*self.N_epoch + len(self.sv_list)]

        amb_dict = {sv: float(amb[i]) for sv, i in self.sv_to_idx.items()}

        # Extract orbit trajectory
        r0 = x_nl[0:3]
        v0 = x_nl[3:6]
        a_emp_segments = x_nl[6:].reshape(self.K_seg, 3)
        r_eci, v_eci, phi = _propagate_and_interp_segments(
            r0, v0, a_emp_segments, self.segment_dt, self.t_epochs,
            self.force_fn, self.mjd_utc_start, self.mjd_tt_start)

        return {
            'r_eci': r_eci,
            'v_eci': v_eci,
            'r0': r0,
            'v0': v0,
            'a_emp': a_rtn,
            'clk': clk,
            'zwd': zwd,
            'amb_dict': amb_dict,
            'amb': amb,
            'sv_list': self.sv_list,
            'cost_history': cost_history,
            'converged': iteration < self.max_iter - 1,
            'iterations': iteration + 1,
            'rms_phase': bls._compute_rms(clk, zwd, amb)[0] if bls else 0,
            'rms_code': bls._compute_rms(clk, zwd, amb)[1] if bls else 0,
        }
