"""Full-arc batch least-squares orbit determination (Phase 11.0 / V2.3.0).

Uses the existing STM integrator for analytical orbit Jacobian:
  Φ = ∂(r,v)/∂(r0,v0)    [6×6]
  S = ∂(r,v)/∂(aR,aT,aN)  [6×3]  from _param_partials_eci

1 integration per GN iteration. Jacobian computed analytically from STM.

Key: both geometry and Jacobian use ECEF range (same as EKF pass 1).
  ρ = |sat_ecef - rcv_ecef| + Sagnac - sat_clk
  ∂ρ/∂r0_eci = -e_ecef · R_ecef2eci · Φ[0:3, 0:3]
For short arcs (<0.5h), R_ecef2eci ≈ I (Earth rot <0.01 rad).
"""

import numpy as np
from src.orbit_integrator import integrate_orbit_eci_with_stm
from src.coordinates import eci_to_ecef
from src.batch_solver import BatchLinearSolver

SEC_PER_DAY = 86400.0
C_LIGHT = 299792458.0
OMEGA_E = 7.2921151467e-5


def _rebuild_geometry_ECEF(geometry_template, r_ep, v_ep,
                            t_epochs, mjd_utc_start):
    """Rebuild _geo_full using ECEF range (matches EKF pass 1).

    geo = |sat_ecef - rcv_ecef| + Sagnac - sat_clk + ZHD*mf_h
    """
    new_geo = []
    for i_ep, ep_list in enumerate(geometry_template):
        mjd_utc = mjd_utc_start + t_epochs[i_ep] / SEC_PER_DAY
        r_ecef, _ = eci_to_ecef(r_ep[i_ep], v_ep[i_ep], mjd_utc)

        lat_rad = np.arcsin(r_ecef[2] / np.linalg.norm(r_ecef))
        h_m = np.linalg.norm(r_ecef) - 6378137.0
        from src.troposphere import saastamoinen_zhd
        zhd = saastamoinen_zhd(lat_rad, h_m)

        new_ep = []
        for d in ep_list:
            dd = dict(d)
            sat_ecef = np.asarray(d['sat_pos'], dtype=float)
            sat_clk = float(d.get('sat_clk', 0))
            el = float(d.get('el', 0.5))

            rho = np.linalg.norm(sat_ecef - r_ecef)
            sag = (OMEGA_E / C_LIGHT) * (sat_ecef[0] * r_ecef[1]
                                          - sat_ecef[1] * r_ecef[0])
            mf_h = 1.001 / np.sqrt(0.002001 + np.sin(el)**2)
            dd['_geo_full'] = rho + sag - sat_clk + zhd * mf_h
            new_ep.append(dd)
        new_geo.append(new_ep)
    return new_geo


def _build_obs_jacobian(r_ep, v_ep, phi_ep, S_ep,
                         geometry, t_epochs, mjd_utc_start,
                         sigma_phase, sigma_code):
    """Build observation Jacobian w.r.t (r0, v0, aR, aT, aN).

    Uses ECEF range derivative: e_ecef = (sat_ecef - rcv_ecef)/ρ
    Chain rule: drho/d(r0_eci) = -e_ecef · ∂(r_ecef)/∂(r0_eci)
    ≈ -e_ecef · Φ[0:3, 0:3]  (R_ecef2eci ≈ I for short arcs)
    """
    N_nl = 9
    N_obs = 0
    for ep in geometry:
        for d in ep:
            if '_obs_code' in d: N_obs += 1
            if '_obs_phase' in d: N_obs += 1

    J = np.zeros((N_obs, N_nl))
    obs_idx = 0

    for i_ep, ep_list in enumerate(geometry):
        mjd_utc = mjd_utc_start + t_epochs[i_ep] / SEC_PER_DAY
        r_ecef, _ = eci_to_ecef(r_ep[i_ep], v_ep[i_ep], mjd_utc)
        phi = phi_ep[i_ep]
        S   = S_ep[i_ep]

        # STM position sub-blocks (ECI→ECEF via approximate identity)
        dr_dr0 = phi[0:3, 0:3]  # ∂r_eci/∂r0_eci ≈ ∂r_ecef/∂r0_eci
        dr_dv0 = phi[0:3, 3:6]
        dr_daR = S[0:3, 0]      # aR at param_names[2]=aR → S col 0
        dr_daT = S[0:3, 1]      # aT at S col 1
        dr_daN = S[0:3, 2]      # aN at S col 2

        for d in ep_list:
            sat_ecef = np.asarray(d['sat_pos'], dtype=float)
            delta = sat_ecef - r_ecef
            rho = np.linalg.norm(delta)
            if rho < 1.0:
                e_vec = np.array([0.0, 0.0, 1.0])
            else:
                e_vec = delta / rho

            j_r0 = -e_vec @ dr_dr0   # (3,) @ (3,3) → (3,)
            j_v0 = -e_vec @ dr_dv0
            j_aR = -e_vec @ dr_daR
            j_aT = -e_vec @ dr_daT
            j_aN = -e_vec @ dr_daN

            if '_obs_code' in d:
                J[obs_idx, 0:3] = j_r0 / sigma_code
                J[obs_idx, 3:6] = j_v0 / sigma_code
                J[obs_idx, 6]   = j_aR / sigma_code
                J[obs_idx, 7]   = j_aT / sigma_code
                J[obs_idx, 8]   = j_aN / sigma_code
                obs_idx += 1

            if '_obs_phase' in d:
                J[obs_idx, 0:3] = j_r0 / sigma_phase
                J[obs_idx, 3:6] = j_v0 / sigma_phase
                J[obs_idx, 6]   = j_aR / sigma_phase
                J[obs_idx, 7]   = j_aT / sigma_phase
                J[obs_idx, 8]   = j_aN / sigma_phase
                obs_idx += 1

    return J


class BatchOrbitLSQv2:
    """Full-arc batch least-squares with analytical STM Jacobian.

    1 integration per GN iteration. STM provides all partials.
    K=1 empirical segment (constant aR, aT, aN over full arc).
    """

    def __init__(self, pass1_geometry, force_fn, t_epochs,
                 mjd_utc_start, mjd_tt_start,
                 sigma_phase=0.20, sigma_code=0.30,
                 max_iter=6, dx_tol_pos=0.002, dx_tol_vel=0.0002,
                 damping=1.0):
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
        self.N_nl = 9
        self._index_svs()
        print(f"  [BatchLSQv2] {self.N_nl} NL + "
              f"{2*self.N_epoch + self.N_sv} lin "
              f"({self.N_epoch} ep, {self.N_sv} SV)")

    def _index_svs(self):
        sv_set = set()
        for ep in self.geometry:
            for d in ep:
                sv_set.add(d['sv'])
        self.sv_list = sorted(sv_set)
        self.N_sv = len(self.sv_list)

    def _propagate(self, x_nl):
        r0 = x_nl[0:3]; v0 = x_nl[3:6]
        aR, aT, aN = float(x_nl[6]), float(x_nl[7]), float(x_nl[8])
        t_end = max(self.t_epochs) + 30.0

        integ = integrate_orbit_eci_with_stm(
            r0, v0, (0, t_end), self.force_fn,
            Cd=2.2, CR=1.3,
            area_drag=0.68, area_srp=3.4, mass=580.0,
            empirical_acc_rtn=np.array([aR, aT, aN]),
            param_names=['aR', 'aT', 'aN'],
            dt=10.0,
            mjd_utc=self.mjd_utc_start,
            mjd_tt=self.mjd_tt_start,
            bodies=['Sun', 'Moon'],
        )

        t_integ = integ['t']
        N = len(self.t_epochs)
        r = np.zeros((N, 3)); v = np.zeros((N, 3))
        phi = np.zeros((N, 6, 6)); S = np.zeros((N, 6, 3))

        for i, t_obs in enumerate(self.t_epochs):
            idx = np.argmin(np.abs(t_integ - t_obs))
            r[i] = integ['r'][idx]; v[i] = integ['v'][idx]
            phi[i] = integ['phi'][idx]
            if 'S' in integ and integ['S'].shape[1] > 0:
                Ns = integ['S'].shape[1]
                for j in range(min(Ns, 3)):
                    S[i, :, j] = integ['S'][idx, :, j]

        return r, v, phi, S

    def solve(self, r0_prior, v0_prior):
        x_nl = np.zeros(self.N_nl)
        x_nl[0:3] = np.asarray(r0_prior, dtype=float)
        x_nl[3:6] = np.asarray(v0_prior, dtype=float)

        for iteration in range(self.max_iter):
            # 1 integration
            r_ep, v_ep, phi_ep, S_ep = self._propagate(x_nl)

            # Update geometry
            geo_new = _rebuild_geometry_ECEF(self.geometry, r_ep, v_ep,
                                              self.t_epochs, self.mjd_utc_start)

            # Linear solver
            bls = BatchLinearSolver(geo_new,
                                     sigma_phase=self.sigma_phase,
                                     sigma_code=self.sigma_code)
            x_lin, res_norm, _ = bls._solve_internal()
            cost = 0.5 * np.sum(res_norm**2)

            # Analytical Jacobian
            J = _build_obs_jacobian(r_ep, v_ep, phi_ep, S_ep,
                                     self.geometry, self.t_epochs,
                                     self.mjd_utc_start,
                                     self.sigma_phase, self.sigma_code)

            # GN update (r0, v0, aR, aT, aN)
            H = J.T @ J
            g = J.T @ res_norm
            if self.damping > 0:
                H += self.damping * np.eye(self.N_nl) * 1e-6

            try:
                dx = -np.linalg.solve(H, g)
            except np.linalg.LinAlgError:
                dx = -np.linalg.lstsq(H, g, rcond=1e-6)[0]

            x_nl += dx
            dr_norm = np.linalg.norm(dx[0:3])
            dv_norm = np.linalg.norm(dx[3:6])
            da_norm = np.linalg.norm(dx[6:9])

            conv = dr_norm < self.dx_tol_pos and dv_norm < self.dx_tol_vel
            tag = " converged" if conv and iteration >= 2 else ""
            print(f"  [GN {iteration}] cost={cost:.1f} dr={dr_norm:.4f}m "
                  f"dv={dv_norm:.6f}m/s da={da_norm:.2e}{tag}")

            if conv and iteration >= 2:
                break
            self.damping = max(1e-4, self.damping * 0.4)

        # Final solution
        r_ep, v_ep, phi_ep, S_ep = self._propagate(x_nl)
        geo_final = _rebuild_geometry_ECEF(self.geometry, r_ep, v_ep,
                                            self.t_epochs, self.mjd_utc_start)
        bls = BatchLinearSolver(geo_final,
                                 sigma_phase=self.sigma_phase,
                                 sigma_code=self.sigma_code)
        x_lin, _, _ = bls._solve_internal()

        clk = x_lin[0:self.N_epoch]
        zwd = x_lin[self.N_epoch:2*self.N_epoch]
        amb = x_lin[2*self.N_epoch:2*self.N_epoch + self.N_sv]
        amb_dict = {sv: float(amb[i]) for i, sv in enumerate(self.sv_list)}
        rms_phase, rms_code = bls._compute_rms(clk, zwd, amb)

        return {
            'r_eci': r_ep, 'v_eci': v_ep,
            'r0': x_nl[0:3], 'v0': x_nl[3:6], 'a_emp': x_nl[6:9],
            'clk': clk, 'zwd': zwd, 'amb_dict': amb_dict,
            'rms_phase': rms_phase, 'rms_code': rms_code,
            'converged': iteration < self.max_iter - 1,
            'iterations': iteration + 1,
            'sv_list': self.sv_list,
        }
