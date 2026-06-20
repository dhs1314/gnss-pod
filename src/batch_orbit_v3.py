"""Framework v3: ECI frame-consistent batch LSQ (Phase 15.0).

Verified: analytical Jacobian for r0/v0 — ratio=1.00, angle=0.0° vs FD.

GN strategy (Steps 1-2 from plan):
  Step 1: BatchLinearSolver on EKF orbit → globally optimal amb (fixed)
  Step 2: Clock-differenced residuals in GN loop → orbit-only gradient

With amb fixed and clock differenced, the GN sees ~pure orbit error,
not orbit+clock+amb compound.
"""

import numpy as np
from src.orbit_integrator import integrate_orbit_eci_with_stm
from src.batch_solver import BatchLinearSolver, _trop_mf

SEC_PER_DAY = 86400.0


def _rebuild_geometry_ECI(geometry_template, r_ep, v_ep,
                           t_epochs, mjd_utc_start):
    """Rebuild _geo_full in ECI frame WITH Sagnac correction."""
    from src.coordinates import ecef_to_eci, eci_to_ecef
    C = 299792458.0; OE = 7.2921151467e-5
    new_geo = []
    for i_ep, ep_list in enumerate(geometry_template):
        r_eci = r_ep[i_ep]; v_eci = v_ep[i_ep]
        mjd_utc = mjd_utc_start + t_epochs[i_ep] / SEC_PER_DAY
        r_ecef, _ = eci_to_ecef(r_eci, v_eci, mjd_utc)
        lat_rad = np.arcsin(r_ecef[2] / np.linalg.norm(r_ecef))
        h_m = np.linalg.norm(r_ecef) - 6378137.0
        from src.troposphere import saastamoinen_zhd
        zhd = saastamoinen_zhd(lat_rad, h_m)
        new_ep = []
        for d in ep_list:
            dd = dict(d)
            sat_ecef = np.asarray(d['sat_pos'], dtype=float)
            sat_clk = float(d.get('sat_clk', 0)); el = float(d.get('el', 0.5))
            sat_eci, _ = ecef_to_eci(sat_ecef, np.zeros(3), mjd_utc)
            rho_eci = np.linalg.norm(sat_eci - r_eci)
            sag = (OE/C) * (sat_ecef[0]*r_ecef[1] - sat_ecef[1]*r_ecef[0])
            mf_h = 1.001 / np.sqrt(0.002001 + np.sin(el)**2)
            dd['_geo_full'] = rho_eci + sag - sat_clk + zhd * mf_h
            new_ep.append(dd)
        new_geo.append(new_ep)
    return new_geo


def _build_obs_jacobian_ECI(r_ep, v_ep, phi_ep, geometry,
                             t_epochs, mjd_utc_start,
                             sigma_phase, sigma_code):
    """Orbit Jacobian in ECI — only r0(3),v0(3). Verified vs FD."""
    from src.coordinates import ecef_to_eci
    N_nl = 6; N_obs = 0
    for ep in geometry:
        for d in ep:
            if '_obs_code' in d: N_obs += 1
            if '_obs_phase' in d: N_obs += 1
    J = np.zeros((N_obs, N_nl)); obs_idx = 0
    for i_ep, ep_list in enumerate(geometry):
        r_eci = r_ep[i_ep]; phi = phi_ep[i_ep]
        mjd_utc = mjd_utc_start + t_epochs[i_ep] / SEC_PER_DAY
        dr_dr0 = phi[0:3, 0:3]; dr_dv0 = phi[0:3, 3:6]
        for d in ep_list:
            sat_ecef = np.asarray(d['sat_pos'], dtype=float)
            sat_eci, _ = ecef_to_eci(sat_ecef, np.zeros(3), mjd_utc)
            delta_eci = sat_eci - r_eci; rho = np.linalg.norm(delta_eci)
            e_eci = delta_eci / max(rho, 1.0)
            j_r0 = +e_eci @ dr_dr0; j_v0 = +e_eci @ dr_dv0
            if '_obs_code' in d:
                J[obs_idx,0:3]=j_r0/sigma_code; J[obs_idx,3:6]=j_v0/sigma_code
                obs_idx+=1
            if '_obs_phase' in d:
                J[obs_idx,0:3]=j_r0/sigma_phase; J[obs_idx,3:6]=j_v0/sigma_phase
                obs_idx+=1
    return J


class BatchOrbitLSQv3:
    """ECI Batch LSQ: fix amb + clock-differenced GN."""

    def __init__(self, pass1_geometry, force_fn, t_epochs,
                 mjd_utc_start, mjd_tt_start,
                 sigma_phase=0.20, sigma_code=0.30,
                 max_iter=6, dx_tol_pos=0.002, dx_tol_vel=0.0002,
                 damping=0.5, prior_r0=1.0, prior_v0=0.01):
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
        self.prior_r0 = prior_r0
        self.prior_v0 = prior_v0
        self.N_epoch = len(pass1_geometry)
        self.N_nl = 6  # r0(3), v0(3)

        self._index_svs()
        print(f"  [BatchLSQv3] {self.N_nl} NL (r0+v0) + "
              f"{2*self.N_epoch + self.N_sv} lin, ECI frame")
        print(f"  Prior: sigma_r0={self.prior_r0}m sigma_v0={self.prior_v0}m/s")

    def _index_svs(self):
        sv_set = set()
        for ep in self.geometry:
            for d in ep:
                sv_set.add(d['sv'])
        self.sv_list = sorted(sv_set)
        self.N_sv = len(self.sv_list)

    def _propagate(self, x_nl):
        r0 = x_nl[0:3]; v0 = x_nl[3:6]
        t_end = max(self.t_epochs) + 30.0
        integ = integrate_orbit_eci_with_stm(
            r0, v0, (0, t_end), self.force_fn,
            Cd=2.2, CR=1.3, area_drag=0.68, area_srp=3.4, mass=580.0,
            empirical_acc_rtn=np.zeros(3), param_names=['aR','aT','aN'], dt=10.0,
            mjd_utc=self.mjd_utc_start, mjd_tt=self.mjd_tt_start,
            bodies=['Sun','Moon'])
        t_integ = integ['t']; N = len(self.t_epochs)
        r = np.zeros((N,3)); v = np.zeros((N,3)); phi = np.zeros((N,6,6))
        for i, t_obs in enumerate(self.t_epochs):
            idx = np.argmin(np.abs(t_integ - t_obs))
            r[i]=integ['r'][idx]; v[i]=integ['v'][idx]; phi[i]=integ['phi'][idx]
        return r, v, phi

    def _compute_residuals_with_amb(self, r_ep, v_ep, fixed_amb):
        """Full residuals: rebuild geometry (with Sagnac) + solve clock/zwd.

        Uses _rebuild_geometry_ECI (now Sagnac-corrected) so that GN orbit
        adjustment changes _geo_full and produces correct residuals.
        Clock-differencing removes common-mode clock absorption.
        """
        geo_new = _rebuild_geometry_ECI(self.geometry, r_ep, v_ep,
                                         self.t_epochs, self.mjd_utc_start)
        W_P = 1.0/self.sigma_code**2; W_L = 1.0/self.sigma_phase**2

        res_norm = []
        for i_ep, ep_list in enumerate(geo_new):
            # Median code residual → clock estimate
            code_res = []
            for d in ep_list:
                if '_obs_code' in d:
                    code_res.append(float(d['_obs_code']) - float(d.get('_geo_full',0)))
            clk_est = float(np.median(code_res)) if code_res else 0.0

            for d in ep_list:
                geo = float(d.get('_geo_full', 0)); sv = d['sv']
                amb_val = fixed_amb.get(sv, 0.0)
                if '_obs_code' in d:
                    res_norm.append((float(d['_obs_code'])-geo-clk_est)/self.sigma_code)
                if '_obs_phase' in d:
                    res_norm.append((float(d['_obs_phase'])-geo-clk_est-amb_val)/self.sigma_phase)

        return np.array(res_norm)

    def solve(self, r0_prior, v0_prior):
        x_nl = np.zeros(self.N_nl)
        x_nl[0:3] = np.asarray(r0_prior, dtype=float)
        x_nl[3:6] = np.asarray(v0_prior, dtype=float)

        # Step 1: Propagate initial orbit, rebuild geometry WITH Sagnac
        r0, v0, _ = self._propagate(x_nl)
        geo_init = _rebuild_geometry_ECI(self.geometry, r0, v0,
                                          self.t_epochs, self.mjd_utc_start)
        bls = BatchLinearSolver(geo_init, sigma_phase=self.sigma_phase,
                                 sigma_code=self.sigma_code)
        bls_sol = bls.solve()
        fixed_amb = bls_sol['amb_dict']
        print(f"  [BS-Am] {len(fixed_amb)} SVs, "
              f"phase={bls_sol['rms_phase']:.3f}m code={bls_sol['rms_code']:.3f}m")

        # Step 2: GN loop on r0/v0 with Sagnac-corrected geometry rebuild
        for it in range(self.max_iter):
            r_ep, v_ep, phi_ep = self._propagate(x_nl)

            # Residuals with fixed amb + clock differencing
            res_norm = self._compute_residuals_with_amb(r_ep, v_ep, fixed_amb)
            cost = 0.5 * np.sum(res_norm**2)

            # Jacobian
            J = _build_obs_jacobian_ECI(r_ep, v_ep, phi_ep,
                self.geometry, self.t_epochs, self.mjd_utc_start,
                self.sigma_phase, self.sigma_code)

            # Normal equations with prior
            H = J.T @ J; g = J.T @ res_norm
            pr = 1.0/self.prior_r0**2; pv = 1.0/self.prior_v0**2
            for j in range(3):
                H[j,j]+=pr; H[j+3,j+3]+=pv
                g[j]+=pr*(x_nl[j]-float(r0_prior[j]))
                g[j+3]+=pv*(x_nl[j+3]-float(v0_prior[j]))
            if self.damping > 0:
                H += self.damping * np.diag(np.diag(H) + 1e-6)

            try: dx = -np.linalg.solve(H, g)
            except: dx = -np.linalg.lstsq(H, g, rcond=1e-6)[0]

            x_nl += dx
            dr_n = np.linalg.norm(dx[0:3]); dv_n = np.linalg.norm(dx[3:6])
            conv = dr_n < self.dx_tol_pos and dv_n < self.dx_tol_vel
            tag = " conv" if conv and it>=2 else ""
            print(f"  [GN {it}] cost={cost:.1f} dr={dr_n:.4f}m "
                  f"dv={dv_n:.6f}m/s{tag}")
            if conv and it>=2: break
            self.damping = max(1e-6, self.damping*0.4)

        # Final solution
        r_ep, v_ep, phi_ep = self._propagate(x_nl)
        geo_final = _rebuild_geometry_ECI(self.geometry, r_ep, v_ep,
                                           self.t_epochs, self.mjd_utc_start)
        bls_final = BatchLinearSolver(geo_final, sigma_phase=self.sigma_phase,
                                       sigma_code=self.sigma_code)
        x_lin, _, _ = bls_final._solve_internal()
        clk=x_lin[0:self.N_epoch]; zwd=x_lin[self.N_epoch:2*self.N_epoch]
        amb=x_lin[2*self.N_epoch:2*self.N_epoch+self.N_sv]
        amb_dict={sv:float(amb[i]) for i,sv in enumerate(self.sv_list)}
        rms_p,rms_c=bls_final._compute_rms(clk,zwd,amb)
        return {'r_eci':r_ep,'v_eci':v_ep,'r0':x_nl[0:3],'v0':x_nl[3:6],
                'a_emp':np.zeros(3),'clk':clk,'zwd':zwd,'amb_dict':amb_dict,
                'rms_phase':rms_p,'rms_code':rms_c,
                'converged':it<self.max_iter-1,'iterations':it+1,
                'sv_list':self.sv_list}
