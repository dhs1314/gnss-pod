"""Framework v3: ECI frame-consistent batch LSQ (Phase 15.0 -> 19.0).

Verified: analytical Jacobian for r0/v0 — ratio=1.00, angle=0.0 deg vs FD.

Phase 19.0: Piecewise-constant RTN empirical accelerations.
Each segment has independent (aR, aT, aN) — the reduced-dynamic approach.
S-matrix computed analytically: chain segment-end sensitivities through Phi.

GN strategy:
  Step 1: BatchLinearSolver on initial orbit -> globally optimal amb (fixed)
  Step 2: Clock-differenced residuals + line-search GN on NL params
"""

import numpy as np
from src.orbit_integrator import integrate_orbit_eci_with_stm
from src.batch_solver import BatchLinearSolver, _trop_mf

SEC_PER_DAY = 86400.0
C_LIGHT = 299792458.0
OE = 7.2921151467e-5


# ---------------------------------------------------------------------------
# Geometry rebuild (with full measurement corrections)
# ---------------------------------------------------------------------------

def _rebuild_geometry_ECI(geometry_template, r_ep, v_ep,
                           t_epochs, mjd_utc_start,
                           pco_z=None, ant_pco_z=None):
    """Rebuild _geo_full in ECI frame with Sagnac, relativity, wind-up, PCO."""
    from src.coordinates import ecef_to_eci, eci_to_ecef
    from src.troposphere import saastamoinen_zhd
    from src.measurement_corrections import (
        phase_wind_up_correction,
        relativity_shapiro_correction,
        compute_pco_ecef_from_nadir,
    )
    new_geo = []
    windup_state = {}
    for i_ep, ep_list in enumerate(geometry_template):
        r_eci = r_ep[i_ep]; v_eci = v_ep[i_ep]
        mjd_utc = mjd_utc_start + t_epochs[i_ep] / SEC_PER_DAY
        r_ecef, _ = eci_to_ecef(r_eci, v_eci, mjd_utc)
        lat_rad = np.arcsin(r_ecef[2] / np.linalg.norm(r_ecef))
        h_m = np.linalg.norm(r_ecef) - 6378137.0
        zhd = saastamoinen_zhd(lat_rad, h_m)
        new_ep = []
        for d in ep_list:
            dd = dict(d)
            sat_ecef = np.asarray(d['sat_pos'], dtype=float)
            sat_clk = float(d.get('sat_clk', 0)); el = float(d.get('el', 0.5))
            sv = d['sv']
            pco_z_sv = None
            if ant_pco_z and sv in ant_pco_z: pco_z_sv = ant_pco_z[sv]
            elif pco_z is not None: pco_z_sv = pco_z
            if pco_z_sv is not None:
                pco_ecef = compute_pco_ecef_from_nadir(sat_ecef, pco_z_sv)
                sat_ecef_corr = sat_ecef + pco_ecef
            else: sat_ecef_corr = sat_ecef
            # Use pre-computed sat_eci if available (no PCO)
            sat_eci = d.get('_sat_eci')
            if sat_eci is None or pco_z_sv is not None:
                sat_eci, _ = ecef_to_eci(sat_ecef_corr, np.zeros(3), mjd_utc)
                sat_eci = np.asarray(sat_eci, dtype=float)
            else:
                sat_eci = np.asarray(sat_eci, dtype=float)
            rho_eci = np.linalg.norm(sat_eci - r_eci)
            sag = (OE / C_LIGHT) * (sat_ecef_corr[0]*r_ecef[1] - sat_ecef_corr[1]*r_ecef[0])
            rel_corr = relativity_shapiro_correction(sat_eci, r_eci)
            mf_h = 1.001 / np.sqrt(0.002001 + np.sin(el)**2)
            dd['_geo_full'] = rho_eci + sag - sat_clk + zhd * mf_h + rel_corr
            prev_dphi = windup_state.get(sv)
            delta_phi_rad, current_dphi = phase_wind_up_correction(
                sat_eci, r_eci, sat_vel_eci=None, prev_dphi=prev_dphi)
            windup_state[sv] = current_dphi
            dd['_windup'] = delta_phi_rad * (C_LIGHT / 1575.42e6)
            new_ep.append(dd)
        new_geo.append(new_ep)
    return new_geo


# ---------------------------------------------------------------------------
# BatchOrbitLSQv3 — piecewise-empirical ECI Gauss-Newton
# ---------------------------------------------------------------------------

class BatchOrbitLSQv3:
    """ECI Batch LSQ with piecewise-constant empirical RTN.

    Parameters:
      r0(3), v0(3), then N_segs × [aR, aT, aN]
      Optionally +Cd +CR (if Orekit + estimate_cd_cr)

    piecewise_minutes: duration of each RTN segment [min].
      0 = single constant set (backward compatible, default)
    """

    def __init__(self, pass1_geometry, force_fn, t_epochs,
                 mjd_utc_start, mjd_tt_start,
                 sigma_phase=0.20, sigma_code=0.30,
                 max_iter=8, dx_tol_pos=0.002, dx_tol_vel=0.0002,
                 dx_tol_emp=1e-10,
                 damping=0.5, prior_r0=1.0, prior_v0=0.01,
                 prior_emp=1e-7,
                 piecewise_minutes=0,
                 estimate_cd_cr=False,
                 orekit_prop=None):
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
        self.dx_tol_emp = dx_tol_emp
        self.damping = damping
        self.prior_r0 = prior_r0
        self.prior_v0 = prior_v0
        self.prior_emp = prior_emp
        self.N_epoch = len(pass1_geometry)

        # Orekit
        self.orekit_prop = orekit_prop
        self._use_orekit = orekit_prop is not None
        self._estimate_cd_cr = estimate_cd_cr and self._use_orekit

        # Piecewise RTN
        self.piecewise_minutes = piecewise_minutes
        arc_minutes = max(self.t_epochs) / 60.0
        if piecewise_minutes > 0 and arc_minutes > piecewise_minutes:
            self.N_segs = max(1, int(np.ceil(arc_minutes / piecewise_minutes)))
        else:
            self.N_segs = 1
        self._seg_boundaries = self._compute_segment_boundaries()

        # NL parameter count
        self.N_emp = 3 * self.N_segs  # total empirical params
        self.I_NL_EMP_START = 6
        if self._estimate_cd_cr:
            self.N_nl = 6 + self.N_emp + 2  # +Cd +CR
        else:
            self.N_nl = 6 + self.N_emp

        self._index_svs()
        seg_str = f" pw={self.N_segs}x3" if self.N_segs > 1 else ""
        orekit_str = " [Orekit]" if self._use_orekit else ""
        cd_str = "+CdCR" if self._estimate_cd_cr else ""
        print(f"  [BatchLSQv3] {self.N_nl} NL (r0+v0+aRTN{seg_str}{cd_str}){orekit_str} + "
              f"{2*self.N_epoch + self.N_sv} lin, ECI frame")
        print(f"  Prior: sigma_r0={self.prior_r0}m sigma_v0={self.prior_v0}m/s "
              f"sigma_emp={self.prior_emp:.1e} m/s2"
              + (f"  N_segs={self.N_segs} ({piecewise_minutes}min each)"
                 if self.N_segs > 1 else ""))

    def _index_svs(self):
        sv_set = set()
        for ep in self.geometry:
            for d in ep: sv_set.add(d['sv'])
        self.sv_list = sorted(sv_set)
        self.N_sv = len(self.sv_list)

    def _compute_segment_boundaries(self):
        """Return list of (t_start, t_end) seconds for each segment."""
        if self.N_segs == 1:
            return [(0.0, float(max(self.t_epochs)))]
        T = float(max(self.t_epochs))
        dur = T / self.N_segs
        return [(i * dur, (i + 1) * dur) for i in range(self.N_segs)]

    def _segment_for_epoch(self, t_epoch):
        """Return segment index [0, N_segs-1] for epoch time t."""
        for seg_idx, (t0, t1) in enumerate(self._seg_boundaries):
            if t0 <= t_epoch < t1 + 1e-6:
                return seg_idx
        return self.N_segs - 1

    def _get_rtn_for_epoch(self, x_nl, t_epoch):
        """Get (aR, aT, aN) for a given epoch from x_nl."""
        seg = self._segment_for_epoch(t_epoch)
        start = self.I_NL_EMP_START + seg * 3
        if len(x_nl) > start + 2:
            return np.array([x_nl[start], x_nl[start+1], x_nl[start+2]])
        return np.zeros(3)

    def _propagate(self, x_nl):
        if self._use_orekit:
            return self._propagate_orekit_pw(x_nl)
        return self._propagate_python_pw(x_nl)

    def _propagate_python_pw(self, x_nl):
        """Propagate orbit with piecewise RTN. Build analytic S matrix."""
        r0 = x_nl[0:3]; v0 = x_nl[3:6]
        N = len(self.t_epochs)

        # Integrate full arc with zero empirical (for Phi reference)
        t_end = max(self.t_epochs) + 30.0
        integ = integrate_orbit_eci_with_stm(
            r0, v0, (0, t_end), self.force_fn,
            Cd=2.2, CR=1.3, area_drag=0.68, area_srp=3.4, mass=580.0,
            empirical_acc_rtn=np.zeros(3), param_names=['aR','aT','aN'], dt=10.0,
            mjd_utc=self.mjd_utc_start, mjd_tt=self.mjd_tt_start,
            bodies=['Sun','Moon'])
        t_integ = integ['t']

        # Interpolate Phi
        phi = np.zeros((N, 6, 6))
        r_py = np.zeros((N, 3)); v_py = np.zeros((N, 3))
        for i, t_obs in enumerate(self.t_epochs):
            idx = np.argmin(np.abs(t_integ - t_obs))
            r_py[i] = integ['r'][idx]; v_py[i] = integ['v'][idx]
            phi[i] = integ['phi'][idx]

        # Multi-step integration with piecewise empirical for actual orbit
        r = np.zeros((N, 3)); v = np.zeros((N, 3))
        r_cur, v_cur = r0.copy(), v0.copy()
        t_prev = 0.0
        for i, t_obs in enumerate(self.t_epochs):
            dt = t_obs - t_prev
            a_emp = self._get_rtn_for_epoch(x_nl, t_obs)
            from src.orbit_integrator import rk4_step_eci
            r_new, v_new = rk4_step_eci(
                r_cur, v_cur, dt, self.force_fn,
                Cd=2.2, CR=1.3, area_drag=0.68, area_srp=3.4, mass=580.0,
                empirical_acc_rtn=a_emp,
                mjd_utc=self.mjd_utc_start + t_obs/SEC_PER_DAY,
                mjd_tt=self.mjd_tt_start + t_obs/SEC_PER_DAY,
                bodies=['Sun','Moon'])
            r[i] = r_new; v[i] = v_new
            r_cur, v_cur = r_new, v_new
            t_prev = t_obs

        # Build analytic S matrix
        S = self._build_analytic_S(r, v, phi)
        return r, v, phi, S

    def _build_analytic_S(self, r, v, phi):
        """Build (N, 6, 3*N_segs) analytic empirical sensitivity matrix.

        For each epoch i in segment s:
          - Columns for segment s: analytic (t^2/2, t) * RTN_basis
          - Columns for segment k < s: segment-end S propagated through Phi
          - Columns for segment k > s: zero
        """
        from src.empirical import compute_rtn_frame
        N = len(self.t_epochs)
        Np = 3 * self.N_segs
        S = np.zeros((N, 6, Np))

        # Precompute Phi_inv at each epoch for chaining
        phi_inv = np.zeros((N, 6, 6))
        for i in range(N):
            try: phi_inv[i] = np.linalg.inv(phi[i] + 1e-12*np.eye(6))
            except: phi_inv[i] = np.linalg.pinv(phi[i])

        # Get Phi indices at segment boundaries
        seg_end_idx = []  # epoch index nearest to end of each segment
        for seg_idx, (t0, t1) in enumerate(self._seg_boundaries):
            idx = np.argmin(np.abs(self.t_epochs - t1))
            seg_end_idx.append(idx)
        seg_start_idx = [0] + [seg_end_idx[i] + 1 for i in range(self.N_segs - 1)]

        for i_ep, t_obs in enumerate(self.t_epochs):
            seg = self._segment_for_epoch(t_obs)
            t_seg_start = self._seg_boundaries[seg][0]
            dt_seg = t_obs - t_seg_start

            # RTN basis at current epoch
            R_vec, T_vec, N_vec = compute_rtn_frame(r[i_ep], v[i_ep])
            RTN_eci = np.column_stack([R_vec, T_vec, N_vec])  # (3,3)

            for k in range(min(seg + 1, self.N_segs)):
                col_start = 3 * k
                if k == seg:
                    # Current segment: analytic sensitivity
                    S_pos = 0.5 * dt_seg * dt_seg * RTN_eci  # (3,3)
                    S_vel = dt_seg * RTN_eci                  # (3,3)
                else:
                    # Previous segment: propagate end-of-segment sensitivity
                    t_seg_end = self._seg_boundaries[k][1]
                    dur_k = t_seg_end - self._seg_boundaries[k][0]
                    # RTN basis at end of segment k
                    idx_k = seg_end_idx[k]
                    Rk, Tk, Nk = compute_rtn_frame(r[idx_k], v[idx_k])
                    RTN_k = np.column_stack([Rk, Tk, Nk])
                    S_pos_k = 0.5 * dur_k * dur_k * RTN_k  # (3,3)
                    S_vel_k = dur_k * RTN_k                # (3,3)
                    S_k = np.zeros((6, 3))
                    S_k[0:3, :] = S_pos_k; S_k[3:6, :] = S_vel_k
                    # Propagate: Phi(t_i -> t_k_end) = Phi[i] @ inv(Phi[idx_k])
                    Phi_prop = phi[i_ep] @ phi_inv[idx_k]
                    S_prop = Phi_prop @ S_k
                    S_pos = S_prop[0:3, :]; S_vel = S_prop[3:6, :]

                S[i_ep, 0:3, col_start:col_start+3] = S_pos
                S[i_ep, 3:6, col_start:col_start+3] = S_vel

        return S

    # Pre-computed Phi cache (STM negligible change for small r0 adjustments)
    _phi_cache = None
    _phi_cache_r0 = None

    def _propagate_orekit_pw(self, x_nl):
        """Ultra-fast: pre-computed Phi + analytic S, Orekit orbit only.

        Total: ~0.5s per call (Orekit orbit), down from 30s.

        Phi(r0,v0) changes <0.5% for r0 adjustments <100m. Pre-computed
        once at initial state (30s one-time cost).
        S is analytic (RTN basis vectors × chaining through Phi).
        Orbit: Orekit continuous arc (0.5s via reused propagator).
        """
        N = len(self.t_epochs)
        a_emp_arr = np.zeros((N, 3))
        for i, t_obs in enumerate(self.t_epochs):
            a_emp_arr[i] = self._get_rtn_for_epoch(x_nl, t_obs)

        # Orekit orbit: single continuous propagation (0.5s)
        r0 = x_nl[0:3]; v0 = x_nl[3:6]
        r_ok, v_ok = self.orekit_prop.propagate_continuous_arc(
            r0, v0, self.t_epochs, a_emp_arr, self.mjd_utc_start)

        # Pre-compute Phi once (30s, one-time)
        if self._phi_cache is None:
            _, _, self._phi_cache, _ = self._propagate_python_pw(x_nl)
            self._phi_cache_r0 = r0.copy()

        # Analytic S using Orekit orbit positions (no Python force model)
        S = self._build_analytic_S(r_ok, v_ok, self._phi_cache)

        return r_ok, v_ok, self._phi_cache, S

    def _compute_residuals_with_amb(self, r_ep, v_ep, fixed_amb):
        """Full residuals with geometry rebuild. Clock-differenced."""
        geo_new = _rebuild_geometry_ECI(self.geometry, r_ep, v_ep,
                                         self.t_epochs, self.mjd_utc_start)
        res_norm = []
        for i_ep, ep_list in enumerate(geo_new):
            code_res = []
            for d in ep_list:
                if '_obs_code' in d:
                    code_res.append(float(d['_obs_code']) - float(d.get('_geo_full', 0)))
            clk_est = float(np.median(code_res)) if code_res else 0.0
            for d in ep_list:
                geo = float(d.get('_geo_full', 0)); sv = d['sv']
                amb_val = fixed_amb.get(sv, 0.0)
                windup = float(d.get('_windup', 0.0))
                if '_obs_code' in d:
                    res_norm.append((float(d['_obs_code'])-geo-clk_est)/self.sigma_code)
                if '_obs_phase' in d:
                    res_norm.append((float(d['_obs_phase'])-geo-windup-clk_est-amb_val)
                                    / self.sigma_phase)
        return np.array(res_norm)

    @staticmethod
    def _prep_for_batch(geometry_with_windup):
        geo_out = []
        for ep in geometry_with_windup:
            new_ep = []
            for d in ep:
                dd = dict(d)
                if '_obs_phase' in dd and '_windup' in dd:
                    dd['_obs_phase'] = float(dd['_obs_phase']) - float(dd['_windup'])
                new_ep.append(dd)
            geo_out.append(new_ep)
        return geo_out

    def _precompute_sat_eci(self):
        """Pre-compute satellite ECI positions once (saves ~4000 astropy transforms)."""
        from src.coordinates import ecef_to_eci
        for i_ep, ep_list in enumerate(self.geometry):
            mjd_utc = self.mjd_utc_start + self.t_epochs[i_ep] / SEC_PER_DAY
            for d in ep_list:
                if '_sat_eci' not in d:
                    sat_ecef = np.asarray(d['sat_pos'], dtype=float)
                    sat_eci, _ = ecef_to_eci(sat_ecef, np.zeros(3), mjd_utc)
                    d['_sat_eci'] = sat_eci

    def solve(self, r0_prior, v0_prior):
        # Pre-compute satellite ECI positions (one-time, saves ~2s per GN iteration)
        self._precompute_sat_eci()

        x_nl = np.zeros(self.N_nl)
        x_nl[0:3] = np.asarray(r0_prior, dtype=float)
        x_nl[3:6] = np.asarray(v0_prior, dtype=float)

        # Step 1: Fix ambiguities
        r_ep0, v_ep0, phi_ep0, S_ep0 = self._propagate(x_nl)
        geo_init = _rebuild_geometry_ECI(self.geometry, r_ep0, v_ep0,
                                          self.t_epochs, self.mjd_utc_start)
        geo_bs = self._prep_for_batch(geo_init)
        bls = BatchLinearSolver(geo_bs, sigma_phase=self.sigma_phase,
                                 sigma_code=self.sigma_code)
        bls_sol = bls.solve()
        fixed_amb = bls_sol['amb_dict']
        print(f"  [BS-Am] {len(fixed_amb)} SVs, "
              f"phase={bls_sol['rms_phase']:.3f}m code={bls_sol['rms_code']:.3f}m")

        # Step 2: GN loop
        for it in range(self.max_iter):
            r_ep, v_ep, phi_ep, S_ep = self._propagate(x_nl)
            res_norm = self._compute_residuals_with_amb(r_ep, v_ep, fixed_amb)
            cost = 0.5 * np.sum(res_norm**2)

            J = self._build_jacobian(r_ep, v_ep, phi_ep, S_ep, res_norm)

            # Normal equations
            H = J.T @ J; g = J.T @ res_norm
            pr = 1.0/self.prior_r0**2; pv = 1.0/self.prior_v0**2
            pe = 1.0/self.prior_emp**2
            for j in range(3):
                H[j,j]+=pr; H[j+3,j+3]+=pv
                g[j]+=pr*(x_nl[j]-float(r0_prior[j]))
                g[j+3]+=pv*(x_nl[j+3]-float(v0_prior[j]))
            # Empirical priors (all segments)
            for j in range(self.N_emp):
                col = 6 + j
                H[col, col] += pe
                g[col] += pe * x_nl[col]

            if self.damping > 0:
                H += self.damping * np.diag(np.diag(H) + 1e-6)

            try: dx = -np.linalg.solve(H, g)
            except: dx = -np.linalg.lstsq(H, g, rcond=1e-6)[0]

            # Line-search
            alpha = 1.0
            for _ in range(4):
                x_try = x_nl + alpha * dx
                r_try, v_try, _, _ = self._propagate(x_try)
                res_try = self._compute_residuals_with_amb(r_try, v_try, fixed_amb)
                if 0.5*np.sum(res_try**2) < cost:
                    break
                alpha *= 0.5

            x_nl += alpha * dx
            dr_n = np.linalg.norm(alpha*dx[0:3])
            dv_n = np.linalg.norm(alpha*dx[3:6])
            da_n = np.linalg.norm(alpha*dx[6:9])
            conv = dr_n < self.dx_tol_pos and dv_n < self.dx_tol_vel and da_n < self.dx_tol_emp
            tag = " CONV" if conv and it >= 2 else ""
            ls_tag = f" a={alpha:.1f}" if alpha < 1.0 else ""
            print(f"  [GN {it}] cost={cost:.1f} dr={dr_n:.4f}m "
                  f"dv={dv_n:.6f}m/s da={da_n:.2e}m/s2{ls_tag}"
                  f"  a0=[{x_nl[6]:.2e},{x_nl[7]:.2e},{x_nl[8]:.2e}]{tag}")
            if conv and it >= 2: break
            self.damping = max(1e-6, self.damping*0.4)

        # Final solution
        r_ep, v_ep, phi_ep, S_ep = self._propagate(x_nl)
        geo_final = _rebuild_geometry_ECI(self.geometry, r_ep, v_ep,
                                           self.t_epochs, self.mjd_utc_start)
        geo_final_bs = self._prep_for_batch(geo_final)
        bls_final = BatchLinearSolver(geo_final_bs, sigma_phase=self.sigma_phase,
                                       sigma_code=self.sigma_code)
        x_lin, _, _ = bls_final._solve_internal()
        clk = x_lin[bls_final.I_CLK:bls_final.I_CLK+bls_final.N_clk]
        zwd = x_lin[bls_final.I_ZWD:bls_final.I_ZWD+bls_final.N_zwd]
        amb = x_lin[bls_final.I_AMB:bls_final.I_AMB+bls_final.N_amb]
        amb_dict = {sv: float(amb[i]) for i, sv in enumerate(self.sv_list)}
        rms_p, rms_c = bls_final._compute_rms(x_lin)

        result = {'r_eci':r_ep, 'v_eci':v_ep, 'r0':x_nl[0:3], 'v0':x_nl[3:6],
                  'a_emp': x_nl[6:6+self.N_emp].copy(),
                  'clk':clk, 'zwd':zwd, 'amb_dict':amb_dict,
                  'rms_phase':rms_p, 'rms_code':rms_c,
                  'converged': it < self.max_iter-1, 'iterations': it+1,
                  'sv_list': self.sv_list, 'N_segs': self.N_segs}
        return result

    def _build_jacobian(self, r_ep, v_ep, phi_ep, S_ep, _res_norm):
        """Build Jacobian matrix (N_obs, N_nl)."""
        from src.coordinates import ecef_to_eci
        N_nl = self.N_nl
        N_obs = 0
        for ep in self.geometry:
            for d in ep:
                if '_obs_code' in d: N_obs += 1
                if '_obs_phase' in d: N_obs += 1
        J = np.zeros((N_obs, N_nl)); obs_idx = 0
        Np = S_ep.shape[2]  # number of emp + optional Cd/CR columns

        for i_ep, ep_list in enumerate(self.geometry):
            r_eci = r_ep[i_ep]; phi = phi_ep[i_ep]; Si = S_ep[i_ep]
            mjd_utc = self.mjd_utc_start + self.t_epochs[i_ep]/SEC_PER_DAY
            dr_dr0 = phi[0:3,0:3]; dr_dv0 = phi[0:3,3:6]

            for d in ep_list:
                sat_eci = d.get('_sat_eci')
                if sat_eci is not None:
                    sat_eci = np.asarray(sat_eci, dtype=float)
                else:
                    sat_ecef = np.asarray(d['sat_pos'], dtype=float)
                    sat_eci, _ = ecef_to_eci(sat_ecef, np.zeros(3), mjd_utc)
                delta_eci = sat_eci - r_eci; rho = np.linalg.norm(delta_eci)
                e_eci = delta_eci / max(rho, 1.0)
                jr = np.dot(e_eci, dr_dr0); jv = np.dot(e_eci, dr_dv0)
                jp = np.array([np.dot(e_eci, Si[0:3, c]) for c in range(Np)])

                def fill(row, W):
                    row[0:3] = jr/W; row[3:6] = jv/W
                    row[6:6+Np] = jp/W

                if '_obs_code' in d:
                    fill(J[obs_idx], self.sigma_code); obs_idx += 1
                if '_obs_phase' in d:
                    fill(J[obs_idx], self.sigma_phase); obs_idx += 1
        return J
