"""Sequential Extended Kalman Filter for reduced-dynamic POD.

Processes epochs one at a time with STM-based state and covariance
propagation. Avoids the S-matrix cross-segment instability that
plagues the batch LSQ approach.

State vector (ECI, dimension = 11 + N_sv):
    [r(3), v(3), aR, aT, aN, zwd, clk, amb_1..amb_Nsv]
"""
import numpy as np
from src.orbit_integrator import integrate_orbit_eci_with_stm
from src.orbit_dynamics import total_acc_eci
from src.coordinates import eci_to_ecef, ecef_to_eci
from src.troposphere import compute_troposphere, ecef_to_geodetic

C_LIGHT = 299792458.0
OMEGA_E = 7.2921151467e-5
SEC_PER_DAY = 86400.0

# State indices (before ambiguity block)
I_R = slice(0, 3)
I_V = slice(3, 6)
I_EMP = slice(6, 9)
I_ZWD = 9
I_CLK = 10
I_AMB_START = 11
N_BASE = 11


class EKFState:
    """Container for EKF state at one epoch."""

    def __init__(self, t, x, P, sv_list):
        self.t = t              # GPS seconds-of-day
        self.x = np.asarray(x, dtype=float)
        self.P = np.asarray(P, dtype=float)
        self.sv_list = list(sv_list)  # ordered list of SV IDs
        self._rebuild_sv_map()

    def _rebuild_sv_map(self):
        self.sv_to_idx = {sv: i for i, sv in enumerate(self.sv_list)}

    @property
    def r_eci(self): return self.x[I_R]

    @property
    def v_eci(self): return self.x[I_V]

    @property
    def a_rtn(self): return self.x[I_EMP]

    @property
    def zwd(self): return float(self.x[I_ZWD])

    @property
    def clk(self): return float(self.x[I_CLK])

    @property
    def n_sv(self): return len(self.sv_list)

    @property
    def n_state(self): return len(self.x)

    def amb(self, sv):
        """Get float ambiguity for SV."""
        if sv not in self.sv_to_idx:
            return 0.0
        return float(self.x[I_AMB_START + self.sv_to_idx[sv]])

    def set_amb(self, sv, value):
        if sv in self.sv_to_idx:
            self.x[I_AMB_START + self.sv_to_idx[sv]] = float(value)

    def amb_array(self):
        return self.x[I_AMB_START:]


class SequentialEKF:
    """Sequential Extended Kalman Filter for reduced-dynamic POD.

    Usage:
        ekf = SequentialEKF(config)
        state = ekf.initialize(r0_eci, v0_eci, mjd_utc, gps_sod)
        for epoch in epochs:
            if not first:
                state = ekf.predict(state, gps_sod, mjd_utc, mjd_tt)
            state = ekf.process_epoch(state, epoch_data, ...)
    """

    def __init__(self, config=None):
        cfg = config or {}
        # Dynamics
        self.Cd = cfg.get('Cd', 2.2)
        self.CR = cfg.get('CR', 1.3)
        self.area_drag = cfg.get('area_drag', 0.68)
        self.area_srp = cfg.get('area_srp', 3.4)
        self.mass = cfg.get('mass', 580.0)
        self.dt_integ = cfg.get('dt_integ', 10.0)
        self.bodies = cfg.get('bodies', ['Sun', 'Moon'])

        # Force model (set separately)
        self.Cnm = cfg.get('Cnm', None)
        self.Snm = cfg.get('Snm', None)
        self.GM_grav = cfg.get('GM_grav', None)
        self.R_grav = cfg.get('R_grav', None)
        self.gravity_nmax = cfg.get('gravity_nmax', 90)

        # Process noise
        self.sigma_acc = cfg.get('sigma_acc_process', 1e-7)    # m/s^2
        self.tau_emp = cfg.get('tau_emp', 600.0)                # s
        self.sigma_emp_ss = cfg.get('sigma_emp_ss', 1e-8)       # m/s^2
        self.sigma_zwd_rw = cfg.get('sigma_zwd_rw', 1e-9)      # m/√s (m^2/s)
        self.sigma_clk = cfg.get('sigma_clk', 1e5)              # m

        # Measurement noise
        self.sigma_phase = cfg.get('sigma_phase', 0.01)         # m
        self.sigma_code = cfg.get('sigma_code', 0.30)           # m
        self.R_phase = self.sigma_phase ** 2
        self.R_code = self.sigma_code ** 2

        # Outlier detection
        self.chi2_threshold = cfg.get('chi2_threshold', 10.828)  # χ²(1, α=0.001)

        # Initial covariance
        self.P0_pos = cfg.get('P0_pos', 100.0)          # m²
        self.P0_vel = cfg.get('P0_vel', 1.0)            # (m/s)²
        self.P0_emp = cfg.get('P0_emp', 1e-12)          # (m/s²)² (σ=1e-6)
        self.P0_zwd = cfg.get('P0_zwd', 0.25)           # m²
        self.P0_clk = cfg.get('P0_clk', 1e10)           # m²
        self.P0_amb = cfg.get('P0_amb', 100.0)          # m²

        # Elevation cutoff
        self.el_min = cfg.get('el_min', 0.087)  # rad (~5 deg)

        # Buffer for on-the-fly code bias estimation (accumulate per-SV)
        self._bias_buf = {}          # sv -> list of P_r values
        self._bias_buf_min = cfg.get('bias_buf_min', 5)  # epochs before commit

    def initialize(self, r0_eci, v0_eci, mjd_utc, gps_sod):
        """Create initial EKF state from GNV1B reference.

        Args:
            r0_eci, v0_eci: ECI position [m] and velocity [m/s]
            mjd_utc: MJD (UTC)
            gps_sod: GPS seconds of day

        Returns:
            EKFState with initial covariance
        """
        n_state = N_BASE  # no SVs yet
        x = np.zeros(n_state)
        x[I_R] = np.asarray(r0_eci)
        x[I_V] = np.asarray(v0_eci)
        x[I_EMP] = 0.0       # aR, aT, aN = 0
        x[I_ZWD] = 0.0       # zwd = 0
        x[I_CLK] = 0.0       # clock = 0

        P = np.zeros((n_state, n_state))
        P[0, 0] = P[1, 1] = P[2, 2] = self.P0_pos
        P[3, 3] = P[4, 4] = P[5, 5] = self.P0_vel
        P[6, 6] = P[7, 7] = P[8, 8] = self.P0_emp
        P[9, 9] = self.P0_zwd
        P[10, 10] = self.P0_clk

        return EKFState(gps_sod, x, P, [])

    def predict(self, state, t_next, mjd_utc_start, mjd_tt_start):
        """EKF time update: integrate orbit + propagate covariance.

        Args:
            state: EKFState at current epoch
            t_next: GPS seconds-of-day at next epoch
            mjd_utc_start: MJD(UTC) at current epoch
            mjd_tt_start: MJD(TT) at current epoch

        Returns:
            New EKFState after prediction
        """
        dt = t_next - state.t
        if dt <= 0:
            return state

        # Build force model closure
        Cnm, Snm = self.Cnm, self.Snm
        GM_grav, R_grav = self.GM_grav, self.R_grav
        Nmax = self.gravity_nmax
        bodies = self.bodies

        def force_model(pos_eci, vel_eci, **kwargs):
            return total_acc_eci(
                pos_eci, vel_eci,
                Cnm=Cnm, Snm=Snm, Nmax=Nmax,
                GM_gravity=GM_grav, R_gravity=R_grav,
                **kwargs,
            )

        # Integrate state with empirical RTN, get STM + S
        integ = integrate_orbit_eci_with_stm(
            state.r_eci, state.v_eci, (0.0, dt), force_model,
            Cd=self.Cd, CR=self.CR,
            area_drag=self.area_drag, area_srp=self.area_srp,
            mass=self.mass,
            empirical_acc_rtn=state.a_rtn,
            param_names=['aR', 'aT', 'aN'],
            dt=self.dt_integ,
            mjd_tt=mjd_tt_start, mjd_utc=mjd_utc_start,
            bodies=bodies,
        )

        r_new = integ['r'][-1].copy()
        v_new = integ['v'][-1].copy()
        Phi_6 = integ['phi'][-1]   # (6,6) STM
        S_6x3 = integ['S'][-1]     # (6,3) empirical sensitivity

        n_state = state.n_state
        n_sv = state.n_sv

        # -- Build augmented state transition matrix --
        Phi_aug = np.eye(n_state)

        # r,v block: STM from integrator
        Phi_aug[0:6, 0:6] = Phi_6

        # Empirical cross-coupling: S maps [aR,aT,aN] → [r,v]
        Phi_aug[0:6, 6:9] = S_6x3

        # Gauss-Markov decay for empirical RTN
        decay = float(np.exp(-dt / self.tau_emp))
        Phi_aug[6, 6] = decay
        Phi_aug[7, 7] = decay
        Phi_aug[8, 8] = decay

        # ZWD: identity (random walk — no deterministic change)
        Phi_aug[9, 9] = 1.0

        # Clock: random walk (continuous, not reset)
        Phi_aug[10, 10] = 1.0

        # Ambiguities: identity (constant)
        # Already eye from np.eye(n_state) initialization

        # -- Build process noise Q --
        Q = np.zeros((n_state, n_state))

        # Position/velocity: unmodeled acceleration mapped via Γ
        # Γ = [0.5*dt²*I, dt*I]^T  (3→6 mapping)
        q_acc = self.sigma_acc ** 2 * dt  # acceleration variance (scalar per axis)
        Gamma_rv = np.zeros((6, 3))
        Gamma_rv[0:3, :] = 0.5 * dt ** 2 * np.eye(3)
        Gamma_rv[3:6, :] = dt * np.eye(3)
        Q[0:6, 0:6] = Gamma_rv @ (q_acc * np.eye(3)) @ Gamma_rv.T

        # Empirical: Gauss-Markov process noise
        q_emp = self.sigma_emp_ss ** 2 * (1.0 - decay ** 2)
        Q[6, 6] = Q[7, 7] = Q[8, 8] = q_emp

        # ZWD: random walk
        Q[9, 9] = self.sigma_zwd_rw * dt

        # Clock: random walk (σ_rw=0.032 m/√s, stable TCXO)
        Q[10, 10] = 0.001 * dt

        # Ambiguities: zero process noise (constants)

        # -- Propagate covariance --
        P_pred = Phi_aug @ state.P @ Phi_aug.T + Q

        # -- Predict state --
        x_pred = np.zeros(n_state)
        x_pred[I_R] = r_new
        x_pred[I_V] = v_new
        x_pred[I_EMP] = state.a_rtn * decay   # GM decay
        x_pred[I_ZWD] = state.zwd              # no deterministic change
        x_pred[I_CLK] = state.clk              # clock: random walk (continuous)
        if n_sv > 0:
            x_pred[I_AMB_START:] = state.amb_array()

        # Clock covariance: propagated via Phi_aug (no reset)

        return EKFState(t_next, x_pred, P_pred, state.sv_list)

    def _ensure_sv(self, state, sv_id, amb_init=None):
        """Add SV to state if not already present.

        Args:
            state: current EKFState
            sv_id: satellite ID to add
            amb_init: initial ambiguity value (default: 0.0)
        """
        if sv_id in state.sv_to_idx:
            return state

        # Expand state vector and covariance
        n_old = state.n_state
        n_new = n_old + 1

        x_new = np.zeros(n_new)
        x_new[:n_old] = state.x
        x_new[-1] = float(amb_init) if amb_init is not None else 0.0

        P_new = np.zeros((n_new, n_new))
        P_new[:n_old, :n_old] = state.P
        P_new[-1, -1] = self.P0_amb

        sv_list_new = state.sv_list + [sv_id]
        return EKFState(state.t, x_new, P_new, sv_list_new)

    def _update_scalar(self, state, sv_id, obs, modeled, H, R_meas):
        """Scalar Kalman update with innovation testing.

        Returns (state, accepted, innovation, chi2_stat) or (state, False, nu, chi2)
        """
        # Innovation
        nu = obs - modeled

        # Innovation covariance S = H·P·H^T + R
        PHt = state.P @ H
        S = float(H @ PHt) + R_meas

        # Chi-square innovation test
        chi2 = nu ** 2 / max(S, 1e-12)

        if chi2 > self.chi2_threshold:
            return state, False, nu, chi2

        # Kalman gain
        K = PHt / S

        # Joseph-form covariance update (guarantees symmetry)
        I_KH = np.eye(state.n_state) - np.outer(K, H)
        P_new = I_KH @ state.P @ I_KH.T + np.outer(K, K) * R_meas

        # State update
        x_new = state.x + K * nu

        state_new = EKFState(state.t, x_new, P_new, state.sv_list)
        return state_new, True, nu, chi2

    def _compute_modeled(self, state, sv_id, sat_pos_eci, sat_pos_ecef,
                         sat_clk, rcv_pos_ecef, zhd, mf_h, mf_w, with_amb):
        """Compute modeled observation (range + clock + tropo + optional amb).

        Geometric range is computed in ECI to match the state frame.
        Sagnac correction uses ECEF positions.
        """
        rho = float(np.linalg.norm(sat_pos_eci - state.r_eci))
        sag = (OMEGA_E / C_LIGHT) * (sat_pos_ecef[0] * rcv_pos_ecef[1]
                                      - sat_pos_ecef[1] * rcv_pos_ecef[0])
        modeled = (rho + sag - sat_clk + state.clk
                   + zhd * mf_h + state.zwd * mf_w)
        if with_amb:
            modeled += state.amb(sv_id)
        return modeled

    def _build_H(self, state, sv_id, sat_pos_eci, mf_w, with_amb):
        """Build observation matrix H for phase or code update."""
        n_state = state.n_state
        e_vec_eci = sat_pos_eci - state.r_eci
        e_vec_eci = e_vec_eci / float(np.linalg.norm(e_vec_eci))

        H = np.zeros(n_state)
        H[I_R] = -e_vec_eci
        H[I_ZWD] = mf_w
        H[I_CLK] = 1.0
        if with_amb and sv_id in state.sv_to_idx:
            H[I_AMB_START + state.sv_to_idx[sv_id]] = 1.0
        return H

    def process_epoch(self, state, epoch_data, sp3, sv_bias, sv_bias_ref,
                      mjd_utc, mjd_tt, doy):
        """Process all observations at one epoch.

        Args:
            state: EKFState before this epoch's updates
            epoch_data: list of per-SV dicts with 'sv', 'sat_pos', 'sat_clk',
                        'L_if_raw', 'P_if_raw', 'el'
            sp3: SP3 data (for future use)
            sv_bias: per-SV code bias dict (zero-meaned, mutated with new SVs)
            sv_bias_ref: reference mean subtracted from biases
            mjd_utc: MJD(UTC) for this epoch
            mjd_tt: MJD(TT) for this epoch
            doy: day of year

        Returns:
            (state_new, stats_dict)
            stats: {'n_phase': int, 'n_code': int, 'n_rej': int,
                    'rms_phase': float, 'rms_code': float}
        """
        rcv_pos_ecef, _ = eci_to_ecef(state.r_eci, state.v_eci, mjd_utc)
        lat_rad, lon_rad, h_m = ecef_to_geodetic(rcv_pos_ecef)

        # Compute ZHD once per epoch (same for all SVs)
        from src.troposphere import saastamoinen_zhd
        zhd = saastamoinen_zhd(lat_rad, h_m)

        n_phase, n_code, n_rej = 0, 0, 0
        sum_phase_nu2, sum_code_nu2 = 0.0, 0.0

        # Sort by elevation (highest first) — process best geometry first
        ed_sorted = sorted(epoch_data, key=lambda d: d.get('el', 0), reverse=True)

        from src.troposphere import gmf_mapping

        for d in ed_sorted:
            sv = d['sv']
            if sv not in sv_bias:
                # On-the-fly code bias: remove receiver clock via filter estimate.
                # P_r = P_if + sat_clk - rho_corr = clk_rcv + B_sat - B_rcv.
                # Subtracting state.clk (≈ clk_rcv - B_rcv + <B_sat>) gives
                # B_sat - <B_sat>, consistent with zero-meaned initial biases.
                if 'P_if_raw' in d:
                    sat_pos_ecef_tmp = d['sat_pos']
                    rho_tmp = float(np.linalg.norm(sat_pos_ecef_tmp - rcv_pos_ecef))
                    sag_tmp = (OMEGA_E / C_LIGHT) * (sat_pos_ecef_tmp[0] * rcv_pos_ecef[1]
                                                      - sat_pos_ecef_tmp[1] * rcv_pos_ecef[0])
                    sv_bias[sv] = float(d['P_if_raw']) + d['sat_clk'] - (rho_tmp + sag_tmp + state.clk)
                else:
                    continue

            sat_pos_ecef = d['sat_pos']   # ECEF (original SP3)
            sat_clk = d['sat_clk']
            el = d['el']
            if el < self.el_min:
                continue

            # Convert satellite position to ECI at reception epoch
            # Using same epoch for both keeps rotation invariant: |R·sat - R·rcv| = |sat - rcv|
            # Sagnac correction then gives the correct ECI path length.
            sat_pos_eci, _ = ecef_to_eci(sat_pos_ecef, np.zeros(3), mjd_utc)

            # Compute mapping functions
            mf_h, mf_w = gmf_mapping(el, lat_rad, lon_rad, h_m, doy)

            # Ensure SV is in the state (code or phase may add it)
            amb_init = None
            if sv not in state.sv_to_idx and 'L_if_raw' in d and 'P_if_raw' in d:
                amb_init = float(d['L_if_raw']) - float(d['P_if_raw'])

            # Code update FIRST (anchors clock, no ambiguity term)
            if 'P_if_raw' in d:
                obs_code = float(d['P_if_raw']) - sv_bias[sv]
                if 1.5e7 <= abs(obs_code) <= 3.5e7:
                    state = self._ensure_sv(state, sv, amb_init=amb_init)
                    modeled = self._compute_modeled(
                        state, sv, sat_pos_eci, sat_pos_ecef, sat_clk,
                        rcv_pos_ecef, zhd, mf_h, mf_w, with_amb=False)
                    H = self._build_H(state, sv, sat_pos_eci, mf_w, with_amb=False)
                    state, accepted, nu, chi2 = self._update_scalar(
                        state, sv, obs_code, modeled, H, self.R_code)
                    if accepted:
                        n_code += 1
                        sum_code_nu2 += nu ** 2
                    else:
                        n_rej += 1

            # Phase update SECOND (ambiguity absorbs phase bias)
            if 'L_if_raw' in d:
                obs_phase = float(d['L_if_raw']) - sv_bias[sv]
                if 1.5e7 <= abs(obs_phase) <= 3.5e7:
                    state = self._ensure_sv(state, sv, amb_init=amb_init)
                    modeled = self._compute_modeled(
                        state, sv, sat_pos_eci, sat_pos_ecef, sat_clk,
                        rcv_pos_ecef, zhd, mf_h, mf_w, with_amb=True)
                    H = self._build_H(state, sv, sat_pos_eci, mf_w, with_amb=True)
                    state, accepted, nu, chi2 = self._update_scalar(
                        state, sv, obs_phase, modeled, H, self.R_phase)
                    if accepted:
                        n_phase += 1
                        sum_phase_nu2 += nu ** 2
                    else:
                        n_rej += 1

        rms_phase = np.sqrt(sum_phase_nu2 / max(n_phase, 1))
        rms_code = np.sqrt(sum_code_nu2 / max(n_code, 1))

        stats = {
            'n_phase': n_phase, 'n_code': n_code, 'n_rej': n_rej,
            'rms_phase': rms_phase, 'rms_code': rms_code,
        }
        return state, stats
