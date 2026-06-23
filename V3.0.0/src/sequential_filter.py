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


class MWBuffer:
    """Accumulate smoothed MW per SV for wide-lane ambiguity fixing.

    The Melbourne-Wübbena combination eliminates geometry, clock, troposphere,
    and ionosphere (first order), leaving:

        MW = N_wl + b_r_wl - b_s_wl  [cycles]

    where b_r_wl is the receiver WL bias (common to all SVs) and b_s_wl is
    the satellite WL bias.  Without external OSB products we use
    self-calibration: estimate b_r_wl as the median fractional part across
    all tracked SVs, then N_w_fixed = round(MW_mean - b_r_wl).
    """

    def __init__(self, min_epochs=10, max_epochs=50):
        self.min_epochs = min_epochs
        self.max_epochs = max_epochs
        self._mw = {}       # sv → list of MW values [cycles]
        self._b_r_wl = None  # receiver WL bias estimate (None = not yet calibrated)

    def add(self, sv, L1_cyc, L2_cyc, P1_m, P2_m):
        from src.ambiguity import compute_mw
        mw = compute_mw(L1_cyc, L2_cyc, P1_m, P2_m)
        if sv not in self._mw:
            self._mw[sv] = []
        self._mw[sv].append(mw)
        if len(self._mw[sv]) > self.max_epochs:
            self._mw[sv].pop(0)

    def _estimate_receiver_bias(self):
        """Estimate b_r_wl from the median fractional part across SVs.

        For each SV with enough data, the fractional part of mean(MW) is
        b_r_wl - b_s_wl.  Satellite WL biases are small (<0.15 cyc for GPS)
        so the median across SVs isolates the receiver term.
        """
        fracs = []
        for sv, vals in self._mw.items():
            if len(vals) >= self.min_epochs:
                mean_mw = float(np.mean(vals))
                frac = mean_mw - round(mean_mw)  # fractional part in [-0.5, 0.5]
                fracs.append(frac)
        if len(fracs) >= 3:
            self._b_r_wl = float(np.median(fracs))
            print(f"    [MW-BIAS] b_r_wl={self._b_r_wl:+.4f} cyc from {len(fracs)} SVs "
                  f"(fracs={[f'{f:+.3f}' for f in sorted(fracs)[:8]]})")
            return self._b_r_wl
        return None

    def try_fix_wl(self, sv):
        """Return N_w_fixed (int) or None if not enough data / too noisy."""
        vals = self._mw.get(sv, [])
        if len(vals) < self.min_epochs:
            return None
        mean_mw = float(np.mean(vals))
        std_mw = float(np.std(vals))

        # Estimate receiver WL bias if not yet done
        if self._b_r_wl is None:
            self._estimate_receiver_bias()

        if self._b_r_wl is not None:
            # Bias-corrected mean
            mean_corrected = mean_mw - self._b_r_wl
            N_w = int(round(mean_corrected))
            residual = mean_corrected - N_w
            if abs(residual) > 0.35:
                return None
            print(f"    [MW-FIX] {sv}: N_w={N_w} (n={len(vals)} "
                  f"mean_corr={mean_corrected:+.3f} std={std_mw:.3f})")
            return N_w

        # No receiver bias estimate yet — fall back to raw rounding
        if std_mw > 0.35:
            return None
        N_w = int(round(mean_mw))
        print(f"    [MW-FIX] {sv}: N_w={N_w} (n={len(vals)} "
              f"mean={mean_mw:+.3f} std={std_mw:.3f}, no b_r)")
        return N_w

    def reset(self, sv):
        self._mw.pop(sv, None)

    def n_epochs(self, sv):
        return len(self._mw.get(sv, []))


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

    def remove_sv(self, sv):
        """Remove an SV from the state, shrinking x and P.

        Returns a new EKFState (immutable update).  No-op if sv not present.
        """
        if sv not in self.sv_to_idx:
            return self
        idx = self.sv_to_idx[sv]
        amb_idx = I_AMB_START + idx
        keep = np.ones(self.n_state, dtype=bool)
        keep[amb_idx] = False
        sv_list_new = [s for s in self.sv_list if s != sv]
        return EKFState(self.t, self.x[keep], self.P[keep][:, keep], sv_list_new)


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
        self.pulse_interval = cfg.get('pulse_interval', 0)       # epochs: 0=off, 12=6min
        self.pulse_amplify = cfg.get('pulse_amplify', 50.0)      # Q multiplier at pulse
        self._epoch_counter = 0

        # CLK1B receiver clock (Phase 16.0, used in batch solver only)
        self._clk1b = cfg.get('clk1b_data', None)
        self._clk1b_sigma = cfg.get('clk1b_sigma', 0.10)  # m
        self.sigma_zwd_rw = cfg.get('sigma_zwd_rw', 1e-9)       # m/√s (m^2/s)
        self.sigma_clk = cfg.get('sigma_clk', 1e5)              # m
        self._clock_rw = cfg.get('clock_rw', 0.001)             # m²/s per s (0.001→σ=0.032m/s)

        # Measurement noise (base values at zenith)
        self.sigma_phase = cfg.get('sigma_phase', 0.01)         # m (zenith)
        self.sigma_code = cfg.get('sigma_code', 0.30)           # m (zenith)
        self._elev_exp_phase = cfg.get('elev_exp_phase', 1.0)   # elev weighting exponent
        self._elev_exp_code = cfg.get('elev_exp_code', 0.70)    # milder for clock
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

        # SV pruning (remove SVs that haven't been observed for prune_timeout)
        self.prune_timeout = cfg.get('prune_timeout', 1800.0)  # s (30 min)
        self._sv_last_seen = {}  # sv -> gps_sod of last observation

        # Buffer for on-the-fly code bias estimation (accumulate per-SV)
        self._bias_buf = {}          # sv -> list of P_r values
        self._bias_buf_min = cfg.get('bias_buf_min', 5)  # epochs before commit

        # ── Measurement corrections (Phase 2.3) ──
        self.use_phase_windup = cfg.get('use_phase_windup', False)
        self.use_relativity = cfg.get('use_relativity', False)
        self.use_cycle_slip = cfg.get('use_cycle_slip', False)
        self.antex_data = cfg.get('antex_data', None)       # parsed ANTEX dict
        self.dcb_data = cfg.get('dcb_data', None)           # parsed DCB dict
        self._dcb_if = {}  # sv → IF DCB correction [m], pre-computed
        if self.dcb_data is not None:
            from src.precision_products import compute_dcb_if_correction
            # Pre-compute IF DCB corrections for all SVs
            for prn in self.dcb_data:
                self._dcb_if[prn] = compute_dcb_if_correction(self.dcb_data, prn)

        # Internal state for cumulative phase wind-up (per SV)
        self._windup_state = {}  # sv → accumulated wind-up angle [rad]

        # TurboEdit cycle slip detector (lazy init)
        self._turbo_edit = None

        # MW buffer for WL ambiguity fixing (Phase 3.0 PPP-AR)
        # Phase 12.0: max_epochs=200 = arc-level (effectively unlimited)
        _mw_max = cfg.get('mw_max_epochs', 200)
        self._mw_buf = MWBuffer(min_epochs=cfg.get('ar_min_epochs', 10),
                                 max_epochs=_mw_max)
        self._wl_fixed = {}  # sv → N_w_fixed (int)
        self._wl_epochs = {}  # sv → count of phase epochs since WL fix (progressive tightening)
        self._osb_wl = cfg.get('osb_wl', {})  # Phase 9.0: {sv: b_wl_sat} satellite WL bias
        self._osb_nl = cfg.get('osb_nl', {})  # Phase 9.0: {sv: b_nl_sat} satellite NL bias
        # Code-phase smoothing buffer (Phase 7.1): per-SV sliding window
        # of L_if-P_if values, used for median-based amb_init instead of
        # single-epoch code-phase.  Reduces σ from ~0.6m to ~0.35m.
        self._cp_buf = {}       # sv → list of float (L_if - P_if) [m], max 5 entries
        self._cp_buf_size = 3   # median of last 3 code-phase values
        self._nl_fixed = {}  # sv → ΔN1_fixed (int, relative to reference SV)
        self._nl_ref_sv = None  # reference SV for SD NL fixing

        # Batch-fixed ambiguities (Phase 6.0): pre-computed IF ambiguity [m] per SV.
        # When set, these override the code-phase initialization and lock P_amb.
        self._amb_batch_fixed = cfg.get('amb_batch_fixed', {})  # {sv: B_if_fixed [m]}
        self._amb_batch_var = cfg.get('amb_batch_var', 0.0001)  # fixed ambiguity variance [m²]

        # ── Dynamics mode (Phase 2.4) ──
        self.dynamics_mode = cfg.get('dynamics_mode', 'simplified')
        self._orekit_prop = None  # lazy-initialized

        if self.dynamics_mode == 'orekit':
            from src.orekit_bridge import create_propagator
            self._orekit_prop, self.dynamics_mode, _warn = create_propagator(
                mode='orekit',
                orekit_data_path=cfg.get('orekit_data_path',
                                         r'd:\prj\gnss_pod\data\orekit'),
                gravity_field=cfg.get('gravity_field',
                                      r'd:\prj\gnss_pod\data\GGM05C.gfc'),
                gravity_degree=cfg.get('orekit_gravity_degree', 150),
                solid_tides=cfg.get('orekit_solid_tides', True),
                ocean_tides=cfg.get('orekit_ocean_tides', True),
                ocean_tide_degree=cfg.get('orekit_ocean_tide_degree', 50),
                relativity=cfg.get('orekit_relativity', True),
                mass=cfg.get('mass', 580.0),
                area_drag=cfg.get('area_drag', 0.68),
                area_srp=cfg.get('area_srp', 3.4),
                CR=cfg.get('CR', 1.3),
                CD=cfg.get('CD', 2.2),
                stm_perturb=cfg.get('stm_perturb', 1.0),
            )

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

        self._epoch_counter += 1  # pulse counter

        # ── Orekit Dynamics Branch (Phase 2.4) ──
        if self.dynamics_mode == 'orekit' and self._orekit_prop is not None:
            try:
                r_new, v_new, Phi_6_orekit, S_6x3_orekit = (
                    self._orekit_prop.propagate(
                        state.r_eci, state.v_eci, state.a_rtn, dt,
                        mjd_utc_start, mjd_tt_start))
                if Phi_6_orekit is None:
                    # Orekit STM not yet implemented; fall through to simplified
                    pass
                else:
                    Phi_6 = Phi_6_orekit
                    S_6x3 = S_6x3_orekit
                    # Skip integration; jump to covariance propagation
                    n_state = state.n_state
                    # Gauss-Markov decay for empirical RTN
                    decay = float(np.exp(-dt / self.tau_emp))
                    a_emp_new = state.a_rtn * decay

                    Phi_aug = np.eye(n_state)
                    Phi_aug[0:6, 0:6] = Phi_6
                    Phi_aug[0:6, 6:9] = S_6x3
                    Phi_aug[6, 6] = decay
                    Phi_aug[7, 7] = decay
                    Phi_aug[8, 8] = decay
                    Phi_aug[9, 9] = 1.0
                    Phi_aug[10, 10] = 1.0
                    # Ambiguities already eye

                    Q = self._build_process_noise(n_state, dt, decay)
                    x_new_orekit = np.zeros(n_state)
                    x_new_orekit[I_R] = r_new
                    x_new_orekit[I_V] = v_new
                    x_new_orekit[I_EMP] = a_emp_new
                    x_new_orekit[I_ZWD] = state.zwd
                    x_new_orekit[I_CLK] = state.clk
                    x_new_orekit[I_AMB_START:] = state.x[I_AMB_START:]

                    P_pred = Phi_aug @ state.P @ Phi_aug.T + Q
                    return EKFState(t_next, x_new_orekit, P_pred, state.sv_list)
            except Exception as e:
                import traceback
                print(f"  [Orekit] Prediction failed: {e}")
                traceback.print_exc()
                print(f"  [Orekit] Falling back to simplified dynamics for this step")

        # ── Simplified Dynamics Branch (V2.2.1) ──
        # Build force model closure
        Cnm, Snm = self.Cnm, self.Snm
        GM_grav, R_grav = self.GM_grav, self.R_grav
        Nmax = self.gravity_nmax
        bodies = self.bodies

        # Compute solid Earth tide corrections once per EKF step (slow sun/moon
        # motion means these are effectively constant over 30s integration).
        tide_corr = {}
        try:
            from src.solid_tides import (compute_solid_tide_corrections,
                                         compute_time_varying_gravity,
                                         merge_tide_corrections)
            if not hasattr(self, '_tvgrav'):
                self._tvgrav = {}
            # Solid tides: recompute each step (changes slowly)
            tide_corr = compute_solid_tide_corrections(mjd_utc_start, mjd_tt_start)
            # Time-varying gravity: compute once, cache
            if not self._tvgrav:
                self._tvgrav = compute_time_varying_gravity(mjd_tt_start)
            tide_corr = merge_tide_corrections(tide_corr, self._tvgrav)
        except Exception:
            pass  # graceful fallback if solid_tides unavailable

        def force_model(pos_eci, vel_eci, **kwargs):
            return total_acc_eci(
                pos_eci, vel_eci,
                Cnm=Cnm, Snm=Snm, Nmax=Nmax,
                GM_gravity=GM_grav, R_gravity=R_grav,
                tide_corrections=tide_corr,
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
        Q = self._build_process_noise(n_state, dt, decay)

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

    def _build_process_noise(self, n_state, dt, decay):
        """Build process noise covariance matrix Q.

        Args:
            n_state: total state dimension
            dt: time step [s]
            decay: empirical acceleration decay factor exp(-dt/tau)

        Returns:
            Q: (n_state, n_state) process noise matrix
        """
        Q = np.zeros((n_state, n_state))

        # Position/velocity: unmodeled acceleration mapped via Γ
        q_acc = self.sigma_acc ** 2 * dt
        Gamma_rv = np.zeros((6, 3))
        Gamma_rv[0:3, :] = 0.5 * dt ** 2 * np.eye(3)
        Gamma_rv[3:6, :] = dt * np.eye(3)
        Q[0:6, 0:6] = Gamma_rv @ (q_acc * np.eye(3)) @ Gamma_rv.T

        # Empirical: Gauss-Markov process noise (amplified at pulse epochs)
        q_emp = self.sigma_emp_ss ** 2 * (1.0 - decay ** 2)
        if self.pulse_interval > 0 and self._epoch_counter % self.pulse_interval == 0:
            q_emp *= self.pulse_amplify
        Q[6, 6] = Q[7, 7] = Q[8, 8] = q_emp

        # ZWD: random walk
        Q[9, 9] = self.sigma_zwd_rw * dt

        # Clock: random walk
        Q[10, 10] = self._clock_rw * dt

        # Ambiguities: zero process noise (constants)
        return Q

    def _ensure_sv(self, state, sv_id, amb_init=None):
        """Add SV to state if not already present.

        If sv_id is in _amb_batch_fixed, uses the pre-computed value with
        very tight variance (Phase 6.0 batch-fixed ambiguity).
        """
        if sv_id in state.sv_to_idx:
            return state

        # Expand state vector and covariance
        n_old = state.n_state
        n_new = n_old + 1

        x_new = np.zeros(n_new)
        x_new[:n_old] = state.x

        P_new = np.zeros((n_new, n_new))
        P_new[:n_old, :n_old] = state.P

        if sv_id in self._amb_batch_fixed:
            # Batch-fixed absolute B_if (Phase 6.0 _amb_batch_fixed path)
            x_new[-1] = float(self._amb_batch_fixed[sv_id])
            P_new[-1, -1] = self._amb_batch_var
        elif sv_id in self._wl_fixed:
            # WL-fixed from batch AR or sequential MW: moderate initial P.
            # P_amb=0.10 (σ=0.32m) allows NL to float while WL is locked.
            # This avoids the 100 m² drift with raw code-phase while giving
            # the clock room to absorb the epoch-0 clock reference offset.
            x_new[-1] = float(amb_init) if amb_init is not None else 0.0
            P_new[-1, -1] = 0.10
        else:
            x_new[-1] = float(amb_init) if amb_init is not None else 0.0
            P_new[-1, -1] = self.P0_amb

        sv_list_new = state.sv_list + [sv_id]
        return EKFState(state.t, x_new, P_new, sv_list_new)

    @staticmethod
    def _meas_noise(sigma_zenith, el_rad, exponent=1.0):
        """Elevation-dependent measurement noise variance.

        σ²(el) = σ²_zenith / sin(el)^exponent
        Phase: exponent=1.0 (full elev dependence, higher precision at zenith)
        Code:  exponent=1.0 (mild, preserves clock constraint from all SVs)
        """
        sin_el = max(np.sin(el_rad), 0.1)
        return sigma_zenith**2 / sin_el ** exponent

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

    def _try_nl_sd_fix(self, state, ep_data):
        """Attempt NL ambiguity fixing via SD constraints.

        Phase 6.0: if _nl_batch is set (batch AR pass 2), applies pre-computed
        SD integers directly (no bootstrapping needed).
        Phase 5.0: falls back to real-time integer bootstrapping.

        Requires ≥3 WL-fixed SVs visible in the current epoch.
        """
        from src.ambiguity import COEFF_W

        LAM_NL = 0.107
        MIN_WL_EPOCHS = 10
        MIN_SD_SVS = 3
        MAX_RESIDUAL = 0.40

        # ── Phase 9.0: Undifferenced NL diagnostics (CODE OSB) ──
        # Log OSB-corrected NL integers without applying constraints.
        # This validates the OSB product by checking if N1 residuals converge.
        if self._osb_nl:
            for sv in self._wl_fixed:
                if sv not in self._osb_nl or sv not in state.sv_to_idx:
                    continue
                N_w = self._wl_fixed[sv]
                b_nl = self._osb_nl[sv]
                amb_idx = I_AMB_START + state.sv_to_idx[sv]
                B_if = float(state.x[amb_idx])
                N1_float = (B_if - COEFF_W * N_w) / LAM_NL - b_nl
                N1_int = int(round(N1_float))
                residual = abs(N1_float - N1_int)
                # Log at first WL fix for this SV
                sv_key = ('_osb_seen', sv)
                if not hasattr(self, '_osb_seen_set'):
                    self._osb_seen_set = set()
                if sv not in self._osb_seen_set:
                    self._osb_seen_set.add(sv)
                    print(f"    [OSB-DIAG] {sv}: N_w={N_w} N1={N1_int} "
                          f"resid={residual:.3f}cyc b_nl={b_nl:+.3f}cyc "
                          f"(wl_ep={self._wl_epochs.get(sv,0)})")
                if residual >= 0.25:
                    continue  # NL not converged yet — constraint would be premature
                # Apply mild OSB NL constraint (Phase 9.0)
                if residual < 0.25 and self._wl_epochs.get(sv, 0) >= 5:
                    target = LAM_NL * (N1_int + b_nl) + COEFF_W * N_w
                    H = np.zeros(state.n_state)
                    H[amb_idx] = 1.0
                    nu = target - B_if
                    PHt = state.P @ H
                    S = float(H @ PHt) + 0.01  # σ=0.1m
                    K = PHt / S
                    state.x = state.x + K * nu
                    state.P = state.P - np.outer(K, PHt)
                    self._nl_fixed[sv] = True
            return  # skip SD path when OSB applied

        # ── Phase 6.0: Batch NL SD constraints (clock-independent) ──
        if hasattr(self, '_nl_batch') and self._nl_batch:
            nl_ref = getattr(self, '_nl_batch_ref', None)
            n_applied = 0
            if nl_ref and nl_ref in self._wl_fixed and nl_ref in state.sv_to_idx:
                N_w_ref = self._wl_fixed[nl_ref]
                amb_idx_ref = I_AMB_START + state.sv_to_idx[nl_ref]
                R_c = 0.01  # moderate: clock-consistent batch SD

                for sv, dN1 in self._nl_batch.items():
                    if sv == nl_ref or sv not in self._wl_fixed:
                        continue
                    if sv not in state.sv_to_idx:
                        continue
                    # Require WL-fixed SV to have ≥3 epochs of phase
                    # (NL float needs basic convergence before constraint)
                    # Batch NL: pre-computed, apply immediately

                    N_w_k = self._wl_fixed[sv]
                    amb_idx_k = I_AMB_START + state.sv_to_idx[sv]

                    # Check SD consistency
                    sd_val = (float(state.x[amb_idx_k]) - float(state.x[amb_idx_ref])
                              - COEFF_W * (N_w_k - N_w_ref))
                    residual = abs(sd_val / LAM_NL - dN1)
                    if residual > 0.40:
                        continue

                    # Apply SD constraint (clock-independent)
                    rhs = LAM_NL * dN1 + COEFF_W * (N_w_k - N_w_ref)
                    H = np.zeros(state.n_state)
                    H[amb_idx_k] = 1.0
                    H[amb_idx_ref] = -1.0
                    nu = rhs - (float(state.x[amb_idx_k]) - float(state.x[amb_idx_ref]))
                    PHt = state.P @ H
                    S = float(H @ PHt) + R_c
                    K = PHt / S
                    state.x = state.x + K * nu
                    state.P = state.P - np.outer(K, PHt)
                    self._nl_fixed[sv] = True
                    n_applied += 1

            if n_applied > 0:
                return  # batch NL applied, skip real-time bootstrapping

        # ── Phase 5.0: Real-time integer bootstrapping ──

        # Find WL-fixed SVs visible in this epoch with sufficient convergence time
        wl_svs = {}
        for sv, N_w in self._wl_fixed.items():
            wl_ep = self._wl_epochs.get(sv, 0)
            if wl_ep < MIN_WL_EPOCHS:
                continue
            # Check SV is visible in this epoch
            for d in ep_data:
                if d.get('sv') == sv and sv in state.sv_to_idx:
                    wl_svs[sv] = d.get('el', 0.5)
                    break

        if len(wl_svs) < MIN_SD_SVS + 1:  # need ref + ≥3 others
            return

        # Choose reference SV (highest elevation, most stable)
        sv_ref = max(wl_svs, key=wl_svs.get)
        N_w_ref = self._wl_fixed[sv_ref]
        amb_idx_ref = I_AMB_START + state.sv_to_idx[sv_ref]

        # Build SD float vector and covariance (non-reference SVs, not yet NL-fixed)
        sd_svs = sorted(sv for sv in wl_svs if sv != sv_ref
                       and sv in state.sv_to_idx
                       and sv not in self._nl_fixed)
        if len(sd_svs) < MIN_SD_SVS:
            return

        n_sd = len(sd_svs)
        z_float = np.zeros(n_sd)  # SD NL float [cycles]
        Q_z = np.zeros((n_sd, n_sd))  # covariance in [cycles²]

        for i, sv in enumerate(sd_svs):
            amb_idx_i = I_AMB_START + state.sv_to_idx[sv]
            N_w_i = self._wl_fixed[sv]
            sd_m = float(state.x[amb_idx_i] - state.x[amb_idx_ref]
                         - COEFF_W * (N_w_i - N_w_ref))
            z_float[i] = sd_m / LAM_NL

            for j, sv2 in enumerate(sd_svs):
                amb_idx_j = I_AMB_START + state.sv_to_idx[sv2]
                Q_z[i, j] = (state.P[amb_idx_i, amb_idx_j]
                            - state.P[amb_idx_i, amb_idx_ref]
                            - state.P[amb_idx_ref, amb_idx_j]
                            + state.P[amb_idx_ref, amb_idx_ref]) / LAM_NL**2

        # Integer bootstrapping: sort by precision, sequentially round
        idx_order = np.argsort(np.diag(Q_z))  # most precise first
        z_fixed = np.zeros(n_sd)
        ok_svs = []

        for rank, idx in enumerate(idx_order):
            # Conditional float given previous fixes
            if rank == 0:
                z_cond = z_float[idx]
                sigma_cond = np.sqrt(Q_z[idx, idx])
            else:
                # Update via Gaussian conditioning on already-fixed components
                fixed_indices = idx_order[:rank]
                Q_ff = Q_z[np.ix_(fixed_indices, fixed_indices)]
                Q_fu = Q_z[np.ix_(fixed_indices, [idx])]
                Q_uf = Q_z[np.ix_([idx], fixed_indices)]
                Q_uu = Q_z[idx, idx]

                z_f = z_float[fixed_indices] - z_fixed[fixed_indices]
                try:
                    Q_ff_inv = np.linalg.inv(Q_ff)
                except np.linalg.LinAlgError:
                    continue

                z_cond = z_float[idx] - (Q_uf @ Q_ff_inv @ z_f)[0]
                sigma_cond = np.sqrt(max(Q_uu - (Q_uf @ Q_ff_inv @ Q_fu)[0, 0], 0.0))

            z_int = int(round(z_cond))
            residual = abs(z_cond - z_int)

            if residual > MAX_RESIDUAL or residual > 3.0 * sigma_cond:
                # This component can't be reliably fixed
                continue

            z_fixed[idx] = z_cond - z_int  # correction applied to float
            ok_svs.append(idx)

        if len(ok_svs) < MIN_SD_SVS:
            return

        # Apply constraints for fixed SD ambiguities
        R_constraint = 0.0001  # m²
        n_applied = 0

        for idx in ok_svs:
            sv = sd_svs[idx]
            dN1_int = int(round(z_float[idx] - z_fixed[idx]))
            N_w_i = self._wl_fixed[sv]
            rhs = LAM_NL * dN1_int + COEFF_W * (N_w_i - N_w_ref)

            amb_idx_i = I_AMB_START + state.sv_to_idx[sv]
            H = np.zeros(state.n_state)
            H[amb_idx_i] = 1.0
            H[amb_idx_ref] = -1.0

            nu = rhs - (float(state.x[amb_idx_i]) - float(state.x[amb_idx_ref]))
            PHt = state.P @ H
            S = float(H @ PHt) + R_constraint
            K = PHt / S
            state.x = state.x + K * nu
            state.P = state.P - np.outer(K, PHt)
            self._nl_fixed[sv] = True  # mark as NL-fixed
            n_applied += 1

        if n_applied > 0:
            n_total = len(self._nl_fixed)
            res_str = ' '.join(f'{r:.2f}' for r in z_float[idx_order[:len(ok_svs)]]
                              - np.round(z_float[idx_order[:len(ok_svs)]]))
            print(f"    [NL-FIX] {n_applied} SVs via SD bootstrapping "
                  f"(ref={sv_ref}, resid={res_str} cyc)")

    def _compute_modeled(self, state, sv_id, sat_pos_eci, sat_pos_ecef,
                         sat_clk, rcv_pos_ecef, zhd, mf_h, mf_w, with_amb,
                         rel_corr=0.0, windup_corr=0.0):
        """Compute modeled observation (range + clock + tropo + optional amb).

        Geometric range is computed in ECI to match the state frame.
        Sagnac correction uses ECEF positions.

        Args:
            rel_corr: relativistic Shapiro correction [m] to ADD
            windup_corr: phase wind-up correction [m] to ADD (phase only)
        """
        rho = float(np.linalg.norm(sat_pos_eci - state.r_eci))
        sag = (OMEGA_E / C_LIGHT) * (sat_pos_ecef[0] * rcv_pos_ecef[1]
                                      - sat_pos_ecef[1] * rcv_pos_ecef[0])
        modeled = (rho + sag - sat_clk + state.clk
                   + zhd * mf_h + state.zwd * mf_w
                   + rel_corr + windup_corr)
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

        Applies optional measurement corrections:
          - Satellite antenna PCO (from ANTEX)
          - Relativistic Shapiro delay
          - Phase wind-up
          - Cycle slip detection + ambiguity reset
        """
        rcv_pos_ecef, rcv_vel_ecef = eci_to_ecef(state.r_eci, state.v_eci, mjd_utc)
        lat_rad, lon_rad, h_m = ecef_to_geodetic(rcv_pos_ecef)

        # Compute ZHD once per epoch (same for all SVs)
        from src.troposphere import saastamoinen_zhd
        zhd = saastamoinen_zhd(lat_rad, h_m)

        # Lazy-init cycle slip detector
        if self.use_cycle_slip and self._turbo_edit is None:
            from src.cycle_slip import TurboEdit
            self._turbo_edit = TurboEdit(gf_threshold=0.08)

        n_phase, n_code, n_rej = 0, 0, 0
        sum_phase_nu2, sum_code_nu2 = 0.0, 0.0

        # Sort by elevation (highest first)
        ed_sorted = sorted(epoch_data, key=lambda d: d.get('el', 0), reverse=True)

        from src.troposphere import gmf_mapping

        for d in ed_sorted:
            sv = d['sv']
            dcb_corr = self._dcb_if.get(sv, 0.0)

            if sv not in sv_bias:
                if 'P_if_raw' in d:
                    sat_pos_ecef_tmp = np.asarray(d['sat_pos'], dtype=float)
                    # Apply PCO correction for consistency with measurement model
                    if self.antex_data is not None:
                        from src.measurement_corrections import compute_pco_ecef_from_nadir
                        from src.precision_products import get_satellite_pco
                        pco_neu_tmp = get_satellite_pco(self.antex_data, sv, 'L1')
                        pco_z_tmp = float(pco_neu_tmp[2])
                        sat_pos_ecef_tmp = sat_pos_ecef_tmp + compute_pco_ecef_from_nadir(
                            sat_pos_ecef_tmp, pco_z_tmp)
                    rho_tmp = float(np.linalg.norm(sat_pos_ecef_tmp - rcv_pos_ecef))
                    sag_tmp = (OMEGA_E / C_LIGHT) * (sat_pos_ecef_tmp[0] * rcv_pos_ecef[1]
                                                      - sat_pos_ecef_tmp[1] * rcv_pos_ecef[0])
                    sv_bias[sv] = float(d['P_if_raw']) + dcb_corr + d['sat_clk'] - (rho_tmp + sag_tmp + state.clk)
                else:
                    continue

            sat_pos_ecef = np.asarray(d['sat_pos'], dtype=float)
            sat_clk = d['sat_clk']
            el = d['el']
            if el < self.el_min:
                continue

            # ── Satellite Antenna PCO Correction ──
            if self.antex_data is not None:
                from src.measurement_corrections import compute_pco_ecef_from_nadir
                from src.precision_products import get_satellite_pco
                pco_neu = get_satellite_pco(self.antex_data, sv, 'L1')
                # PCO Z component (up = away from Earth = toward satellite -Z_body)
                pco_z = float(pco_neu[2])  # UP component in ANTEX
                pco_ecef = compute_pco_ecef_from_nadir(sat_pos_ecef, pco_z)
                sat_pos_ecef_corr = sat_pos_ecef + pco_ecef
            else:
                sat_pos_ecef_corr = sat_pos_ecef

            # Convert satellite position to ECI
            sat_pos_eci, sat_vel_eci = ecef_to_eci(sat_pos_ecef_corr, np.zeros(3), mjd_utc)

            # Mapping functions
            mf_h, mf_w = gmf_mapping(el, lat_rad, lon_rad, h_m, doy)

            # ── Relativistic Shapiro Correction ──
            rel_corr = 0.0
            if self.use_relativity:
                from src.measurement_corrections import relativity_shapiro_correction
                rel_corr = relativity_shapiro_correction(sat_pos_eci, state.r_eci)

            # ── Phase Wind-Up Correction ──
            windup_corr = 0.0
            if self.use_phase_windup:
                from src.measurement_corrections import phase_wind_up_correction
                prev_windup = self._windup_state.get(sv)
                delta_phi, current_phi = phase_wind_up_correction(
                    sat_pos_eci, state.r_eci, sat_vel_eci, prev_windup)
                self._windup_state[sv] = current_phi
                # Convert radians to metres (L1 wavelength)
                lambda_l1 = C_LIGHT / 1575.42e6
                windup_corr = delta_phi * lambda_l1

            # ── Cycle Slip Detection ──
            if self.use_cycle_slip and self._turbo_edit is not None:
                # Use raw L1/L2/P1/P2 for proper MW+GF detection.
                # TurboEdit requires raw dual-frequency data, not IF combinations.
                L1_raw = float(d.get('L1_raw', d.get('L_if_raw', 0)))
                L2_raw = float(d.get('L2_raw', d.get('L_if_raw', 0)))
                P1_raw = float(d.get('P1_raw', d.get('P_if_raw', 0)))
                P2_raw = float(d.get('P2_raw', d.get('P_if_raw', 0)))
                has_slip_lli = bool(d.get('slip_lli', False))
                has_slip, dN1, dN2, slip_mw = self._turbo_edit.detect_slip_L1(
                    sv, L1_raw, L2_raw, P1_raw, P2_raw)
                if (has_slip or has_slip_lli) and sv in state.sv_to_idx:
                    # MW slip → WL ambiguity changed, clear the fix
                    if slip_mw and sv in self._wl_fixed:
                        print(f"    [SLIP-MW] {sv}: WL fix cleared! dN1={dN1} dN2={dN2}")
                        self._wl_fixed.pop(sv, None)
                        self._wl_epochs.pop(sv, None)
                        self._nl_fixed.pop(sv, None)
                        if self._nl_ref_sv == sv:
                            self._nl_ref_sv = None
                            self._nl_fixed.clear()
                        self._mw_buf.reset(sv)
                        self._cp_buf.pop(sv, None)  # Phase 7.1: clear code-phase buffer
                        self._turbo_edit.reset_sv(sv)
                    elif not slip_mw and sv in self._wl_fixed:
                        # GF-only or LLI: N_w unchanged, keep WL fix.
                        # Only log substantial GF jumps (> 0.08 m L1-L2 change).
                        pass

                    # Reset ambiguity with WL-constrained re-init if WL is fixed
                    amb_new = float(d['L_if_raw']) - (float(d['P_if_raw']) + dcb_corr)
                    if sv in self._wl_fixed:
                        from src.ambiguity import compute_mw, COEFF_W
                        L1c = float(d.get('L1_cyc', 0))
                        L2c = float(d.get('L2_cyc', 0))
                        P1 = float(d.get('P1_raw', d.get('P_if_raw', 0)))
                        P2 = float(d.get('P2_raw', d.get('P_if_raw', 0)))
                        mw_float = compute_mw(L1c, L2c, P1, P2)
                        N_w_float = round(mw_float)
                        amb_new = amb_new + COEFF_W * (self._wl_fixed[sv] - N_w_float)
                    amb_idx = I_AMB_START + state.sv_to_idx[sv]
                    state.x[amb_idx] = amb_new
                    state.P[amb_idx, :] = 0.0
                    state.P[:, amb_idx] = 0.0
                    state.P[amb_idx, amb_idx] = self.P0_amb

            # Ensure SV is in the state
            amb_init = None
            if sv not in state.sv_to_idx and 'L_if_raw' in d and 'P_if_raw' in d:
                # ── Multi-epoch code-phase smoothed amb_base (Phase 7.1) ──
                # Use median of buffered L_if-P_if instead of single-epoch.
                # Reduces initialization σ from ~0.6m to ~0.35m.
                cp_val = float(d['L_if_raw']) - (float(d['P_if_raw']) + dcb_corr)
                buf = self._cp_buf.setdefault(sv, [])
                buf.append(cp_val)
                if len(buf) > self._cp_buf_size:
                    buf.pop(0)

                if len(buf) >= min(3, self._cp_buf_size):
                    amb_base = float(np.median(buf))
                else:
                    amb_base = cp_val  # fallback: single-epoch
                # If WL already fixed from earlier epochs of this SV pass, use
                # the WL-constrained initialization.  This anchors the wide-lane
                # component and lets NL float absorb only NL noise.
                if sv in self._wl_fixed:
                    from src.ambiguity import compute_mw, COEFF_W
                    L1c = float(d.get('L1_cyc', 0))
                    L2c = float(d.get('L2_cyc', 0))
                    P1 = float(d.get('P1_raw', d.get('P_if_raw', 0)))
                    P2 = float(d.get('P2_raw', d.get('P_if_raw', 0)))
                    mw_float = compute_mw(L1c, L2c, P1, P2)
                    N_w_float = round(mw_float)
                    # Use WL-fixed to correct the code-phase initialization
                    amb_init = amb_base + COEFF_W * (self._wl_fixed[sv] - N_w_float)
                else:
                    amb_init = amb_base

            # ── Code update FIRST (anchors clock, no ambiguity term) ──
            if 'P_if_raw' in d:
                obs_code = float(d['P_if_raw']) + dcb_corr - sv_bias[sv]
                if 1.5e7 <= abs(obs_code) <= 3.5e7:
                    state = self._ensure_sv(state, sv, amb_init=amb_init)
                    if sv not in self._sv_last_seen:
                        self._sv_last_seen[sv] = state.t
                    modeled = self._compute_modeled(
                        state, sv, sat_pos_eci, sat_pos_ecef_corr, sat_clk,
                        rcv_pos_ecef, zhd, mf_h, mf_w, with_amb=False,
                        rel_corr=rel_corr)
                    H = self._build_H(state, sv, sat_pos_eci, mf_w, with_amb=False)
                    R_code_el = self._meas_noise(self.sigma_code, el, exponent=self._elev_exp_code)
                    state, accepted, nu, chi2 = self._update_scalar(
                        state, sv, obs_code, modeled, H, R_code_el)
                    if accepted:
                        n_code += 1
                        sum_code_nu2 += nu ** 2
                        self._sv_last_seen[sv] = state.t
                    else:
                        n_rej += 1

            # ── Phase update SECOND (ambiguity absorbs phase bias) ──
            if 'L_if_raw' in d:
                obs_phase = float(d['L_if_raw']) - sv_bias[sv]
                if 1.5e7 <= abs(obs_phase) <= 3.5e7:
                    state = self._ensure_sv(state, sv, amb_init=amb_init)
                    if sv not in self._sv_last_seen:
                        self._sv_last_seen[sv] = state.t
                    modeled = self._compute_modeled(
                        state, sv, sat_pos_eci, sat_pos_ecef_corr, sat_clk,
                        rcv_pos_ecef, zhd, mf_h, mf_w, with_amb=True,
                        rel_corr=rel_corr, windup_corr=windup_corr)
                    H = self._build_H(state, sv, sat_pos_eci, mf_w, with_amb=True)
                    R_phase_el = self._meas_noise(self.sigma_phase, el, exponent=self._elev_exp_phase)
                    state, accepted, nu, chi2 = self._update_scalar(
                        state, sv, obs_phase, modeled, H, R_phase_el)
                    if accepted:
                        n_phase += 1
                        sum_phase_nu2 += nu ** 2
                        self._sv_last_seen[sv] = state.t
                        # Progressive WL tightening: reduce P_amb as phase confirms fix
                        if sv in self._wl_fixed and sv in state.sv_to_idx:
                            wl_ep = self._wl_epochs.get(sv, 0) + 1
                            self._wl_epochs[sv] = wl_ep
                            # Phase residuals < 0.5m → WL fix is consistent
                            if abs(nu) < 0.5:
                                amb_idx = I_AMB_START + state.sv_to_idx[sv]
                                if wl_ep <= 3:
                                    state.P[amb_idx, amb_idx] = min(
                                        state.P[amb_idx, amb_idx], 0.10)
                                elif wl_ep <= 7:
                                    state.P[amb_idx, amb_idx] = min(
                                        state.P[amb_idx, amb_idx], 0.04)
                                else:
                                    state.P[amb_idx, amb_idx] = min(
                                        state.P[amb_idx, amb_idx], 0.01)
                    else:
                        n_rej += 1

            # ── MW Accumulation + WL Fixing (Phase 3.0 PPP-AR) ──
            L1c = float(d.get('L1_cyc', 0))
            L2c = float(d.get('L2_cyc', 0))
            if L1c != 0 and L2c != 0:
                P1 = float(d.get('P1_raw', d.get('P_if_raw', 0)))
                P2 = float(d.get('P2_raw', d.get('P_if_raw', 0)))
                self._mw_buf.add(sv, L1c, L2c, P1, P2)

                # Attempt WL fix after sufficient accumulation
                if sv not in self._wl_fixed and sv in state.sv_to_idx:
                    N_w = self._mw_buf.try_fix_wl(sv)
                    if N_w is not None:
                        self._wl_fixed[sv] = N_w
                        # Correct the ambiguity state to be consistent with WL-fixed N_w.
                        # B_if = λ_nl*N1 + coeff_w*N_w  where N1 is still float (unknown).
                        # The filter currently has B_if_float.  We shift it so that
                        # the WL component matches N_w_fixed.
                        from src.ambiguity import COEFF_W
                        vals = self._mw_buf._mw.get(sv, [])
                        N_w_float_mw = float(np.mean(vals))
                        if self._mw_buf._b_r_wl is not None:
                            N_w_float_mw -= self._mw_buf._b_r_wl
                        amb_correction = COEFF_W * (N_w - N_w_float_mw)
                        amb_idx = I_AMB_START + state.sv_to_idx[sv]
                        state.x[amb_idx] += amb_correction
                        # Tighten ambiguity variance: WL fixed → only NL unknown
                        # σ_B_if = λ_nl * σ_N1 ≈ 0.107m * 3 ≈ 0.32m (NL still float)
                        # P_amb = (3*0.107)^2 ≈ 0.10 m^2
                        state.P[amb_idx, amb_idx] = min(
                            state.P[amb_idx, amb_idx], 0.10)
                        self._wl_epochs[sv] = 0  # start progressive tightening counter

        # ── NL Ambiguity Fixing via SD (Phase 5.0) ──
        # When ≥4 SVs have WL fixed, attempt between-satellite single-difference
        # NL fixing.  SD eliminates the receiver NL bias, leaving integer ΔN1.
        # Called once per epoch after all SVs have been processed.
        self._try_nl_sd_fix(state, ed_sorted)

        # ── Prune dead SVs ──
        # Remove SVs that haven't been observed for prune_timeout seconds.
        # Dead ambiguities from set SVs poison the covariance for multi-hour arcs.
        pruned = []
        for sv in list(state.sv_list):
            last_seen = self._sv_last_seen.get(sv, state.t)
            if state.t - last_seen > self.prune_timeout:
                state = state.remove_sv(sv)
                self._mw_buf.reset(sv)
                self._cp_buf.pop(sv, None)  # Phase 7.1
                self._wl_fixed.pop(sv, None)
                self._wl_epochs.pop(sv, None)
                self._nl_fixed.pop(sv, None)
                if self._nl_ref_sv == sv:
                    self._nl_ref_sv = None
                    self._nl_fixed.clear()
                self._sv_last_seen.pop(sv, None)
                pruned.append(sv)
        if pruned:
            n_pruned = len(pruned)
            # Only log if more than a couple
            if n_pruned > 2:
                print(f"  [PRUNE] {n_pruned} SVs pruned: {', '.join(sorted(pruned)[:6])}"
                      f"{'...' if n_pruned > 6 else ''}")

        rms_phase = np.sqrt(sum_phase_nu2 / max(n_phase, 1))
        rms_code = np.sqrt(sum_code_nu2 / max(n_code, 1))

        stats = {
            'n_phase': n_phase, 'n_code': n_code, 'n_rej': n_rej,
            'rms_phase': rms_phase, 'rms_code': rms_code,
            'n_wl_fixed': len(self._wl_fixed),
        }
        return state, stats
