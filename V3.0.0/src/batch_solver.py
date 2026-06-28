"""Batch linear solver for clock + ZWD + ambiguity estimation (Phase 6.0).

Given a fixed orbit trajectory from the EKF, jointly estimates:
  - Receiver clock clk_i [m] per epoch i
  - Zenith wet delay zwd_i [m] per epoch i
  - IF carrier-phase ambiguity amb_k [m] per SV k

All parameters are solved from a single over-determined linear system
A.x = b.  This eliminates the EKF's sequential initialization noise and
cross-clock-reference-frame ambiguity inconsistency.

The solved ambiguities are SELF-CONSISTENT (same clock reference) and
can be used directly in the EKF pass 2 via _amb_batch_fixed.

Observations:
  P_if = rho_orbit + clk_i + zwd_i*mf + noise            (code)
  L_if = rho_orbit + clk_i + zwd_i*mf + amb_k + noise    (phase)

Kinematic mode (Phase C):
  Adds per-epoch ECEF position corrections dr_i(3) as free parameters.
  Model: |sat - (rcv0+dr)| = |sat-rcv0| + e_vec·dr + O(|dr|²)/R
  where e_vec = (rcv0 - sat)/|rcv0 - sat| (sat->rcv LOS in ECEF).
  Position prior sigma = 10m keeps the system well-conditioned
  without forcing dynamics constraints.

Optional CLK1B prior (Phase 16.0):
  clk_i - clk_{i-1} ~ clk1b_i - clk1b_{i-1}   (sigma ~ 0.50m)

Optional OSB NL fixing (Phase 9.0):
  After float solve, use CODE satellite biases to compute
  absolute N1 integers: N1 = round((B_if - COEFF_W*(N_w+b_wl))/LAM_NL - b_nl)
  Then re-solve with fixed ambiguities as strong priors.
"""

import numpy as np
from src.ambiguity import COEFF_W, F1, F2
C_LIGHT = 299792458.0
LAM_NL = C_LIGHT / (F1 + F2)  # ~0.1070 m
LAM_WL = C_LIGHT / (F1 - F2)  # ~0.8619 m

SEC_PER_DAY = 86400.0


def _trop_mf(el_rad):
    """GMP wet mapping function."""
    if el_rad < 0.05:
        return 20.0
    return 1.001 / np.sqrt(0.002001 + np.sin(el_rad)**2)


def _e_vec_ecef(rcv_ecef, sat_ecef):
    """Unit LOS vector FROM satellite TOWARD receiver in ECEF.

    e = (rcv - sat) / |rcv - sat|

    Range = |sat - rcv|.  For small dr:
      |sat - (rcv0+dr)| = |sat - rcv0| + e·dr + O(|dr|²/R)
    """
    d = np.asarray(rcv_ecef, dtype=float) - np.asarray(sat_ecef, dtype=float)
    r = np.linalg.norm(d)
    if r < 1.0:
        return np.array([0.0, 0.0, 1.0])
    return d / r


def _add_kinematic_block(ATA, ATb, W, e_ecef, res_m,
                         i_drx, i_clk, i_zwd, i_amb, mf):
    """Add per-epoch position (dr) contributions to normal equations.

    Model:  model = geo + e·dr + clk + zwd*mf + amb
    Jacobian: J = [+e, +1, +mf, +1(phase)]
    Normal equations: J^T @ W @ J  and  J^T @ W @ (obs - geo)
    """
    for j in range(3):
        jj = i_drx + j
        # Diagonal position block: e[j] * e[j]
        ATA[jj, jj] += W * e_ecef[j] * e_ecef[j]
        # Cross with clock: e[j] * 1
        ATA[jj, i_clk] += W * e_ecef[j]
        ATA[i_clk, jj] += W * e_ecef[j]
        # Cross with zwd: e[j] * mf
        ATA[jj, i_zwd] += W * e_ecef[j] * mf
        ATA[i_zwd, jj] += W * e_ecef[j] * mf
        # Cross with amb (phase only): e[j] * 1
        if i_amb is not None:
            ATA[jj, i_amb] += W * e_ecef[j]
            ATA[i_amb, jj] += W * e_ecef[j]
        # RHS: e * (obs - geo)
        ATb[jj] += W * e_ecef[j] * res_m
    # Off-diagonal position blocks: e[j] * e[k] for j != k
    for j in range(3):
        jj = i_drx + j
        for k in range(j + 1, 3):
            kk = i_drx + k
            val = W * e_ecef[j] * e_ecef[k]
            ATA[jj, kk] += val
            ATA[kk, jj] += val


class BatchLinearSolver:
    """Joint clock + ZWD + ambiguity (+position in kinematic mode) estimation.

    Usage:
        # EKF-orbit-fixed (existing)
        solver = BatchLinearSolver(epoch_geometry,
                                    sigma_phase=0.20, sigma_code=0.30)
        solution = solver.solve()

        # Kinematic: per-epoch position corrections
        solver = BatchLinearSolver(epoch_geometry,
                                    sigma_phase=0.20, sigma_code=0.30,
                                    kinematic=True,
                                    rcv_ecef_per_epoch=rcv_ecef_list)
        solution = solver.solve()
        # solution['r_ecef'] = (N_epoch, 3)  -- corrected positions
    """

    def __init__(self, epoch_geometry,
                 sigma_phase=0.20, sigma_code=0.30,
                 clk1b_data=None, clk1b_sigma=0.03,
                 epoch_gps_sod=None,
                 wl_fixed=None, osb_wl=None, osb_nl=None,
                 kinematic=False, rcv_ecef_per_epoch=None,
                 sigma_pos_prior=3.0, sigma_pos_smooth=0.5):
        """
        Args:
            epoch_geometry: list of lists, epoch_geometry[i_epoch] = [ep_dict, ...]
                where each ep_dict has keys: sv, el, sat_pos, _obs_code, _obs_phase, _geo_full
            sigma_phase, sigma_code: measurement noise [m]
            clk1b_data: dict {gps_sod: clk_m} -- CLK1B USO clock prior [m]
            clk1b_sigma: CLK1B prior standard deviation [m]
            epoch_gps_sod: list of GPS seconds-of-day per epoch
            wl_fixed: {sv: N_w} -- WL integers from EKF MW
            osb_wl, osb_nl: {sv: bias_cycles} -- CODE OSB satellite biases
            kinematic: if True, add per-epoch ECEF dr parameters
            rcv_ecef_per_epoch: list or (N_epoch,3) array of ECEF positions [m]
                (required when kinematic=True)
            sigma_pos_prior: position prior sigma [m] (default 10m)
        """
        self.epoch_geometry = epoch_geometry
        self.sigma_phase = sigma_phase
        self.sigma_code = sigma_code
        self.clk1b_data = clk1b_data if clk1b_data else {}
        self.clk1b_sigma = clk1b_sigma
        self.epoch_gps_sod = epoch_gps_sod
        self._has_clk1b = bool(self.clk1b_data) and self.epoch_gps_sod is not None

        # OSB NL fixing
        self.wl_fixed = wl_fixed if wl_fixed else {}
        self.osb_wl = osb_wl if osb_wl else {}
        self.osb_nl = osb_nl if osb_nl else {}
        self._has_osb = bool(self.wl_fixed) and bool(self.osb_nl)
        self._amb_fixed = {}
        self._sigma_amb_fixed = 0.01  # 1cm prior

        # Kinematic mode
        self.kinematic = kinematic
        self.sigma_pos_prior = sigma_pos_prior
        self.sigma_pos_smooth = sigma_pos_smooth
        if kinematic and rcv_ecef_per_epoch is not None:
            self.rcv_ecef0 = np.asarray(rcv_ecef_per_epoch, dtype=float)
        else:
            self.rcv_ecef0 = None

        self._index()

    def _index(self):
        """Index SVs, count observations, set up parameter layout."""
        sv_set = set()
        n_code, n_phase = 0, 0
        for ep_list in self.epoch_geometry:
            for d in ep_list:
                if '_obs_code' in d or 'P_if_raw' in d:
                    sv_set.add(d['sv'])
                    n_code += 1
                if '_obs_phase' in d or 'L_if_raw' in d:
                    sv_set.add(d['sv'])
                    n_phase += 1

        self.sv_list = sorted(sv_set)
        self.sv_to_idx = {sv: i for i, sv in enumerate(self.sv_list)}
        self.N_sv = len(self.sv_list)
        self.N_epoch = len(self.epoch_geometry)
        self.N_obs = n_code + n_phase

        # Parameter layout — shifts when kinematic adds position params
        if self.kinematic:
            self.N_pos = self.N_epoch * 3
            self.I_CLK = self.N_pos
        else:
            self.N_pos = 0
            self.I_CLK = 0
        self.N_clk = self.N_epoch
        self.N_zwd = self.N_epoch
        self.N_amb = self.N_sv
        self.I_ZWD = self.I_CLK + self.N_epoch
        self.I_AMB = self.I_ZWD + self.N_epoch
        self.N_total = self.N_pos + self.N_clk + self.N_zwd + self.N_amb

        kin_str = f" +pos({self.N_pos})" if self.kinematic else ""
        clk1b_str = f" +CLK1B({len(self.clk1b_data)}ep)" if self._has_clk1b else ""
        osb_str = f" +OSB({len(self.osb_nl)}SV)" if self._has_osb else ""
        print(f"  [BatchSolver] {self.N_epoch} epochs x {self.N_sv} SVs = {self.N_obs} obs")
        print(f"  [BatchSolver] {self.N_total} params: "
              f"{'dr='+str(self.N_pos)+'+' if self.kinematic else ''}"
              f"clk={self.N_clk} + zwd={self.N_zwd} + amb={self.N_amb}"
              f"{clk1b_str}{osb_str}{kin_str}")

    # ------------------------------------------------------------------
    # Core linear system assembly
    # ------------------------------------------------------------------

    def _solve_internal(self):
        """Build and solve normal equations, return raw (x, res_norm, cost)."""
        W_P = 1.0 / self.sigma_code**2
        W_L = 1.0 / self.sigma_phase**2

        N_param = self.N_total
        ATA = np.zeros((N_param, N_param))
        ATb = np.zeros(N_param)
        res_norm_list = []

        for i_ep, ep_list in enumerate(self.epoch_geometry):
            i_clk = self.I_CLK + i_ep
            i_zwd = self.I_ZWD + i_ep
            i_drx = i_ep * 3  # position index (only used in kinematic mode)

            for d in ep_list:
                sv = d['sv']
                if sv not in self.sv_to_idx:
                    continue
                i_amb = self.I_AMB + self.sv_to_idx[sv]
                geo = float(d.get('_geo_full', 0))
                el = float(d.get('el', 0.5))
                mf = _trop_mf(el)

                # LOS vector: satellite -> receiver in ECEF
                e_ecef = None
                if self.kinematic and self.rcv_ecef0 is not None:
                    sat_ecef = np.asarray(d.get('sat_pos', [0, 0, 0]), dtype=float)
                    e_ecef = _e_vec_ecef(self.rcv_ecef0[i_ep], sat_ecef)

                if '_obs_code' in d:
                    res_m = float(d['_obs_code']) - geo
                    res_norm_list.append(res_m / self.sigma_code)
                    ATA[i_clk, i_clk] += W_P
                    ATA[i_clk, i_zwd] += W_P * mf
                    ATA[i_zwd, i_clk] += W_P * mf
                    ATA[i_zwd, i_zwd] += W_P * mf * mf
                    ATb[i_clk] += W_P * res_m
                    ATb[i_zwd] += W_P * res_m * mf
                    if self.kinematic and e_ecef is not None:
                        _add_kinematic_block(ATA, ATb, W_P, e_ecef, res_m,
                                            i_drx, i_clk, i_zwd, None, mf)

                if '_obs_phase' in d:
                    res_m = float(d['_obs_phase']) - geo
                    res_norm_list.append(res_m / self.sigma_phase)
                    ATA[i_clk, i_clk] += W_L
                    ATA[i_clk, i_zwd] += W_L * mf
                    ATA[i_zwd, i_clk] += W_L * mf
                    ATA[i_zwd, i_zwd] += W_L * mf * mf
                    ATA[i_clk, i_amb] += W_L
                    ATA[i_amb, i_clk] += W_L
                    ATA[i_zwd, i_amb] += W_L * mf
                    ATA[i_amb, i_zwd] += W_L * mf
                    ATA[i_amb, i_amb] += W_L
                    ATb[i_clk] += W_L * res_m
                    ATb[i_zwd] += W_L * res_m * mf
                    ATb[i_amb] += W_L * res_m
                    if self.kinematic and e_ecef is not None:
                        _add_kinematic_block(ATA, ATb, W_L, e_ecef, res_m,
                                            i_drx, i_clk, i_zwd, i_amb, mf)

        # ---- Position priors + smoothness (kinematic mode) ----
        if self.kinematic:
            # Absolute prior: dr ~ N(0, sigma_pos)
            W_pos = 1.0 / self.sigma_pos_prior**2
            for i_ep in range(self.N_epoch):
                for j in range(3):
                    jj = i_ep * 3 + j
                    ATA[jj, jj] += W_pos

            # Epoch-to-epoch smoothness: dr_i - dr_{i-1} ~ N(0, sigma_sm)
            sigma_smooth = self.sigma_pos_smooth  # allows 0.017m/s drift
            W_sm = 1.0 / sigma_smooth**2
            for i_ep in range(1, self.N_epoch):
                for j in range(3):
                    jj_p = (i_ep - 1) * 3 + j
                    jj_c = i_ep * 3 + j
                    ATA[jj_p, jj_p] += W_sm
                    ATA[jj_p, jj_c] -= W_sm
                    ATA[jj_c, jj_p] -= W_sm
                    ATA[jj_c, jj_c] += W_sm

        # ---- CLK1B prior constraints ----
        if self._has_clk1b:
            W_clk1b = 1.0 / self.clk1b_sigma**2
            n_prior = 0
            for i_ep in range(1, self.N_epoch):
                gps_prev = self.epoch_gps_sod[i_ep - 1]
                gps_curr = self.epoch_gps_sod[i_ep]
                clk_prev = clk_curr = None
                if gps_prev in self.clk1b_data:
                    clk_prev = self.clk1b_data[gps_prev]
                else:
                    for clk_sod in self.clk1b_data:
                        if abs(clk_sod - gps_prev) <= 5.0:
                            clk_prev = self.clk1b_data[clk_sod]; break
                if gps_curr in self.clk1b_data:
                    clk_curr = self.clk1b_data[gps_curr]
                else:
                    for clk_sod in self.clk1b_data:
                        if abs(clk_sod - gps_curr) <= 5.0:
                            clk_curr = self.clk1b_data[clk_sod]; break
                if clk_prev is not None and clk_curr is not None:
                    raw_diff = clk_curr - clk_prev
                    clk1b_diff = -raw_diff
                    i_p = self.I_CLK + i_ep - 1
                    i_c = self.I_CLK + i_ep
                    ATA[i_p, i_p] += W_clk1b
                    ATA[i_p, i_c] -= W_clk1b
                    ATA[i_c, i_p] -= W_clk1b
                    ATA[i_c, i_c] += W_clk1b
                    ATb[i_p] -= W_clk1b * clk1b_diff
                    ATb[i_c] += W_clk1b * clk1b_diff
                    n_prior += 1
            if n_prior > 0:
                print(f"  [BS-CLK1B] {n_prior} clock-diff priors "
                      f"(sigma={self.clk1b_sigma:.3f}m)")

        # ---- Fixed ambiguity priors (OSB NL fixing) ----
        if self._amb_fixed:
            W_amb = 1.0 / self._sigma_amb_fixed**2
            for sv, B_if_fixed in self._amb_fixed.items():
                if sv in self.sv_to_idx:
                    i_amb = self.I_AMB + self.sv_to_idx[sv]
                    ATA[i_amb, i_amb] += W_amb
                    ATb[i_amb] += W_amb * B_if_fixed

        # ---- Regularization (rank deficiency removal) ----
        for i in range(self.N_epoch):
            ATA[self.I_ZWD + i, self.I_ZWD + i] += 1e-4
        for i in range(self.N_sv):
            ATA[self.I_AMB + i, self.I_AMB + i] += 1e-4

        try:
            x = np.linalg.solve(ATA, ATb)
        except np.linalg.LinAlgError:
            x = np.linalg.lstsq(ATA, ATb, rcond=1e-8)[0]

        res_norm = np.array(res_norm_list)
        cost = 0.5 * np.sum(res_norm**2)
        return x, res_norm, cost

    # ------------------------------------------------------------------
    # RMS computation and solving
    # ------------------------------------------------------------------

    def _compute_rms(self, x):
        """Compute separate phase/code RMS from solved parameters."""
        n_phase, n_code = 0, 0
        sum_phase, sum_code = 0.0, 0.0
        for i_ep, ep_list in enumerate(self.epoch_geometry):
            clk = x[self.I_CLK + i_ep]
            zwd = x[self.I_ZWD + i_ep]

            # Position correction in kinematic mode
            dr = np.zeros(3)
            if self.kinematic:
                dr = x[i_ep * 3 : i_ep * 3 + 3]

            # LOS vector reuse
            rcv0 = self.rcv_ecef0[i_ep] if self.rcv_ecef0 is not None else np.zeros(3)

            for d in ep_list:
                sv = d['sv']
                if sv not in self.sv_to_idx:
                    continue
                geo = float(d.get('_geo_full', 0))
                el = float(d.get('el', 0.5))
                mf = _trop_mf(el)
                amb = x[self.I_AMB + self.sv_to_idx[sv]]

                # Geometry correction from position update
                geo_corr = 0.0
                if self.kinematic and np.any(dr):
                    sat_ecef = np.asarray(d.get('sat_pos', [0, 0, 0]), dtype=float)
                    e_ecef = _e_vec_ecef(rcv0, sat_ecef)
                    geo_corr = float(np.dot(e_ecef, dr))

                if '_obs_phase' in d:
                    model = geo + geo_corr + clk + zwd * mf + amb
                    err = float(d['_obs_phase']) - model
                    sum_phase += err**2
                    n_phase += 1
                if '_obs_code' in d:
                    model = geo + geo_corr + clk + zwd * mf
                    err = float(d['_obs_code']) - model
                    sum_code += err**2
                    n_code += 1
        rms_phase = np.sqrt(sum_phase / max(n_phase, 1))
        rms_code = np.sqrt(sum_code / max(n_code, 1))
        return rms_phase, rms_code

    def _fix_nl_with_osb(self, amb_float):
        """Compute absolute N1 integers from float B_if + OSB + WL."""
        nl_fixed = {}
        amb_fixed = {}
        for sv, N_w in self.wl_fixed.items():
            if sv not in amb_float or sv not in self.osb_nl:
                continue
            b_wl = self.osb_wl.get(sv, 0.0)
            b_nl = self.osb_nl.get(sv, 0.0)
            B_if = float(amb_float[sv])
            N1_float = (B_if - COEFF_W * (N_w + b_wl)) / LAM_NL - b_nl
            N1_int = int(round(N1_float))
            resid = abs(N1_float - N1_int)
            if resid < 0.30:
                nl_fixed[sv] = N1_int
                B_if_fixed = LAM_NL * (N1_int + b_nl) + COEFF_W * (N_w + b_wl)
                amb_fixed[sv] = B_if_fixed
        return nl_fixed, amb_fixed

    def solve(self):
        """Public API: two-pass float then OSB-fixed solution.

        In kinematic mode, also returns corrected ECEF positions.
        """
        self._amb_fixed = {}

        # ---- Pass 1: Float solution ----
        x_float, res_norm_float, cost_float = self._solve_internal()
        clk_f = x_float[self.I_CLK : self.I_CLK + self.N_clk]
        zwd_f = x_float[self.I_ZWD : self.I_ZWD + self.N_zwd]
        amb_arr_f = x_float[self.I_AMB : self.I_AMB + self.N_amb]
        amb_float = {sv: float(amb_arr_f[i]) for sv, i in self.sv_to_idx.items()}
        rms_p_float, rms_c_float = self._compute_rms(x_float)

        # ---- Pass 2: OSB NL fixing ----
        nl_fixed = {}
        amb_dict = dict(amb_float)
        x_final = x_float
        rms_p, rms_c = rms_p_float, rms_c_float

        if self._has_osb:
            nl_fixed, self._amb_fixed = self._fix_nl_with_osb(amb_float)
            n_fixed = len(self._amb_fixed)
            if n_fixed > 0:
                print(f"  [BS-OSB] {n_fixed}/{len(self.wl_fixed)} SVs NL-fixed "
                      f"(sigma_amb={self._sigma_amb_fixed:.3f}m)")
                x_final, _, _ = self._solve_internal()
                rms_p, rms_c = self._compute_rms(x_final)
                # Update dicts from re-solved values
                amb_arr_new = x_final[self.I_AMB : self.I_AMB + self.N_amb]
                amb_dict = {sv: float(amb_arr_new[i]) for sv, i in self.sv_to_idx.items()}
                clk_f = x_final[self.I_CLK : self.I_CLK + self.N_clk]
                zwd_f = x_final[self.I_ZWD : self.I_ZWD + self.N_zwd]
            else:
                print(f"  [BS-OSB] No SVs NL-fixed "
                      f"(OSB={len(self.osb_nl)}, WL={len(self.wl_fixed)})")

        # ---- Extract kinematic positions ----
        r_ecef_out = None
        if self.kinematic and self.rcv_ecef0 is not None:
            dr_out = np.zeros((self.N_epoch, 3))
            for i_ep in range(self.N_epoch):
                dr_out[i_ep] = x_final[i_ep * 3 : i_ep * 3 + 3]
            r_ecef_out = self.rcv_ecef0 + dr_out

        rms_total = np.sqrt(rms_p**2 * 0.5 + rms_c**2 * 0.5)

        return {
            'clk': clk_f, 'zwd': zwd_f,
            'amb': np.array([amb_dict.get(sv, 0.0) for sv in self.sv_list]),
            'amb_dict': amb_dict,
            'amb_float': amb_float,
            'nl_fixed': nl_fixed,
            'rms_phase': rms_p, 'rms_code': rms_c,
            'rms_phase_float': rms_p_float,
            'rms_total': rms_total, 'cost': cost_float,
            'sv_list': self.sv_list,
            'r_ecef': r_ecef_out,
            'dr_ecef': (x_final[:self.N_pos].reshape(-1, 3)
                         if self.kinematic else None),
        }
