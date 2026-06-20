"""Batch linear solver for clock + ZWD + ambiguity estimation (Phase 6.0).

Given a fixed orbit trajectory from the EKF, jointly estimates:
  - Receiver clock clk_i [m] per epoch i
  - Zenith wet delay zwd_i [m] per epoch i
  - IF carrier-phase ambiguity amb_k [m] per SV k

All parameters are solved from a single over-determined linear system
A·x = b.  This eliminates the EKF's sequential initialization noise and
cross-clock-reference-frame ambiguity inconsistency.

The solved ambiguities are SELF-CONSISTENT (same clock reference) and
can be used directly in the EKF pass 2 via _amb_batch_fixed.

Observations:
  P_if = rho_orbit + clk_i + zwd_i*mf + noise            (code)
  L_if = rho_orbit + clk_i + zwd_i*mf + amb_k + noise    (phase)

Auxiliary constraints remove rank deficiency:
  mean(zwd) = 0, mean(amb) = 0
"""

import numpy as np

SEC_PER_DAY = 86400.0
C_LIGHT = 299792458.0


def _trop_mf(el_rad):
    """GMP wet mapping function."""
    if el_rad < 0.05:
        return 20.0
    return 1.001 / np.sqrt(0.002001 + np.sin(el_rad)**2)


class BatchLinearSolver:
    """Joint clock + ZWD + ambiguity estimation from all epochs.

    Uses pre-computed geometric range base (_geo_base) from the EKF pass 1.
    This includes ECI range + Sagnac - satellite clock, so the batch
    solver only needs to estimate clock, ZWD, and ambiguity.

    Usage:
        solver = BatchLinearSolver(epoch_geometry,
                                    sigma_phase=0.20, sigma_code=0.30)
        solution = solver.solve()
        # solution['amb_dict'] = {sv: B_if [m]}  ← pass to _amb_batch_fixed
    """

    def __init__(self, epoch_geometry,
                 sigma_phase=0.20, sigma_code=0.30):
        """
        Args:
            epoch_geometry: list of lists, epoch_geometry[i_epoch] = [ep_dict, ...]
                where each ep_dict has keys: sv, el, P_if_raw, L_if_raw, _geo_base
            sigma_phase, sigma_code: measurement noise [m]
        """
        self.epoch_geometry = epoch_geometry
        self.sigma_phase = sigma_phase
        self.sigma_code = sigma_code

        self._index()

    def _index(self):
        """Index SVs and count observations."""
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

        # Parameter layout
        self.N_clk = self.N_epoch
        self.N_zwd = self.N_epoch
        self.N_amb = self.N_sv
        self.N_total = self.N_clk + self.N_zwd + self.N_amb
        self.I_ZWD = self.N_epoch
        self.I_AMB = 2 * self.N_epoch

        print(f"  [BatchSolver] {self.N_epoch} epochs × {self.N_sv} SVs = {self.N_obs} obs")
        print(f"  [BatchSolver] {self.N_total} params: clk={self.N_clk} + "
              f"zwd={self.N_zwd} + amb={self.N_amb}")

    def _solve_internal(self):
        """Build and solve normal equations, return raw (x, res_norm, cost)."""
        W_P = 1.0 / self.sigma_code**2
        W_L = 1.0 / self.sigma_phase**2

        N_param = self.N_total
        ATA = np.zeros((N_param, N_param))
        ATb = np.zeros(N_param)
        res_norm_list = []

        for i_ep, ep_list in enumerate(self.epoch_geometry):
            i_clk = i_ep
            i_zwd = self.I_ZWD + i_ep

            for d in ep_list:
                sv = d['sv']
                if sv not in self.sv_to_idx:
                    continue
                i_amb = self.I_AMB + self.sv_to_idx[sv]
                geo = float(d.get('_geo_full', 0))
                el = float(d.get('el', 0.5))
                mf = _trop_mf(el)

                if '_obs_code' in d:
                    res_m = float(d['_obs_code']) - geo
                    res_norm_list.append(res_m / self.sigma_code)
                    ATA[i_clk, i_clk] += W_P
                    ATA[i_clk, i_zwd] += W_P * mf
                    ATA[i_zwd, i_clk] += W_P * mf
                    ATA[i_zwd, i_zwd] += W_P * mf * mf
                    ATb[i_clk] += W_P * res_m
                    ATb[i_zwd] += W_P * res_m * mf

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

        # Regularization
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

    def _compute_rms(self, clk, zwd, amb):
        """Compute separate phase/code RMS from solved parameters."""
        n_phase, n_code = 0, 0
        sum_phase, sum_code = 0.0, 0.0
        for i_ep, ep_list in enumerate(self.epoch_geometry):
            for d in ep_list:
                sv = d['sv']
                if sv not in self.sv_to_idx:
                    continue
                geo = float(d.get('_geo_full', 0))
                el = float(d.get('el', 0.5))
                mf = _trop_mf(el)

                if '_obs_phase' in d:
                    model = geo + clk[i_ep] + zwd[i_ep] * mf + amb[self.sv_to_idx[sv]]
                    err = float(d['_obs_phase']) - model
                    sum_phase += err**2
                    n_phase += 1
                if '_obs_code' in d:
                    model = geo + clk[i_ep] + zwd[i_ep] * mf
                    err = float(d['_obs_code']) - model
                    sum_code += err**2
                    n_code += 1
        rms_phase = np.sqrt(sum_phase / max(n_phase, 1))
        rms_code = np.sqrt(sum_code / max(n_code, 1))
        return rms_phase, rms_code

    def solve(self):
        """Public API: returns solution dict with clk, zwd, amb_dict etc."""
        x, res_norm, cost = self._solve_internal()

        clk = x[0:self.N_clk]
        zwd = x[self.I_ZWD:self.I_ZWD + self.N_zwd]
        amb = x[self.I_AMB:self.I_AMB + self.N_amb]
        amb_dict = {sv: float(amb[i]) for sv, i in self.sv_to_idx.items()}

        rms_phase, rms_code = self._compute_rms(clk, zwd, amb)
        rms_total = np.sqrt(rms_phase**2 * 0.5 + rms_code**2 * 0.5)

        return {
            'clk': clk, 'zwd': zwd, 'amb': amb,
            'amb_dict': amb_dict,
            'rms_phase': rms_phase, 'rms_code': rms_code,
            'rms_total': rms_total, 'cost': cost,
            'sv_list': self.sv_list,
        }
