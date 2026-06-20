"""Batch ambiguity resolution for POD (Phase 6.0).

Multi-pass strategy:
  1. First EKF pass (float) — collect ambiguity states
  2. Full-arc MW smoothing → WL integer per SV pass
  3. NL float extraction from EKF ambiguity estimates
  4. Integer bootstrapping with proper NL covariance
  5. Second EKF pass with fixed ambiguities (P_amb = 0.001 m²)

This achieves the batch AR benefit without rebuilding the POD pipeline.
"""
import numpy as np

# GPS constants
F1, F2 = 1575.42e6, 1227.60e6
C_LIGHT = 299792458.0
WL_WAVELENGTH = C_LIGHT / (F1 - F2)  # ~0.862 m
LAM_NL = C_LIGHT / (F1 + F2)  # ~0.107 m (narrow-lane wavelength)
COEFF_W = F2 / (F1**2 - F2**2) * C_LIGHT  # ~0.3776 m/cycle (WL→IF factor)


def mw_combination(L1_cyc, L2_cyc, P1_m, P2_m):
    """Melbourne-Wübbena wide-lane [cycles]."""
    nl_pr = (F1 * P1_m + F2 * P2_m) / (F1 + F2)
    wl_phase = C_LIGHT * (L1_cyc / F1 - L2_cyc / F2) / (F1 - F2) * (F1 * F2 / C_LIGHT)
    # Simplified: MW = (L1 - L2) - (f1*P1 + f2*P2)/(f1+f2) / wl_wavelength
    wl_phase_m = (F1 * L1_cyc * C_LIGHT / F1 - F2 * L2_cyc * C_LIGHT / F2) / (F1 - F2)
    return (wl_phase_m - nl_pr) / WL_WAVELENGTH


def compute_mw(L1_cyc, L2_cyc, P1_m, P2_m):
    """MW from L1/L2 in cycles, P1/P2 in meters."""
    # Wide-lane phase in meters: (f1*L1 - f2*L2) / (f1 - f2)
    L1_m = L1_cyc * C_LIGHT / F1
    L2_m = L2_cyc * C_LIGHT / F2
    wl_phase = (F1 * L1_m - F2 * L2_m) / (F1 - F2)
    # Narrow-lane code in meters: (f1*P1 + f2*P2) / (f1 + f2)
    nl_code = (F1 * P1_m + F2 * P2_m) / (F1 + F2)
    # MW in cycles: (wl_phase - nl_code) / wl_wavelength
    return (wl_phase - nl_code) / WL_WAVELENGTH


class BatchAmbiguityResolver:
    """Full-arc batch ambiguity resolver.

    Processes the complete arc of observations to produce WL and NL
    integer ambiguity estimates using multi-epoch smoothing rather
    than the EKF's single-epoch initialization.

    Usage:
        resolver = BatchAmbiguityResolver(epoch_data_list, ekf_float_results)
        amb_fixed = resolver.resolve()
        # amb_fixed['N_w']['G05'] = -2 (WL integer)
        # amb_fixed['N1']['G05'] = 12345 (NL integer)
    """

    def __init__(self, epoch_data_list, sv_to_idx_map=None,
                 min_epochs_wl=5, max_wl_std=0.30,
                 max_nl_residual=0.30, min_sv_for_nl=4):
        """
        Args:
            epoch_data_list: list of per-epoch dicts with keys:
                sv, L1_cyc, L2_cyc, P1_raw, P2_raw, L_if_raw, P_if_raw,
                sat_pos, sat_clk, rho_corr, el
            sv_to_idx_map: optional {sv: index} from EKF state
            min_epochs_wl: minimum epochs for WL arc-average
            max_wl_std: maximum WL std [cycles] for fixing
            max_nl_residual: maximum NL rounding residual [cycles]
            min_sv_for_nl: minimum SVs for NL SD fixing
        """
        self.ep_data = epoch_data_list
        self.min_epochs_wl = min_epochs_wl
        self.max_wl_std = max_wl_std
        self.max_nl_residual = max_nl_residual
        self.min_sv_for_nl = min_sv_for_nl
        self.sv_to_idx = sv_to_idx_map or {}

        self._collect_sv_passes()

    def _collect_sv_passes(self):
        """Group observations by SV for full-arc processing."""
        sv_epochs = {}
        for d in self.ep_data:
            sv = d['sv']
            sv_epochs.setdefault(sv, []).append(d)

        self.sv_passes = {}
        for sv, epochs in sv_epochs.items():
            if len(epochs) >= self.min_epochs_wl:
                self.sv_passes[sv] = epochs

    def resolve_wl(self):
        """Full-arc MW smoothing → WL integer per SV pass.

        Uses the full arc's MW observations to estimate N_w with
        arc-averaging precision: σ_MW ≈ 0.3/√N → 0.04 cyc for N=60.

        Returns:
            dict: {sv: N_w_fixed} for successfully fixed SVs
            dict: {sv: std_mw} MW standard deviation per SV
        """
        wl_fixed = {}
        wl_stats = {}

        for sv, epochs in self.sv_passes.items():
            mw_vals = []
            for d in epochs:
                L1c = float(d.get('L1_cyc', 0))
                L2c = float(d.get('L2_cyc', 0))
                P1 = float(d.get('P1_raw', 0))
                P2 = float(d.get('P2_raw', 0))
                if L1c == 0 or L2c == 0:
                    continue
                mw = compute_mw(L1c, L2c, P1, P2)
                mw_vals.append(mw)

            if len(mw_vals) < self.min_epochs_wl:
                continue

            mean_mw = float(np.mean(mw_vals))
            std_mw = float(np.std(mw_vals))
            wl_stats[sv] = std_mw

            if std_mw > self.max_wl_std:
                continue

            N_w = int(round(mean_mw))
            residual = abs(mean_mw - N_w)
            if residual > 0.35:
                continue

            wl_fixed[sv] = N_w

        # Estimate receiver WL bias from median fractional part
        if len(wl_fixed) >= 3:
            fracs = []
            for sv in wl_fixed:
                mw_vals = []
                for d in self.sv_passes.get(sv, []):
                    L1c = float(d.get('L1_cyc', 0))
                    L2c = float(d.get('L2_cyc', 0))
                    P1 = float(d.get('P1_raw', 0))
                    P2 = float(d.get('P2_raw', 0))
                    if L1c and L2c:
                        mw_vals.append(compute_mw(L1c, L2c, P1, P2))
                if mw_vals:
                    mean_mw = float(np.mean(mw_vals))
                    frac = mean_mw - round(mean_mw)
                    fracs.append(frac)

            self.b_r_wl = float(np.median(fracs))

            # Re-fix with bias correction
            wl_fixed_corrected = {}
            for sv, N_w in wl_fixed.items():
                mw_vals = []
                for d in self.sv_passes.get(sv, []):
                    L1c = float(d.get('L1_cyc', 0))
                    L2c = float(d.get('L2_cyc', 0))
                    P1 = float(d.get('P1_raw', 0))
                    P2 = float(d.get('P2_raw', 0))
                    if L1c and L2c:
                        mw_vals.append(compute_mw(L1c, L2c, P1, P2))
                mean_corrected = float(np.mean(mw_vals)) - self.b_r_wl
                N_w_corr = int(round(mean_corrected))
                if abs(mean_corrected - N_w_corr) < 0.35:
                    wl_fixed_corrected[sv] = N_w_corr

            print(f"    [BATCH-WL] {len(wl_fixed_corrected)}/{len(self.sv_passes)} SVs fixed "
                  f"(b_r_wl={self.b_r_wl:+.4f} cyc)")
            return wl_fixed_corrected, wl_stats

        print(f"    [BATCH-WL] {len(wl_fixed)}/{len(self.sv_passes)} SVs fixed "
              f"(no receiver bias)")
        return wl_fixed, wl_stats

    def resolve_nl(self, wl_fixed, ekf_amb, ekf_amb_cov=None):
        """NL ambiguity resolution using SD integer bootstrapping.

        For WL-fixed SVs, extracts N1_float = (B_if - coeff_w*N_w) / λ_nl.
        Forms between-satellite single differences to eliminate receiver NL bias.
        Uses integer bootstrapping with conditional variance ordering.

        Args:
            wl_fixed: {sv: N_w} from resolve_wl()
            ekf_amb: {sv: B_if_float} from EKF state [meters]
            ekf_amb_cov: optional covariance sub-matrix for ambiguities

        Returns:
            nl_fixed: {sv: N1_fixed} NL integer per SV (relative to reference)
            nl_ref_sv: reference SV for SD
        """
        if len(wl_fixed) < self.min_sv_for_nl:
            return {}, None

        # Get ambiguity estimates for WL-fixed SVs
        svs = sorted(sv for sv in wl_fixed if sv in ekf_amb)
        if len(svs) < self.min_sv_for_nl:
            return {}, None

        # Select reference SV (one with most epochs = longest arc)
        sv_ref = max(svs, key=lambda sv: len(self.sv_passes.get(sv, [])))
        N_w_ref = wl_fixed[sv_ref]
        B_if_ref = ekf_amb[sv_ref]

        sd_svs = [sv for sv in svs if sv != sv_ref]
        n_sd = len(sd_svs)

        # SD NL float estimates
        z_float = np.zeros(n_sd)
        Q_z = np.eye(n_sd) * 0.01  # default covariance [cycles²]

        for i, sv in enumerate(sd_svs):
            N_w_k = wl_fixed[sv]
            B_if_k = ekf_amb[sv]
            sd_m = B_if_k - B_if_ref - COEFF_W * (N_w_k - N_w_ref)
            z_float[i] = sd_m / LAM_NL

        # Build SD covariance from EKF if available
        if ekf_amb_cov is not None:
            for i, sv_i in enumerate(sd_svs):
                for j, sv_j in enumerate(sd_svs):
                    idx_i = self.sv_to_idx.get(sv_i)
                    idx_j = self.sv_to_idx.get(sv_j)
                    idx_ref = self.sv_to_idx.get(sv_ref)
                    if all(x is not None for x in [idx_i, idx_j, idx_ref]):
                        Q_z[i, j] = (ekf_amb_cov[idx_i, idx_j]
                                     - ekf_amb_cov[idx_i, idx_ref]
                                     - ekf_amb_cov[idx_ref, idx_j]
                                     + ekf_amb_cov[idx_ref, idx_ref]) / LAM_NL**2

        # Integer bootstrapping
        idx_order = np.argsort(np.diag(Q_z))
        z_fixed = np.zeros(n_sd)
        nl_fixed = {}
        ok_svs = []

        for rank, idx in enumerate(idx_order):
            if rank == 0:
                z_cond = z_float[idx]
                sigma_cond = np.sqrt(Q_z[idx, idx])
            else:
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

            if residual > self.max_nl_residual:
                continue

            z_fixed[idx] = z_cond - z_int
            ok_svs.append(idx)
            nl_fixed[sd_svs[idx]] = z_int

        if len(nl_fixed) >= 3:
            res_str = ' '.join(f'{r:.2f}' for r in (z_float[ok_svs]
                               - np.round(z_float[ok_svs])))
            print(f"    [BATCH-NL] {len(nl_fixed)} SVs via SD bootstrapping "
                  f"(ref={sv_ref}, resid={res_str} cyc)")
        else:
            print(f"    [BATCH-NL] failed: only {len(nl_fixed)} SVs fixed "
                  f"(need {self.min_sv_for_nl})")
            return {}, None

        return nl_fixed, sv_ref

    def resolve(self, ekf_amb, ekf_amb_cov=None):
        """Full batch ambiguity resolution.

        Returns:
            dict with keys:
                'wl_fixed': {sv: N_w}
                'nl_fixed': {sv: dN1} SD NL integers
                'nl_ref': reference SV
                'b_r_wl': receiver WL bias
                'wl_stats': {sv: std_mw}
                'amb_if_fixed': {sv: B_if_fixed} full IF ambiguity [m]
        """
        wl_fixed, wl_stats = self.resolve_wl()

        if not wl_fixed:
            return {
                'wl_fixed': {}, 'nl_fixed': {}, 'nl_ref': None,
                'b_r_wl': getattr(self, 'b_r_wl', 0.0),
                'wl_stats': wl_stats,
                'amb_if_fixed': {},
            }

        nl_fixed, nl_ref = self.resolve_nl(wl_fixed, ekf_amb, ekf_amb_cov)
        amb_if_fixed = self._build_absolute_if(wl_fixed, nl_fixed, nl_ref, ekf_amb)

        return {
            'wl_fixed': wl_fixed,
            'nl_fixed': nl_fixed,
            'nl_ref': nl_ref,
            'b_r_wl': getattr(self, 'b_r_wl', 0.0),
            'wl_stats': wl_stats,
            'amb_if_fixed': amb_if_fixed,
        }

    def _build_absolute_if(self, wl_fixed, nl_fixed, nl_ref, ekf_amb):
        """Convert WL+NL integers to absolute IF ambiguities [m].

        Uses the EKF's float IF ambiguity for the reference SV to determine
        the absolute N1 reference.  SD ΔN1 values from bootstrapping are
        clock-independent (clock cancels in between-satellite differences).

        Strategy:
          N1_ref = round((B_if_ekf[ref] - coeff_w*N_w[ref]) / λ_nl)
          B_if[k] = λ_nl * (N1_ref + ΔN1[k]) + coeff_w * N_w[k]

        The absolute B_if values are self-consistent (same N1 reference)
        and compatible with the EKF's final clock state.  For Pass 2,
        the common clock offset is absorbed during the first epoch's
        code update.

        Only NL-fixed SVs are included — WL-only SVs use sequential init.
        """
        amb_if = {}
        if not nl_fixed:
            return amb_if

        # Anchor N1_ref from EKF float B_if of the reference SV.
        # Use median of multiple candidates if reference not available.
        if nl_ref and nl_ref in ekf_amb and nl_ref in wl_fixed:
            B_if_ref = float(ekf_amb[nl_ref])
            N_w_ref = wl_fixed[nl_ref]
            N1_ref = int(round((B_if_ref - COEFF_W * N_w_ref) / LAM_NL))
        else:
            # Fallback: use median N1 across all NL-fixed SVs
            n1_candidates = []
            for sv in nl_fixed:
                if sv in ekf_amb and sv in wl_fixed:
                    B_if = float(ekf_amb[sv])
                    N_w = wl_fixed[sv]
                    dN1 = nl_fixed[sv]
                    n1_abs = round((B_if - COEFF_W * N_w) / LAM_NL)
                    n1_candidates.append(n1_abs - dN1)  # back-compute N1_ref
            if n1_candidates:
                N1_ref = int(round(float(np.median(n1_candidates))))
            else:
                N1_ref = 0

        for sv in nl_fixed:
            if sv not in wl_fixed:
                continue
            N_w = wl_fixed[sv]
            dN1 = nl_fixed[sv]
            amb_if[sv] = LAM_NL * (N1_ref + dN1) + COEFF_W * N_w

        # Reference SV (ΔN1 = 0 by definition), if not already in nl_fixed
        if nl_ref and nl_ref in wl_fixed and nl_ref not in amb_if:
            N_w_ref = wl_fixed[nl_ref]
            amb_if[nl_ref] = LAM_NL * N1_ref + COEFF_W * N_w_ref

        # WL-only SVs (not NL-fixed): include with EKF float B_if rounded
        for sv, N_w in wl_fixed.items():
            if sv not in amb_if and sv in ekf_amb:
                B_if = float(ekf_amb[sv])
                N1_wl = int(round((B_if - COEFF_W * N_w) / LAM_NL))
                amb_if[sv] = LAM_NL * N1_wl + COEFF_W * N_w

        return amb_if


# ═══════════════════════════════════════════════════════════════
# Phase 9.0: CODE OSB integer clock support (undifferenced AR)
# ═══════════════════════════════════════════════════════════════

def read_code_osb(osb_path):
    """Read CODE OSB (SINEX BIAS) file and extract per-satellite WL/NL biases.

    SINEX BIAS v1.00 format (CODE):
      OSB  SITE_CODE  SAT_ID  OBS_TYPE  START  END  UNIT  BIAS  STD

    Example:
      OSB  G049 G01  C1W  2024:120:00000 2024:121:00000 ns  -6.0657  0.0001
      OSB  G049 G01  L1W  2024:120:00000 2024:121:00000 ns  -0.48122 0.00000

    Download:
      http://ftp.aiub.unibe.ch/CODE/2024/COD0OPSFIN_20241200000_01D_01D_OSB.BIA.gz

    Returns:
        osb_wl: {sv: b_wl_sat} in WL cycles
        osb_nl: {sv: b_nl_sat} in NL cycles
    """
    import gzip
    print(f"  [OSB] Reading {osb_path}...")

    sv_data = {}  # sv → {obs_type: bias_value_meters}

    opener = gzip.open if osb_path.endswith('.gz') else open
    with opener(osb_path, 'rt' if osb_path.endswith('.gz') else 'r') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('OSB'):
                continue
            parts = line.split()
            if len(parts) < 9:
                continue
            # parts[0]='OSB', parts[1]=site, parts[2]=SV, parts[3]=obs_type,
            # parts[4-5]=time range, parts[6]=unit, parts[7]=value, parts[8]=std
            sv = parts[2]  # e.g., G01
            obs_type = parts[3]  # e.g., C1W, L1W
            unit = parts[6]  # e.g., ns
            try:
                val_ns = float(parts[7])  # bias value in nanoseconds
            except ValueError:
                continue
            # Convert ns to meters: c * ns * 1e-9 = ns * 0.299792458
            val_m = val_ns * 0.299792458
            sv_data.setdefault(sv, {})[obs_type] = val_m

    # Compute per-SV WL and NL biases from OSB values
    # WL satellite bias:
    #   b_wl_sat = WL_phase_bias - WL_code_bias  [cycles]
    #   WL_code  = (f1*C1W + f2*C2W) / (f1+f2)   [m]
    #   WL_phase = (f1*L1W - f2*L2W) / (f1-f2)   [m]
    #   b_wl_sat = (WL_phase - WL_code) / λ_wl   [cycles]
    #
    # NL satellite bias:
    #   b_nl_sat = NL_phase_bias / λ_nl  [cycles]
    #   NL_phase = (f1*L1W + f2*L2W) / (f1+f2)   [m]
    osb_wl = {}
    osb_nl = {}

    for sv, b in sv_data.items():
        if not sv.startswith('G'):
            continue
        if not ('C1W' in b and 'C2W' in b and 'L1W' in b and 'L2W' in b):
            continue

        # Wide-lane code and phase combinations
        wl_code_m = (F1 * b['C1W'] + F2 * b['C2W']) / (F1 + F2)
        wl_phase_m = (F1 * b['L1W'] - F2 * b['L2W']) / (F1 - F2)
        b_wl_sat = (wl_phase_m - wl_code_m) / WL_WAVELENGTH
        osb_wl[sv] = b_wl_sat

        # Narrow-lane phase combination
        nl_phase_m = (F1 * b['L1W'] + F2 * b['L2W']) / (F1 + F2)
        b_nl_sat = nl_phase_m / LAM_NL
        osb_nl[sv] = b_nl_sat

    print(f"  [OSB] {len(osb_wl)} GPS SVs loaded "
          f"(b_wl range=[{min(osb_wl.values()):+.3f}, {max(osb_wl.values()):+.3f}] cyc, "
          f"b_nl range=[{min(osb_nl.values()):+.3f}, {max(osb_nl.values()):+.3f}] cyc)")
    return osb_wl, osb_nl


def compute_osb_abs_nl(wl_fixed, osb_wl, osb_nl, ekf_amb):
    """With CODE OSB satellite biases, compute absolute N1 per SV.

    Unlike SD (which needs a reference SV), each SV's N1 is independently
    integer-fixable when the satellite bias is known.

    B_if[k] = λ_nl * (N1_true[k] + b_nl_sat[k]) + coeff_w * (N_w_true[k] + b_wl_sat[k])
    → N1[k] = round((B_if[k] - coeff_w * (N_w_fixed[k] + b_wl_sat[k])) / λ_nl - b_nl_sat[k])

    Args:
        wl_fixed: {sv: N_w} WL integers
        osb_wl: {sv: b_wl_sat} satellite WL biases [cyc]
        osb_nl: {sv: b_nl_sat} satellite NL biases [cyc]
        ekf_amb: {sv: B_if} EKF float IF ambiguity [m]

    Returns:
        nl_fixed: {sv: N1} absolute NL integers
    """
    LAM_NL = C_LIGHT / (F1 + F2)
    nl_fixed = {}

    for sv, N_w in wl_fixed.items():
        if sv not in osb_nl or sv not in ekf_amb:
            continue
        b_wl = osb_wl.get(sv, 0.0)
        b_nl = osb_nl.get(sv, 0.0)
        B_if = float(ekf_amb[sv])

        # Absolute N1 (undifferenced, directly from OSB-corrected B_if)
        N1_float = (B_if - COEFF_W * (N_w + b_wl)) / LAM_NL - b_nl
        N1_int = int(round(N1_float))
        residual = abs(N1_float - N1_int)

        if residual < 0.30:
            nl_fixed[sv] = N1_int

    return nl_fixed

    def build_absolute_if_epoch0(self, wl_fixed, nl_fixed, nl_ref, epoch0_data):
        """Build absolute IF ambiguities from epoch-0 code-phase.

        Uses L_if - P_if at epoch 0 (same clock as Pass 2 epoch 0) to
        determine N1_ref, then propagates via SD integers.  This eliminates
        the clock drift between EKF final state and Pass 2 initial epoch.

        Args:
            wl_fixed: {sv: N_w} from resolve_wl()
            nl_fixed: {sv: dN1} SD NL integers
            nl_ref: reference SV ID
            epoch0_data: list of epoch dicts from compute_epoch_geometry
                         at the first epoch (gps_sod_start)

        Returns:
            {sv: B_if_fixed} absolute IF ambiguities [m]
        """
        amb_if = {}
        if not nl_fixed:
            return amb_if

        # Build epoch-0 SV lookup
        ep0_svs = {}
        for d in epoch0_data:
            sv = d['sv']
            if 'L_if_raw' in d and 'P_if_raw' in d:
                # code-phase B_if = L_if - P_if at epoch 0
                ep0_svs[sv] = float(d['L_if_raw']) - float(d['P_if_raw'])

        # Anchor N1_ref from epoch-0 code-phase of reference SV
        N1_ref = 0
        if nl_ref and nl_ref in ep0_svs and nl_ref in wl_fixed:
            B_if_ref = ep0_svs[nl_ref]
            N_w_ref = wl_fixed[nl_ref]
            N1_ref = int(round((B_if_ref - COEFF_W * N_w_ref) / LAM_NL))
        else:
            # Reference not visible at epoch 0: use first visible NL-fixed SV
            for sv in nl_fixed:
                if sv in ep0_svs and sv in wl_fixed:
                    B_if_sv = ep0_svs[sv]
                    N_w_sv = wl_fixed[sv]
                    dN1_sv = nl_fixed[sv]
                    # N1_abs = N1_ref + dN1_sv → N1_ref = N1_abs - dN1_sv
                    N1_ref = int(round((B_if_sv - COEFF_W * N_w_sv) / LAM_NL)) - dN1_sv
                    break

        # Build absolute B_if from N1_ref + SD integers
        for sv in nl_fixed:
            if sv not in wl_fixed:
                continue
            N_w = wl_fixed[sv]
            dN1 = nl_fixed[sv]
            amb_if[sv] = LAM_NL * (N1_ref + dN1) + COEFF_W * N_w

        # Reference SV
        if nl_ref and nl_ref in wl_fixed and nl_ref not in amb_if:
            N_w_ref = wl_fixed[nl_ref]
            amb_if[nl_ref] = LAM_NL * N1_ref + COEFF_W * N_w_ref

        # WL-only SVs visible at epoch 0: use their own code-phase B_if
        for sv, N_w in wl_fixed.items():
            if sv not in amb_if and sv in ep0_svs:
                B_if_sv = ep0_svs[sv]
                N1_sv = int(round((B_if_sv - COEFF_W * N_w) / LAM_NL))
                amb_if[sv] = LAM_NL * N1_sv + COEFF_W * N_w

        return amb_if
