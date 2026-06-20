"""Arc-based ambiguity resolution (Phase 12.0 / V2.3.0).

Instead of per-epoch MW sliding window, identifies continuous tracking arcs
(free of cycle slips) and estimates one WL + one NL integer per arc.

This follows the literature approach:
  - Identify slip-free arcs using MW+GF TurboEdit
  - Arc-averaged MW → N_w (σ ∝ 1/√N_epoch)
  - Arc-averaged B_if (corrected by OSB) → N1 (σ ∝ 1/√N_epoch)
  - Fixed ambiguities applied as virtual observations in EKF

Benefits over per-epoch fixing:
  - More epochs → better WL precision (σ=0.03 cyc vs 0.08)
  - Full arc B_if → NL converges naturally
  - No "cold start" problem — arc's ambiguity is fixed retroactively
"""

import numpy as np

# GPS constants
F1, F2 = 1575.42e6, 1227.60e6
C_LIGHT = 299792458.0
WL_WAVELENGTH = C_LIGHT / (F1 - F2)  # ~0.862 m
LAM_NL = C_LIGHT / (F1 + F2)         # ~0.107 m
COEFF_W = F2 / (F1**2 - F2**2) * C_LIGHT  # ~0.3776


def compute_mw(L1_cyc, L2_cyc, P1_m, P2_m):
    """Melbourne-Wübbena wide-lane combination [cycles]."""
    L1_m = L1_cyc * C_LIGHT / F1
    L2_m = L2_cyc * C_LIGHT / F2
    wl_phase = (F1 * L1_m - F2 * L2_m) / (F1 - F2)
    nl_code = (F1 * P1_m + F2 * P2_m) / (F1 + F2)
    return (wl_phase - nl_code) / WL_WAVELENGTH


class ArcTracker:
    """Track continuous observation arcs per SV using TurboEdit slip flags.

    An "arc" is a contiguous sequence of epochs where no MW-detected
    cycle slip occurred.  GF-only detections do NOT break the arc
    (N_w is unchanged by equal-frequency slips).

    Usage:
        tracker = ArcTracker()
        for each epoch:
            for each SV:
                slip_mw = tracker.add_epoch(sv, has_slip_mw)
                if slip_mw:
                    # new arc started
    """

    def __init__(self):
        self._arcs = {}       # sv → list of ArcData
        self._current = {}    # sv → ArcData (currently open)

    def add_epoch(self, sv, L1_cyc, L2_cyc, P1_m, P2_m,
                   B_if, mjd_utc, has_slip_mw):
        """Record one epoch of data.  Returns True if a new arc started."""
        mw = compute_mw(L1_cyc, L2_cyc, P1_m, P2_m)

        if has_slip_mw or sv not in self._current:
            # Close previous arc
            self._close_arc(sv)
            # Start new arc
            self._current[sv] = ArcData(sv, mjd_utc)

        arc = self._current[sv]
        arc.mw_vals.append(mw)
        arc.bif_vals.append(B_if)
        arc.mjd_end = mjd_utc
        arc.n_epochs += 1
        return has_slip_mw

    def _close_arc(self, sv):
        if sv in self._current:
            arc = self._current.pop(sv)
            if arc.n_epochs >= 5:
                self._arcs.setdefault(sv, []).append(arc)

    def finalize(self):
        """Close all open arcs and return arc-level ambiguity estimates."""
        for sv in list(self._current.keys()):
            self._close_arc(sv)
        return self._arcs

    def get_arcs(self, sv):
        return self._arcs.get(sv, [])


class ArcData:
    """Data for one continuous tracking arc of a single SV."""
    def __init__(self, sv, mjd_start):
        self.sv = sv
        self.mjd_start = mjd_start
        self.mjd_end = mjd_start
        self.mw_vals = []     # MW per epoch [cycles]
        self.bif_vals = []    # B_if per epoch [m]
        self.n_epochs = 0
        self.N_w = None       # fixed WL integer
        self.N1 = None        # fixed NL integer
        self.B_if_fixed = None  # fixed IF ambiguity [m]


class ArcAmbiguityResolver:
    """Arc-based ambiguity resolution.

    For each slip-free arc, estimates WL integer from arc-averaged MW
    and NL integer from arc-averaged B_if (optionally OSB-corrected).
    """

    def __init__(self, wl_bias=None, nl_bias=None,
                 min_epochs=5, max_wl_std=0.30, max_nl_resid=0.30):
        """
        Args:
            wl_bias: {sv: b_wl_sat} from CODE OSB product (optional)
            nl_bias: {sv: b_nl_sat} from CODE OSB product (optional)
            min_epochs: minimum epochs in arc for fixing
            max_wl_std: maximum MW std for WL fixing [cycles]
            max_nl_resid: maximum NL rounding residual [cycles]
        """
        self.wl_bias = wl_bias or {}
        self.nl_bias = nl_bias or {}
        self.min_epochs = min_epochs
        self.max_wl_std = max_wl_std
        self.max_nl_resid = max_nl_resid

    def resolve(self, arc_dict, pass1_amb=None):
        """Resolve WL+NL integers for all arcs.

        Args:
            arc_dict: {sv: [ArcData, ...]} from ArcTracker.finalize()
            pass1_amb: {sv: B_if_float} from EKF for arcs with <2 epochs

        Returns:
            wl_fixed: {sv: N_w} WL integer per SV
            nl_fixed: {sv: N1}  NL integer per SV
            amb_fixed: {sv: B_if_fixed} IF ambiguity [m] per SV
        """
        wl_fixed = {}
        nl_fixed = {}

        for sv, arcs in arc_dict.items():
            # Use longest arc for this SV
            best_arc = max(arcs, key=lambda a: a.n_epochs)
            if best_arc.n_epochs < self.min_epochs:
                continue

            # ── WL fixing ──
            mw_arr = np.array(best_arc.mw_vals)
            mw_mean = float(np.mean(mw_arr))
            mw_std = float(np.std(mw_arr))

            if mw_std > self.max_wl_std:
                continue

            # Apply OSB WL bias if available
            if sv in self.wl_bias:
                mw_mean_corr = mw_mean - self.wl_bias[sv]
            else:
                mw_mean_corr = mw_mean

            N_w = int(round(mw_mean_corr))
            if abs(mw_mean_corr - N_w) > 0.35:
                continue

            wl_fixed[sv] = N_w

            # ── NL fixing ──
            bif_arr = np.array(best_arc.bif_vals)
            bif_mean = float(np.mean(bif_arr))
            # B_if = λ_nl * (N1 + b_nl) + coeff_w * N_w
            # → N1 = (B_if - coeff_w * N_w) / λ_nl - b_nl
            b_nl = self.nl_bias.get(sv, 0.0)
            N1_float = (bif_mean - COEFF_W * N_w) / LAM_NL - b_nl
            N1_int = int(round(N1_float))
            N1_resid = abs(N1_float - N1_int)

            if N1_resid <= self.max_nl_resid:
                nl_fixed[sv] = N1_int
                best_arc.N_w = N_w
                best_arc.N1 = N1_int
                best_arc.B_if_fixed = LAM_NL * (N1_int + b_nl) + COEFF_W * N_w

        # Build IF ambiguity dict for EKF pass 2
        amb_fixed = {}
        for sv in wl_fixed:
            N_w = wl_fixed[sv]
            if sv in nl_fixed:
                b_nl = self.nl_bias.get(sv, 0.0)
                amb_fixed[sv] = LAM_NL * (nl_fixed[sv] + b_nl) + COEFF_W * N_w
            elif sv in (pass1_amb or {}):
                # WL-only: estimate N1 from EKF float
                bif_float = pass1_amb[sv]
                N1_est = int(round((bif_float - COEFF_W * N_w) / LAM_NL))
                amb_fixed[sv] = LAM_NL * N1_est + COEFF_W * N_w

        return wl_fixed, nl_fixed, amb_fixed
