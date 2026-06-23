"""Cycle slip detection using TurboEdit algorithm (Blewitt, 1990).

Melbourne-Wübbena (MW) wide-lane and Geometry-Free (GF) combinations
for detecting cycle slips in dual-frequency GPS data.
"""
import numpy as np

# GPS frequencies
F1 = 1575.42e6   # L1 [Hz]
F2 = 1227.60e6   # L2 [Hz]
WL_WAVELENGTH = 299792458.0 / (F1 - F2)  # ~0.862 m (wide-lane)


def mw_combination(L1, L2, P1, P2):
    """Compute Melbourne-Wübbena wide-lane combination [cycles].

    MW = (L1 - L2) - (f1*P1 + f2*P2) / (f1 + f2) * (f1 - f2) / c

    The MW combination eliminates geometry, clock, troposphere, and
    ionosphere (first order), leaving only wide-lane ambiguity and
    multipath/noise.

    Args:
        L1, L2: carrier phase [m]
        P1, P2: pseudorange [m]

    Returns:
        MW value in wide-lane cycles
    """
    # Narrow-lane pseudorange in metres
    nl_pr = (F1 * P1 + F2 * P2) / (F1 + F2)
    # Wide-lane phase in metres
    wl_phase = (F1 * L1 - F2 * L2) / (F1 - F2)
    # MW = (wl_phase - nl_pr) / wl_wavelength
    return (wl_phase - nl_pr) / WL_WAVELENGTH


def gf_combination(L1, L2):
    """Compute Geometry-Free combination [m].

    GF = L1 - L2

    Eliminates geometry, clock, and troposphere. Changes indicate
    either ionospheric variation or cycle slips. Time-differencing
    removes the ionospheric trend, leaving cycle slip jumps.

    Args:
        L1, L2: carrier phase [m]

    Returns:
        GF value in metres
    """
    return L1 - L2


class TurboEdit:
    """TurboEdit cycle slip detector with sliding-window statistics.

    Usage:
        te = TurboEdit(window_size=50, mw_threshold=4.0, gf_threshold=0.02)

        for each epoch:
            for each SV:
                mw = te.mw_combination(L1, L2, P1, P2)
                gf = te.gf_combination(L1, L2)
                has_slip = te.detect(sv, mw, gf)
                if has_slip:
                    # reset ambiguity for this SV
    """

    def __init__(self, window_size=50, mw_threshold=4.0, gf_threshold=0.02):
        """
        Args:
            window_size: number of epochs for sliding median/std window
            mw_threshold: MW jump threshold [wide-lane cycles]
                          (4.0 * sigma_MW ≈ 4.0 * 0.5 = 2.0 cycles typical)
            gf_threshold: GF time-difference threshold [m]
                          (0.02 m/s typical for quiet ionosphere)
        """
        self.window_size = window_size
        self.mw_threshold = mw_threshold
        self.gf_threshold = gf_threshold

        self._mw_buf = {}    # sv → deque of MW values
        self._gf_buf = {}    # sv → deque of GF values
        self._mw_mean = {}   # sv → running MW mean
        self._mw_std = {}    # sv → running MW std

    def detect(self, sv, L1, L2, P1, P2):
        """Detect cycle slip for one SV at one epoch.

        Returns True if a cycle slip is detected.
        """
        mw = mw_combination(L1, L2, P1, P2)
        gf = gf_combination(L1, L2)
        return self._detect_from_values(sv, mw, gf)

    def _detect_from_values(self, sv, mw, gf):
        """Core detection using pre-computed MW and GF values."""
        slip_mw = False
        slip_gf = False

        # ── MW detection ──
        if sv not in self._mw_buf:
            self._mw_buf[sv] = []
            self._gf_buf[sv] = []
            self._mw_buf[sv].append(mw)
            self._gf_buf[sv].append(gf)
            return False  # need at least 2 epochs

        buf_mw = self._mw_buf[sv]
        buf_gf = self._gf_buf[sv]

        if len(buf_mw) >= 5:
            # Use median as robust estimator of MW ambiguity
            mw_median = float(np.median(buf_mw))
            mw_std = max(float(np.std(buf_mw)), 0.25)  # floor at 0.25 cycles
            if abs(mw - mw_median) > self.mw_threshold * mw_std:
                slip_mw = True

        # ── GF time-difference detection ──
        if len(buf_gf) >= 2:
            gf_prev = buf_gf[-1]
            gf_diff = abs(gf - gf_prev)
            if gf_diff > self.gf_threshold:
                slip_gf = True

        # Update buffers
        buf_mw.append(mw)
        buf_gf.append(gf)
        if len(buf_mw) > self.window_size:
            buf_mw.pop(0)
        if len(buf_gf) > self.window_size:
            buf_gf.pop(0)

        return slip_mw or slip_gf

    def detect_slip_L1(self, sv, L1, L2, P1, P2):
        """Convenience: detect cycle slip, return True/False.

        Also returns estimated slip magnitude on L1 and L2 [cycles] if
        slip detected (for ambiguity reset). Otherwise (0, 0).
        """
        mw = mw_combination(L1, L2, P1, P2)
        gf = gf_combination(L1, L2)
        has_slip = self._detect_from_values(sv, mw, gf)

        if has_slip and sv in self._mw_buf and len(self._mw_buf[sv]) >= 2:
            # Estimate jump sizes
            # MW jump = wide-lane slip: dN_wl = dN1 - dN2
            # GF jump / lambda_1 ≈ dN1 - (lambda_2/lambda_1) * dN2
            lambda_1 = 299792458.0 / F1
            lambda_2 = 299792458.0 / F2
            mw_buf = self._mw_buf[sv]
            gf_buf = self._gf_buf[sv]
            dmw = mw - mw_buf[-2]
            dgf = gf - gf_buf[-2]

            # Solve: dMW = dN1 - dN2,  dGF ≈ lambda_1*dN1 - lambda_2*dN2
            # [1, -1; lambda_1, -lambda_2] * [dN1, dN2] = [dMW, dGF/lambda_1]
            dN1 = (dmw * lambda_2 - dgf) / (lambda_2 - lambda_1)
            dN2 = dN1 - dmw
            return True, round(dN1), round(dN2)

        return has_slip, 0, 0

    def reset_sv(self, sv):
        """Clear buffers for a satellite (e.g., after ambiguity reset)."""
        self._mw_buf.pop(sv, None)
        self._gf_buf.pop(sv, None)
