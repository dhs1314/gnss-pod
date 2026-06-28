"""
Integer ambiguity resolution: MW wide-lane + narrow-lane

GPS1A phase is in raw cycles (integer ambiguity preserved).
MW gives N_w (N1-N2). Arc-averaged B_if gives N1.

Reference: Teunissen (1995) The LAMBDA method
"""
import numpy as np
from collections import defaultdict

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
F1_SQ, F2_SQ = F1 * F1, F2 * F2
LAM1 = C / F1
LAM2 = C / F2
LAM_W = C / (F1 - F2)    # ~0.862 m  (wide-lane)
LAM_N = C / (F1 + F2)    # ~0.107 m  (narrow-lane)
ALPHA = F1_SQ / (F1_SQ - F2_SQ)   # ~2.546
BETA = -F2_SQ / (F1_SQ - F2_SQ)   # ~-1.546
# Coefficient: B_if = lam_n * N1 + coeff_w * N_w
COEFF_W = C * F2 / (F1_SQ - F2_SQ)  # ~0.3776 m/cycle


def compute_mw(L1_cyc, L2_cyc, P1_m, P2_m):
    """Melbourne-Wubbena wide-lane in cycles.

    N_w = L1_cyc - L2_cyc - (f1*P1 + f2*P2)/((f1+f2)*lambda_w)
    """
    wl_phase_cyc = L1_cyc - L2_cyc
    nl_code_m = (F1 * P1_m + F2 * P2_m) / (F1 + F2)
    return wl_phase_cyc - nl_code_m / LAM_W


def compute_if_m(L1_cyc, L2_cyc, P1_m, P2_m):
    """Ionosphere-free combinations in meters."""
    L1_m = L1_cyc * LAM1
    L2_m = L2_cyc * LAM2
    L_if = ALPHA * L1_m + BETA * L2_m
    P_if = ALPHA * P1_m + BETA * P2_m
    B_if = L_if - P_if
    return L_if, P_if, B_if


def estimate_wide_lane(arc_epochs):
    """Estimate N_w integer for one arc."""
    mw_vals = [compute_mw(d['L1_cyc'], d['L2_cyc'], d['P1'], d['P2'])
               for d in arc_epochs]
    mw_arr = np.array(mw_vals)
    N_w_float = float(np.mean(mw_arr))
    sigma = float(np.std(mw_arr))
    n_obs = len(mw_arr)

    if sigma > 0.4 or n_obs < 5:
        return None, N_w_float, sigma, n_obs

    N_w_fixed = int(round(N_w_float))
    return N_w_fixed, N_w_float, sigma, n_obs


def narrow_lane_from_if(B_if_mean, B_if_std, N_w_fixed):
    """Derive narrow-lane N1 from IF ambiguity and fixed wide-lane.

    B_if_mean: mean IF ambiguity over arc (meters)
    B_if_std:  std of IF ambiguity over arc (meters)
    N_w_fixed: fixed wide-lane integer (cycles)

    B_if = lam_n * N1 + coeff_w * N_w
    → N1 = (B_if - coeff_w * N_w) / lam_n
    """
    N1_float = (B_if_mean - COEFF_W * N_w_fixed) / LAM_N
    N1_std = B_if_std / LAM_N   # Uncertainty in cycles

    if N1_std > 0.4:
        N1_fixed = None
    else:
        N1_fixed = int(round(N1_float))
        if abs(N1_float - N1_fixed) > 0.35:
            N1_fixed = None

    return N1_float, N1_fixed, N1_std


def b_if_from_ints(N_w, N1):
    """Compute IF ambiguity in meters from N_w and N1 integers."""
    return LAM_N * N1 + COEFF_W * N_w


def segment_arcs(sv_obs_sorted, max_gap_s=15.0):
    """Segment SV observations into continuous tracking arcs."""
    arcs = []
    current_arc = []

    for obs in sv_obs_sorted:
        is_break = False
        if current_arc:
            gap = obs['gps_sod'] - current_arc[-1]['gps_sod']
            if gap > max_gap_s:
                is_break = True
        if obs.get('slip_L1') or obs.get('slip_L2'):
            is_break = True

        if is_break and current_arc:
            arcs.append(current_arc)
            current_arc = []
        current_arc.append(obs)

    if current_arc:
        arcs.append(current_arc)
    return arcs
