"""Standard POD I/O abstraction (Phase 13.0).

Decouples satellite-specific data loading from the core EKF algorithm.
Supports:
  - GRACE-FO (GPS1B .rnx/.pkl + GNV1B reference)
  - Generic RINEX observation files (any LEO satellite)
  - CODE/IGS precision products (SP3, CLK, DCB, ANTEX, IERS)
"""

import os, pickle
from pathlib import Path
from datetime import datetime


# ═══════════════════════════════════════════════════════════════
# Standard data paths — auto-detection by mission
# ═══════════════════════════════════════════════════════════════

def find_rinex_obs(data_root, mission, sat_id, date_str):
    """Find RINEX observation file for a given satellite and date.

    Priority: .pkl (pre-processed) > .rnx (raw RINEX)

    Returns:
        Path to observation file, or None if not found.
    """
    dp = Path(data_root)

    # GRACE-FO naming: GPS1B_YYYY-MM-DD_X_04.rnx
    if mission.upper() in ('GRACE-FO', 'GRACE'):
        patterns = [
            dp / f"GPS1B_{date_str}_{sat_id}_04.pkl",
            dp / f"GPS1B_{date_str}_{sat_id}_04.rnx",
            dp / 'gracefo' / date_str[:4] / date_str / f"GPS1B_{date_str}_{sat_id}_04.pkl",
        ]
        for p in patterns:
            if p.exists():
                return p

    # Generic RINEX naming: SITE*_YYYYDDD*.rnx or SITE*_YYYYDDD*.YYo
    y, m, d = [int(x) for x in date_str.split('-')]
    doy = (datetime(y, m, d) - datetime(y, 1, 1)).days + 1
    generic_patterns = [
        dp / f"{sat_id}*{y}{doy:03d}*.rnx",
        dp / f"*{y}{doy:03d}*_{sat_id}_*.rnx",
        dp / f"{mission}*_{date_str}*.rnx",
    ]
    import glob
    for pat in generic_patterns:
        matches = glob.glob(str(pat))
        if matches:
            return Path(matches[0])

    return None


def find_reference_orbit(data_root, mission, sat_id, date_str):
    """Find reference/precise orbit for a satellite.

    GRACE-FO: GNV1B within gracefo/YYYY/YYYY-MM-DD/
    SWARM:    SP3 precise orbit from ESA
    Generic:  external SP3 file

    Returns:
        (file_path, orbit_type) where orbit_type is 'gnv1b' or 'sp3'
    """
    dp = Path(data_root)
    y, m, d = [int(x) for x in date_str.split('-')]

    if mission.upper() in ('GRACE-FO', 'GRACE'):
        gnv_path = dp / 'gracefo' / str(y) / date_str / f'GNV1B_{date_str}_{sat_id}_04.txt'
        if gnv_path.exists():
            return gnv_path, 'gnv1b'
        gnv_path2 = dp / f'GNV1B_{date_str}_{sat_id}_04.txt'
        if gnv_path2.exists():
            return gnv_path2, 'gnv1b'

    # Generic SP3
    sp3_pat = dp / f"*{sat_id}*{y}{doy:03d}*.sp3"
    matches = glob.glob(str(sp3_pat)) if 'glob' in dir() else []
    if matches:
        return Path(matches[0]), 'sp3'

    return None, None


def load_precision_products(data_root, year, doy):
    """Load CODE precision products (SP3, CLK, DCB, ANTEX, IERS).

    Returns dict with keys: sp3, clk_data, antex_data, dcb_data, iers_ok
    """
    dp = Path(data_root)
    products = {}

    # SP3
    sp3_pkl = dp / 'CODE' / str(year) / f"COD0OPSFIN_{year}{doy:03d}0000_01D_05M_ORB.pkl"
    sp3_txt = dp / 'CODE' / str(year) / f"COD0OPSFIN_{year}{doy:03d}0000_01D_05M_ORB.SP3"
    if sp3_pkl.exists():
        products['sp3'] = pickle.load(open(str(sp3_pkl), 'rb'))
    elif sp3_txt.exists():
        from src.sp3_loader import parse_sp3_text
        with open(str(sp3_txt), 'r') as f:
            epochs_dict, ts_list = parse_sp3_text(f.read())
        products['sp3'] = {'ts': ts_list, 'epochs': epochs_dict, 'source': 'CODE', 'product': 'FIN'}
    else:
        products['sp3'] = None

    # CLK
    clk_txt = dp / 'CODE' / str(year) / f"COD0OPSFIN_{year}{doy:03d}0000_01D_30S_CLK.CLK"
    products['clk_data'] = None
    if clk_txt.exists():
        from src.precision_products import read_rinex_clk
        products['clk_data'] = read_rinex_clk(str(clk_txt))

    # ANTEX
    antex_path = dp / 'igs14.atx'
    products['antex_data'] = None
    if antex_path.exists():
        from src.precision_products import read_antex
        products['antex_data'] = read_antex(str(antex_path))

    # DCB
    products['dcb_data'] = None
    dcb_p1p2 = dp / 'CODE' / str(year) / f"P1P2{doy:03d}.DCB"
    dcb_p1c1 = dp / 'CODE' / str(year) / f"P1C1{doy:03d}.DCB"
    if dcb_p1p2.exists() and dcb_p1c1.exists():
        from src.precision_products import load_code_dcb_pair
        products['dcb_data'] = load_code_dcb_pair(str(dcb_p1c1), str(dcb_p1p2))

    # IERS
    iers_path = dp / 'IERS' / 'eopc04_IAU2000.txt'
    if iers_path.exists():
        from src.precision_products import setup_iers_from_c04
        setup_iers_from_c04(str(iers_path))
        products['iers_ok'] = True
    else:
        products['iers_ok'] = False

    return products


# ═══════════════════════════════════════════════════════════════
# Standard metrics
# ═══════════════════════════════════════════════════════════════

def compute_metrics(ecef_positions, ref_positions, epochs):
    """Compute standard POD accuracy metrics.

    Returns dict with:
      rms_3d, mean_3d, max_3d, std_3d
      rms_r, rms_a, rms_c (radial/along-track/cross-track)
      n_epochs, n_gaps
    """
    import numpy as np

    dr_vals = []
    for i in range(len(ecef_positions)):
        if ref_positions[i] is not None:
            dr_vals.append(np.linalg.norm(
                np.array(ecef_positions[i]) - np.array(ref_positions[i])))

    if not dr_vals:
        return {'rms_3d': 0, 'n_epochs': 0, 'error': 'no valid epochs'}

    dr = np.array(dr_vals)
    return {
        'rms_3d': float(np.sqrt(np.mean(dr**2))),
        'mean_3d': float(np.mean(dr)),
        'max_3d': float(np.max(dr)),
        'std_3d': float(np.std(dr)),
        'median_3d': float(np.median(dr)),
        'n_epochs': len(dr),
        'p50': float(np.percentile(dr, 50)),
        'p95': float(np.percentile(dr, 95)),
    }


# ═══════════════════════════════════════════════════════════════
# Summary report generation
# ═══════════════════════════════════════════════════════════════

def print_summary(results_dict):
    """Print a formatted summary table from multi-date results."""
    print(f"\n{'='*80}")
    print(f"POD Accuracy Summary")
    print(f"{'='*80}")
    print(f"{'Date':<12} {'Arc':>6} {'RMS_3D':>8} {'Mean':>8} {'Max':>8} "
          f"{'<0.5m':>6} {'n_ep':>5} {'SVs':>4}")
    print(f"{'-'*80}")

    all_rms = []
    for date_str, arcs in sorted(results_dict.items()):
        for arc_h, metrics in sorted(arcs.items()):
            rms = metrics.get('rms_3d', 0)
            all_rms.append(rms)
            n_lt05 = metrics.get('lt_05m_windows', 0)
            n_ep = metrics.get('n_epochs', 0)
            n_sv = metrics.get('n_sv', 0)
            print(f"{date_str:<12} {arc_h:>5.2f}h {rms:>8.3f} "
                  f"{metrics.get('mean_3d', 0):>8.3f} "
                  f"{metrics.get('max_3d', 0):>8.3f} "
                  f"{n_lt05:>6} {n_ep:>5} {n_sv:>4}")

    if all_rms:
        print(f"{'-'*80}")
        print(f"{'OVERALL':<12} {'':>6} {np.mean(all_rms):>8.3f} "
              f"{'':>8} {'':>8} {'':>6} {'':>5} {'':>4}")
        print(f"  n={len(all_rms)}  median={np.median(all_rms):.3f}  "
              f"best={np.min(all_rms):.3f}  worst={np.max(all_rms):.3f}")
