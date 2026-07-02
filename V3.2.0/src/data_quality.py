"""Automated data quality monitoring and anomaly detection (Phase 24.0).

Layered QC architecture:
  Layer 1: Raw GPS1B / product integrity (pre-EKF)
  Layer 2: Per-SV MW/MP statistics (EKF→Batch bridge)
  Layer 3: Solution quality (post-Batch/GN)
  Layer 4: Aggregate quality score + auto-diagnosis

All functions are stateless — they take data, return diagnostics.
"""

import numpy as np
from collections import defaultdict

# ── Constants ──
F1, F2, C_L = 1575.42e6, 1227.60e6, 299792458.0
LAM1 = C_L / F1  # ~0.190 m
LAM2 = C_L / F2  # ~0.244 m
GAMMA = (F1 / F2) ** 2  # (77/60)^2


# ══════════════════════════════════════════════════════════════════════
# Layer 1: Raw data integrity (pre-EKF)
# ══════════════════════════════════════════════════════════════════════

def check_gps1b_epoch_continuity(gps1b, max_gap_s=120.0):
    """Check GPS1B epoch spacing for gaps.

    Args:
        gps1b: dict {gps_sec: {sv: {obs...}}}
        max_gap_s: maximum allowed gap in seconds before flagging

    Returns:
        (is_continuous, gap_epochs, gap_durations)
        is_continuous: True if no gaps > max_gap_s
        gap_epochs: list of (start_gps_sec, end_gps_sec) for each gap
        gap_durations: list of gap durations in seconds
    """
    epochs = sorted(gps1b.keys())
    if len(epochs) < 2:
        return True, [], []

    diffs = np.diff(epochs)
    gap_idx = np.where(diffs > max_gap_s)[0]
    gap_epochs = [(epochs[i], epochs[i + 1]) for i in gap_idx]
    gap_durations = [float(diffs[i]) for i in gap_idx]
    is_continuous = len(gap_idx) == 0

    if gap_epochs:
        for (s, e), d in zip(gap_epochs, gap_durations):
            print(f"  [QC-Gap] epoch gap {d:.0f}s at GPS {s:.0f} → {e:.0f}")

    return is_continuous, gap_epochs, gap_durations


def check_sv_coverage(gps1b, epochs, min_coverage_pct=0.50,
                       exclude_svs=None):
    """Per-SV epoch coverage fraction across the arc.

    SVs appearing in fewer than min_coverage_pct of epochs are flagged.
    These SVs have too few observations for reliable ambiguity estimation.

    Args:
        gps1b: dict {gps_sec: {sv: {obs...}}}
        epochs: list of GPS seconds in the arc
        min_coverage_pct: minimum fraction of epochs a SV must be present
        exclude_svs: set of SV names to skip (e.g. already segmented)

    Returns:
        {sv: {'n_epoch': int, 'total': int, 'pct': float, 'usable': bool}}
    """
    exclude_svs = exclude_svs or set()
    n_total = len(epochs)
    results = {}
    for sv in _all_svs(gps1b, epochs):
        if sv in exclude_svs:
            continue
        n_present = sum(1 for ep in epochs if sv in gps1b.get(int(ep), gps1b.get(ep, {})))
        pct = n_present / max(n_total, 1)
        usable = pct >= min_coverage_pct
        results[sv] = {
            'n_epoch': n_present, 'total': n_total,
            'pct': pct, 'usable': usable,
        }
    n_unusable = sum(1 for v in results.values() if not v['usable'])
    if n_unusable > 0:
        bad_svs = [sv for sv, v in results.items() if not v['usable']]
        print(f"  [QC-Coverage] {n_unusable}/{len(results)} SVs < {min_coverage_pct*100:.0f}%: {bad_svs}")
    return results


def check_snr(gps1b, epochs, snr_l1_min=30.0, snr_l2_min=25.0,
               flag_pct_threshold=0.30):
    """Per-SV SNR screening.

    Flags SVs where >30% of observations have SNR below threshold.
    Low SNR indicates multipath, obstruction, or low elevation tracking.

    Args:
        gps1b, epochs: as above
        snr_l1_min, snr_l2_min: SNR thresholds in dB-Hz
        flag_pct_threshold: fraction of low-SNR epochs that triggers a flag

    Returns:
        {sv: {'snr_l1_mean': float, 'snr_l2_mean': float,
              'pct_low_l1': float, 'pct_low_l2': float, 'flagged': bool}}
    """
    results = {}
    for sv in _all_svs(gps1b, epochs):
        l1_snr_vals, l2_snr_vals = [], []
        for ep in epochs:
            rec = gps1b.get(int(ep), gps1b.get(ep, {}))
            if sv not in rec:
                continue
            if 'L1_SNR' in rec[sv]:
                l1_snr_vals.append(float(rec[sv]['L1_SNR']))
            if 'L2_SNR' in rec[sv]:
                l2_snr_vals.append(float(rec[sv]['L2_SNR']))

        l1_arr = np.array(l1_snr_vals) if l1_snr_vals else np.array([999])
        l2_arr = np.array(l2_snr_vals) if l2_snr_vals else np.array([999])
        pct_low_l1 = np.mean(l1_arr < snr_l1_min)
        pct_low_l2 = np.mean(l2_arr < snr_l2_min)

        flagged = (pct_low_l1 > flag_pct_threshold or
                    pct_low_l2 > flag_pct_threshold)
        results[sv] = {
            'snr_l1_mean': float(np.mean(l1_arr)),
            'snr_l2_mean': float(np.mean(l2_arr)),
            'pct_low_l1': float(pct_low_l1),
            'pct_low_l2': float(pct_low_l2),
            'flagged': flagged,
        }
    n_flagged = sum(1 for v in results.values() if v['flagged'])
    if n_flagged > 0:
        flagged_svs = [sv for sv, v in results.items() if v['flagged']]
        print(f"  [QC-SNR] {n_flagged} SV(s) with low SNR: {flagged_svs}")
    return results


# ══════════════════════════════════════════════════════════════════════
# Layer 2: Per-SV statistics (EKF → Batch bridge)
# ══════════════════════════════════════════════════════════════════════

def compute_code_multipath(L1, L2, P1, P2):
    """Compute MP1 multipath combination.

    MP1 = P1 - L1 - 2*(L1-L2)/(gamma-1)
    where gamma = (f1/f2)^2

    For dual-frequency GPS: MP1 absorbs code multipath on P1.
    RMS > 0.5m suggests moderate multipath; > 2.0m is severe.

    Returns:
        mp1 value in meters, or NaN if any input is zero/None.
    """
    try:
        L1_m = float(L1)
        L2_m = float(L2)
        P1_m = float(P1)
        P2_m = float(P2)
    except (TypeError, ValueError):
        return np.nan

    if any(x == 0 for x in (L1_m, L2_m, P1_m, P2_m)):
        return np.nan
    # MP1 = P1 - L1 - 2*(f1^2)*(L1-L2)/(f1^2-f2^2)
    # Actually: MP1 = P1 - L1 - 2 * f2^2/(f1^2-f2^2) * (L1-L2)
    # Simpler: MP1 = P1 - L1 - (2/(gamma-1)) * (L1-L2)
    coeff = 2.0 / (GAMMA - 1.0)
    mp1 = P1_m - L1_m - coeff * (L1_m - L2_m)
    return float(mp1)


def compute_per_sv_multipath(gps1b, epochs, mp_threshold=3.5):
    """Compute MP1 RMS per SV across the arc.

    Returns:
        {sv: {'mp1_rms': float, 'mp1_max': float, 'n_obs': int, 'flagged': bool}}
    """
    results = {}
    for sv in _all_svs(gps1b, epochs):
        mp_vals = []
        for ep in epochs:
            rec = gps1b.get(int(ep), gps1b.get(ep, {}))
            if sv not in rec:
                continue
            L1 = rec[sv].get('L1', rec[sv].get('L1_phase'))
            L2 = rec[sv].get('L2', rec[sv].get('L2_phase'))
            P1 = rec[sv].get('P1', rec[sv].get('CA_range'))
            P2 = rec[sv].get('P2')
            if None in (L1, L2, P1, P2):
                continue
            mp1 = compute_code_multipath(L1, L2, P1, P2)
            if not np.isnan(mp1):
                mp_vals.append(mp1)

        if len(mp_vals) < 3:
            results[sv] = {'mp1_rms': 0.0, 'mp1_max': 0.0,
                           'n_obs': len(mp_vals), 'flagged': False}
            continue

        mp_arr = np.array(mp_vals)
        mp_rms = np.sqrt(np.mean(mp_arr**2))
        mp_max = np.max(np.abs(mp_arr))
        results[sv] = {
            'mp1_rms': float(mp_rms),
            'mp1_max': float(mp_max),
            'n_obs': len(mp_vals),
            'flagged': mp_rms > mp_threshold,
        }
    n_flagged = sum(1 for v in results.values() if v['flagged'])
    if n_flagged > 0:
        flagged_svs = [sv for sv, v in results.items() if v['flagged']]
        print(f"  [QC-MP] {n_flagged} SV(s) with MP1 RMS > {mp_threshold}m: {flagged_svs}")
    return results


def compute_mw_stability(mw_list):
    """Analyze MW stability for one SV.

    Args:
        mw_list: list of MW float values across epochs

    Returns:
        (mean, std, is_stable, n_jumps, jump_magnitudes)
        is_stable: True if std < 0.5 cyc and no jumps > 5 cyc
        n_jumps: count of consecutive differences > 5 cyc
    """
    if len(mw_list) < 3:
        return {'mean': 0.0, 'std': 0.0, 'n_points': len(mw_list),
                'is_stable': False, 'is_corrupted': False,
                'n_jumps': 0, 'jump_mags': []}

    arr = np.array(mw_list)
    mean = float(np.mean(arr))
    std = float(np.std(arr))

    n_jumps = 0
    jump_mags = []
    if len(arr) >= 2:
        diffs = np.abs(np.diff(arr))
        jump_idx = np.where(diffs > 5.0)[0]  # >5 cycle jump = gap/corruption
        n_jumps = len(jump_idx)
        jump_mags = [float(diffs[i]) for i in jump_idx]

    # Stable: low std AND no large jumps
    is_stable = (std <= 0.50) and (n_jumps == 0)

    # Corrupted: very large std (>10 cyc) or jumps indicate data gap
    is_corrupted = (std > 10.0) or (n_jumps >= 2)

    return {
        'mean': mean, 'std': std, 'n_points': len(mw_list),
        'is_stable': is_stable, 'is_corrupted': is_corrupted,
        'n_jumps': n_jumps, 'jump_mags': jump_mags,
    }


def screen_sv_for_batch(sv, cov_stats, mw_stats, mp_stats, snr_stats):
    """Decide whether a SV is usable for batch ambiguity estimation.

    Decision logic (priority order):
      1. Coverage < 50% → REJECT (too few epochs)
      2. MW corrupted (std > 10 cyc or >= 2 jumps) → REJECT (gap/cycle slip)
      3. MW unstable (std > 0.30 cyc) → REJECT_LOW_QUALITY (noisy, discard for AR)
      4. MW jump detected (1 jump > 5 cyc) → SEGMENT (split at gap, if not done)
      5. MP1 RMS > 2.0m → FLAG_MULTIPATH (keep but note)
      6. SNR flagged → FLAG_SNR (keep but note)
      7. Otherwise → KEEP

    Returns:
        {'decision': 'KEEP'|'REJECT'|'SEGMENT'|'FLAG_*',
         'reason': str, 'wl_eligible': bool}
    """
    # Coverage check
    cov = cov_stats.get(sv, {})
    if not cov.get('usable', True):
        return {'decision': 'REJECT',
                'reason': f'coverage={cov.get("pct", 0)*100:.0f}%',
                'wl_eligible': False}

    # MW checks
    mw = mw_stats.get(sv, {})
    if mw.get('is_corrupted', False):
        return {'decision': 'REJECT',
                'reason': f'MW corrupted (std={mw.get("std", 0):.1f}cyc jumps={mw.get("n_jumps", 0)})',
                'wl_eligible': False}

    n_jumps = mw.get('n_jumps', 0)
    if n_jumps == 1 and not mw.get('is_corrupted', False):
        return {'decision': 'SEGMENT',
                'reason': f'MW single jump ({mw.get("jump_mags", [0])[0]:.0f}cyc) → segment',
                'wl_eligible': True}

    if not mw.get('is_stable', False):
        mw_std = mw.get('std', 99)
        if mw_std > 0.35:
            return {'decision': 'REJECT',
                    'reason': f'MW unstable (std={mw_std:.2f}cyc)',
                    'wl_eligible': False}
        else:
            return {'decision': 'REJECT',
                    'reason': f'MW marginal (std={mw_std:.2f}cyc)',
                    'wl_eligible': False}

    # Multipath check (flag only — keep the SV)
    mp = mp_stats.get(sv, {})
    if mp.get('flagged', False):
        return {'decision': 'KEEP',
                'reason': f'MP1={mp.get("mp1_rms", 0):.2f}m (flagged)',
                'wl_eligible': True}

    # SNR check (flag only — keep the SV)
    snr = snr_stats.get(sv, {})
    if snr.get('flagged', False):
        return {'decision': 'KEEP',
                'reason': f'SNR L1={snr.get("snr_l1_mean", 0):.0f} dB-Hz (low)',
                'wl_eligible': True}

    return {'decision': 'KEEP', 'reason': 'OK', 'wl_eligible': True}


# ══════════════════════════════════════════════════════════════════════
# Layer 3: Product integrity
# ══════════════════════════════════════════════════════════════════════

def check_sp3_integrity(sp3_data):
    """Validate SP3 product completeness.

    Expects the canonical pickle format: {'epochs': {datetime: {sv: [x,y,z,clk]}}, 'ts': [...]}
    or the legacy format with epochs dict keyed by datetime.

    Returns:
        {'n_epoch': int, 'n_sv': int, 'sv_list': list,
         'is_valid': bool, 'issues': list}
    """
    issues = []
    epochs = sp3_data.get('epochs', sp3_data)
    if not isinstance(epochs, dict):
        return {'n_epoch': 0, 'n_sv': 0, 'sv_list': [],
                'is_valid': False, 'issues': ['SP3 format unrecognized']}

    epoch_keys = sorted(epochs.keys())
    n_epoch = len(epoch_keys)
    sv_set = set()
    for ek in epoch_keys:
        sv_set.update(epochs[ek].keys())

    if n_epoch < 280:
        issues.append(f'SP3 epochs={n_epoch} (<280 expected)')
    gps_svs = sorted([s for s in sv_set if s.startswith('G')])
    if len(gps_svs) < 28:
        issues.append(f'GPS SVs={len(gps_svs)} (<28 expected)')

    return {
        'n_epoch': n_epoch, 'n_sv': len(sv_set),
        'sv_list': sorted(sv_set),
        'is_valid': len(issues) == 0,
        'issues': issues,
    }


def check_clk_integrity(clk_data):
    """Validate CLK product completeness.

    clk_data is the output of read_rinex_clk(): dict keyed by SV name (e.g. 'G01'),
    each value is {'epochs': [...], 'clk': [...], 'epochs_sod': [...]}.

    Returns:
        {'n_sv': int, 'n_gps_sv': int, 'n_epochs_per_sv': int,
         'is_valid': bool, 'issues': list}
    """
    issues = []
    if not isinstance(clk_data, dict) or len(clk_data) == 0:
        return {'n_sv': 0, 'n_gps_sv': 0, 'n_epochs_per_sv': 0,
                'is_valid': False, 'issues': ['CLK data empty']}

    gps_svs = sorted([s for s in clk_data.keys() if s.startswith('G')])
    n_gps_sv = len(gps_svs)

    # Get epoch count from first GPS SV
    n_epochs_per_sv = 0
    for sv in gps_svs:
        if isinstance(clk_data[sv], dict) and 'epochs' in clk_data[sv]:
            n_epochs_per_sv = max(n_epochs_per_sv, len(clk_data[sv]['epochs']))

    if n_gps_sv < 28:
        issues.append(f'GPS SVs={n_gps_sv} (<28 expected)')
    if n_epochs_per_sv < 2800:
        issues.append(f'CLK epochs/SV={n_epochs_per_sv} (<2800 expected)')

    return {
        'n_sv': len(clk_data),
        'n_gps_sv': n_gps_sv, 'n_epochs_per_sv': n_epochs_per_sv,
        'is_valid': len(issues) == 0,
        'issues': issues,
    }


# ══════════════════════════════════════════════════════════════════════
# Layer 4: Aggregate quality score + report
# ══════════════════════════════════════════════════════════════════════

def compute_quality_score(cov_stats, mw_stats, snr_stats, mp_stats,
                           sp3_ok, clk_ok, batch_phase_rms, n_gap_svs,
                           n_sv_total):
    """Compute 0-1 quality score for an arc.

    Weights reflect impact on 3D RMS:
      - MW stability (0.30): strongest predictor of AR quality
      - SV coverage (0.20): data completeness
      - Gap freedom (0.15): phase continuity
      - Product integrity (0.15): SP3+CLK valid
      - SNR quality (0.10): observation noise
      - Batch residual (0.10): solution consistency

    Returns float 0.0 (worst) to 1.0 (best).
    """
    score = 0.0
    details = {}

    # MW stability score
    n_mw_total = max(len(mw_stats), 1)
    n_stable = sum(1 for v in mw_stats.values()
                   if v.get('is_stable', False))
    mw_frac = n_stable / n_mw_total
    score += 0.30 * mw_frac
    details['mw_stable'] = f'{n_stable}/{n_mw_total}'

    # Coverage score
    if cov_stats:
        n_usable = sum(1 for v in cov_stats.values() if v.get('usable', False))
        cov_frac = n_usable / max(len(cov_stats), 1)
    else:
        cov_frac = 1.0  # no stats → assume OK
    score += 0.20 * cov_frac
    details['coverage'] = f'{cov_frac:.0%}'

    # Gap freedom score
    n_clean_svs = max(n_sv_total - n_gap_svs, 0)
    gap_score = n_clean_svs / max(n_sv_total, 1)
    score += 0.15 * gap_score
    details['gaps'] = f'{n_gap_svs} gap SV(s)'

    # Product integrity
    prod_score = (1.0 if sp3_ok.get('is_valid', False) else 0.5) * 0.5
    prod_score += (1.0 if clk_ok.get('is_valid', False) else 0.5) * 0.5
    score += 0.15 * prod_score
    details['products'] = f'SP3={sp3_ok.get("is_valid","?")} CLK={clk_ok.get("is_valid","?")}'

    # SNR score
    if snr_stats:
        n_good_snr = sum(1 for v in snr_stats.values() if not v.get('flagged', False))
        snr_frac = n_good_snr / max(len(snr_stats), 1)
    else:
        snr_frac = 1.0
    score += 0.10 * snr_frac
    details['snr'] = f'{snr_frac:.0%}'

    # Batch residual score: <0.5m = good, >5m = bad
    if np.isfinite(batch_phase_rms):
        if batch_phase_rms < 0.3:
            resid_score = 1.0
        elif batch_phase_rms < 0.8:
            resid_score = 0.7
        elif batch_phase_rms < 2.0:
            resid_score = 0.4
        elif batch_phase_rms < 10.0:
            resid_score = 0.2
        else:
            resid_score = 0.0
    else:
        resid_score = 0.0
    score += 0.10 * resid_score
    details['batch_resid'] = f'{batch_phase_rms:.3f}m'

    return float(np.clip(score, 0.0, 1.0)), details


def generate_qc_report(checks):
    """Generate a structured QC report from accumulated check results.

    Args:
        checks: dict with keys from each QC layer:
            'cov_stats', 'snr_stats', 'mp_stats', 'mw_stats',
            'sp3_ok', 'clk_ok', 'batch_phase_rms',
            'n_gap_svs', 'n_sv_total', 'sv_decision'

    Returns:
        {'score': float, 'grade': str, 'flags': list, 'actions': list,
         'details': dict}
    """
    score, details = compute_quality_score(
        cov_stats=checks.get('cov_stats', {}),
        mw_stats=checks.get('mw_stats', {}),
        snr_stats=checks.get('snr_stats', {}),
        mp_stats=checks.get('mp_stats', {}),
        sp3_ok=checks.get('sp3_ok', {'is_valid': True}),
        clk_ok=checks.get('clk_ok', {'is_valid': True}),
        batch_phase_rms=checks.get('batch_phase_rms', 0.0),
        n_gap_svs=checks.get('n_gap_svs', 0),
        n_sv_total=checks.get('n_sv_total', 1),
    )

    # Grade
    if score >= 0.85:
        grade = 'A'
    elif score >= 0.70:
        grade = 'B'
    elif score >= 0.50:
        grade = 'C'
    elif score >= 0.30:
        grade = 'D'
    else:
        grade = 'F'

    # Collect flags
    flags = []
    if checks.get('n_gap_svs', 0) > 0:
        flags.append(f'{checks["n_gap_svs"]} gap-SV(s)')
    n_rejected = sum(1 for v in checks.get('sv_decision', {}).values()
                     if v.get('decision') == 'REJECT')
    if n_rejected > 0:
        flags.append(f'{n_rejected} rejected SV(s)')
    if not checks.get('sp3_ok', {}).get('is_valid', True):
        flags.append('SP3!')
    if not checks.get('clk_ok', {}).get('is_valid', True):
        flags.append('CLK!')

    # Actions taken
    actions = []
    if checks.get('n_gap_svs', 0) > 0:
        actions.append(f'{checks["n_gap_svs"]} SV(s) segmented (gap split)')
    if n_rejected > 0:
        actions.append(f'{n_rejected} SV(s) rejected (MW/coverage)')
    if not actions:
        actions.append('none')

    return {
        'score': score, 'grade': grade,
        'flags': flags, 'actions': actions,
        'details': details,
    }


def print_qc_line(report):
    """Print a one-line QC summary."""
    flag_str = ', '.join(report['flags']) if report['flags'] else 'clean'
    print(f"  [QC] score={report['score']:.2f}({report['grade']}) "
          f"| flags: {flag_str} "
          f"| actions: {', '.join(report['actions'])}")


# ══════════════════════════════════════════════════════════════════════
# Helpers
# ══════════════════════════════════════════════════════════════════════

def _all_svs(gps1b, epochs):
    """Get set of all unique SV names in the arc."""
    svs = set()
    for ep in epochs:
        recs = gps1b.get(int(ep), gps1b.get(ep, {}))
        svs.update(recs.keys())
    return sorted(svs)
