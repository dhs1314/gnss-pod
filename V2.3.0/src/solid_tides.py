"""Solid Earth tide corrections to geopotential coefficients (IERS 2010).

The tide-generating potential from Sun and Moon deforms the Earth, changing
the Stokes coefficients Cnm, Snm. The frequency-independent (elastic) model
uses nominal Love numbers knm.

References:
  - IERS Conventions (2010), Section 6.2
  - McCarthy & Petit (2004), IERS Technical Note 32
"""
import numpy as np

# Earth constants
GM_EARTH = 3.986004415e14
R_EARTH = 6378136.6

# IERS 2010 nominal Love numbers (frequency-independent)
# Table 6.1, elastic Earth
K20 = 0.30190
K21 = 0.29830
K22 = 0.30102
K30 = 0.093
K31 = 0.093
K32 = 0.093
K33 = 0.094

# Secular rate of C20 from post-glacial rebound (IERS 2010, Eq 6.12)
# For GGM05C (tide-free), the zero-frequency term must be added
C20_DOT = 1.163e-11  # per year
C30_DOT = 4.9e-12   # per year (approximate)
C40_DOT = 5.4e-12   # per year (approximate)


def _legendre_n2(sin_phi):
    """Associated Legendre functions P_2m(sin_phi)."""
    cos_phi = np.sqrt(1.0 - sin_phi ** 2)
    p20 = 0.5 * (3.0 * sin_phi ** 2 - 1.0)          # P20
    p21 = 3.0 * sin_phi * cos_phi                    # P21
    p22 = 3.0 * cos_phi ** 2                          # P22
    return np.sqrt(5.0) * p20, np.sqrt(5.0 / 3.0) * p21, np.sqrt(5.0 / 12.0) * p22


def _legendre_n3(sin_phi):
    """Associated Legendre functions P_3m(sin_phi)."""
    cos_phi = np.sqrt(1.0 - sin_phi ** 2)
    p30 = 0.5 * sin_phi * (5.0 * sin_phi ** 2 - 3.0)
    p31 = 1.5 * cos_phi * (5.0 * sin_phi ** 2 - 1.0)
    p32 = 15.0 * sin_phi * cos_phi ** 2
    p33 = 15.0 * cos_phi ** 3
    return p30, p31, p32, p33


def compute_solid_tide_corrections(mjd_utc, mjd_tt):
    """Compute solid Earth tide corrections to low-degree Stokes coefficients.

    Computes ΔCnm, ΔSnm for n=2 using frequency-independent IERS 2010 model.
    Sun and Moon positions are obtained via astropy (GCRS → ITRF).

    Args:
        mjd_utc: MJD in UTC (for ECEF rotation)
        mjd_tt: MJD in TT (for sun/moon position)

    Returns:
        dict with keys 'C20','C21','S21','C22','S22','C30','C31','S31',...
             values in fully-normalized units
    """
    from src.third_body import get_body_position_eci
    from src.coordinates import ecef_to_eci, eci_to_ecef

    corrections = {}

    for body_name in ['Moon', 'Sun']:
        if body_name == 'Moon':
            GMj = 4.90280031e12
        else:
            GMj = 1.32712440018e20

        # Get body position in ECI (GCRS)
        pos_body_eci = get_body_position_eci(body_name, mjd_tt)
        # Convert to ECEF for spherical coordinates relative to Earth
        pos_body_ecef, _ = eci_to_ecef(pos_body_eci, np.zeros(3), mjd_utc)

        rj = np.linalg.norm(pos_body_ecef)
        if rj < 1.0:
            continue

        sin_phi = pos_body_ecef[2] / rj
        lon = np.arctan2(pos_body_ecef[1], pos_body_ecef[0])

        # Scaling factor: (GMj/GM) * (Re/rj)^3 for n=2, (Re/rj)^4 for n=3
        g_factor = GMj / GM_EARTH
        r_factor_2 = (R_EARTH / rj) ** 3
        r_factor_3 = (R_EARTH / rj) ** 4

        # ── Degree 2 ──
        p2 = _legendre_n2(sin_phi)
        # Normalized Legendre functions: P̄_nm = N_nm * P_nm
        # where N_nm = sqrt((2-δ0m)*(2n+1)*(n-m)!/(n+m)!)
        # We use the already-normalized forms
        scaling_20 = 1.0 / 5.0       # k2m / (2n+1) = k20/5
        scaling_21 = 1.0 / 5.0       # k21/5
        scaling_22 = 1.0 / 5.0       # k22/5

        factor = g_factor * r_factor_2
        # P̄₂₀ = √5 * P₂₀, so ΔC̄₂₀ = k20/5 * factor * P̄₂₀
        # Using complex formulation: ΔCnm = knm/(2n+1) * (GMj/GM) * (Re/r)^(n+1) * Pnm * cos(mλ)
        corrections['C20'] = corrections.get('C20', 0.0) + K20 * factor * p2[0] / (2 * 2 + 1)
        corrections['C21'] = corrections.get('C21', 0.0) + K21 * factor * p2[1] * np.cos(lon) / (2 * 2 + 1)
        corrections['S21'] = corrections.get('S21', 0.0) + K21 * factor * p2[1] * np.sin(lon) / (2 * 2 + 1)
        corrections['C22'] = corrections.get('C22', 0.0) + K22 * factor * p2[2] * np.cos(2 * lon) / (2 * 2 + 1)
        corrections['S22'] = corrections.get('S22', 0.0) + K22 * factor * p2[2] * np.sin(2 * lon) / (2 * 2 + 1)

        # ── Degree 3 ──
        p3 = _legendre_n3(sin_phi)
        factor3 = g_factor * r_factor_3
        corrections['C30'] = corrections.get('C30', 0.0) + K30 * factor3 * p3[0] / (2 * 3 + 1)
        corrections['C31'] = corrections.get('C31', 0.0) + K31 * factor3 * p3[1] * np.cos(lon) / (2 * 3 + 1)
        corrections['S31'] = corrections.get('S31', 0.0) + K31 * factor3 * p3[1] * np.sin(lon) / (2 * 3 + 1)
        corrections['C32'] = corrections.get('C32', 0.0) + K32 * factor3 * p3[2] * np.cos(2 * lon) / (2 * 3 + 1)
        corrections['S32'] = corrections.get('S32', 0.0) + K32 * factor3 * p3[2] * np.sin(2 * lon) / (2 * 3 + 1)
        corrections['C33'] = corrections.get('C33', 0.0) + K33 * factor3 * p3[3] * np.cos(3 * lon) / (2 * 3 + 1)
        corrections['S33'] = corrections.get('S33', 0.0) + K33 * factor3 * p3[3] * np.sin(3 * lon) / (2 * 3 + 1)

    # ── Zero-frequency (permanent tide) correction ──
    # GGM05C is tide-free: permanent tide removed. For POD we need the
    # "zero-tide" system. The correction is ΔC20_perm ≈ 4.1736e-9
    # for a tide-free→zero-tide conversion (IERS 2010, Eq 6.15).
    # The frequency-independent model above gives the time-varying part
    # relative to the zero-tide mean. So we add the permanent offset back.
    corrections['C20'] = corrections.get('C20', 0.0) + 4.1736e-9

    return corrections


def compute_time_varying_gravity(mjd_tt, ref_mjd=51544.5):
    """Compute secular time-varying corrections to C20/C30/C40.

    Applies linear secular rates from post-glacial rebound and other
    long-term mass redistribution. Reference epoch is J2000.0.

    Args:
        mjd_tt: Current MJD in TT
        ref_mjd: Reference epoch MJD (default J2000.0 = 51544.5)

    Returns:
        dict with keys 'C20','C30','C40' (secular change only)
    """
    dT_years = (mjd_tt - ref_mjd) / 365.25
    return {
        'C20': C20_DOT * dT_years,
        'C30': C30_DOT * dT_years,
        'C40': C40_DOT * dT_years,
    }


def merge_tide_corrections(*corrections):
    """Merge multiple tide correction dicts, summing overlapping keys."""
    merged = {}
    for corr in corrections:
        if corr:
            for key, val in corr.items():
                merged[key] = merged.get(key, 0.0) + val
    return merged
