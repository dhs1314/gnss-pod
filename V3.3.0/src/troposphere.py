"""Troposphere delay model: Saastamoinen ZHD + GMF mapping.

Replaces the trivial compute_tropo(el) with a physically motivated model.
ZHD is computed from the Saastamoinen (1972) formula. The mapping function
uses a simplified GMF (Boehm et al., 2006) with latitude-band tabulated
coefficients and annual cosine variation.
"""
import numpy as np

# GMF mean coefficients by latitude band (15-degree steps, centre of band).
# ah: hydrostatic a coefficient (continued-fraction form)
# bh, ch: hydrostatic b, c coefficients
# aw, bw, cw: wet coefficients
# Annual amplitude for ah and aw (cosine of day-of-year term).
# Values derived from Boehm et al. (2006) GMF spherical harmonic evaluation
# at band centres with zero longitude.

_GMF_TABLE = {
    # lat_deg: (ah_mean, bh, ch, aw_mean, bw, cw, ah_amp, aw_amp)
    90: (1.277e-3, 2.915e-3, 6.246e-2, 5.762e-4, 1.491e-3, 5.500e-2, 0.0, 0.0),
    75: (1.270e-3, 2.915e-3, 6.246e-2, 5.768e-4, 1.495e-3, 5.510e-2, 5.0e-6, 2.0e-6),
    60: (1.243e-3, 2.917e-3, 6.248e-2, 5.789e-4, 1.508e-3, 5.540e-2, 1.2e-5, 5.0e-6),
    45: (1.201e-3, 2.920e-3, 6.252e-2, 5.825e-4, 1.527e-3, 5.586e-2, 1.8e-5, 8.0e-6),
    30: (1.153e-3, 2.924e-3, 6.257e-2, 5.875e-4, 1.551e-3, 5.648e-2, 2.2e-5, 1.0e-5),
    15: (1.113e-3, 2.927e-3, 6.262e-2, 5.928e-4, 1.572e-3, 5.716e-2, 2.4e-5, 1.2e-5),
    0: (1.097e-3, 2.928e-3, 6.264e-2, 5.957e-4, 1.582e-3, 5.754e-2, 2.5e-5, 1.3e-5),
    -15: (1.113e-3, 2.927e-3, 6.262e-2, 5.928e-4, 1.572e-3, 5.716e-2, 2.4e-5, 1.2e-5),
    -30: (1.153e-3, 2.924e-3, 6.257e-2, 5.875e-4, 1.551e-3, 5.648e-2, 2.2e-5, 1.0e-5),
    -45: (1.201e-3, 2.920e-3, 6.252e-2, 5.825e-4, 1.527e-3, 5.586e-2, 1.8e-5, 8.0e-6),
    -60: (1.243e-3, 2.917e-3, 6.248e-2, 5.789e-4, 1.508e-3, 5.540e-2, 1.2e-5, 5.0e-6),
    -75: (1.270e-3, 2.915e-3, 6.246e-2, 5.768e-4, 1.495e-3, 5.510e-2, 5.0e-6, 2.0e-6),
    -90: (1.277e-3, 2.915e-3, 6.246e-2, 5.762e-4, 1.491e-3, 5.500e-2, 0.0, 0.0),
}

_LAT_BANDS = np.array(sorted(_GMF_TABLE.keys()))


def _interp_gmf_coeffs(lat_deg, doy):
    """Interpolate GMF coefficients to given latitude and day-of-year."""
    # Find bracketing latitude bands
    lat_deg = float(np.clip(lat_deg, -90.0, 90.0))
    idx = int(np.searchsorted(_LAT_BANDS, lat_deg))
    if idx == 0:
        idx = 1
    if idx >= len(_LAT_BANDS):
        idx = len(_LAT_BANDS) - 1

    lat_lo = _LAT_BANDS[idx - 1]
    lat_hi = _LAT_BANDS[idx]
    frac = (lat_deg - lat_lo) / (lat_hi - lat_lo) if lat_hi != lat_lo else 0.0

    row_lo = _GMF_TABLE[lat_lo]
    row_hi = _GMF_TABLE[lat_hi]

    # Annual variation: cos(2*pi*(doy - 28)/365.25), peak at doy=28 (late January)
    annual = np.cos(2.0 * np.pi * (doy - 28.0) / 365.25)

    ah_lo = row_lo[0] + row_lo[6] * annual
    ah_hi = row_hi[0] + row_hi[6] * annual
    aw_lo = row_lo[3] + row_lo[7] * annual
    aw_hi = row_hi[3] + row_hi[7] * annual

    ah = ah_lo + frac * (ah_hi - ah_lo)
    aw = aw_lo + frac * (aw_hi - aw_lo)
    # b and c coefficients vary slowly with latitude; interpolate directly
    bh = row_lo[1] + frac * (row_hi[1] - row_lo[1])
    ch = row_lo[2] + frac * (row_hi[2] - row_lo[2])
    bw = row_lo[4] + frac * (row_hi[4] - row_lo[4])
    cw = row_lo[5] + frac * (row_hi[5] - row_lo[5])

    return ah, bh, ch, aw, bw, cw


def _mf_continued_fraction(el_rad, a, b, c):
    """Evaluate Niell continued-fraction mapping function."""
    sin_el = np.sin(el_rad)
    denom = sin_el + a / (sin_el + b / (sin_el + c))
    return 1.0 / denom


def saastamoinen_zhd(lat_rad, h_m):
    """Compute Saastamoinen zenith hydrostatic delay [m].

    For spaceborne receivers above the neutral atmosphere (h > 80 km),
    ZHD is zero — the hydrostatic delay comes from the ray path through
    the lower atmosphere, not from the receiver zenith.

    Args:
        lat_rad: geodetic latitude [rad]
        h_m: ellipsoidal height [m]

    Returns:
        ZHD in metres.
    """
    if h_m > 80000.0:
        return 0.0
    # Surface pressure from standard exponential model
    # Valid only for h < ~10 km; capped for stratospheric heights
    h_clip = min(h_m, 10000.0)
    P = 1013.25 * (1.0 - 2.25577e-5 * h_clip) ** 5.25588
    # Saastamoinen (1972) / Davis et al. (1985)
    cos2lat = np.cos(2.0 * lat_rad)
    zhd = 0.0022768 * P / (1.0 - 0.00266 * cos2lat - 2.8e-7 * h_clip)
    return float(zhd)


def gmf_mapping(el_rad, lat_rad, lon_rad, h_m, doy):
    """Compute GMF hydrostatic and wet mapping functions.

    Args:
        el_rad: elevation angle [rad]
        lat_rad: geodetic latitude [rad]
        lon_rad: geodetic longitude [rad] (reserved for future enhancement)
        h_m: ellipsoidal height [m]
        doy: day of year (1-366)

    Returns:
        (mf_h, mf_w): hydrostatic and wet mapping function values.
    """
    lat_deg = float(np.rad2deg(lat_rad))
    ah, bh, ch, aw, bw, cw = _interp_gmf_coeffs(lat_deg, int(doy))
    mf_h = _mf_continued_fraction(el_rad, ah, bh, ch)
    mf_w = _mf_continued_fraction(el_rad, aw, bw, cw)
    return float(mf_h), float(mf_w)


def compute_troposphere(el_rad, lat_rad, lon_rad, h_m, doy, zwd=0.1):
    """Compute total tropospheric delay.

    Args:
        el_rad: elevation angle [rad]
        lat_rad, lon_rad: receiver position [rad]
        h_m: ellipsoidal height [m]
        doy: day of year (1-366)
        zwd: zenith wet delay [m] (estimated by filter, default 0.1)

    Returns:
        (total_tropo, zhd, mf_h, mf_w)
    """
    zhd = saastamoinen_zhd(lat_rad, h_m)
    mf_h, mf_w = gmf_mapping(el_rad, lat_rad, lon_rad, h_m, doy)
    total = zhd * mf_h + zwd * mf_w
    return total, zhd, mf_h, mf_w


def ecef_to_geodetic(pos_ecef):
    """Convert ECEF position to geodetic lat/lon/height.

    Args:
        pos_ecef: (3,) ECEF position [m]

    Returns:
        (lat_rad, lon_rad, h_m)
    """
    x, y, z = float(pos_ecef[0]), float(pos_ecef[1]), float(pos_ecef[2])
    # WGS84
    a = 6378137.0
    f = 1.0 / 298.257223563
    e2 = 2.0 * f - f * f

    lon = np.arctan2(y, x)
    p = np.sqrt(x * x + y * y)
    # Iterative solution for latitude
    lat = np.arctan2(z, p * (1.0 - e2))
    for _ in range(5):
        N = a / np.sqrt(1.0 - e2 * np.sin(lat) ** 2)
        h = p / np.cos(lat) - N
        lat = np.arctan2(z, p * (1.0 - e2 * N / (N + h)))

    N = a / np.sqrt(1.0 - e2 * np.sin(lat) ** 2)
    h = p / np.cos(lat) - N
    return float(lat), float(lon), float(h)
