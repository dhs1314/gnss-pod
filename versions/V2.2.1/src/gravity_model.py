"""Spherical harmonic gravity field model with Clenshaw summation.

Models the Earth's gravitational potential V(r,theta,lambda) at a given position:

    V(r,theta,lambda) = (GM/r) * sum(n=0..Nmax) sum(m=0..n) (R/r)^n *
               P_nm(sin_theta) * [C_nm*cos(m*lambda) + S_nm*sin(m*lambda)]

where:
    P_nm: Fully normalized associated Legendre functions
    C_nm, S_nm: Fully normalized Stokes coefficients
    R: Earth equatorial radius
    GM: Earth gravitational constant

References:
  - Montenbruck & Gill (2000), "Satellite Orbits", Chapter 3
  - Holmes & Featherstone (2002), "A unified approach to the Clenshaw summation..."
  - IERS Conventions (2010), Chapter 6
"""
import numpy as np

GM_EARTH = 3.986004415e14
R_EARTH = 6378136.6
DEG2RAD = np.pi / 180.0


def read_icgem_gfc(filepath):
    """Read gravity field coefficients from ICGEM .gfc format.

    Returns:
        Cnm: (Nmax+1, Nmax+1) array of C coefficients (fully normalized)
        Snm: (Nmax+1, Nmax+1) array of S coefficients (fully normalized)
        Nmax: Maximum degree
        GM: Gravitational constant from file (falls back to GM_EARTH)
        R: Equatorial radius from file (falls back to R_EARTH)
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    Nmax = 0
    for line in lines:
        if line.strip().startswith('max_degree') or line.strip().startswith('degree'):
            try:
                Nmax = int(line.split()[-1])
            except (ValueError, IndexError):
                pass

    if Nmax == 0:
        for line in lines:
            if line.strip().startswith('gfc') or (
                len(line.split()) >= 4 and not line.startswith(('#', '%', '!'))
            ):
                try:
                    parts = line.split()
                    n = int(parts[1])
                    Nmax = max(Nmax, n)
                except (ValueError, IndexError):
                    continue

    Cnm = np.zeros((Nmax + 1, Nmax + 1))
    Snm = np.zeros((Nmax + 1, Nmax + 1))
    GM = GM_EARTH
    R = R_EARTH

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith(('#', '%', '!')):
            continue

        if stripped.startswith('earth_gravity_constant'):
            try:
                GM = float(line.split()[-1])
            except (ValueError, IndexError):
                pass
            continue

        if stripped.startswith('radius'):
            try:
                R = float(line.split()[-1])
            except (ValueError, IndexError):
                pass
            continue

        parts = stripped.split()
        if len(parts) >= 4:
            try:
                keyword = parts[0]
                n = int(parts[1])
                m = int(parts[2])
                # Handle Fortran D-format (e.g. 1.0D+00 -> 1.0E+00)
                C_val = float(parts[3].replace('D', 'E'))
                S_val = float(parts[4].replace('D', 'E')) if len(parts) >= 5 else 0.0

                if keyword == 'gfc' and n <= Nmax and m <= n:
                    Cnm[n, m] = C_val
                    Snm[n, m] = S_val
            except (ValueError, IndexError):
                continue

    return Cnm, Snm, Nmax, GM, R


def compute_gravity_acceleration(pos_ecef, Cnm, Snm, Nmax=120, tide_corrections=None,
                                 GM=None, R=None):
    """Compute gravitational acceleration using Clenshaw summation.

    Args:
        pos_ecef: Position vector in ECEF [m], shape (3,)
        Cnm, Snm: Gravity field coefficients (fully normalized)
        Nmax: Maximum degree to use
        tide_corrections: Optional dict of delta_C20, delta_C21, delta_S21, etc.
        GM: Gravitational constant (defaults to GM_EARTH)
        R: Equatorial radius (defaults to R_EARTH)

    Returns:
        acc_ecef: Acceleration vector in ECEF [m/s^2], shape (3,)
    """
    _GM = GM if GM is not None else GM_EARTH
    _R = R if R is not None else R_EARTH

    x, y, z = pos_ecef[0], pos_ecef[1], pos_ecef[2]
    r = np.sqrt(x**2 + y**2 + z**2)

    if r < _R * 0.5:
        return np.zeros(3)

    lat = np.arcsin(z / r)
    lon = np.arctan2(y, x)

    # Apply tide corrections to Cnm/Snm if provided
    Cnm_t = Cnm
    Snm_t = Snm
    if tide_corrections:
        Cnm_t = Cnm.copy()
        Snm_t = Snm.copy()
        for key, val in tide_corrections.items():
            if key == 'C20' and Nmax >= 2:
                Cnm_t[2, 0] += val
            elif key == 'C21' and Nmax >= 2:
                Cnm_t[2, 1] += val
            elif key == 'S21' and Nmax >= 2:
                Snm_t[2, 1] += val
            elif key == 'C22' and Nmax >= 2:
                Cnm_t[2, 2] += val
            elif key == 'S22' and Nmax >= 2:
                Snm_t[2, 2] += val

    sin_lat = np.sin(lat)
    cos_lat = np.cos(lat)

    # Pre-compute (R/r)^n factors
    ratio = _R / r
    r_pow = np.ones(Nmax + 2)
    for n in range(1, Nmax + 2):
        r_pow[n] = r_pow[n - 1] * ratio

    # Pre-compute cos(m*lon) and sin(m*lon) for all m
    cos_m_lon = np.ones(Nmax + 1)
    sin_m_lon = np.zeros(Nmax + 1)
    for m in range(1, Nmax + 1):
        cos_m_lon[m] = cos_m_lon[m - 1] * np.cos(lon) - sin_m_lon[m - 1] * np.sin(lon)
        sin_m_lon[m] = sin_m_lon[m - 1] * np.cos(lon) + cos_m_lon[m - 1] * np.sin(lon)

    # Legendre polynomials via recurrence
    Pnm = np.zeros((Nmax + 2, Nmax + 2))
    dPnm = np.zeros((Nmax + 2, Nmax + 2))

    Pnm[0, 0] = 1.0
    if Nmax >= 1:
        Pnm[1, 0] = np.sqrt(3.0) * sin_lat
        Pnm[1, 1] = np.sqrt(3.0) * cos_lat
        dPnm[1, 0] = np.sqrt(3.0) * cos_lat
        dPnm[1, 1] = -np.sqrt(3.0) * sin_lat

    for n in range(2, Nmax + 1):
        # Sectoral (m = n)
        Pnm[n, n] = np.sqrt((2.0 * n + 1.0) / (2.0 * n)) * cos_lat * Pnm[n - 1, n - 1]

        # Near-sectoral (m = n-1)
        Pnm[n, n - 1] = np.sqrt(2.0 * n + 1.0) * sin_lat * Pnm[n - 1, n - 1]

        # Full recurrence
        for m in range(min(n - 2, Nmax) + 1):
            if m <= n - 2:
                a = np.sqrt(((2.0 * n - 1.0) * (2.0 * n + 1.0)) /
                           ((n - m) * (n + m)))
                b = np.sqrt(((2.0 * n + 1.0) * (n + m - 1.0) * (n - m - 1.0)) /
                           ((n - m) * (n + m) * (2.0 * n - 3.0)))
                Pnm[n, m] = a * sin_lat * Pnm[n - 1, m] - b * Pnm[n - 2, m]

    # Compute dPnm/d(lat)
    for n in range(1, Nmax + 1):
        for m in range(min(n, Nmax) + 1):
            if m == n:
                dPnm[n, m] = -n * sin_lat / max(cos_lat, 1e-12) * Pnm[n, m]
            elif m <= n:
                f = np.sqrt(float(n**2 - m**2) * (2.0 * n + 1.0) / (2.0 * n - 1.0))
                dPnm[n, m] = (-n * sin_lat * Pnm[n, m] + f * Pnm[n - 1, m]) / max(cos_lat, 1e-12)

    # Accumulate potential derivatives
    GM_over_r = _GM / r

    dU_dr = 1.0  # central term (n=0, m=0): point-mass contribution
    dU_dlat = 0.0
    dU_dlon = 0.0

    for n in range(1, Nmax + 1):
        rn_factor = r_pow[n]
        for m in range(min(n, Nmax) + 1):
            C = Cnm_t[n, m]
            S = Snm_t[n, m]

            if abs(C) < 1e-30 and abs(S) < 1e-30:
                continue

            term_cos = cos_m_lon[m] * C + sin_m_lon[m] * S
            term_sin = sin_m_lon[m] * C - cos_m_lon[m] * S

            dU_dr += (n + 1.0) * rn_factor * Pnm[n, m] * term_cos
            dU_dlat += rn_factor * dPnm[n, m] * term_cos
            dU_dlon += rn_factor * Pnm[n, m] * m * term_sin

    # Scale by GM/r^2
    dU_dr *= -_GM / r**2
    dU_dlat *= GM_over_r
    dU_dlon *= GM_over_r

    # Rotate from local (radial, lat, lon) frame to ECEF
    cos_lon = np.cos(lon)
    sin_lon = np.sin(lon)

    dU_dlat_a = dU_dlat / r
    dU_dlon_a = dU_dlon / (r * cos_lat) if abs(cos_lat) > 1e-12 else 0.0

    ax = cos_lat * cos_lon * dU_dr - sin_lat * cos_lon * dU_dlat_a - sin_lon * dU_dlon_a
    ay = cos_lat * sin_lon * dU_dr - sin_lat * sin_lon * dU_dlat_a + cos_lon * dU_dlon_a
    az = sin_lat * dU_dr + cos_lat * dU_dlat_a

    return np.array([ax, ay, az])
