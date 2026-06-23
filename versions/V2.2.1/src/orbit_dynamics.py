"""
Force models for LEO orbit propagation.

Phase 2 (ECEF): two-body, J2, exponential drag, Coriolis + centrifugal
Phase 3 (ECI):  full spherical harmonic gravity (GGM05C), third-body
                perturbations (Sun/Moon/planets), cannonball SRP with
                conical shadow, empirical RTN accelerations.

The ECI formulation avoids fictitious Coriolis/centrifugal forces.
Gravity and drag are computed in ECEF then rotated to ECI; third-body
and SRP are computed natively in ECI.

Reference: Montenbruck & Gill, "Satellite Orbits", Springer 2000.
"""
import numpy as np

# Earth constants
GM = 3.986004415e14
J2 = 1.082626683e-3
R_EQUATOR = 6378137.0
F_EARTH = 1.0 / 298.257223563
OMEGA_E = 7.2921151467e-5

# Drag: exponential atmosphere (GRACE-FO altitude ~490 km)
RHO0 = 1.0e-13
H0_REF = 500000.0
H_SCALE = 60000.0


# ---------------------------------------------------------------------------
# Phase 2 ECEF functions (kept for backward compatibility)
# ---------------------------------------------------------------------------

def earth_radius_lat(lat_geodetic):
    sin_lat = np.sin(lat_geodetic)
    return R_EQUATOR / np.sqrt(1.0 - (2.0 * F_EARTH - F_EARTH**2) * sin_lat**2)


def altitude_from_r(r_ecef):
    r_mag = np.linalg.norm(r_ecef)
    lat = np.arcsin(r_ecef[2] / r_mag)
    return r_mag - earth_radius_lat(lat)


def two_body_acc(r_ecef):
    r = np.linalg.norm(r_ecef)
    if r < 1.0:
        return np.zeros(3)
    return -GM * r_ecef / r**3


def j2_acc(r_ecef):
    x, y, z = r_ecef
    r = np.linalg.norm(r_ecef)
    if r < 1.0:
        return np.zeros(3)
    r2 = r * r
    z2_r2 = (z / r)**2
    coeff = -1.5 * J2 * GM * R_EQUATOR**2 / (r2 * r2 * r)
    a_x = coeff * x * (1.0 - 5.0 * z2_r2)
    a_y = coeff * y * (1.0 - 5.0 * z2_r2)
    a_z = coeff * z * (3.0 - 5.0 * z2_r2)
    return np.array([a_x, a_y, a_z])


def drag_acc(r_ecef, v_ecef, Cd, area_to_mass=0.002):
    alt = altitude_from_r(r_ecef)
    if alt > 2000000.0:
        return np.zeros(3)
    rho = RHO0 * np.exp(-(alt - H0_REF) / H_SCALE)
    omega_vec = np.array([0.0, 0.0, OMEGA_E])
    v_rel = v_ecef - np.cross(omega_vec, r_ecef)
    v_rel_mag = np.linalg.norm(v_rel)
    if v_rel_mag < 1e-9:
        return np.zeros(3)
    return -0.5 * Cd * area_to_mass * rho * v_rel_mag * v_rel


def total_acc(r_ecef, v_ecef, Cd=2.2, area_to_mass=0.002):
    """Total acceleration in ECEF frame (m/s^2). Phase 2.

    Includes two-body, J2, drag, Coriolis, and centrifugal.
    """
    a = two_body_acc(r_ecef)
    a += j2_acc(r_ecef)
    a += drag_acc(r_ecef, v_ecef, Cd, area_to_mass)
    a[0] += 2.0 * OMEGA_E * v_ecef[1]
    a[1] += -2.0 * OMEGA_E * v_ecef[0]
    a[0] += OMEGA_E * OMEGA_E * r_ecef[0]
    a[1] += OMEGA_E * OMEGA_E * r_ecef[1]
    return a


# ---------------------------------------------------------------------------
# Phase 3 ECI force model
# ---------------------------------------------------------------------------

def total_acc_eci(pos_eci, vel_eci, mjd_tt, mjd_utc,
                  Cnm, Snm, Nmax,
                  CD=2.2, CR=1.3,
                  area_drag=0.68, area_srp=3.4, mass=580.0,
                  empirical_acc_rtn=None,
                  tide_corrections=None,
                  bodies=None,
                  GM_gravity=None, R_gravity=None):
    """Total acceleration in ECI frame (m/s^2). Phase 3.

    Combines all force models:
      1. Earth gravity (spherical harmonic, computed in ECEF)
      2. Third-body perturbations (Sun, Moon, planets)
      3. Atmospheric drag (computed in ECEF)
      4. Solar radiation pressure (cannonball + conical shadow)
      5. Empirical RTN accelerations

    Args:
        pos_eci: ECI position [m], shape (3,)
        vel_eci: ECI velocity [m/s], shape (3,)
        mjd_tt: MJD in Terrestrial Time (for third-body/SRP)
        mjd_utc: MJD in UTC (for ECEF transforms)
        Cnm, Snm: Gravity field coefficients
        Nmax: Maximum gravity degree
        CD: Drag coefficient
        CR: SRP coefficient
        area_drag: Drag cross-section [m^2]
        area_srp: SRP cross-section [m^2]
        mass: Satellite mass [kg]
        empirical_acc_rtn: (3,) empirical RTN acceleration [m/s^2] or None
        tide_corrections: Optional dict of delta_C20, etc.
        bodies: List of third body names, default ['Sun','Moon','Jupiter','Venus']
        GM_gravity: Gravitational constant for gravity model
        R_gravity: Equatorial radius for gravity model

    Returns:
        total_acc: Total acceleration in ECI [m/s^2], shape (3,)
    """
    from src.coordinates import eci_to_ecef, ecef_to_eci
    from src.gravity_model import compute_gravity_acceleration
    from src.third_body import compute_total_third_body

    # Convert to ECEF for gravity and drag computation
    pos_ecef, vel_ecef = eci_to_ecef(pos_eci, vel_eci, mjd_utc)

    # 1. Earth gravity in ECEF, rotate to ECI
    acc_grav_ecef = compute_gravity_acceleration(
        pos_ecef, Cnm, Snm, Nmax, tide_corrections,
        GM=GM_gravity, R=R_gravity,
    )
    acc_grav_eci, _ = ecef_to_eci(acc_grav_ecef, np.zeros(3), mjd_utc)

    # 2. Third-body perturbations (computed in ECI)
    acc_third = compute_total_third_body(pos_eci, mjd_tt, bodies)

    # 3. Atmospheric drag in ECEF, rotate to ECI
    acc_drag_ecef = drag_acc(pos_ecef, vel_ecef, CD, area_drag / mass)
    acc_drag_eci, _ = ecef_to_eci(acc_drag_ecef, np.zeros(3), mjd_utc)

    # 4. Solar radiation pressure in ECI
    from src.srp import compute_srp_acceleration
    acc_srp = compute_srp_acceleration(pos_eci, mjd_tt, CR, area_srp, mass)

    # 5. Empirical RTN accelerations
    total = acc_grav_eci + acc_third + acc_drag_eci + acc_srp

    if empirical_acc_rtn is not None:
        from src.empirical import rtn_to_eci
        acc_emp = rtn_to_eci(empirical_acc_rtn, pos_eci, vel_eci)
        total += acc_emp

    return total
