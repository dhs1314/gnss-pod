"""Solar radiation pressure (SRP) acceleration.

For a cannonball model (spherical satellite approximation):
    a_srp = CR * (A/m) * (S0/c) * (1 AU / r_sun)^2 * e_sun

where:
    CR = radiation pressure coefficient (~1.0-1.4)
    A = cross-sectional area [m^2]
    m = satellite mass [kg]
    S0 = solar constant = 1361 W/m^2
    c = speed of light
    AU = astronomical unit
    r_sun = distance from satellite to Sun
    e_sun = unit vector from satellite to Sun

Includes conical Earth shadow model for eclipse handling.

For GRACE-FO at ~500 km, SRP acceleration is ~1-5 x 10^-7 m/s^2.
"""
import numpy as np
from astropy.time import Time
from astropy.coordinates import get_sun
from astropy import units as u

C_LIGHT = 299792458.0
AU = 149597870700.0
R_EARTH = 6378136.6
SOLAR_CONSTANT = 1361.0  # W/m^2 at 1 AU


def get_sun_position_eci(mjd_tt):
    """Get Sun position in GCRS (ECI) at given time.

    Args:
        mjd_tt: MJD in Terrestrial Time

    Returns:
        pos_sun_eci: (3,) Sun position in GCRS [m]
    """
    t = Time(mjd_tt, format='mjd', scale='tt')
    sun = get_sun(t)
    pos = sun.cartesian.xyz.to(u.m).value
    return np.array([pos[0], pos[1], pos[2]])


def conical_earth_shadow(sat_pos_eci, sun_pos_eci):
    """Compute shadow factor using conical Earth shadow model.

    Shadow factor nu:
        nu = 0  -> full shadow (umbra)
        nu = 1  -> full sunlight
        0 < nu < 1 -> partial shadow (penumbra)

    Args:
        sat_pos_eci: Satellite position in ECI [m], shape (3,)
        sun_pos_eci: Sun position in ECI [m], shape (3,)

    Returns:
        nu: Shadow factor in [0, 1]
    """
    r_sat = np.linalg.norm(sat_pos_eci)
    r_sun = np.linalg.norm(sun_pos_eci)

    s_hat = -sun_pos_eci / r_sun

    R_SUN = 696340e3

    a_sun = np.arcsin(np.clip(R_SUN / r_sun, -1.0, 1.0))
    cos_psi = np.dot(-sat_pos_eci, s_hat) / r_sat
    a_earth = np.arcsin(np.clip(R_EARTH / r_sat, -1.0, 1.0))

    if cos_psi >= 1.0:
        psi = 0.0
    elif cos_psi <= -1.0:
        psi = np.pi
    else:
        psi = np.arccos(cos_psi)

    y = a_sun + a_earth  # penumbra edge

    if psi > y:
        return 1.0

    x = abs(a_earth - a_sun)  # umbra edge

    if psi < x and a_earth > a_sun:
        return 0.0

    if a_earth > a_sun:
        A = (psi - x) / (y - x)
    else:
        A = 1.0 - (y - psi) / (y - x)

    return float(np.clip(A, 0.0, 1.0))


def compute_srp_acceleration(sat_pos_eci, mjd_tt, CR=1.3, area=3.4, mass=580.0):
    """Compute solar radiation pressure acceleration (cannonball model).

    Args:
        sat_pos_eci: Satellite position in ECI [m], shape (3,)
        mjd_tt: MJD in TT
        CR: Solar radiation pressure coefficient
        area: Cross-sectional area [m^2]
        mass: Satellite mass [kg]

    Returns:
        acc_srp_eci: SRP acceleration in ECI [m/s^2], shape (3,)
    """
    sun_pos_eci = get_sun_position_eci(mjd_tt)

    d = sun_pos_eci - sat_pos_eci
    d_norm = np.linalg.norm(d)
    sun_hat = d / d_norm

    flux_ratio = (AU / d_norm) ** 2
    P0 = SOLAR_CONSTANT / C_LIGHT

    nu = conical_earth_shadow(sat_pos_eci, sun_pos_eci)

    acc_mag = nu * CR * (area / mass) * P0 * flux_ratio
    return acc_mag * sun_hat
