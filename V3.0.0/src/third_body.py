"""Third-body gravitational perturbations from Sun, Moon, and planets.

The perturbation acceleration on the LEO satellite from a third body j
with mass Mj is:

    a_j = GMj * [(r_j - r)/(|r_j - r|^3) - r_j/|r_j|^3]

where:
    r_j = position of body j relative to Earth center
    r = position of the satellite relative to Earth center

The first term is the direct attraction from body j.
The second term accounts for the fact that the Earth is also
accelerated by body j (indirect effect).
"""
import numpy as np
from astropy.time import Time
from astropy.coordinates import get_body, get_sun
from astropy import units as u

C_LIGHT = 299792458.0

BODY_GM = {
    'Sun': 1.32712440018e20,
    'Moon': 4.90280031e12,
    'Jupiter': 1.26686534e17,
    'Venus': 3.24858599e14,
    'Mars': 4.282837e13,
    'Saturn': 3.7931187e16,
}


def get_body_position_eci(body_name, mjd_tt):
    """Get position of a celestial body in GCRS (ECI) using astropy/DE440.

    Args:
        body_name: Name of the body ('Sun', 'Moon', 'Jupiter', etc.)
        mjd_tt: MJD in Terrestrial Time

    Returns:
        pos_eci: (3,) position in GCRS [m]
    """
    t = Time(mjd_tt, format='mjd', scale='tt')

    if body_name.lower() == 'sun':
        body = get_sun(t)
    else:
        body = get_body(body_name.lower(), t)

    pos = body.cartesian.xyz.to(u.m).value
    return np.array([pos[0], pos[1], pos[2]])


def compute_third_body_acceleration(sat_pos_eci, body_name, mjd_tt):
    """Compute third-body perturbation acceleration.

    Args:
        sat_pos_eci: Satellite position in GCRS [m], shape (3,)
        body_name: Name of the perturbing body
        mjd_tt: MJD in TT

    Returns:
        acc_eci: Acceleration in GCRS [m/s^2], shape (3,)
    """
    gm = BODY_GM.get(body_name)
    if gm is None:
        return np.zeros(3)

    r_body = get_body_position_eci(body_name, mjd_tt)
    r_sat = sat_pos_eci

    d = r_body - r_sat
    d_norm = np.linalg.norm(d)
    r_body_norm = np.linalg.norm(r_body)

    if d_norm < 1e-10 or r_body_norm < 1e-10:
        return np.zeros(3)

    # Direct acceleration from body on satellite
    acc_direct = gm * d / d_norm**3

    # Indirect acceleration on Earth from body (acts oppositely on satellite)
    acc_indirect = gm * r_body / r_body_norm**3

    return acc_direct - acc_indirect


def compute_total_third_body(sat_pos_eci, mjd_tt, bodies=None):
    """Total third-body acceleration from multiple bodies.

    Args:
        sat_pos_eci: Satellite position in ECI [m], shape (3,)
        mjd_tt: MJD in TT
        bodies: List of body names. Default: ['Sun', 'Moon', 'Jupiter', 'Venus']

    Returns:
        total_acc_eci: Total third-body acceleration in ECI [m/s^2], shape (3,)
    """
    if bodies is None:
        bodies = ['Sun', 'Moon', 'Jupiter', 'Venus']

    total_acc = np.zeros(3)
    for body in bodies:
        total_acc += compute_third_body_acceleration(sat_pos_eci, body, mjd_tt)

    return total_acc
