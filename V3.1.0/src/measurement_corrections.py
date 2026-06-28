"""Measurement model corrections for precision orbit determination.

Phase wind-up (Wu et al. 1993), relativistic Shapiro delay, and satellite
antenna phase center offset (PCO) corrections.
"""
import numpy as np

C_LIGHT = 299792458.0
MU_EARTH = 3.986004415e14


# ---------------------------------------------------------------------------
# Phase Wind-Up (Wu et al., 1993)
# ---------------------------------------------------------------------------

def phase_wind_up_correction(sat_pos_eci, leo_pos_eci, sat_vel_eci=None,
                             prev_dphi=None):
    """Compute phase wind-up correction for a GPS-LEO link [radians].

    The GPS transmit antenna rotates to keep its +z axis pointing to Earth
    centre. The LEO receive antenna orientation changes with spacecraft
    attitude. This rotation introduces a fractional-cycle phase error that
    accumulates over time.

    Simplified formulation: assumes GPS antenna +z points toward geocentre
    (valid for Block II/III satellites under nominal yaw-steering), and
    LEO antenna +z points toward local zenith.

    Args:
        sat_pos_eci: GPS satellite ECI position (3,) [m]
        leo_pos_eci: LEO receiver ECI position (3,) [m]
        sat_vel_eci: GPS satellite ECI velocity (3,) [m] (optional, for
                     computing satellite body axes)
        prev_dphi: accumulated wind-up from previous epoch [radians]

    Returns:
        delta_phi: phase correction to ADD to modeled phase [radians]
        current_dphi: accumulated wind-up angle for next epoch [radians]
    """
    # Line-of-sight unit vector (receiver → transmitter)
    los = sat_pos_eci - leo_pos_eci
    rho = float(np.linalg.norm(los))
    k = los / rho  # unit vector from LEO to GPS

    # GPS satellite effective dipole
    # Nominal yaw-steering: +z body axis points to Earth centre,
    # +y body axis = cross(sun_dir, z_body) normalized
    r_sat = float(np.linalg.norm(sat_pos_eci))
    z_gps = -sat_pos_eci / r_sat  # toward Earth centre

    # Sun direction (approximate — good to ~0.5° which is sufficient
    # for the dipole cross-product geometry)
    sun_dir = sat_pos_eci - leo_pos_eci  # reuse as rough sun direction
    # Better: use the actual satellite velocity for y-body
    if sat_vel_eci is not None:
        y_gps_n = np.cross(z_gps, sat_vel_eci)
        y_gps_norm = float(np.linalg.norm(y_gps_n))
        if y_gps_norm > 1e-12:
            y_gps = y_gps_n / y_gps_norm
        else:
            y_gps = np.array([0.0, 1.0, 0.0])
    else:
        y_gps = np.array([0.0, 1.0, 0.0])
        y_gps = y_gps - np.dot(y_gps, z_gps) * z_gps
        y_gps = y_gps / float(np.linalg.norm(y_gps))

    x_gps = np.cross(y_gps, z_gps)

    # GPS effective dipole
    D_gps = x_gps - k * np.dot(k, x_gps) - np.cross(k, y_gps)

    # LEO effective dipole (zenith-pointing, +z toward local zenith)
    r_leo = float(np.linalg.norm(leo_pos_eci))
    z_leo = leo_pos_eci / r_leo  # toward local zenith
    # LEO y-axis: cross(z_leo, velocity) — use GPS approximation
    y_leo = np.array([0.0, 1.0, 0.0])
    y_leo = y_leo - np.dot(y_leo, z_leo) * z_leo
    y_leo_norm = float(np.linalg.norm(y_leo))
    if y_leo_norm > 1e-12:
        y_leo = y_leo / y_leo_norm
    x_leo = np.cross(y_leo, z_leo)

    # LEO effective dipole
    D_leo = x_leo - k * np.dot(k, x_leo) + np.cross(k, y_leo)

    # Wind-up angle
    cos_phi = np.dot(D_leo, D_gps)
    sin_phi = np.dot(np.cross(k, D_leo), D_gps)
    phi = float(np.arctan2(sin_phi, cos_phi))

    # Accumulate fractional turns
    if prev_dphi is None:
        current_dphi = phi
        delta_phi = 0.0
    else:
        delta = phi - prev_dphi
        # Bring delta into [-pi, pi]
        delta = (delta + np.pi) % (2.0 * np.pi) - np.pi
        current_dphi = prev_dphi + delta
        delta_phi = delta

    return delta_phi, current_dphi


# ---------------------------------------------------------------------------
# Relativistic Shapiro Delay
# ---------------------------------------------------------------------------

def relativity_shapiro_correction(sat_pos_eci, leo_pos_eci):
    """Compute relativistic Shapiro signal delay for GPS-LEO link [m].

    General relativistic correction due to Earth's gravitational field
    bending the signal path. Magnitude ~1-2 cm for GPS→LEO.

    Args:
        sat_pos_eci: GPS satellite ECI position (3,) [m]
        leo_pos_eci: LEO receiver ECI position (3,) [m]

    Returns:
        Range correction to ADD to modeled observation [m]
    """
    R_sat = float(np.linalg.norm(sat_pos_eci))
    R_leo = float(np.linalg.norm(leo_pos_eci))
    rho = float(np.linalg.norm(sat_pos_eci - leo_pos_eci))

    # Shapiro formula
    num = R_sat + R_leo + rho
    den = R_sat + R_leo - rho
    if den < 1e-9:
        return 0.0

    dt_rel = 2.0 * MU_EARTH / C_LIGHT**3 * np.log(num / den)
    return float(C_LIGHT * dt_rel)


# ---------------------------------------------------------------------------
# Satellite Antenna Phase Centre Offset (PCO)
# ---------------------------------------------------------------------------

def apply_satellite_pco(sat_pos_ecef, leo_pos_ecef, pco_ecef):
    """Apply GPS satellite antenna PCO to satellite position.

    The phase centre offset shifts the effective measurement point from
    the satellite centre of mass (SP3 orbit) to the antenna phase centre.

    Args:
        sat_pos_ecef: satellite CoM ECEF position (3,) [m]
        leo_pos_ecef: receiver ECEF position (3,) [m]
        pco_ecef: PCO vector in ECEF (3,) [m] — typically ~[0, 0, 0.10] m
                  for GPS Block IIR/IIR-M (Z offset toward Earth)

    Returns:
        corrected_sat_pos: satellite antenna phase centre ECEF position (3,) [m]
    """
    return sat_pos_ecef + np.asarray(pco_ecef)


def compute_pco_ecef_from_nadir(sat_pos_ecef, pco_z):
    """Compute PCO vector in ECEF for a nadir-pointing GPS antenna.

    The GPS satellite antenna's +Z body axis nominally points toward the
    Earth centre under yaw-steering. The PCO Z component (~0.1 m) is
    along this direction away from the satellite CoM.

    Args:
        sat_pos_ecef: satellite CoM ECEF position (3,) [m]
        pco_z: PCO Z component [m] (positive = away from Earth centre,
               toward antenna boresight)

    Returns:
        pco_ecef: PCO vector in ECEF (3,) [m]
    """
    r_sat = float(np.linalg.norm(sat_pos_ecef))
    if r_sat < 1e-6:
        return np.zeros(3)
    nadir = -sat_pos_ecef / r_sat  # toward Earth centre
    # PCO positive Z = away from Earth = -nadir
    return pco_z * (-nadir)
