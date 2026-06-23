"""Coordinate system transformations: ECEF (ITRF) <-> ECI (GCRS/J2000).

Uses astropy's ERFA-based transformations implementing the
IAU 2006/2000A precession-nutation model and IERS polar motion,
including the full velocity transformation.

All ``mjd`` parameters are in UTC timescale.
"""
import numpy as np
from astropy.time import Time
from astropy.coordinates import (
    ITRS, GCRS, CartesianRepresentation, CartesianDifferential,
)
from astropy import units as u
from astropy.utils import iers

iers.conf.auto_download = False

SEC_PER_DAY = 86400.0


def ecef_to_eci(pos_ecef, vel_ecef, mjd_utc):
    """Convert position and velocity from ECEF (ITRF) to ECI (GCRS).

    Args:
        pos_ecef: Position vector(s) in ITRF [m], shape (3,) or (N,3)
        vel_ecef: Velocity vector(s) in ITRF [m/s], shape (3,) or (N,3)
        mjd_utc: Modified Julian Date (UTC)

    Returns:
        pos_eci: Position in GCRS [m]
        vel_eci: Velocity in GCRS [m/s]
    """
    obs_time = Time(mjd_utc, format='mjd', scale='utc')

    rep = CartesianRepresentation(
        pos_ecef * u.m,
        differentials={'s': CartesianDifferential(vel_ecef * u.m / u.s)}
    )
    itrs = ITRS(rep, obstime=obs_time)
    gcrs = itrs.transform_to(GCRS(obstime=obs_time))

    pos_eci = gcrs.cartesian.xyz.to(u.m).value.T
    vel_eci = gcrs.velocity.d_xyz.to(u.m / u.s).value.T

    return pos_eci, vel_eci


def eci_to_ecef(pos_eci, vel_eci, mjd_utc):
    """Convert position and velocity from ECI (GCRS) to ECEF (ITRF).

    Args:
        pos_eci: Position vector(s) in GCRS [m]
        vel_eci: Velocity vector(s) in GCRS [m/s]
        mjd_utc: Modified Julian Date (UTC)

    Returns:
        pos_ecef: Position in ITRF [m]
        vel_ecef: Velocity in ITRF [m/s]
    """
    obs_time = Time(mjd_utc, format='mjd', scale='utc')

    rep = CartesianRepresentation(
        pos_eci * u.m,
        differentials={'s': CartesianDifferential(vel_eci * u.m / u.s)}
    )
    gcrs = GCRS(rep, obstime=obs_time)
    itrs = gcrs.transform_to(ITRS(obstime=obs_time))

    pos_ecef = itrs.cartesian.xyz.to(u.m).value.T
    vel_ecef = itrs.velocity.d_xyz.to(u.m / u.s).value.T

    return pos_ecef, vel_ecef
