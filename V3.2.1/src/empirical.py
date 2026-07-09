"""Empirical acceleration model (RTN piecewise-constant).

In reduced-dynamic POD, empirical accelerations absorb the residual
force modeling errors. They are parameterized as piecewise-constant
accelerations in the RTN frame:

    R = Radial (away from Earth center)
    T = Transverse / Along-track (velocity direction)
    N = Normal (cross-track, completes right-handed system)

The intervals are typically 15 minutes, giving 96 x 3 = 288 parameters
per 24-hour arc.
"""
import numpy as np


def compute_rtn_frame(pos_eci, vel_eci):
    """Compute RTN (Radial, Transverse, Normal) unit vectors.

    R = pos / |pos|                          (radial)
    T = (vel - radial_component) / norm      (along-track)
    N = R x T                                 (cross-track)

    Args:
        pos_eci: Position in ECI [m], shape (3,) or (N,3)
        vel_eci: Velocity in ECI [m/s], shape (3,) or (N,3)

    Returns:
        R, T, N: Unit vectors in ECI, each shape (3,) or (N,3)
    """
    squeeze = pos_eci.ndim == 1
    if squeeze:
        pos_eci = pos_eci.reshape(1, 3)
        vel_eci = vel_eci.reshape(1, 3)

    pos_mag = np.linalg.norm(pos_eci, axis=-1, keepdims=True)
    R = pos_eci / pos_mag

    vel_radial = np.sum(vel_eci * R, axis=-1, keepdims=True) * R
    vel_tangential = vel_eci - vel_radial

    T = vel_tangential / (np.linalg.norm(vel_tangential, axis=-1, keepdims=True) + 1e-15)
    N = np.cross(R, T)
    N = N / (np.linalg.norm(N, axis=-1, keepdims=True) + 1e-15)

    if squeeze:
        return R[0], T[0], N[0]
    return R, T, N


def rtn_to_eci(acc_rtn, pos_eci, vel_eci):
    """Convert RTN-frame accelerations to ECI frame.

    Args:
        acc_rtn: (a_R, a_T, a_N) in RTN frame [m/s^2], shape (3,)
        pos_eci: Position in ECI [m], shape (3,)
        vel_eci: Velocity in ECI [m/s], shape (3,)

    Returns:
        acc_eci: Acceleration in ECI frame [m/s^2], shape (3,)
    """
    R, T, N = compute_rtn_frame(pos_eci, vel_eci)
    return acc_rtn[0] * R + acc_rtn[1] * T + acc_rtn[2] * N


def eci_to_rtn(acc_eci, pos_eci, vel_eci):
    """Convert ECI acceleration to RTN frame.

    Args:
        acc_eci: Acceleration in ECI frame [m/s^2], shape (3,)
        pos_eci: Position in ECI [m], shape (3,)
        vel_eci: Velocity in ECI [m/s], shape (3,)

    Returns:
        acc_rtn: (a_R, a_T, a_N) in RTN frame [m/s^2], shape (3,)
    """
    R, T, N = compute_rtn_frame(pos_eci, vel_eci)
    return np.array([np.dot(acc_eci, R), np.dot(acc_eci, T), np.dot(acc_eci, N)])


def build_empirical_intervals(arc_duration_seconds, interval_minutes=15.0):
    """Build time boundaries for piecewise-constant empirical parameters.

    Args:
        arc_duration_seconds: Total arc duration [s]
        interval_minutes: Interval length [minutes]

    Returns:
        boundaries: Array of interval boundary times [s] from arc start
    """
    interval_seconds = interval_minutes * 60.0
    n_intervals = int(np.ceil(arc_duration_seconds / interval_seconds))
    boundaries = np.linspace(0, n_intervals * interval_seconds, n_intervals + 1)
    return boundaries


def get_empirical_interval_index(t, boundaries):
    """Find which empirical interval a given time belongs to.

    Args:
        t: Time from arc start [s]
        boundaries: Interval boundaries [s]

    Returns:
        index: Interval index (0-based)
    """
    return int(np.searchsorted(boundaries, t) - 1)
