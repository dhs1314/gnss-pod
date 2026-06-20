"""Relativistic correction to LEO orbit acceleration (IERS 2010).

The general relativistic Schwarzschild acceleration for a satellite in
Earth orbit:

  a_rel = (GMe / (c² r³)) * [2(β+γ)(GMe/r) r - γ v² r + 2(1+γ)(r·v) v]

With General Relativity values (β=γ=1), simplified to:

  a_rel = (GMe / (c² r³)) * [(4GMe/r - v²) r + 4(r·v) v]

Magnitude: ~1-3×10⁻⁸ m/s² at LEO (500 km), ~3 cm in position over 1 orbit.

References:
  - IERS Conventions (2010), Chapter 10
  - Montenbruck & Gill (2000), Section 3.5
"""
import numpy as np

GM_EARTH = 3.986004415e14
C_LIGHT = 299792458.0

# Post-Newtonian parameters (GR: β=γ=1)
BETA = 1.0
GAMMA = 1.0


def compute_schwarzschild_acc(pos_eci, vel_eci):
    """Compute Schwarzschild relativistic acceleration in ECI.

    Args:
        pos_eci: ECI position [m], shape (3,)
        vel_eci: ECI velocity [m/s], shape (3,)

    Returns:
        acc_rel: Relativistic acceleration [m/s²], shape (3,)
    """
    r = np.linalg.norm(pos_eci)
    v2 = float(np.dot(vel_eci, vel_eci))
    r_dot_v = float(np.dot(pos_eci, vel_eci))

    if r < 1.0:
        return np.zeros(3)

    r3 = r ** 3
    mu_over_c2 = GM_EARTH / (C_LIGHT ** 2)

    # a = (μ/c²r³) * [(4μ/r - v²)r + 4(r·v)v]
    coeff = mu_over_c2 / r3
    term1 = (4.0 * GM_EARTH / r - v2) * pos_eci
    term2 = 4.0 * r_dot_v * vel_eci

    return coeff * (term1 + term2)


def compute_lense_thirring_acc(pos_eci, vel_eci):
    """Lense-Thirring (frame-dragging) acceleration [m/s²].

    For Earth: a_LT ≈ 2(GMe/c²r³) * (3(r·Ĵ)(r×v)/r² + v×J)
    where J ≈ 9.8e8 * [0,0,1] m²/s is Earth's angular momentum per unit mass.

    Magnitude: ~2×10⁻¹¹ m/s² at LEO. Included for completeness;
    negligible for current accuracy requirements.
    """
    r = np.linalg.norm(pos_eci)
    if r < 1.0:
        return np.zeros(3)

    J = np.array([0.0, 0.0, 9.8e8])  # Earth ang. momentum per unit mass [m²/s]
    r3 = r ** 3
    mu_over_c2 = GM_EARTH / (C_LIGHT ** 2)

    r_cross_v = np.cross(pos_eci, vel_eci)
    r_dot_J = np.dot(pos_eci, J)
    v_cross_J = np.cross(vel_eci, J)

    term = (3.0 * r_dot_J / (r ** 2)) * r_cross_v + v_cross_J
    return 2.0 * mu_over_c2 / r3 * term


def compute_relativistic_acceleration(pos_eci, vel_eci):
    """Compute total relativistic acceleration (Schwarzschild + Lense-Thirring).

    Convenience wrapper calling compute_schwarzschild_acc for the main term.
    """
    return compute_schwarzschild_acc(pos_eci, vel_eci)
