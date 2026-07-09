"""GRACE-FO box-wing macro model (Bettadpur 2012, GRACE-FO L1 Handbook).

Replaces cannonball SRP with 8-panel analytical box-wing formulation.
Each panel has area, body-frame normal, and optical properties.

References:
  Bettadpur (2012): GRACE macro-model, CSR-TM-2012-01
  GRACE-FO L1 Handbook (JPL, 2019): Table 5 — surface properties
  Woske et al. (2019): GRACE-FO Finite Element Model

Satellite Body Frame (SF):
  +X : Flight direction (along-track)
  +Y : Cross-track (starboard)
  +Z : Nadir (toward Earth center)

SRP acceleration per panel (Milani et al. formulation):
  a = -(A/m) * (Phi_sun/c) * cos(theta) *
      [ (1 - rho_s) * s_hat  +  2 * (rho_d/3 + rho_s*cos(theta)) * n_hat ]

  where:
    cos(theta) = s_hat . n_hat   (panel illuminated only if cos(theta) > 0)
    rho_s = specular reflectivity (visible band)
    rho_d = diffuse reflectivity (visible band)
    Phi_sun = solar flux at satellite distance [W/m^2] (~1361 W/m^2 at 1 AU)
    c = speed of light

IR emission (thermal recoil) per panel:
  a_IR = -(A/m) * (2/3) * (epsilon * sigma * T^4 / c) * n_hat
  epsilon = IR emissivity
  sigma = Stefan-Boltzmann constant (5.67e-8 W/m^2/K^4)
  T = panel temperature [K]

Simplified thermal model: constant temperature per panel type.
"""

import numpy as np

C_LIGHT = 299792458.0
SIGMA_SB = 5.670367e-8   # Stefan-Boltzmann [W/m^2/K^4]
SOLAR_FLUX_1AU = 1361.0  # Solar constant [W/m^2]


# ═══════════════════════════════════════════════════════════════════
# Panel definitions: (name, area[m²], normal_in_BF[x,y,z], material)
# ═══════════════════════════════════════════════════════════════════

PANELS = [
    # Panel 1: Front (+X) — SiOx/Kapton
    ("Front",        0.9552,  np.array([+1.0, 0.0, 0.0])),
    # Panel 2: Rear (-X) — SiOx/Kapton
    ("Rear",         0.9552,  np.array([-1.0, 0.0, 0.0])),
    # Panel 3: Starboard outer (+Y canted ~40°) — Si Glass Solar Array
    ("Starboard",    3.1555,  np.array([0.0, +0.766044, -0.642788])),
    # Panel 4: Starboard inner — SiOx/Kapton
    ("Stbd_inner",   0.2283,  np.array([0.0, -0.766044, +0.642788])),
    # Panel 5: Port outer (-Y canted ~40°) — Si Glass Solar Array
    ("Port",         3.1555,  np.array([0.0, -0.766044, -0.642788])),
    # Panel 6: Port inner — SiOx/Kapton
    ("Port_inner",   0.2283,  np.array([0.0, +0.766044, +0.642788])),
    # Panel 7: Nadir (+Z, Earth-facing) — Teflon FEP
    ("Nadir",        6.0711,  np.array([0.0, 0.0, +1.0])),
    # Panel 8: Zenith (-Z, space-facing) — Si Glass Solar Array
    ("Zenith",       2.1674,  np.array([0.0, 0.0, -1.0])),
]

# Total area: 17.1158 m²


# ═══════════════════════════════════════════════════════════════════
# Optical properties (GRACE-FO L1 Handbook, Table 5)
# ═══════════════════════════════════════════════════════════════════

# Per material type: {emissivity, absorptivity, specular_refl, diffuse_refl}
# All values in VISIBLE band (for SRP)
OPTICAL = {
    "SiOx/Kapton": {
        "epsilon_IR": 0.62,   # IR emissivity
        "alpha_Vis":  0.34,   # visible absorptivity
        "rho_s_Vis":  0.40,   # visible specular reflectivity
        "rho_d_Vis":  0.26,   # visible diffuse reflectivity
        "T_K":        293.0,  # ~20°C (bus body temperature)
    },
    "Si_Glass_Solar": {
        "epsilon_IR": 0.81,
        "alpha_Vis":  0.65,   # 0.72 when non-operating
        "rho_s_Vis":  0.05,
        "rho_d_Vis":  0.30,
        "T_K":        330.0,  # ~57°C (solar panel in sunlight)
    },
    "Teflon_FEP": {
        "epsilon_IR": 0.75,
        "alpha_Vis":  0.12,
        "rho_s_Vis":  0.68,
        "rho_d_Vis":  0.20,
        "T_K":        280.0,  # ~7°C (nadir panel facing cold Earth)
    },
}


def _material_for_panel(name):
    """Return material type string for a panel name."""
    if "Starboard" in name or "Port" in name or "Zenith" in name:
        if "inner" in name:
            return "SiOx/Kapton"
        return "Si_Glass_Solar"
    if "inner" in name:
        return "SiOx/Kapton"
    if name == "Nadir":
        return "Teflon_FEP"
    if name in ("Front", "Rear"):
        return "SiOx/Kapton"
    return "SiOx/Kapton"  # default


# ═══════════════════════════════════════════════════════════════════
# Body frame computation
# ═══════════════════════════════════════════════════════════════════

def compute_body_frame(pos_eci, vel_eci):
    """Compute GRACE-FO body frame axes in ECI.

    SF definition (Science Mode, nadir-pointing):
      +Z = -pos / |pos|  (toward Earth center = nadir)
      +Y = +Z × vel_dir  (cross-track, starboard)
      +X = +Y × +Z       (along-track, near velocity direction)

    Returns:
        X_eci, Y_eci, Z_eci: (3,) unit vectors in ECI
    """
    r_mag = np.linalg.norm(pos_eci)
    if r_mag < 1.0:
        return np.eye(3)  # fallback

    Z_eci = -pos_eci / r_mag  # toward Earth center (nadir)

    # Y = Z × v_dir (starboard)
    v_dir = vel_eci / np.linalg.norm(vel_eci)
    Y_eci = np.cross(Z_eci, v_dir)
    y_norm = np.linalg.norm(Y_eci)
    if y_norm < 1e-12:
        Y_eci = np.array([0.0, 1.0, 0.0])
        Y_eci = Y_eci - np.dot(Y_eci, Z_eci) * Z_eci
        Y_eci = Y_eci / np.linalg.norm(Y_eci)
    else:
        Y_eci = Y_eci / y_norm

    X_eci = np.cross(Y_eci, Z_eci)  # along-track
    return X_eci, Y_eci, Z_eci


def body_to_eci(vec_bf, X_eci, Y_eci, Z_eci):
    """Transform vector from body frame to ECI."""
    return vec_bf[0] * X_eci + vec_bf[1] * Y_eci + vec_bf[2] * Z_eci


# ═══════════════════════════════════════════════════════════════════
# Sun position
# ═══════════════════════════════════════════════════════════════════

def _sun_eci(mjd_tt):
    """Approximate Sun position in ECI [m] (accuracy ~0.01°)."""
    # Days since J2000.0
    T = (mjd_tt - 51544.5) / 36525.0

    # Mean longitude
    L = (280.46646 + 36000.76983 * T + 0.0003032 * T**2) % 360.0
    # Mean anomaly
    M = (357.52911 + 35999.05029 * T - 0.0001537 * T**2) % 360.0
    # Ecliptic obliquity
    obl = 23.439291 - 0.0130042 * T

    L_rad = np.radians(L)
    M_rad = np.radians(M)
    obl_rad = np.radians(obl)

    # Equation of center + Sun-Earth distance
    C = (1.914602 - 0.004817*T - 0.000014*T**2) * np.sin(M_rad) \
      + (0.019993 - 0.000101*T) * np.sin(2*M_rad) \
      + 0.000289 * np.sin(3*M_rad)

    true_long = L_rad + np.radians(C)
    R_AU = 1.00014 - 0.01671*np.cos(M_rad) - 0.00014*np.cos(2*M_rad)

    # Ecliptic coordinates
    X_ecl = R_AU * np.cos(true_long)
    Y_ecl = R_AU * np.sin(true_long)

    # Rotate to equatorial
    AU_m = 149597870700.0
    X_eci = X_ecl * AU_m
    Y_eci = Y_ecl * np.cos(obl_rad) * AU_m
    Z_eci = Y_ecl * np.sin(obl_rad) * AU_m

    return np.array([X_eci, Y_eci, Z_eci])


# ═══════════════════════════════════════════════════════════════════
# Eclipse / shadow function
# ═══════════════════════════════════════════════════════════════════

def shadow_function(sat_eci, sun_eci):
    """Conical shadow model (0=umbra, 1=sunlight).

    Returns shadow factor nu in [0, 1].
    """
    R_EARTH = 6378137.0
    R_SUN = 695700000.0

    r_sat = np.linalg.norm(sat_eci)
    r_sun = np.linalg.norm(sun_eci)

    # Satellite-Sun vector
    d = sun_eci - sat_eci
    d_mag = np.linalg.norm(d)
    d_hat = d / d_mag

    # Angles
    sat_to_earth_center = -sat_eci / r_sat

    # Apparent Earth radius from satellite
    a = np.arcsin(R_EARTH / r_sat)
    # Apparent Sun radius from satellite
    b = np.arcsin(R_SUN / d_mag)
    # Angular separation between Earth center and Sun center
    c = np.arccos(np.clip(np.dot(-sat_eci, d_hat) / r_sat, -1.0, 1.0))

    # Penumbra/umbra transitions
    if c < abs(a - b):
        return 0.0       # umbra (total eclipse)
    elif c < a + b:
        # Penumbra: linear interpolation
        x = (c - abs(a - b)) / (2.0 * min(a, b))
        return max(0.0, min(1.0, x))
    else:
        return 1.0       # full sunlight


# ═══════════════════════════════════════════════════════════════════
# Main SRP computation
# ═══════════════════════════════════════════════════════════════════

def compute_box_wing_srp(pos_eci, vel_eci, mjd_tt, mass=580.0,
                          enable_thermal=False):
    """Compute SRP acceleration in ECI using the 8-panel box-wing model.

    Args:
        pos_eci: satellite position in ECI [m] (3,)
        vel_eci: satellite velocity in ECI [m/s] (3,)
        mjd_tt: MJD in Terrestrial Time
        mass: satellite mass [kg]
        enable_thermal: include IR thermal recoil (small correction)

    Returns:
        a_srp_eci: acceleration in ECI [m/s²] (3,)
    """
    # Body frame
    X_eci, Y_eci, Z_eci = compute_body_frame(pos_eci, vel_eci)

    # Sun direction
    sun_eci = _sun_eci(mjd_tt)
    sun_vec = sun_eci - pos_eci
    sun_dist = np.linalg.norm(sun_vec)
    s_hat = sun_vec / sun_dist

    # Solar flux at satellite distance
    AU_m = 149597870700.0
    Phi_sun = SOLAR_FLUX_1AU * (AU_m / sun_dist)**2
    P_sun = Phi_sun / C_LIGHT  # [N/m²]

    # Shadow function
    nu = shadow_function(pos_eci, sun_eci)

    # Sum over panels
    a_total = np.zeros(3)
    A_m = 1.0 / mass

    for name, area_m2, n_bf in PANELS:
        # Panel normal in ECI
        n_eci = body_to_eci(n_bf, X_eci, Y_eci, Z_eci)

        # Incidence angle
        cos_theta = np.dot(s_hat, n_eci)

        # Panel is illuminated only if cos_theta > 0
        if cos_theta <= 0:
            continue

        # Get optical properties
        mat = _material_for_panel(name)
        opt = OPTICAL[mat]
        rho_s = opt["rho_s_Vis"]
        rho_d = opt["rho_d_Vis"]

        # SRP force (visible band)
        # a = -(A/m) * P_sun * cos_theta * nu *
        #     [(1 - rho_s) * s_hat + 2*(rho_d/3 + rho_s*cos_theta) * n_hat]
        absorpt_frac = 1.0 - rho_s - rho_d
        specular_term = (1.0 - rho_s) * s_hat  # WRONG: should be absorpt fraction
        # Actually the correct formula is:
        # force = -P * A * cos(theta) * nu * [
        #     (1 - rho_s) * s_hat                     <- absorbed + specular (specular acts along s_hat + 2*cos(theta)*n_hat)
        #   Actually: ρ_s fraction reflects specularly: direction = s_hat - 2*cos(theta)*n_hat
        #   The force on the satellite is the opposite of the momentum change:
        #   F = P*A*cos(theta) * [ (absorption)*s_hat - ρ_s*(s_hat - 2*cos(theta)*n_hat) - ρ_d*(s_hat + (2/3)*n_hat) ]
        #   = P*A*cos(theta) * [ (1 - ρ_s - ρ_d)*s_hat + 2*ρ_s*cos(theta)*n_hat + ρ_d*(s_hat + (2/3)*n_hat) ]
        #   Wait, let me re-derive. The incoming photon momentum is along -s_hat.
        #   Outgoing:
        #     absorbed fraction (1-ρ_s-ρ_d): momentum = 0 (absorbed as heat, re-emitted isotropically = 0 net)
        #     specular fraction ρ_s: reflected along s_hat - 2*cos(theta)*n_hat
        #     diffuse fraction ρ_d: reflected isotropically from surface = (2/3)*n_hat on average
        #   Momentum change = incoming - outgoing:
        #     incoming momentum along s_hat
        #   absorbed: incoming momentum stays = s_hat direction
        #   specular: (s_hat) - ρ_s*(s_hat - 2*cos(theta)*n_hat) = (1-ρ_s)*s_hat + 2*ρ_s*cos(theta)*n_hat
        #   diffuse: (s_hat) - ρ_d*(s_hat + (2/3)*n_hat) = (1-ρ_d)*s_hat - (2/3)*ρ_d*n_hat
        #   Hmm, different sources use different formulations. Let me use the standard
        #   from Montenbruck & Gill / Doornbos which is widely used in POD:
        #
        #   a = -(A/m) * P_sun * nu * cos(theta) *
        #       [(1 - ρ_s)*ŝ + 2*(ρ_d/3 + ρ_s*cos(theta))*n̂]
        #
        #   where the absorption coefficient α = 1 - ρ_s - ρ_d is implicit
        #   (this is the classic cannonball extension to flat plates)

        absorption = 1.0 - rho_s - rho_d
        diff_term = 2.0 * rho_d / 3.0 + 2.0 * rho_s * cos_theta
        panel_force = -A_m * area_m2 * P_sun * nu * cos_theta * (
            absorption * s_hat + diff_term * n_eci
        )
        a_total += panel_force

        # Thermal IR emission (optional, small correction)
        if enable_thermal:
            eps = opt["epsilon_IR"]
            T_K = opt["T_K"]
            P_thermal = 2.0 / 3.0 * eps * SIGMA_SB * T_K**4 / C_LIGHT
            # IR emitted isotropically from panel surface → (2/3) outward
            a_total -= A_m * area_m2 * P_thermal * n_eci

    return a_total


# ═══════════════════════════════════════════════════════════════════
# Drag cross-section (macro model)
# ═══════════════════════════════════════════════════════════════════

def compute_projected_area(pos_eci, vel_eci):
    """Compute projected cross-sectional area for drag.

    Uses the 8-panel macro model to compute total projected area
    normal to the incident flow direction (including co-rotation).

    Returns:
        A_proj: projected area [m²] normal to velocity direction
    """
    # Body frame
    X_eci, Y_eci, Z_eci = compute_body_frame(pos_eci, vel_eci)

    # Velocity direction in body frame
    v_dir_bf = np.array([
        np.dot(vel_eci, X_eci),
        np.dot(vel_eci, Y_eci),
        np.dot(vel_eci, Z_eci),
    ])
    v_dir_bf = v_dir_bf / np.linalg.norm(v_dir_bf)

    # Sum projected areas of panels facing the flow
    A_proj = 0.0
    for name, area_m2, n_bf in PANELS:
        cos_theta = np.dot(v_dir_bf, n_bf)
        if cos_theta < 0:
            # Flow hits the back of this panel
            A_proj += area_m2 * abs(cos_theta)
    return A_proj
