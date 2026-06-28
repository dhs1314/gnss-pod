"""Satellite configuration database for multi-satellite POD support (Phase 13.0).

Provides standard physical parameters for LEO satellites with onboard
GNSS receivers.  Each satellite entry includes:
  - Mass [kg], cross-sectional area [m²] (drag, SRP)
  - Surface coefficients: CD (drag), CR (SRP)
  - GNSS receiver type and supported systems
  - Orbit altitude [km] and inclination [deg]

Sources:
  - GRACE-FO: JPL L1B Product Handbook (GRACE-FO mission)
  - SWARM: ESA Swarm Product Data Handbook
  - GRACE: JPL L1B Product Handbook
  - FY-3C/D: CMA/NSMC documentation
  - Jason-3: CNES/NASA POD standards
  - COSMIC-2: UCAR/COSMIC documentation
  - Sentinel-3: ESA Sentinel-3 POD service

Usage:
    from satellite_config import get_config
    cfg = get_config('GRACE-FO', 'C')
    ekf_cfg.update(cfg['ekf_params'])
"""

SATELLITE_DB = {
    # ═══ GRACE-FO (2018-present) ═══
    'GRACE-FO': {
        'C': {
            'mass_kg': 580.0,
            'area_drag_m2': 0.68,
            'area_srp_m2': 3.4,
            'CD': 2.2,
            'CR': 1.3,
            'altitude_km': 490,
            'inclination_deg': 89.0,
            'gnss_systems': ['GPS'],
            'receiver': 'BlackJack (modified)',
            'description': 'GRACE Follow-On, satellite C (trailing)',
        },
        'D': {
            'mass_kg': 580.0,
            'area_drag_m2': 0.68,
            'area_srp_m2': 3.4,
            'CD': 2.2,
            'CR': 1.3,
            'altitude_km': 490,
            'inclination_deg': 89.0,
            'gnss_systems': ['GPS'],
            'receiver': 'BlackJack (modified)',
            'description': 'GRACE Follow-On, satellite D (leading)',
        },
    },

    # ═══ GRACE (2002-2017) ═══
    'GRACE': {
        'A': {
            'mass_kg': 487.0,
            'area_drag_m2': 0.50,
            'area_srp_m2': 3.1,
            'CD': 2.3,
            'CR': 1.3,
            'altitude_km': 490,
            'inclination_deg': 89.0,
            'gnss_systems': ['GPS'],
            'receiver': 'BlackJack',
            'description': 'GRACE, satellite A (trailing)',
        },
        'B': {
            'mass_kg': 487.0,
            'area_drag_m2': 0.50,
            'area_srp_m2': 3.1,
            'CD': 2.3,
            'CR': 1.3,
            'altitude_km': 490,
            'inclination_deg': 89.0,
            'gnss_systems': ['GPS'],
            'receiver': 'BlackJack',
            'description': 'GRACE, satellite B (leading)',
        },
    },

    # ═══ SWARM (2013-present) ═══
    'SWARM': {
        'A': {
            'mass_kg': 473.0,
            'area_drag_m2': 1.0,
            'area_srp_m2': 4.5,
            'CD': 2.4,
            'CR': 1.2,
            'altitude_km': 460,
            'inclination_deg': 87.4,
            'gnss_systems': ['GPS'],
            'receiver': 'RUAG 8-channel GPS',
            'description': 'SWARM Alpha (lower pair, 460 km)',
        },
        'B': {
            'mass_kg': 473.0,
            'area_drag_m2': 1.0,
            'area_srp_m2': 4.5,
            'CD': 2.4,
            'CR': 1.2,
            'altitude_km': 460,
            'inclination_deg': 87.4,
            'gnss_systems': ['GPS'],
            'receiver': 'RUAG 8-channel GPS',
            'description': 'SWARM Bravo (lower pair, 460 km)',
        },
        'C': {
            'mass_kg': 473.0,
            'area_drag_m2': 1.0,
            'area_srp_m2': 4.5,
            'CD': 2.4,
            'CR': 1.2,
            'altitude_km': 530,
            'inclination_deg': 87.4,
            'gnss_systems': ['GPS'],
            'receiver': 'RUAG 8-channel GPS',
            'description': 'SWARM Charlie (upper, 530 km)',
        },
    },

    # ═══ FengYun-3C/3D (Chinese meteorological, BDS+GPS) ═══
    'FY-3': {
        'C': {
            'mass_kg': 2250.0,
            'area_drag_m2': 4.5,
            'area_srp_m2': 12.0,
            'CD': 2.2,
            'CR': 1.3,
            'altitude_km': 836,
            'inclination_deg': 98.8,
            'gnss_systems': ['GPS', 'BDS'],
            'receiver': 'GNOS (BDS+GPS dual-system)',
            'description': 'FengYun-3C, sun-sync orbit, BDS+GPS receiver',
        },
        'D': {
            'mass_kg': 2250.0,
            'area_drag_m2': 4.5,
            'area_srp_m2': 12.0,
            'CD': 2.2,
            'CR': 1.3,
            'altitude_km': 836,
            'inclination_deg': 98.8,
            'gnss_systems': ['GPS', 'BDS'],
            'receiver': 'GNOS (BDS+GPS dual-system)',
            'description': 'FengYun-3D, sun-sync orbit, BDS+GPS receiver',
        },
    },

    # ═══ COSMIC-2 (2019-present, 6 satellites) ═══
    'COSMIC-2': {
        'E1': {
            'mass_kg': 278.0,
            'area_drag_m2': 0.8,
            'area_srp_m2': 2.0,
            'CD': 2.2,
            'CR': 1.2,
            'altitude_km': 520,
            'inclination_deg': 24.0,
            'gnss_systems': ['GPS', 'GLONASS'],
            'receiver': 'TGRS (TriG GNSS Receiver System)',
            'description': 'COSMIC-2 FM1, low-inclination, GPS+GLO',
        },
    },

    # ═══ Jason-3 (2016-present, altimetry reference) ═══
    'Jason-3': {
        'A': {
            'mass_kg': 510.0,
            'area_drag_m2': 1.2,
            'area_srp_m2': 5.0,
            'CD': 2.3,
            'CR': 1.4,
            'altitude_km': 1336,
            'inclination_deg': 66.0,
            'gnss_systems': ['GPS'],
            'receiver': 'GPSP (TurboRogue)',
            'description': 'Jason-3, 1336 km altimetry reference orbit',
        },
    },
}


def get_config(mission, sat_id=None):
    """Get satellite configuration for a mission and satellite ID.

    Args:
        mission: 'GRACE-FO', 'GRACE', 'SWARM', 'FY-3', 'COSMIC-2', 'Jason-3'
        sat_id: 'C', 'D', 'A', 'B', etc. (default: first satellite)

    Returns:
        dict with: mass_kg, area_drag_m2, area_srp_m2, CD, CR,
                   altitude_km, inclination_deg, gnss_systems, receiver, description
    """
    if mission not in SATELLITE_DB:
        raise ValueError(f"Unknown mission: {mission}. "
                         f"Available: {list(SATELLITE_DB.keys())}")

    mission_sats = SATELLITE_DB[mission]
    if sat_id is None:
        sat_id = sorted(mission_sats.keys())[0]

    if sat_id not in mission_sats:
        raise ValueError(f"Unknown satellite {mission}-{sat_id}. "
                         f"Available: {list(mission_sats.keys())}")

    return dict(mission_sats[sat_id])


def list_satellites():
    """List all available satellite configurations."""
    for mission, sats in SATELLITE_DB.items():
        for sid, cfg in sats.items():
            print(f"  {mission}-{sid}: {cfg['description']} "
                  f"({cfg['altitude_km']}km, {', '.join(cfg['gnss_systems'])})")


def get_ekf_params(mission, sat_id=None):
    """Get satellite-specific EKF initialization parameters."""
    cfg = get_config(mission, sat_id)
    return {
        'mass': cfg['mass_kg'],
        'area_drag': cfg['area_drag_m2'],
        'area_srp': cfg['area_srp_m2'],
        'CD': cfg['CD'],
        'CR': cfg['CR'],
    }
