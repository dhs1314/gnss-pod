"""SWARM data adapter: load SWARM GPS1B pickle into pipeline format.

Converts SWARM GPS1B dict (from convert_swarm_rnx.py) to eval_5day_orekit.py
input format. SWARM GPS data uses the same dict structure as GRACE-FO:
  {gps_sod_int: {sv: {L_if, P_if, L1, L2, P1, P2, L1_cyc, L2_cyc, ...}}}

This module handles satellite-specific config loading from satellite_config.py.
"""
import sys, os, pickle
from pathlib import Path


def load_swarm_gps1b(date_str, satellite='A', data_root=None):
    """Load SWARM GPS1B pickle for a given date and satellite.

    Args:
        date_str: 'YYYY-MM-DD'
        satellite: 'A', 'B', or 'C'
        data_root: optional Path override for data directory

    Returns:
        dict {gps_sod: {sv: {fields...}}} or None if not found
    """
    if data_root is None:
        data_root = Path('data')

    y = date_str[:4]
    pkl_path = data_root / 'swarm' / y / date_str / f'GPS1B_{date_str}_SWARM_{satellite}.pkl'

    if not pkl_path.exists():
        print(f"[SWARM] GPS1B not found: {pkl_path}")
        print(f"[SWARM] Run: python scripts/convert_swarm_rnx.py --date {date_str} --sat {satellite}")
        return None

    with open(str(pkl_path), 'rb') as f:
        return pickle.load(f)


def load_swarm_ref_orbit(date_str, satellite='A', data_root=None):
    """Load SWARM precise reference orbit (SP3 format) for validation.

    CODE produces combined LEO SP3 files containing SWARM-A/B/C orbits.
    Returns dict {gps_sod: np.array([x,y,z])} in ECEF.

    Args:
        date_str: 'YYYY-MM-DD'
        satellite: 'A', 'B', or 'C'
        data_root: optional Path override

    Returns:
        (pos_dict, vel_dict) or (None, None)
    """
    import numpy as np
    from datetime import datetime, timedelta

    if data_root is None:
        data_root = Path('data')

    y = date_str[:4]
    dt = datetime.strptime(date_str, '%Y-%m-%d')
    doy = dt.strftime('%j')

    # SWARM satellite IDs in CODE SP3: 'SA', 'SB', 'SC'
    sp3_sv_name = f'S{satellite}'

    # Check multiple possible locations for SWARM precise orbit
    candidates = [
        data_root / 'CODE' / y / f'COD0OPSFIN_{y}{doy}0000_01D_05M_LEO_ORB.SP3',
        data_root / 'swarm' / y / date_str / f'SWARM_{satellite}_{date_str}_ORB.SP3',
    ]

    for sp3_path in candidates:
        if not sp3_path.exists():
            continue

        pos_dict = {}
        vel_dict = {}
        gps_origin = datetime(1980, 1, 6)

        with open(str(sp3_path)) as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            if line.startswith('* '):
                parts = line.split()
                if len(parts) >= 6:
                    try:
                        ep_dt = datetime(int(parts[1]), int(parts[2]), int(parts[3]),
                                        int(parts[4]), int(parts[5]),
                                        int(float(parts[6])))
                        gps_sec = (ep_dt - gps_origin).total_seconds()
                    except: continue
            elif line.startswith('P'):
                sv = line[1:4].strip()
                if sv == sp3_sv_name:
                    try:
                        x = float(line[4:18]) * 1000  # km -> m
                        y = float(line[18:32]) * 1000
                        z = float(line[32:46]) * 1000
                        clk = float(line[46:60])  # microseconds
                        pos_dict[gps_sec] = np.array([x, y, z])
                    except: pass
            elif line.startswith('V'):
                sv = line[1:4].strip()
                if sv == sp3_sv_name:
                    try:
                        vx = float(line[4:18]) * 0.1  # dm/s -> m/s
                        vy = float(line[18:32]) * 0.1
                        vz = float(line[32:46]) * 0.1
                        vel_dict[gps_sec] = np.array([vx, vy, vz])
                    except: pass

        if pos_dict:
            print(f"[SWARM] Reference orbit loaded: {len(pos_dict)} epochs from {sp3_path.name}")
            return pos_dict, vel_dict

    print(f"[SWARM] No reference orbit found for SWARM-{satellite} on {date_str}")
    return None, None


def get_swarm_config(satellite='A'):
    """Get SWARM satellite dynamics parameters from satellite_config.py."""
    from src.satellite_config import get_config

    cfg = get_config('SWARM', satellite)
    return {
        'mass': cfg['mass_kg'],
        'area_drag': cfg['area_drag_m2'],
        'area_srp': cfg['area_srp_m2'],
        'CD': cfg['CD'],
        'CR': cfg['CR'],
    }
