"""Convert SWARM RINEX 3.00 GPS observations to standard GPS1B pickle format.

Input:  SWARM RINEX 3.00 observation file (.rnx)
Output: GPS1B_YYYY-MM-DD_SWARM_A.pkl

RINEX 3 format specifics (SWARM GPSR receiver):
  - Header: SYS / # / OBS TYPES = G 8 C1C L1C S1C C1W S1W C2W L2W S2W
  - Epoch lines start with '>'
  - SV PRNs are numeric (1-32) without 'G' prefix
  - 8 values per SV per epoch

Usage:
  python scripts/convert_swarm_rnx3.py --date 2014-04-29
  python scripts/convert_swarm_rnx3.py --all
"""
import sys, os, pickle, re, argparse, numpy as np
from pathlib import Path
from datetime import datetime, timedelta

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / 'src'))

C_LIGHT = 299792458.0
F1 = 1575.42e6
F2 = 1227.60e6
GAMMA = (F1 / F2) ** 2
GPS_ORIGIN = datetime(1980, 1, 6)

parser = argparse.ArgumentParser(description='Convert SWARM RINEX 3 to GPS1B pkl')
parser.add_argument('--date', help='Date YYYY-MM-DD')
parser.add_argument('--all', action='store_true', help='Process all available dates')
parser.add_argument('--sat', default='A', choices=['A','B','C'])
args = parser.parse_args()

SWARM_DIR = ROOT / 'data' / 'swarm'

def get_dates():
    if args.date:
        return [args.date]
    dates = []
    for d in sorted(SWARM_DIR.glob('2014/2014-*')):
        if (d / f'SWARM_A_{d.name}.rnx').exists():
            dates.append(d.name)
    return dates

DATES = get_dates()
for date_str in DATES:
    dt = datetime.strptime(date_str, '%Y-%m-%d')
    y = dt.year
    in_dir = SWARM_DIR / str(y) / date_str
    rnx_path = in_dir / f'SWARM_{args.sat}_{date_str}.rnx'
    if not rnx_path.exists():
        print(f'{date_str}: RINEX not found at {rnx_path}')
        continue

    print(f'{date_str}: reading {rnx_path.name} ({rnx_path.stat().st_size//1024}KB)...')

    with open(str(rnx_path), encoding='utf-8', errors='replace') as f:
        lines = f.readlines()

    # Parse header
    obs_map = {}
    h_end = 0
    for i, line in enumerate(lines):
        tag = line[60:].strip()
        if 'END OF HEADER' in tag:
            h_end = i
            break
        if 'SYS / # / OBS TYPES' in tag:
            # Format: "G    8 C1C L1C S1C C1W S1W C2W L2W S2W"
            parts = line[:60].split()
            if len(parts) >= 3 and parts[0] == 'G':
                n_obs = int(parts[1])
                obs_types = parts[2:2+n_obs]
                for idx, name in enumerate(obs_types):
                    if name in ('C1C', 'C1W'): obs_map['P1'] = idx
                    if name == 'C2W': obs_map['P2'] = idx
                    if name == 'L1C': obs_map['L1'] = idx
                    if name == 'L2W': obs_map['L2'] = idx
                    if name in ('S1C', 'S1W'): obs_map.setdefault('S1', idx)
                    if name == 'S2W': obs_map.setdefault('S2', idx)

    print(f'  Obs types: {obs_types}, Map: P1={obs_map.get("P1")} P2={obs_map.get("P2")} L1={obs_map.get("L1")} L2={obs_map.get("L2")}')

    # Estimate GPS week from date
    days_since = (dt - GPS_ORIGIN).days
    gps_week = days_since // 7
    sow_base = (days_since % 7) * 86400

    # Parse body
    gps_obs = {}
    epoch_count = 0
    i = h_end + 1

    while i < len(lines):
        line = lines[i].rstrip()
        if not line or line.startswith('COMMENT'):
            i += 1; continue

        if line.startswith('>'):
            # Parse epoch header: "> 2014 04 29 00 00 19.0000000  0  7"
            parts = line[1:].split()
            if len(parts) < 3: i += 1; continue
            try:
                yr = int(parts[0]); mo = int(parts[1]); dy = int(parts[2])
                hr = int(parts[3]); mn = int(parts[4]); sc = float(parts[5])
                flag = int(parts[6]) if len(parts) > 6 else 0
                n_svs = int(parts[7]) if len(parts) > 7 else 0

                if yr < 80: yr += 2000

                # GPS seconds of week
                epoch_dt = datetime(yr, mo, dy, hr, mn, int(sc))
                days_since = (epoch_dt - GPS_ORIGIN).days
                gps_sec = days_since * 86400 + hr * 3600 + mn * 60 + int(sc)
                gps_sec = int(gps_sec)

                i += 1
                sv_count = 0
                while i < len(lines) and sv_count < n_svs:
                    data_line = lines[i].rstrip()
                    if data_line.startswith('>'):
                        break
                    if not data_line or data_line.startswith('COMMENT'):
                        i += 1; continue

                    # Parse G04 field: first char is 'G', then PRN number
                    if len(data_line) > 3 and data_line[0] == 'G':
                        sv_num = data_line[1:3].strip()
                        sv = f'G{int(sv_num):02d}' if sv_num.isdigit() else data_line[:3]
                    else:
                        i += 1; continue

                    # Parse values using RINEX 3.00 fixed-width format.
                    # Each observation: F14.3, I1(LLI), I1(SSI) = 16 chars.
                    # Whitespace-split mixes LLI/SSI digits (e.g. "6","4") into tokens.
                    vals = []
                    for oi in range(len(obs_types)):
                        start = 3 + oi * 16
                        block = data_line[start:start + 16] if start < len(data_line) else ''
                        val_str = block[:14].strip()
                        try:
                            vals.append(float(val_str))
                        except ValueError:
                            vals.append(0.0)

                    if len(vals) < len(obs_types):
                        i += 1; sv_count += 1; continue

                    p1_idx = obs_map.get('P1'); p2_idx = obs_map.get('P2')
                    l1_idx = obs_map.get('L1'); l2_idx = obs_map.get('L2')
                    s1_idx = obs_map.get('S1'); s2_idx = obs_map.get('S2')

                    if None in (p1_idx, p2_idx, l1_idx, l2_idx):
                        i += 1; sv_count += 1; continue

                    p1_m = float(vals[p1_idx])
                    p2_m = float(vals[p2_idx])
                    l1_cyc = float(vals[l1_idx])
                    l2_cyc = float(vals[l2_idx])
                    s1 = float(vals[s1_idx]) if s1_idx is not None and s1_idx < len(vals) else 999
                    s2 = float(vals[s2_idx]) if s2_idx is not None and s2_idx < len(vals) else 999

                    l1_m = l1_cyc * C_LIGHT / F1
                    l2_m = l2_cyc * C_LIGHT / F2

                    p_if = (GAMMA * p2_m - p1_m) / (GAMMA - 1)
                    l_if = (GAMMA * l2_m - l1_m) / (GAMMA - 1)

                    if gps_sec not in gps_obs:
                        gps_obs[gps_sec] = {}
                    gps_obs[gps_sec][sv] = {
                        'L1': l1_m, 'L2': l2_m,
                        'P1': p1_m, 'P2': p2_m,
                        'L_if': float(l_if), 'P_if': float(p_if),
                        'L1_cyc': float(l1_cyc), 'L2_cyc': float(l2_cyc),
                        'L1_SNR': float(s1), 'L2_SNR': float(s2),
                    }

                    i += 1; sv_count += 1
                epoch_count += 1
            except (ValueError, IndexError) as e:
                i += 1; continue
        else:
            i += 1

    # Save
    out_path = in_dir / f'GPS1B_{date_str}_SWARM_{args.sat}.pkl'
    with open(str(out_path), 'wb') as f:
        pickle.dump(gps_obs, f)

    epochs = sorted(gps_obs.keys())
    n_obs = sum(len(v) for v in gps_obs.values())
    print(f'  Saved: {epoch_count} epochs, {n_obs} obs -> {out_path}')

    if epochs:
        e0 = epochs[0]
        svs = sorted(gps_obs[e0].keys())
        print(f'  First epoch GPS {e0}: {len(svs)} SVs [{svs[0]}..{svs[-1]}]')
