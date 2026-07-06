"""Convert SWARM RINEX GPS observations to standard GPS1B pickle format.

Input:  SWARM RINEX 2.11 observation file (.rnx or .yyO)
Output: GPS1B_YYYY-MM-DD_SWARM_X.pkl (same dict format as GRACE-FO)

The converter handles:
  - RINEX 2.11 header parsing (obs types, interval, marker)
  - Epoch records with multi-SV observations
  - IF combination computation (L_if, P_if from L1/L2/P1/P2)
  - Cycle count extraction for MW computation

Usage:
  python scripts/convert_swarm_rnx.py --date 2024-04-29 --sat A
  python scripts/convert_swarm_rnx.py --date 2024-04-29 --sat A --rnx swarm.rnx
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

parser = argparse.ArgumentParser(description='Convert SWARM RINEX to GPS1B pkl')
parser.add_argument('--date', required=True, help='Date YYYY-MM-DD')
parser.add_argument('--sat', choices=['A','B','C'], required=True)
parser.add_argument('--rnx', help='Override RINEX file path')
args = parser.parse_args()

dt = datetime.strptime(args.date, '%Y-%m-%d')
y = dt.year

SWARM_DIR = ROOT / 'data' / 'swarm' / str(y) / args.date
SWARM_DIR.mkdir(parents=True, exist_ok=True)

# ── Find RINEX file ──
if args.rnx:
    rnx_path = Path(args.rnx)
else:
    candidates = sorted(SWARM_DIR.glob('*.rnx')) + sorted(SWARM_DIR.glob('*.yyO')) + \
                 sorted(SWARM_DIR.glob('*RINEX*')) + sorted(SWARM_DIR.glob('*GPS*'))
    if not candidates:
        print(f"ERROR: No RINEX files found in {SWARM_DIR}")
        print("Download SWARM data first: python scripts/download_swarm_data.py")
        sys.exit(1)
    rnx_path = candidates[0]

print(f"Loading: {rnx_path}")

# ── Parse RINEX 2.11 header ──
with open(str(rnx_path), encoding='utf-8', errors='replace') as f:
    lines = f.readlines()

obs_types = []
header_end = 0
for i, line in enumerate(lines):
    tag = line[60:].strip()
    if 'END OF HEADER' in tag:
        header_end = i
        break
    if '# / TYPES OF OBSERV' in tag:
        try:
            n_obs = int(line[:6].strip())
            toks = line[6:60].split()
            obs_types.extend([t.strip() for t in toks if t.strip()])
        except: pass
    elif obs_types and len(obs_types) < 20:
        toks = line[6:60].split()
        obs_types.extend([t.strip() for t in toks if t.strip()])

print(f"  Observation types: {obs_types}")

# Map RINEX obs codes to our field names
# Standard RINEX 2.11 GPS obs: C1, P1, P2, L1, L2, D1, D2, S1, S2
# SWARM may use different codes based on firmware version
obs_map = {}
for idx, name in enumerate(obs_types):
    name_upper = name.upper()
    if name_upper in ('C1', 'P1', 'CA', 'C1C', 'C1W'):
        obs_map['P1'] = idx
    if name_upper in ('P2', 'C2', 'C2W'):
        obs_map['P2'] = idx
    if name_upper in ('L1', 'L1C', 'L1W'):
        obs_map['L1'] = idx
    if name_upper in ('L2', 'L2C', 'L2W'):
        obs_map['L2'] = idx
    if name_upper in ('S1', 'SN1', 'S1C'):
        obs_map['S1'] = idx
    if name_upper in ('S2', 'SN2', 'S2C'):
        obs_map['S2'] = idx

required = ['P1', 'P2', 'L1', 'L2']
missing = [k for k in required if k not in obs_map]
if missing:
    print(f"  WARNING: Missing obs types: {missing}")
    print(f"  Available types: {obs_types}")
    print("  Trying to auto-map legacy RINEX codes...")
    # Try alternate mapping for older GPS receivers
    for idx, name in enumerate(obs_types):
        if name.upper() == 'C1' and 'P1' not in obs_map:
            obs_map['P1'] = idx; missing.remove('P1') if 'P1' in missing else None
        if name.upper() == 'C2' and 'P2' not in obs_map:
            obs_map['P2'] = idx; missing.remove('P2') if 'P2' in missing else None
    if missing:
        print(f"  Still missing: {missing}")
        print(f"  Observation types: {obs_types}")
        sys.exit(1)

print(f"  Mapped: P1={obs_map.get('P1')} P2={obs_map.get('P2')} L1={obs_map.get('L1')} L2={obs_map.get('L2')}")

# ── Parse body ──
gps_obs = {}
epoch_count = 0
i = header_end + 1

# Find GPS week from first epoch
gps_week = None
sample_line = None
while i < len(lines) and gps_week is None:
    line = lines[i].strip()
    if line and len(line) > 25 and not line.startswith('COMMENT'):
        sample_line = line
        i += 1; continue
    i += 1

# Estimate GPS week from date
days_since = (dt - GPS_ORIGIN).days
gps_week = days_since // 7
sow_base = (days_since % 7) * 86400

i = header_end + 1
while i < len(lines):
    line = lines[i].rstrip()
    if len(line) < 25 or line.startswith('COMMENT'):
        i += 1; continue

    # Parse epoch header
    try:
        yr = int(line[0:3].strip() or '0')
        mo = int(line[3:6].strip() or '0')
        dy = int(line[6:9].strip() or '0')
        hr = int(line[9:12].strip() or '0')
        mn = int(line[12:15].strip() or '0')
        sc = float(line[15:26].strip() or '0')
        epoch_flag = int(line[26:29].strip() or '0')
        n_svs = int(line[29:32].strip() or '0')

        if yr < 80: yr += 2000
        epoch_dt = datetime(yr, mo, dy, hr, mn, int(sc))
        epoch_sec = (epoch_dt - timedelta(days=epoch_dt.day-1)).total_seconds()

        # GPS time: seconds of week
        epoch_gps_sec = gps_week * 604800 + sow_base + hr * 3600 + mn * 60 + sc
        epoch_gps_sec = int(epoch_gps_sec)

        # Read SV PRNs in the line (up to 12 per line, 3 chars each)
        sv_line = line[32:68] if len(line) > 68 else line[32:]
        svs = []
        for j in range(0, min(len(sv_line), 36), 3):
            prn = sv_line[j:j+3].strip()
            if prn and prn.startswith('G'):
                svs.append(prn)
            elif prn and prn.isdigit():
                svs.append(f'G{int(prn):02d}')

        # Advance i
        svs_read = svs[:n_svs]
        if not svs_read:
            i += 1; continue

        num_data_lines = (n_svs * len(obs_types) + 4) // 5  # 5 values per line

        for sv_idx, sv in enumerate(svs_read[:n_svs]):
            # Collect values for this SV across the data block
            vals = []
            for dl in range(num_data_lines):
                if i + 1 + dl < len(lines):
                    data_line = lines[i + 1 + dl]
                    # Extract 5 fields of 16 chars each (RINEX 2.11 F14.3 format)
                    for vp in range(0, min(len(data_line), 80), 16):
                        try:
                            v = float(data_line[vp:vp+16].strip())
                            vals.append(v)
                        except: pass

            if len(vals) < len(obs_types):
                i += num_data_lines + 1; continue  # not enough data

            # Extract our fields
            p1_val = vals[obs_map['P1']] if 'P1' in obs_map else 0.0
            p2_val = vals[obs_map['P2']] if 'P2' in obs_map else 0.0
            l1_val = vals[obs_map['L1']] if 'L1' in obs_map else 0.0
            l2_val = vals[obs_map['L2']] if 'L2' in obs_map else 0.0
            s1_val = vals[obs_map['S1']] if 'S1' in obs_map else 999.0
            s2_val = vals[obs_map['S2']] if 'S2' in obs_map else 999.0

            # Skip if no valid data
            if abs(p1_val) < 1e5 or abs(l1_val) < 1e5:
                continue

            # Convert L1/L2 from cycles to meters (RINEX stores phase in cycles)
            # Actually SWARM RINEX phase is in cycles
            l1_m = float(l1_val) * C_LIGHT / F1
            l2_m = float(l2_val) * C_LIGHT / F2

            # Compute IF combinations
            # L_if = (gamma * L2 - L1) / (gamma - 1)
            p1_m = float(p1_val)
            p2_m = float(p2_val)
            l_if = (GAMMA * l2_m - l1_m) / (GAMMA - 1)
            p_if = (GAMMA * p2_m - p1_m) / (GAMMA - 1)

            # Store
            if epoch_gps_sec not in gps_obs:
                gps_obs[epoch_gps_sec] = {}
            gps_obs[epoch_gps_sec][sv] = {
                'L1': l1_m, 'L2': l2_m,
                'P1': p1_m, 'P2': p2_m,
                'L_if': float(l_if),
                'P_if': float(p_if),
                'L1_cyc': float(l1_val),
                'L2_cyc': float(l2_val),
                'L1_SNR': float(s1_val),
                'L2_SNR': float(s2_val),
            }

        i += num_data_lines + 1
        epoch_count += 1
    except (ValueError, IndexError, KeyError) as e:
        i += 1; continue

print(f"  Parsed {epoch_count} epochs, {sum(len(v) for v in gps_obs.values())} total observations")

# ── Save as pickle ──
out_path = SWARM_DIR / f'GPS1B_{args.date}_SWARM_{args.sat}.pkl'
with open(str(out_path), 'wb') as f:
    pickle.dump(gps_obs, f)

print(f"  Saved: {out_path} ({out_path.stat().st_size/1e6:.1f}MB)")
print()

# ── Quick data check ──
epochs = sorted(gps_obs.keys())
if epochs:
    e0 = epochs[0]
    sv_list = sorted(gps_obs[e0].keys())
    print(f"  First epoch: GPS {e0}, {len(sv_list)} SVs: {sv_list}")
    sample = gps_obs[e0][sv_list[0]]
    print(f"  Sample {sv_list[0]}: P1={sample['P1']:.1f} L1={sample['L1']:.1f} "
          f"L_if={sample['L_if']:.1f} SNR={sample['L1_SNR']:.0f}")
print(f"\nSWARM-{args.sat} GPS data ready for POD processing.")
