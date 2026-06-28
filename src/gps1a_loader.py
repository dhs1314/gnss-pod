"""
GRACE-FO GPS Level-1A ASCII data loader

Data source: ISDC GFZ
URL: https://isdc-data.gfz.de/grace-fo/Level-1A/JPL/INSTRUMENT/RL04/{year}/
File: gracefo_1A_{date}_RL04.ascii.noLRI.tgz → GPS1A_{date}_C_04.txt

Key differences from GPS1B:
  - Phase at 1 Hz, pseudorange at 0.1 Hz (only every 10th epoch has valid P1/P2)
  - Phase preserved with integer ambiguity (negative values, internal receiver ref)
  - qualflg has explicit cycle-slip and phase-break bits
  - Phase units: meters (same as GPS1B), but NOT code-aligned
"""
from pathlib import Path
from datetime import datetime, timedelta
import numpy as np
import ssl, urllib.request, tarfile, pickle

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
F1_SQ, F2_SQ = F1 * F1, F2 * F2
LAM1 = C / F1
LAM2 = C / F2
GPS_ORIGIN = datetime(2000, 1, 1, 12, 0, 0)

N_HEADER = 7  # rcvtime_intg, rcvtime_frac, GRACEFO_id, prn_id, ant_id, prod_flag, qualflg

# Product fields in order (same as GPS1B)
PROD_FIELDS = [
    'CA_range', 'L1_range', 'L2_range', 'CA_phase', 'L1_phase', 'L2_phase',
    'CA_SNR', 'L1_SNR', 'L2_SNR', 'CA_chan', 'L1_chan', 'L2_chan',
    'L2_raw', 'Ka_phase', 'K_SNR', 'Ka_SNR',
]

ISDC_TGZ_URL = (
    "https://isdc-data.gfz.de/grace-fo/"
    "Level-1A/JPL/INSTRUMENT/RL04/{year}/"
    "gracefo_1A_{date}_RL04.ascii.noLRI.tgz"
)

# GPS time → UTC
_LEAP_TABLE = [
    (datetime(2000, 1, 1), 32), (datetime(2006, 1, 1), 33),
    (datetime(2009, 1, 1), 34), (datetime(2012, 7, 1), 35),
    (datetime(2015, 7, 1), 36), (datetime(2017, 1, 1), 37),
    (datetime(2024, 1, 1), 18),
]

def gps_sod_to_utc(gps_sod):
    gps_dt = GPS_ORIGIN + timedelta(seconds=gps_sod)
    for d, leap in reversed(_LEAP_TABLE):
        if gps_dt >= d:
            return gps_dt - timedelta(seconds=leap)
    return gps_dt


def parse_prod_flag(s):
    return int(s.strip(), 2)


def parse_qual_flag(s):
    return int(s.strip(), 2)


def parse_gps1a_record(parts, grace_filter='C'):
    """Parse single GPS1A record. Returns dict or None."""
    if len(parts) < N_HEADER:
        return None

    try:
        gps_sod = int(parts[0]) + int(parts[1]) * 1e-6
        grace_id = parts[2].strip()
        prn = int(parts[3])
        prod_flag_int = parse_prod_flag(parts[5])
        qual_flag_int = parse_qual_flag(parts[6])
    except (ValueError, IndexError):
        return None

    if prn < 1 or prn > 32:
        return None
    if grace_id != grace_filter:
        return None

    prod_values = parts[N_HEADER:]
    rec = {
        'sv': f"G{prn:02d}",
        'gps_sod': gps_sod,
        'qualflg': qual_flag_int,
    }

    # Extract product fields
    for bit_pos, field_name in enumerate(PROD_FIELDS):
        if (prod_flag_int >> bit_pos) & 1:
            if bit_pos < len(prod_values):
                try:
                    rec[field_name] = float(prod_values[bit_pos])
                except ValueError:
                    rec[field_name] = None
            else:
                rec[field_name] = None

    # Need L1_phase and L2_phase (in meters, raw)
    L1_ph_m = rec.get('L1_phase')
    L2_ph_m = rec.get('L2_phase')
    if L1_ph_m is None or L2_ph_m is None:
        return None

    # Need pseudorange (L1_range or CA_range)
    P1_m = rec.get('L1_range') or rec.get('CA_range')
    P2_m = rec.get('L2_range')
    if P1_m is None or P2_m is None or P1_m == 0 or P2_m == 0:
        return None  # Pseudorange-only epochs (0.1 Hz data)

    # SNR filter
    snr1 = rec.get('L1_SNR')
    if snr1 is not None and snr1 < 5:
        return None

    # Convert phase from meters to cycles
    rec['L1_cyc'] = L1_ph_m / LAM1
    rec['L2_cyc'] = L2_ph_m / LAM2
    rec['P1'] = P1_m
    rec['P2'] = P2_m
    rec['P_if'] = (F1_SQ * P1_m - F2_SQ * P2_m) / (F1_SQ - F2_SQ)

    # Cycle-slip flags from qualflg
    rec['slip_L1'] = bool(qual_flag_int & 1) or bool(qual_flag_int & 4)
    rec['slip_L2'] = bool(qual_flag_int & 2) or bool(qual_flag_int & 8)

    return rec


def parse_gps1a_text(text, grace_filter='C'):
    """Parse GPS1A ASCII text, returning {gps_sod: {sv: rec}}."""
    lines = text.split('\n')
    header_end = 0
    for i, l in enumerate(lines):
        if l.strip() == '# End of YAML header':
            header_end = i + 1
            break

    gps_obs = {}
    n_total = n_valid = n_no_pr = 0

    for line in lines[header_end:]:
        if not line.strip() or line.startswith('#'):
            continue
        n_total += 1
        parts = line.split()
        rec = parse_gps1a_record(parts, grace_filter)
        if rec is None:
            n_no_pr += 1
            continue
        gps_sod = rec['gps_sod']
        sv = rec['sv']
        if gps_sod not in gps_obs:
            gps_obs[gps_sod] = {}
        gps_obs[gps_sod][sv] = rec
        n_valid += 1

    n_epochs = len(gps_obs)
    print(f"  GPS1A: {n_total} rows → {n_valid} valid ({n_no_pr} no-PR), "
          f"{n_epochs} epochs")
    return gps_obs


def download_gps1a(year, month, day, data_dir="./data", grace_filter='C'):
    """Download and parse GRACE-FO GPS1A data with .pkl cache."""
    date_str = f"{year:04d}-{month:02d}-{day:02d}"
    cache = Path(data_dir) / "gracefo" / str(year) / date_str / f"gps1a_{grace_filter}.pkl"
    cache.parent.mkdir(parents=True, exist_ok=True)

    if cache.exists():
        try:
            gps_obs = pickle.load(open(cache, 'rb'))
            n_obs = sum(len(v) for v in gps_obs.values())
            print(f"  [GPS1A cache] {date_str} ({grace_filter}): "
                  f"{len(gps_obs)} epochs, {n_obs} obs")
            return gps_obs
        except Exception:
            pass

    # Local tgz
    tgz_path = Path(data_dir) / f"gracefo_1A_{date_str}_RL04.ascii.noLRI.tgz"
    if not tgz_path.exists():
        tgz_path = Path(data_dir) / "gracefo" / str(year) / f"gracefo_1A_{date_str}_RL04.ascii.noLRI.tgz"

    if not tgz_path.exists():
        tgz_url = ISDC_TGZ_URL.format(year=year, date=date_str)
        print(f"  [GPS1A download] {tgz_url.split('/')[-1]}")
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        try:
            req = urllib.request.Request(tgz_url, headers={
                'User-Agent': 'Mozilla/5.0', 'Accept-Encoding': 'gzip'
            })
            with urllib.request.urlopen(req, timeout=600, context=ctx) as r:
                data = r.read()
            tgz_path.parent.mkdir(parents=True, exist_ok=True)
            tgz_path.write_bytes(data)
            print(f"  [done] {len(data)//1024//1024} MB")
        except Exception as e:
            print(f"  [download failed] {e}")
            return None

    # Extract GPS1A text
    try:
        tar = tarfile.open(tgz_path)
        names = tar.getnames()
        gps_member = [n for n in names if 'GPS1A' in n and n.endswith('.txt') and f'_{grace_filter}_' in n]
        if not gps_member:
            print(f"  [error] no GPS1A .txt for {grace_filter} in tgz")
            tar.close()
            return None
        gps_member = gps_member[0]
        print(f"  [GPS1A extract] {gps_member}")
        f = tar.extractfile(gps_member)
        text = f.read().decode('ascii', errors='replace')
        f.close()
        tar.close()
    except Exception as e:
        print(f"  [extract failed] {e}")
        return None

    gps_obs = parse_gps1a_text(text, grace_filter=grace_filter)
    if gps_obs:
        pickle.dump(gps_obs, open(cache, 'wb'))
        print(f"  [cached] {cache}")
    return gps_obs
