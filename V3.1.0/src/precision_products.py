"""Precision GNSS product loaders: RINEX CLK, IGS ANTEX, CODE DCB.

Provides access to CODE/IGS final products for cm-level POD.
"""
import os
import re
import gzip
import urllib.request
import ssl
from pathlib import Path
from datetime import datetime, timedelta
import numpy as np

C_LIGHT = 299792458.0
GPS_ORIGIN = datetime(1980, 1, 6)

# GPS L1/L2 frequencies for ANTEX
F1 = 1575.42e6
F2 = 1227.60e6


# ===========================================================================
# RINEX Clock File Reader (.clk)
# ===========================================================================

def _gps_week_seconds(dt):
    """Return (gps_week, seconds_of_week) for a datetime."""
    delta = dt - GPS_ORIGIN
    week = int(delta.total_seconds() // 604800)
    sow = delta.total_seconds() - week * 604800
    return week, sow


def read_rinex_clk(filepath):
    """Read a RINEX 3.04 clock file.

    CODE final clock files have 30s interval, one record per satellite per
    epoch. Format:
      AS G01 2024 04 29 00 00  0.000000  0   -0.123456789012e-04

    Returns:
        dict: sv_id → {'epochs': [datetime, ...], 'clk': [m, ...]}
              Clock values are in metres (c * dt).
    """
    data = {}  # sv → {'epochs': [], 'clk': []}
    current_epoch = None

    with open(filepath, 'r') as f:
        for line in f:
            line = line.rstrip()
            if not line or line.startswith(' ') or len(line) < 40:
                continue
            if line.startswith('AS'):
                # AS G01 2024 04 29 00 00  0.000000  0   -0.123456789012e-04
                parts = line.split()
                if len(parts) < 8:
                    continue
                sv = parts[1]
                try:
                    yr, mo, dy = int(parts[2]), int(parts[3]), int(parts[4])
                    hr, mi = int(parts[5]), int(parts[6])
                    sec = float(parts[7])
                    epoch = datetime(yr, mo, dy, hr, mi, int(sec),
                                     int((sec % 1) * 1e6))
                except (ValueError, IndexError):
                    continue

                clk_val = float(parts[8])  # seconds
                clk_m = clk_val * C_LIGHT  # metres

                if sv not in data:
                    data[sv] = {'epochs': [], 'clk': []}
                data[sv]['epochs'].append(epoch)
                data[sv]['clk'].append(clk_m)

    # Convert to numpy arrays for fast interpolation
    for sv in data:
        data[sv]['clk'] = np.array(data[sv]['clk'], dtype=float)
        data[sv]['epochs_sod'] = np.array(
            [(t - t.replace(hour=0, minute=0, second=0, microsecond=0)
              ).total_seconds() for t in data[sv]['epochs']],
            dtype=float
        )

    return data


def get_clock_from_rinex_clk(clk_data, sv, epoch_dt):
    """Get interpolated clock correction from RINEX CLK data.

    Uses linear interpolation between 30s clock epochs.

    Args:
        clk_data: dict from read_rinex_clk()
        sv: satellite ID (e.g. 'G01')
        epoch_dt: query epoch as datetime

    Returns:
        clock_correction in metres, or 0.0 if SV not found
    """
    if sv not in clk_data:
        return 0.0

    d = clk_data[sv]
    sod = (epoch_dt - epoch_dt.replace(hour=0, minute=0, second=0,
                                       microsecond=0)).total_seconds()

    epochs_sod = d['epochs_sod']
    clk = d['clk']

    # Find bracketing epochs
    idx = int(np.searchsorted(epochs_sod, sod))
    if idx == 0:
        return float(clk[0])
    if idx >= len(epochs_sod):
        return float(clk[-1])

    # Linear interpolation
    t0, t1 = epochs_sod[idx - 1], epochs_sod[idx]
    c0, c1 = clk[idx - 1], clk[idx]
    if t1 - t0 < 1e-9:
        return float(c0)
    frac = (sod - t0) / (t1 - t0)
    return float(c0 + frac * (c1 - c0))


# ===========================================================================
# IGS ANTEX File Reader (.atx)
# ===========================================================================

def _parse_antex_frequency_section(lines, start_idx):
    """Parse frequency-specific PCO and PCV values from ANTEX section."""
    result = {'pco': np.zeros(3), 'pcv_zenith': {}, 'pcv_azimuth': {}}

    i = start_idx
    while i < len(lines):
        line = lines[i]
        if line.strip().startswith('END OF ANTENNA') or 'TYPE / SERIAL NO' in line:
            break
        if 'ZEN1 / ZEN2 / DZEN' in line:
            try:
                parts = line.split()
                if len(parts) >= 3:
                    result['pcv_grid'] = {
                        'zen1': float(parts[0]),
                        'zen2': float(parts[1]),
                        'dzen': float(parts[2]),
                    }
            except (ValueError, IndexError):
                pass
            # Parse PCV values
            n_zen = len(result.get('pcv_grid', {}))
            if n_zen:
                result['pcv_zenith'] = []
                result['pcv_azimuth'] = []
                i += 1
                while i < len(lines):
                    line2 = lines[i].strip()
                    if line2.startswith('END OF ANTENNA') or 'TYPE / SERIAL NO' in line2:
                        i -= 1
                        break
                    if 'NOAZI' in line2:
                        parts = line2.split()
                        try:
                            result['pcv_noazi'] = [float(x) for x in parts[:-2]]
                        except (ValueError, IndexError):
                            pass
                    else:
                        # Azimuth-dependent values
                        parts = line2.split()
                        if len(parts) >= 2 and parts[-1].startswith('S'):
                            result['pcv_azimuth'].append(
                                [float(x) for x in parts[:-1]])
                        else:
                            result['pcv_zenith'].append(
                                [float(x) for x in parts[:-1]])
                    i += 1
        elif 'NORTH / EAST / UP' in line:
            p = line.split()
            try:
                result['pco'] = np.array(
                    [float(p[0]), float(p[1]), float(p[2])])
            except (ValueError, IndexError):
                pass
        i += 1

    return result


def read_antex(filepath):
    """Read an IGS ANTEX file (e.g., igs14.atx).

    Returns a nested dict:
        data['G'][prn]['L1'|'L2'] = {
            'pco': np.array([north, east, up]) [m],
            'pcv_noazi': list of PCV values [mm],
            'pcv_grid': {'zen1': ..., 'zen2': ..., 'dzen': ...},
        }

    ANTEX uses fixed-width columns:
      Cols  1-20: ANTENNA TYPE (e.g. "BLOCK IIA", "GLONASS-M")
      Cols 21-40: SERIAL NUMBER / SVN (e.g. "G01", "R01") — this IS the PRN
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    data = {}
    current_system = None
    current_prn = None
    current_freq = None

    i = 0
    while i < len(lines):
        line = lines[i]
        # Detect satellite system — use fixed-width columns
        if 'TYPE / SERIAL NO' in line:
            # Fixed-width: TYPE(1:20), SERIAL_NO(21:40)
            ant_type = line[0:20].strip()
            ant_sn = line[20:40].strip() if len(line) >= 40 else ''
            # Satellite PRN is in the SERIAL NO field (G01, R01, E01, C01, J01)
            prn = ant_sn
            if prn and len(prn) == 3 and prn[0] in 'GRECJI':
                system = prn[0]
                if system not in data:
                    data[system] = {}
                if prn not in data[system]:
                    data[system][prn] = {}
                current_system = system
                current_prn = prn

        # Detect frequency section
        elif 'START OF FREQUENCY' in line:
            # Frequency code in fixed-width columns (cols ~3-6)
            freq_code = line[0:6].strip()
            if freq_code and freq_code[0] == 'G':
                # GPS frequency codes: G01=L1, G02=L2, G05=L5
                if '01' in freq_code or '1' in freq_code[1:]:
                    current_freq = 'L1'
                elif '02' in freq_code or '2' in freq_code[1:]:
                    current_freq = 'L2'
                elif '05' in freq_code or '5' in freq_code[1:]:
                    current_freq = 'L5'
                if current_system and current_prn and current_freq:
                    freq_data = _parse_antex_frequency_section(lines, i + 1)
                    data[current_system][current_prn][current_freq] = freq_data

        i += 1

    return data


def get_satellite_pco(antex_data, prn, freq='L1'):
    """Extract satellite antenna PCO from ANTEX data.

    Returns (north, east, up) PCO in **metres** (converted from ANTEX mm).
    UP component is positive away from Earth centre (toward antenna boresight).
    Typical values: Block IIA ~1.5m, Block IIF ~0.8m.
    """
    try:
        pco_mm = antex_data['G'][prn][freq]['pco'].copy()
        return pco_mm / 1000.0  # mm → m
    except (KeyError, IndexError):
        return np.zeros(3)


def get_satellite_pcv(antex_data, prn, nadir_deg, freq='L1'):
    """Get satellite antenna PCV at a given nadir angle.

    Args:
        antex_data: parsed ANTEX data
        prn: satellite PRN (e.g. 'G05')
        nadir_deg: nadir angle in degrees (0 = boresight)
        freq: 'L1' or 'L2'

    Returns:
        PCV in millimetres
    """
    try:
        sat_data = antex_data['G'][prn][freq]
    except KeyError:
        return 0.0

    if 'pcv_noazi' not in sat_data or not sat_data['pcv_noazi']:
        return 0.0

    pcv_vals = sat_data['pcv_noazi']
    grid = sat_data.get('pcv_grid', {})
    zen1 = grid.get('zen1', 0.0)
    dzen = grid.get('dzen', 1.0)

    # Linear interpolation in zenith
    idx = (nadir_deg - zen1) / dzen
    i0 = int(np.floor(idx))
    i1 = i0 + 1
    if i0 < 0:
        return float(pcv_vals[0])
    if i1 >= len(pcv_vals):
        return float(pcv_vals[-1])
    frac = idx - i0
    return float(pcv_vals[i0] + frac * (pcv_vals[i1] - pcv_vals[i0]))


# ===========================================================================
# CODE DCB Reader
# ===========================================================================

def read_code_dcb(filepath):
    """Read CODE monthly DCB file.

    CODE monthly DCB format:
      Header: "CODE'S MONTHLY ... P1-C1 DCB SOLUTION ..."
      Then: "DIFFERENTIAL (P1-C1) CODE BIASES ..."
      Data:  PRN  VALUE_NS  RMS_NS
      G01                           0.620       0.213

    Auto-detects DCB type from header line.

    Returns:
        tuple: (dcb_dict, dcb_type_str)
        dcb_dict: prn → value_ns
        dcb_type: 'P1C1' or 'P1P2'
    """
    dcb = {}
    dcb_type = None

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Detect DCB type from header
            if dcb_type is None:
                if 'P1-C1' in line:
                    dcb_type = 'P1C1'
                elif 'P1-P2' in line:
                    dcb_type = 'P1P2'
                # Skip header lines without data
                if '--------' in line or 'PRN / STATION' in line or '****' in line:
                    continue
                if line.startswith('CODE') or line.startswith('DIFFERENTIAL'):
                    continue

            parts = line.split()
            if len(parts) < 2:
                continue
            prn = parts[0]
            if not prn.startswith('G'):
                continue
            try:
                value_ns = float(parts[1])
            except ValueError:
                continue
            dcb[prn] = value_ns

    return dcb, dcb_type


def load_code_dcb_pair(p1c1_path, p1p2_path):
    """Load both CODE DCB files and return combined per-SV corrections.

    Args:
        p1c1_path: path to P1C1*.DCB file
        p1p2_path: path to P1P2*.DCB file

    Returns:
        dict: prn → {'P1C1_ns': float, 'P1P2_ns': float}
    """
    dcb = {}
    for path in (p1c1_path, p1p2_path):
        if path is None:
            continue
        p = Path(path)
        if not p.exists():
            continue
        data, dtype = read_code_dcb(str(p))
        for prn, val_ns in data.items():
            if prn not in dcb:
                dcb[prn] = {}
            dcb[prn][dtype + '_ns'] = val_ns
    return dcb


# IF combination coefficients for GPS L1/L2
_GPS_F1 = 1575.42e6
_GPS_F2 = 1227.60e6
_GPS_ALPHA = _GPS_F1**2 / (_GPS_F1**2 - _GPS_F2**2)  # ≈ 2.5457
_GPS_BETA  = _GPS_F2**2 / (_GPS_F1**2 - _GPS_F2**2)   # ≈ 1.5457


def compute_dcb_if_correction(dcb_pair, prn):
    """Compute iono-free DCB correction for a GPS satellite.

    GPS1B uses C1 (C/A code on L1) and P2 (Z-tracking on L2) for the
    iono-free combination. The true IF uses P1 and P2.  Since we already
    observe P2 directly, only the C1→P1 correction is needed:

        P_if_corrected = P_if_raw + α * DCB(P1-C1)

    where α = f1^2/(f1^2-f2^2) ≈ 2.5457 is the IF coefficient for L1.
    The P1-P2 DCB is NOT applied because P2 is already the P2 observable.

    Args:
        dcb_pair: dict from load_code_dcb_pair()
        prn: satellite ID (e.g. 'G01')

    Returns:
        correction in metres to add to P_if_raw
    """
    if prn not in dcb_pair:
        return 0.0
    p1c1_ns = dcb_pair[prn].get('P1C1_ns', 0.0)
    return _GPS_ALPHA * p1c1_ns * 1e-9 * C_LIGHT


# ===========================================================================
# CODE Product Downloader
# ===========================================================================

CODE_FTP = 'ftp://ftp.aiub.unibe.ch/CODE/'

# GPS week → YYYY directory mapping for CODE products
# CODE organizes by GPS week in YYYY/ subdirectories


def _download_file(url, dest_path, timeout=60):
    """Download a file with SSL context."""
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE

    req = urllib.request.Request(url, headers={
        'User-Agent': 'gnss-pod/2.3',
    })

    try:
        with urllib.request.urlopen(req, timeout=timeout, context=ctx) as r:
            data = r.read()
        with open(dest_path, 'wb') as f:
            f.write(data)
        return True
    except Exception as e:
        print(f"  [Download] Failed: {url} → {e}")
        return False


def download_code_products(year, doy, data_dir='./data/CODE',
                           products=('SP3', 'CLK')):
    """Download CODE final products for a given date.

    CODE final products are available ~12-14 days after the observation date.
    File naming:
      COD0OPSFIN_YYYYDDDHHMM_01D_05M_ORB.SP3.gz   (5-min orbit)
      COD0OPSFIN_YYYYDDDHHMM_01D_30S_CLK.CLK.gz   (30s clock)

    Note: CODE uses the same YYYYDDD naming as IGS, where DDD = GPS DOY.

    Args:
        year: GPS year (e.g. 2024)
        doy: day of year (1-366)
        data_dir: local directory to store downloads
        products: list of product types ('SP3', 'CLK')

    Returns:
        dict: product_type → local_filepath or None
    """
    gps_week, _ = _gps_week_seconds(datetime(year, 1, 1) + timedelta(days=doy - 1))
    week_dir = str(gps_week)

    dest_dir = Path(data_dir) / str(year)
    dest_dir.mkdir(parents=True, exist_ok=True)

    gps_week_str = f'{gps_week:04d}'
    doy_str = f'{doy:01d}'

    result = {}

    # CODE SP3: COD0OPSFIN_YYYYDDDHHMM_01D_05M_ORB.SP3.gz
    if 'SP3' in products:
        # DDD in CODE filename is day-of-week (0-6) concatenated with GPS week?
        # Actually CODE uses YYYY + DDD where DDD is day-of-year
        # File: COD0OPSFIN_20241200000_01D_05M_ORB.SP3.gz
        sp3_name = f'COD0OPSFIN_{year}{doy:03d}0000_01D_05M_ORB.SP3.gz'
        sp3_url = f'{CODE_FTP}{year}/COD0OPSFIN_{year}{doy:03d}0000_01D_05M_ORB.SP3.gz'
        sp3_path = dest_dir / sp3_name
        if sp3_path.exists():
            result['SP3'] = str(sp3_path)
        else:
            print(f"  Downloading CODE SP3: {sp3_name}...")
            if _download_file(sp3_url, str(sp3_path)):
                result['SP3'] = str(sp3_path)

    # CODE CLK: COD0OPSFIN_YYYYDDDHHMM_01D_30S_CLK.CLK.gz
    if 'CLK' in products:
        clk_name = f'COD0OPSFIN_{year}{doy:03d}0000_01D_30S_CLK.CLK.gz'
        clk_url = f'{CODE_FTP}{year}/{clk_name}'
        clk_path = dest_dir / clk_name
        if clk_path.exists():
            result['CLK'] = str(clk_path)
        else:
            print(f"  Downloading CODE CLK: {clk_name}...")
            if _download_file(clk_url, str(clk_path)):
                result['CLK'] = str(clk_path)

    return result


# ===========================================================================
# IERS C04 ERP Reader (Earth Orientation Parameters)
# ===========================================================================

def read_iers_c04(filepath):
    """Read IERS C04 Earth Orientation Parameters file (old format).

    Format:
        Date      MJD      x          y        UT1-UTC       LOD         dX        dY
                   "          "           s           s          "         "
       (0h UTC)
    1962   1   1  37665  -0.012700   0.213000   0.0326338   0.0017230   0.000000   0.000000

    Returns:
        dict: {
            'mjd': np.array (N,),
            'x_arcsec': np.array (N,),     # polar motion x [arcsec]
            'y_arcsec': np.array (N,),     # polar motion y [arcsec]
            'ut1_utc': np.array (N,),      # UT1-UTC [s]
            'lod': np.array (N,),          # length of day excess [s]
            'dx_arcsec': np.array (N,),    # celestial pole offset dX [arcsec]
            'dy_arcsec': np.array (N,),    # celestial pole offset dY [arcsec]
        }
    """
    mjd_list, x_list, y_list, ut1_list, lod_list, dx_list, dy_list = [], [], [], [], [], [], []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.rstrip()
            if not line or line[0] in ('#', ' ', '\n'):
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            try:
                yr, mo, dy = int(parts[0]), int(parts[1]), int(parts[2])
                mjd = int(parts[3])
                x = float(parts[4])
                y = float(parts[5])
                ut1 = float(parts[6])
                lod = float(parts[7]) if len(parts) > 7 else 0.0
                dx = float(parts[8]) if len(parts) > 8 else 0.0
                dy2 = float(parts[9]) if len(parts) > 9 else 0.0
            except (ValueError, IndexError):
                continue

            mjd_list.append(mjd)
            x_list.append(x)
            y_list.append(y)
            ut1_list.append(ut1)
            lod_list.append(lod)
            dx_list.append(dx)
            dy_list.append(dy2)

    return {
        'mjd': np.array(mjd_list, dtype=float),
        'x_arcsec': np.array(x_list, dtype=float),
        'y_arcsec': np.array(y_list, dtype=float),
        'ut1_utc': np.array(ut1_list, dtype=float),
        'lod': np.array(lod_list, dtype=float),
        'dx_arcsec': np.array(dx_list, dtype=float),
        'dy_arcsec': np.array(dy_list, dtype=float),
    }


def get_eop_from_c04(c04_data, mjd):
    """Linearly interpolate EOP from C04 data at given MJD.

    Args:
        c04_data: dict from read_iers_c04()
        mjd: Modified Julian Date (UTC)

    Returns:
        (x_pole_rad, y_pole_rad, ut1_utc_s)
    """
    mjd_arr = c04_data['mjd']
    idx = np.searchsorted(mjd_arr, mjd)
    if idx <= 0:
        x = c04_data['x_arcsec'][0]
        y = c04_data['y_arcsec'][0]
        ut1 = c04_data['ut1_utc'][0]
    elif idx >= len(mjd_arr):
        x = c04_data['x_arcsec'][-1]
        y = c04_data['y_arcsec'][-1]
        ut1 = c04_data['ut1_utc'][-1]
    else:
        mjd0, mjd1 = mjd_arr[idx - 1], mjd_arr[idx]
        frac = (mjd - mjd0) / (mjd1 - mjd0)
        x = c04_data['x_arcsec'][idx - 1] + frac * (c04_data['x_arcsec'][idx] - c04_data['x_arcsec'][idx - 1])
        y = c04_data['y_arcsec'][idx - 1] + frac * (c04_data['y_arcsec'][idx] - c04_data['y_arcsec'][idx - 1])
        ut1 = c04_data['ut1_utc'][idx - 1] + frac * (c04_data['ut1_utc'][idx] - c04_data['ut1_utc'][idx - 1])

    ARC_SEC_TO_RAD = np.pi / (180.0 * 3600.0)
    return x * ARC_SEC_TO_RAD, y * ARC_SEC_TO_RAD, ut1


def setup_iers_from_c04(c04_path):
    """Load IERS C04 data and replace astropy's default EOP table.

    Creates a compatible IERS table from C04 data so astropy's
    ITRS <-> GCRS transformations use precise Earth orientation.

    Args:
        c04_path: path to eopc04_IAU2000.txt

    Returns:
        True on success, False on failure
    """
    try:
        from astropy.utils.iers import IERS_Auto
        from astropy.table import QTable
        import astropy.units as u

        c04 = read_iers_c04(c04_path)
        n = len(c04['mjd'])

        ARC_SEC_TO_RAD = np.pi / (180.0 * 3600.0)

        table = QTable()
        table['MJD'] = c04['mjd']
        table['UT1_UTC'] = c04['ut1_utc'] * u.s
        table['PM_x'] = c04['x_arcsec'] * ARC_SEC_TO_RAD * u.rad
        table['PM_y'] = c04['y_arcsec'] * ARC_SEC_TO_RAD * u.rad
        table['dX_2000A'] = c04['dx_arcsec'] * ARC_SEC_TO_RAD * u.rad
        table['dY_2000A'] = c04['dy_arcsec'] * ARC_SEC_TO_RAD * u.rad
        table['UT1Flag'] = np.full(n, 'B')
        table['PolarFlag'] = np.full(n, 'B')
        table['NutFlag'] = np.full(n, 'B')
        table['LOD'] = c04['lod'] * u.s
        table['e_UT1_UTC'] = np.full(n, 1e-6) * u.s
        table['e_PM_x'] = np.full(n, 1e-9) * u.rad
        table['e_PM_y'] = np.full(n, 1e-9) * u.rad
        table['e_dX_2000A'] = np.full(n, 1e-9) * u.rad
        table['e_dY_2000A'] = np.full(n, 1e-9) * u.rad
        table['Year'] = np.zeros(n, dtype=int)

        IERS_Auto.iers_table = table
        print(f"  [IERS] C04 loaded: {n} records, MJD {c04['mjd'][0]:.0f} to {c04['mjd'][-1]:.0f}")
        return True
    except Exception as e:
        print(f"  [IERS] WARNING: Could not load C04 ERP: {e}")
        return False
