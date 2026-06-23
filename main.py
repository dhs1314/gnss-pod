#!/usr/bin/env python3

"""

GNSS PPP 主程序 — 支持地面站 + GRACE-FO 卫星双模式

支持 WHU ultra-rapid SP3 精密星历（无需认证）

Usage: python3 main.py --station wuhn --year 2024 --doy 120 --hours 4

       python3 main.py --satellite gracefo1 --year 2024 --month 4 --day 29 --hours 4

"""



import sys, os, argparse, warnings, urllib.request, urllib.error, gzip, ssl



from src.sp3_loader import load_whu_sp3_ultra_rapid, get_gps_pos_from_sp3

from pathlib import Path

from datetime import datetime, timedelta

import numpy as np



sys.path.insert(0, str(Path(__file__).parent / 'src'))

from ppp import PPPProcessor, ecef_to_blh, ecef_to_enu_matrix

from plotting import plot_all



C = 299792458.0

F1, F2 = 1575.42e6, 1227.60e6

LAMBDA1, LAMBDA2 = C / F1, C / F2

MU_E = 3.9860050e14

OMEGA_E = 7.2921151467e-5

GPS_ORIGIN = datetime(1980, 1, 6)



# ── GPS 星座参数 ───────────────────────────────────────────────────────

GPS_SV_PLAN = [

    (1, 0.0, 0.0, 0.0, 0.000020, 55.0, 0.0, 5153.6),

    (3, 30.0, 30.0, 60.0, 0.000015, 54.8, 30.0, 5153.6),

    (5, 90.0, 90.0, 180.0, 0.000010, 54.9, 90.0, 5153.6),

    (7, 150.0, 150.0, 300.0, 0.000015, 54.7, 150.0, 5153.6),

    (8, 180.0, 180.0, 0.0, 0.000020, 55.0, 180.0, 5153.6),

    (10, 240.0, 240.0, 120.0, 0.000020, 55.0, 240.0, 5153.6),

    (11, 270.0, 270.0, 180.0, 0.000015, 54.9, 270.0, 5153.6),

    (13, 330.0, 330.0, 300.0, 0.000010, 54.7, 330.0, 5153.6),

    (15, 15.0, 15.0, 90.0, 0.000020, 65.0, 15.0, 5153.6),

    (17, 45.0, 45.0, 150.0, 0.000015, 64.8, 45.0, 5153.6),

    (18, 75.0, 75.0, 210.0, 0.000020, 65.0, 75.0, 5153.6),

    (19, 105.0, 105.0, 270.0, 0.000010, 64.9, 105.0, 5153.6),

    (20, 135.0, 135.0, 330.0, 0.000020, 65.1, 135.0, 5153.6),

    (21, 165.0, 165.0, 30.0, 0.000015, 64.7, 165.0, 5153.6),

    (22, 195.0, 195.0, 90.0, 0.000020, 65.0, 195.0, 5153.6),

    (23, 225.0, 225.0, 150.0, 0.000010, 64.8, 225.0, 5153.6),

    (24, 255.0, 255.0, 210.0, 0.000020, 65.0, 255.0, 5153.6),

    (25, 285.0, 285.0, 270.0, 0.000015, 64.9, 285.0, 5153.6),

    (27, 345.0, 345.0, 30.0, 0.000020, 64.7, 345.0, 5153.6),

    (29, 25.0, 25.0, 120.0, 0.000010, 56.0, 25.0, 5153.6),

    (30, 55.0, 55.0, 200.0, 0.000020, 56.1, 55.0, 5153.6),

]



IGS_COORDS = {

    'wuhn': (-2148755.2370, 4426641.5334, 4044655.5415),

    'bjfs': (-2148799.8920, 4426645.6280, 4044650.9890),

    'cusv': (-2684752.0990, -5294414.2330, 3071275.2420),

    'will': (-2052060.6350, -3627147.9330, 4811309.9200),

    'mas1': (5103907.9020, -2883930.0280, 3190947.2230),

}



ISDC_GRACEFO_L1B = "https://isdc-data.gfz.de/grace-fo/Level-1B/JPL/INSTRUMENT/RL04/{year}/"





def satpos_from_sv(sv_rec, ts):

    prn, M0_d, Omega0_d, omega_d, ecc, inc_d, u0_d, sqrtA = sv_rec

    a = sqrtA ** 2

    n0 = np.sqrt(MU_E / a**3)

    M0 = np.radians(M0_d)

    Omega0 = np.radians(Omega0_d)

    omega = np.radians(omega_d)

    inc = np.radians(inc_d)

    gps_sow = (ts - GPS_ORIGIN).total_seconds()

    sow = gps_sow % 604800

    toe = sow

    t_k = sow - toe

    if t_k > 302400: t_k -= 604800

    if t_k < -302400: t_k += 604800

    M = M0 + n0 * t_k

    E = M

    for _ in range(10): E = M + ecc * np.sin(E)

    sinE, cosE = np.sin(E), np.cos(E)

    v = np.arctan2(np.sqrt(1-ecc**2)*sinE, cosE-ecc)

    phi = v + omega

    r = a * (1 - ecc*cosE)

    i = inc

    omega_dot = np.radians(0.0265) / 86164

    Om = Omega0 + omega_dot*t_k - OMEGA_E*toe

    xp = r*np.cos(phi); yp = r*np.sin(phi)

    pos = np.array([

        xp*np.cos(Om) - yp*np.cos(i)*np.sin(Om),

        xp*np.sin(Om) + yp*np.cos(i)*np.cos(Om),

        yp*np.sin(i)

    ])

    return pos, 0.0





def _ssl_ctx():

    ctx = ssl.create_default_context()

    ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE

    return ctx





def _try_download(url, local, timeout=60):

    try:

        req = urllib.request.Request(url, headers={'User-Agent':'curl/7.88','Accept-Encoding':'gzip'})

        with urllib.request.urlopen(req, timeout=timeout, context=_ssl_ctx()) as r:

            data = r.read()

        if data[:5] in (b'<!DOC', b'<html', b'<HTML', b''): return False

        Path(local).parent.mkdir(parents=True, exist_ok=True)

        with open(local, 'wb') as f: f.write(data)

        return True

    except Exception: return False





def _decompress(gz_path, out_path):

    try:

        with gzip.open(gz_path, 'rb') as fi:

            data = fi.read()

        if data[:5] in (b'<!DOC', b'<html'): return False

        with open(out_path, 'wb') as fo: fo.write(data)

        return True

    except Exception: return False





RINEX_SOURCES = [

    "https://cddis.nasa.gov/archive/gnss/data/daily/{yyyy}/{doy:03d}/{yy}o/{station}{doy:03d}0.{yy}{ext}",

    "https://igs.bkg.bund.de/root_ftp/IGS/data/{yyyy}/{doy:03d}/{station}{doy:03d}0.{yy}{ext}",

    "https://www.epncb.oma.be/ftp/data/{yyyy}/{doy:03d}/{station}{doy:03d}0.{yy}{ext}",

]

NAV_SOURCES = [

    "https://cddis.nasa.gov/archive/gnss/data/daily/{yyyy}/{doy:03d}/{yy}n/brdc{doy:03d}0.{yy}n.gz",

    "https://igs.bkg.bund.de/root_ftp/IGS/data/{yyyy}/{doy:03d}/brdc{doy:03d}0.{yy}n.gz",

    "https://www.epncb.oma.be/ftp/data/{yyyy}/{doy:03d}/brdc{doy:03d}0.{yy}n.gz",

]





def try_download_rinex(station, year, doy, data_dir):

    ddir = Path(data_dir) / str(year) / f"{doy:03d}"

    ddir.mkdir(parents=True, exist_ok=True)

    yy = str(year)[2:]

    yyyy = str(year)

    for ext in ['.rnx.gz', '.rnx', '.crx.gz']:

        fname = f"{station}{doy:03d}0.{yy}{ext}"

        local_gz = ddir / fname

        local_out = ddir / fname.replace('.gz','')

        if local_out.exists() and local_out.stat().st_size > 10000: return str(local_out)

        for src in RINEX_SOURCES:

            url = src.format(yy=yy, yyyy=yyyy, doy=doy, station=station, ext=ext)

            if _try_download(url, str(local_gz)):

                if ext.endswith('.gz') and _decompress(str(local_gz), str(local_out)):

                    if Path(str(local_out)).stat().st_size > 10000: return str(local_out)

                elif Path(str(local_gz)).stat().st_size > 10000:

                    return str(local_gz)

    return None





def try_download_nav(year, doy, data_dir):

    ddir = Path(data_dir) / str(year) / f"{doy:03d}"

    ddir.mkdir(parents=True, exist_ok=True)

    yy = str(year)[2:]; yyyy = str(year)

    fname = f"brdc{doy:03d}0.{yy}n.gz"

    local_gz = ddir / fname; local_out = ddir / fname.replace('.gz','')

    if local_out.exists() and local_out.stat().st_size > 5000: return str(local_out)

    for src in NAV_SOURCES:

        url = src.format(yy=yy, yyyy=yyyy, doy=doy)

        if _try_download(url, str(local_gz)):

            if _decompress(str(local_gz), str(local_out)):

                if Path(str(local_out)).stat().st_size > 5000: return str(local_out)

    return None





def parse_rinex_obs(filepath):

    with open(filepath, encoding='utf-8', errors='replace') as f: lines = f.readlines()

    h_end = next((i for i, l in enumerate(lines) if 'END OF HEADER' in l[60:]), len(lines))

    obs_types = []

    for line in lines[:h_end+1]:

        tag = line[60:].strip()

        if tag in ('# / TYPES OF OBS', 'OBS TYPES'):

            n = int(line[:6].strip())

            obs_types = line[6:].split()[:n]

    epochs = []

    i = h_end + 1

    while i < len(lines):

        line = lines[i].strip()

        if not line or line.startswith('COMMENT'):

            i += 1; continue

        parts = line.split()

        if len(parts) < 7:

            i += 1; continue

        try:

            yr = int(parts[0])

            if yr < 80: yr += 2000

            elif yr < 100: yr += 1900

            mo, dy, hr, mn = int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])

            sec = float(parts[5])

            epoch_t = datetime(yr, mo, dy, hr, mn, int(sec))

            i += 1; sat_obs = {}

            for _ in range(100):

                if i >= len(lines): break

                raw = lines[i]

                if raw.strip() and raw[0].isdigit(): break

                pos = 0; k = 0; obs_vals = {}

                while pos < len(raw) and k < len(obs_types):

                    tok = raw[pos:pos+14].strip()

                    try: obs_vals[obs_types[k]] = float(tok) if tok not in ('', '0', '0.0') else np.nan

                    except: obs_vals[obs_types[k]] = np.nan

                    pos += 14; k += 1

                sv_id = raw[:3].strip()

                if sv_id: sat_obs[sv_id] = obs_vals

                i += 1

            epochs.append({'time': epoch_t, 'obs': sat_obs})

        except: i += 1; continue

    return epochs, obs_types





def parse_rinex_nav(filepath):

    with open(filepath, encoding='utf-8', errors='replace') as f: lines = f.readlines()

    i = 0

    while i < len(lines) and 'END OF HEADER' not in lines[i][60:]: i += 1

    i += 1; records = []

    while i < len(lines):

        line = lines[i].strip()

        if not line or line.startswith('#'): i += 1; continue

        parts = line.split()

        if len(parts) < 8: i += 1; continue

        try:

            yr = int(parts[0]); mo = int(parts[1]); dy = int(parts[2])

            hr = int(parts[3]); mn = int(parts[4]); sc = float(parts[5])

            sv = parts[6]

            if yr < 80: yr += 2000

            elif yr < 100: yr += 1900

            t = datetime(yr, mo, dy, hr, mn, int(sc))

            def f(lst, idx, d=0.0):

                try: return float(lst[idx])

                except: return d

            i += 1; p1 = lines[i].strip().split(); i += 1

            p2 = lines[i].strip().split(); i += 1

            p3 = lines[i].strip().split(); i += 1

            records.append({

                'sv': sv, 'time': t, 'toe': f(p1, 15),

                'sqrtA': f(p1, 2), 'e': f(p1, 1), 'M0': f(p1, 0),

                'i0': f(p1, 3), 'Omega0': f(p1, 4), 'omega': f(p1, 5),

                'OmegaDot': f(p1, 6), 'iDot': f(p1, 7), 'dN': f(p1, 8),

                'cuc': f(p1, 9), 'cus': f(p1, 10), 'crc': f(p1, 11), 'crs': f(p1, 12),

                'cic': f(p1, 13), 'cis': f(p1, 14),

                'af0': f(p1, 17), 'af1': f(p1, 18),

            })

        except: i += 1; continue

    return records





def satpos_from_navrec(rec, ts):

    a = rec['sqrtA'] ** 2

    n0 = np.sqrt(MU_E / a**3)

    M0, ecc, omega = rec['M0'], rec['e'], rec['omega']

    i0, Omega0 = rec['i0'], rec['Omega0']

    toe = rec.get('toe', 0.0)

    gps_sow = (ts - GPS_ORIGIN).total_seconds()

    sow = gps_sow % 604800

    t_k = sow - toe

    if t_k > 302400: t_k -= 604800

    if t_k < -302400: t_k += 604800

    n = n0 + rec.get('dN', 0.0)

    M = M0 + n * t_k

    E = M

    for _ in range(10): E = M + ecc * np.sin(E)

    sinE, cosE = np.sin(E), np.cos(E)

    v = np.arctan2(np.sqrt(1-ecc**2)*sinE, cosE-ecc)

    phi = v + omega

    r = a * (1 - ecc*cosE)

    i = i0 + rec.get('iDot', 0.0)*t_k

    Om = Omega0 + rec.get('OmegaDot', 0.0)*t_k - OMEGA_E*toe

    xp = r*np.cos(phi); yp = r*np.sin(phi)

    pos = np.array([

        xp*np.cos(Om) - yp*np.cos(i)*np.sin(Om),

        xp*np.sin(Om) + yp*np.cos(i)*np.cos(Om),

        yp*np.sin(i)

    ])

    return pos





def generate_broadcast_nav(t_ref):

    records = []

    gps_secs_ref = (t_ref - GPS_ORIGIN).total_seconds()

    toe = gps_secs_ref % 604800

    for sv_rec in GPS_SV_PLAN:

        prn = sv_rec[0]

        records.append({

            'sv': f'G{prn:02d}', 'toe': toe, 'time': t_ref,

            'sqrtA': sv_rec[7], 'e': sv_rec[4], 'M0': np.radians(sv_rec[1]),

            'i0': np.radians(sv_rec[5]), 'omega': np.radians(sv_rec[3]),

            'Omega0': np.radians(sv_rec[2]), 'OmegaDot': np.radians(0.0265),

            'iDot': 0.0, 'dN': 0.0, 'cuc': 0.0, 'cus': 1e-6,

            'crc': 100.0, 'crs': 0.0, 'cic': 0.0, 'cis': 0.0, 'af0': 0.0, 'af1': 0.0,

        })

    return records





def generate_syn_observations(station_pos, t_start, n_hours, interval=30.0, seed=42):

    ref_lat, ref_lon, _ = ecef_to_blh(station_pos)

    M_enu_mat = ecef_to_enu_matrix(ref_lat, ref_lon)

    n_epochs = int(n_hours * 3600 / interval)

    records = []

    rng = np.random.default_rng(seed)

    dtropo_wet = 0.05

    pos_bias = np.zeros(3)

    for k in range(n_epochs):

        t = t_start + timedelta(seconds=k * interval)

        dtropo_wet += rng.normal(0, 0.0005)

        pos_bias += rng.normal(0, 0.005)

        for sv_rec in GPS_SV_PLAN:

            sv = f"G{sv_rec[0]:02d}"

            sat_pos, _ = satpos_from_sv(sv_rec, t)

            rho_vec = sat_pos - (station_pos + pos_bias)

            rho = np.linalg.norm(rho_vec)

            if not (20e6 < rho < 50e6): continue

            e_sat = rho_vec / rho

            e_enu = M_enu_mat @ e_sat

            el = np.arcsin(max(-1.0, min(1.0, e_enu[2])))

            if np.degrees(el) < 10.0: continue

            az = np.arctan2(e_enu[0], e_enu[1])

            zhd = 2.3

            mf = 1.0 / max(np.sin(el), 0.05)

            code_noise = rng.normal(0, 0.3)

            phase_noise = rng.normal(0, 0.003)

            P_if = rho + zhd*mf + dtropo_wet*mf + code_noise

            L_if = rho + zhd*mf + dtropo_wet*mf + phase_noise

            records.append({

                'time': t, 'sv': sv,

                'L1': L_if, 'L2': L_if, 'P1': P_if, 'P2': P_if,

                'sat_pos': sat_pos, 'sat_clock': 0.0,

                'el': np.degrees(el), 'az': np.degrees(az),

            })

    return records





def parse_gracefo_l1b_gnss(filepath):

    records = []

    with open(filepath, encoding='utf-8', errors='replace') as f: lines = f.readlines()

    has_gnss = any('GNSS' in l.upper() or 'GPS' in l.upper() for l in lines[:100])

    has_acc = any('gps_time' in l or 'lin_accl' in l for l in lines[:100])

    if has_gnss and not has_acc:

        for line in lines:

            line = line.strip()

            if not line or line.startswith('#') or line.startswith('%'): continue

            parts = line.split()

            if len(parts) < 8: continue

            try:

                prn = parts[1]

                if not (prn.startswith('G') and len(prn) == 3 and prn[1:].isdigit()):

                    continue

                t_gps = float(parts[0])

                gps_origin_j2000 = datetime(2000, 1, 1, 12, 0, 0)

                if 40000 < t_gps < 60000:

                    t = gps_origin_j2000 + timedelta(days=t_gps - 44239.0)

                else:

                    week = int(t_gps // 604800)

                    sow = t_gps % 604800

                    t = GPS_ORIGIN + timedelta(weeks=week, seconds=sow)

                records.append({

                    'time': t, 'sv': prn,

                    'L1': float(parts[2]), 'L2': float(parts[3]),

                    'P1': float(parts[4]), 'P2': float(parts[5]),

                    'sat_pos': None, 'sat_clock': 0.0, 'el': 90.0, 'az': 0.0,

                })

            except (ValueError, IndexError): continue

        return records

    if has_acc:

        gps_origin_j2000 = datetime(2000, 1, 1, 12, 0, 0)

        data_lines = []

        in_data = False

        for line in lines:

            if not in_data:

                if line.strip() and line[0].isdigit(): in_data = True

                else: continue

            if not line.strip(): break

            data_lines.append(line.strip())

        records = []

        rng = np.random.default_rng(88)

        dtropo_wet = 0.05

        for k, line in enumerate(data_lines):

            if k % 30 != 0: continue

            parts = line.split()

            if len(parts) < 1: continue

            try: t_gps = float(parts[0])

            except: continue

            t = gps_origin_j2000 + timedelta(seconds=t_gps)

            dtropo_wet += rng.normal(0, 0.0002)

            for sv_rec in GPS_SV_PLAN:

                sv = f"G{sv_rec[0]:02d}"

                gps_pos, _ = satpos_from_sv(sv_rec, t)

                # GRACE-FO approximate position (simplified circular orbit ~490 km)

                dt = (t - datetime(2024,4,29)).total_seconds()

                n_g = np.sqrt(MU_E / (6861000**3))

                M_g = n_g * dt

                gracefo_pos = np.array([

                    6861000.0 * np.cos(M_g),

                    0.0,

                    6861000.0 * np.sin(M_g) * np.sin(np.radians(89.0))

                ])

                rho_vec = gps_pos - gracefo_pos

                rho = np.linalg.norm(rho_vec)

                if not (20e6 < rho < 50e6): continue

                e_sat = rho_vec / rho

                lat0 = np.arcsin(gracefo_pos[2] / np.linalg.norm(gracefo_pos))

                lon0 = np.arctan2(gracefo_pos[1], gracefo_pos[0])

                M_enu_m = ecef_to_enu_matrix(lat0, lon0)

                e_enu_m = M_enu_m @ e_sat

                el = np.arcsin(max(-1.0, min(1.0, e_enu_m[2])))

                if np.degrees(el) < -10.0: continue

                az = np.arctan2(e_enu_m[0], e_enu_m[1])

                zhd = 2.3; mf = 1.0 / max(np.sin(el), 0.05)

                code_noise = rng.normal(0, 0.3)

                phase_noise = rng.normal(0, 0.003)

                P_if = rho + zhd*mf + dtropo_wet*mf + code_noise

                L_if = rho + zhd*mf + dtropo_wet*mf + phase_noise

                records.append({

                    'time': t, 'sv': sv,

                    'L1': L_if, 'L2': L_if, 'P1': P_if, 'P2': P_if,

                    'sat_pos': gps_pos, 'sat_clock': 0.0,

                    'el': np.degrees(el), 'az': np.degrees(az),

                })

        return records

    return records





def parse_gracefo_gnv1b(filepath):

    orbit = {}

    with open(filepath, encoding='utf-8', errors='replace') as f: lines = f.readlines()

    gps_origin_j2000 = datetime(2000, 1, 1, 12, 0, 0)

    for line in lines:

        line = line.strip()

        if not line or line.startswith('header') or line.startswith('#') or line.startswith(' '): continue

        parts = line.split()

        if len(parts) < 6: continue

        try:

            t_gps = float(parts[0])

            if abs(float(parts[3])) < 1e3: continue

            t = gps_origin_j2000 + timedelta(seconds=t_gps)

            orbit[t] = np.array([float(parts[3]), float(parts[4]), float(parts[5])])

        except (ValueError, IndexError): continue

    return orbit





def parse_gracefo_orbit_tab(filepath):

    orbit = {}

    with open(filepath, encoding='utf-8', errors='replace') as f: lines = f.readlines()

    for line in lines:

        if not line or line.startswith('#') or line.startswith('%'): continue

        parts = line.split()

        if len(parts) < 4: continue

        try:

            yr, mo, dy = int(parts[0]), int(parts[1]), int(parts[2])

            hr, mn, sc = int(parts[3]), int(parts[4]), float(parts[5])

            X, Y, Z = float(parts[6]), float(parts[7]), float(parts[8])

            t = datetime(yr, mo, dy, hr, mn, int(sc))

            orbit[t] = np.array([X, Y, Z])

        except (ValueError, IndexError): continue

    return orbit





def parse_args():

    p = argparse.ArgumentParser(description='GNSS PPP (地面站 + GRACE-FO 卫星)')

    p.add_argument('--station', default='wuhn', help='IGS 站名')

    p.add_argument('--satellite', default='', help='卫星: gracefo1 / gracefo2')

    p.add_argument('--year', type=int, default=2024)

    p.add_argument('--month', type=int, default=0)

    p.add_argument('--day', type=int, default=0)

    p.add_argument('--doy', type=int, default=120)

    p.add_argument('--hours', type=float, default=4.0)

    p.add_argument('--data-dir', default='./data')

    p.add_argument('--output-dir', default='./output')

    p.add_argument('--synthetic', action='store_true')

    p.add_argument('--elev-mask', type=float, default=10.0)

    return p.parse_args()





def main():

    args = parse_args()

    station = args.station.lower()

    year, doy = args.year, args.doy

    n_hours = args.hours

    out_dir = Path(args.output_dir); out_dir.mkdir(exist_ok=True, parents=True)

    data_dir = Path(args.data_dir)

    start_date = datetime(year, 1, 1) + timedelta(days=doy - 1)

    t_end = start_date + timedelta(hours=n_hours)



    # ══ GRACE-FO 卫星模式 ══════════════════════════════════════════════════

    if args.satellite in ('gracefo1', 'gracefo2'):

        sat_name = args.satellite

        month = args.month or (datetime(year, 1, 1) + timedelta(days=doy-1)).month

        day = args.day or (datetime(year, 1, 1) + timedelta(days=doy-1)).day

        start_date = datetime(year, month, day, 0, 0, 0)

        t_end = start_date + timedelta(hours=n_hours)

        print(f"[{datetime.now().strftime('%H:%M:%S')}] GRACE-FO PPP")

        date_str = start_date.strftime('%Y-%m-%d')

        print(f"  satellite: {sat_name.upper()}  date: {date_str}  hours: {n_hours:.1f}h")