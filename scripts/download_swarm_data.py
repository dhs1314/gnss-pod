"""Download SWARM GPS RINEX observations and precise orbit products.

Data sources:
  1. SWARM GPS RINEX: ESA via https://swarm-diss.eo.esa.int
     (free registration required at https://swarm-diss.eo.esa.int)
  2. SWARM precise orbits (SP3): CODE/AIUB LEO products

Usage:
  python scripts/download_swarm_data.py --date 2024-04-29 --sat A
  python scripts/download_swarm_data.py --date 2024-04-29 --sat A --auto

Data storage:
  data/swarm/{date}/           -- GPS RINEX obs and converted pkl
  data/CODE/{year}/             -- LEO SP3 precise orbits (shared with GRACE)
"""
import sys, os, argparse, urllib.request, gzip, ftplib
from pathlib import Path
from datetime import datetime

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / 'src'))

# ── ESA SWARM FTP credentials ──
# Register at https://swarm-diss.eo.esa.int for free access
ESA_HOST = 'swarm-diss.eo.esa.int'
ESA_USER = 'anonymous'
ESA_PASS = 'guest@example.com'

parser = argparse.ArgumentParser(description='Download SWARM GPS data')
parser.add_argument('--date', required=True, help='Date YYYY-MM-DD')
parser.add_argument('--sat', choices=['A','B','C'], default='A')
parser.add_argument('--auto', action='store_true', help='Use ESA FTP (requires registration)')
args = parser.parse_args()

dt = datetime.strptime(args.date, '%Y-%m-%d')
y = dt.year; doy = dt.strftime('%j'); yy = dt.strftime('%y')
SAT_ID = f'Sat_{args.sat}'
SWARM_DIR = ROOT / 'data' / 'swarm' / dt.strftime('%Y') / args.date
SWARM_DIR.mkdir(parents=True, exist_ok=True)

# ── 1. SWARM GPS RINEX download ──
print(f"=== SWARM-{args.sat} GPS data for {args.date} (DOY {doy}) ===")

if args.auto:
    try:
        print("Connecting to ESA Swarm FTP...")
        ftp = ftplib.FTP(ESA_HOST, timeout=60)
        ftp.login(ESA_USER, ESA_PASS)

        rinex_path = f'/Level1b/Latest_baselines/RINEX/{SAT_ID}/RINEX'
        ftp.cwd(rinex_path)

        files = ftp.nlst()
        matches = [f for f in files if doy in f and ('GPS' in f.upper() or 'RINEX' in f.upper())]

        if not matches:
            print(f"No files found for DOY {doy}")
            ftp.quit()
            sys.exit(1)

        for fname in matches[:3]:
            local_gz = SWARM_DIR / fname
            local_out = SWARM_DIR / fname.replace('.Z', '').replace('.gz', '')

            if local_out.exists() and local_out.stat().st_size > 1000:
                print(f"  {fname}: already exists")
                continue

            print(f"  Downloading {fname}...", end=' ', flush=True)
            with open(str(local_gz), 'wb') as f:
                ftp.retrbinary(f'RETR {fname}', f.write)

            if fname.endswith('.Z') or fname.endswith('.gz'):
                import gzip
                with gzip.open(str(local_gz)) as gf, open(str(local_out), 'wb') as of:
                    of.write(gf.read())

            sz = local_out.stat().st_size
            print(f"OK ({sz/1e6:.1f}MB)")

        ftp.quit()
    except Exception as e:
        print(f"FTP failed: {e}")
        print("Manual download: visit https://swarm-diss.eo.esa.int")
        print(f"  Navigate to: Level1b > Latest_baselines > RINEX > {SAT_ID} > RINEX")
        print(f"  Download files containing '{doy}' with GPS/RINEX")
        sys.exit(1)

# ── 2. CODE SWARM precise orbits (LEO SP3) ──
print()
print("=== SWARM precise orbit (CODE LEO SP3) ===")

CODE_LEO_DIR = ROOT / 'data' / 'CODE' / str(y)

# CODE combined LEO SP3: contains multi-mission orbits
# Naming: COD0OPSFIN_{year}{doy}0000_01D_05M_LEO_ORB.SP3.gz
for fname_gz in [
    f'COD0OPSFIN_{y}{doy}0000_01D_05M_LEO_ORB.SP3.gz',
    f'COD0MGXFIN_{y}{doy}0000_01D_05M_LEO.SP3.gz',
]:
    local_gz = CODE_LEO_DIR / fname_gz
    local_out = CODE_LEO_DIR / fname_gz.replace('.gz', '')

    if local_out.exists() and local_out.stat().st_size > 1000:
        print(f"  {local_out.name}: already exists ({local_out.stat().st_size/1e3:.0f}KB)")
        continue

    url = f'http://ftp.aiub.unibe.ch/CODE/{y}/{fname_gz}'
    try:
        print(f"  Downloading {fname_gz}...", end=' ', flush=True)
        urllib.request.urlretrieve(url, str(local_gz))
        with gzip.open(str(local_gz)) as gf, open(str(local_out), 'wb') as of:
            of.write(gf.read())
        print(f"OK ({local_out.stat().st_size/1e3:.0f}KB)")
    except Exception as e:
        print(f"FAILED: {e}")

print()
print(f"Data stored in: {SWARM_DIR}")
print(f"CODE products in: {CODE_LEO_DIR}")
print()
print("Next: run SWARM conversion: python scripts/convert_swarm_rnx.py --date {args.date} --sat {args.sat}")
