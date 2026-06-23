"""Download CODE products for dates with GPS1B+GNV1B data."""
import subprocess, sys, os

dates = sys.argv[1:] if len(sys.argv) > 1 else []
if not dates:
    print("Usage: py download_code_products.py 2024-04-30 2024-05-01 ...")
    sys.exit(1)

# CDDIS base URL (use AIUB mirror)
BASE = "http://ftp.aiub.unibe.ch/CODE/2024"
OUT_DIR = r"d:\prj\gnss_pod\data\CODE\2024"
os.makedirs(OUT_DIR, exist_ok=True)

from datetime import datetime

for d in dates:
    dt = datetime.strptime(d, "%Y-%m-%d")
    doy = dt.strftime("%j")
    y = dt.year

    files_to_get = [
        f"COD0OPSFIN_{y}{doy}0000_01D_05M_ORB.SP3.gz",
        f"COD0OPSFIN_{y}{doy}0000_01D_30S_CLK.CLK.gz",
    ]

    for fn in files_to_get:
        local_path = os.path.join(OUT_DIR, fn)
        if os.path.exists(local_path):
            print(f"  {fn}: already exists")
            continue
        url = f"{BASE}/{fn}"
        print(f"  Downloading {fn}...", end=" ", flush=True)
        try:
            import urllib.request
            urllib.request.urlretrieve(url, local_path)
            print("OK")
        except Exception as e:
            print(f"FAILED: {e}")

print("\nDone. Run gunzip on .gz files if needed.")
