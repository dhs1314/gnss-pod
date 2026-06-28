"""Download essential Orekit data files from fast sources."""
import urllib.request
import os

dest_dir = r"d:\prj\gnss_pod\data\orekit"
os.makedirs(dest_dir, exist_ok=True)

# 1. EOP: finals2000A.all from IERS (~5 MB)
eop_urls = [
    "https://datacenter.iers.org/data/latestVersion/finals2000A.all",
    "https://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now",
]
eop_dest = os.path.join(dest_dir, "finals2000A.all")
print("Downloading EOP data...")
for url in eop_urls:
    try:
        with urllib.request.urlopen(url, timeout=60) as response, open(eop_dest, 'wb') as f:
            import shutil
            shutil.copyfileobj(response, f)
        size_kb = os.path.getsize(eop_dest) / 1024
        print(f"  Downloaded {size_kb:.0f} KB from {url[:50]}...")
        break
    except Exception as e:
        print(f"  Failed: {url[:50]}...: {e}")
else:
    print("  WARNING: Could not download EOP file. Orekit frame transforms may fail.")

# 2. Also try to get leap seconds from IERS
leap_dest = os.path.join(dest_dir, "tai-utc.dat")
leap_url = "https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat"
print("Downloading leap seconds...")
try:
    with urllib.request.urlopen(leap_url, timeout=30) as response, open(leap_dest, 'wb') as f:
        import shutil
        shutil.copyfileobj(response, f)
    print(f"  Downloaded leap seconds file")
except Exception as e:
    print(f"  Failed to download leap seconds: {e}")

print(f"\nFiles in {dest_dir}:")
for f in os.listdir(dest_dir):
    size = os.path.getsize(os.path.join(dest_dir, f))
    print(f"  {f}: {size} bytes")
