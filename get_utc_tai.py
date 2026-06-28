"""Download the official Orekit UTC-TAI.history file from gitlab raw."""
import urllib.request
import ssl

url = "https://gitlab.orekit.org/orekit/orekit-data/-/raw/main/UTC-TAI.history"
dest = r"d:\prj\gnss_pod\data\orekit\UTC-TAI.history"

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

print(f"Downloading {url}...")
try:
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=60, context=ctx) as response:
        content = response.read()
    with open(dest, 'wb') as f:
        f.write(content)
    print(f"Downloaded {len(content)} bytes")
    print(f"First 300 chars:\n{content[:300].decode('utf-8', errors='replace')}")
except Exception as e:
    print(f"Download failed: {e}")
    # Try alternative: use the IERS official file
    print("\nTrying IERS official file...")
    url2 = "https://maia.usno.navy.mil/ser7/tai-utc.dat"
    try:
        req = urllib.request.Request(url2, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=60, context=ctx) as response:
            content = response.read()
        with open(dest, 'wb') as f:
            f.write(content)
        print(f"Downloaded {len(content)} bytes")
        print(f"First 300 chars:\n{content[:300].decode('utf-8', errors='replace')}")
    except Exception as e2:
        print(f"Alternative also failed: {e2}")
