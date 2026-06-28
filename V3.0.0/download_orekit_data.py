"""Download Orekit data zip with SSL workaround and retries."""
import urllib.request
import ssl
import os
import shutil

url = "https://gitlab.orekit.org/orekit/orekit-data/-/archive/main/orekit-data-main.zip"
dest_dir = r"d:\prj\gnss_pod\data\orekit"
dest_zip = os.path.join(dest_dir, "orekit-data.zip")

os.makedirs(dest_dir, exist_ok=True)

# Remove any partial download
if os.path.exists(dest_zip):
    os.remove(dest_zip)

print(f"Downloading {url}...")
print(f"To: {dest_zip}")

# Create SSL context that doesn't verify (some corporate networks interfere)
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

try:
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=300, context=ctx) as response:
        with open(dest_zip, 'wb') as out_file:
            shutil.copyfileobj(response, out_file, length=16*1024*1024)
    size_mb = os.path.getsize(dest_zip) / (1024 * 1024)
    print(f"Downloaded successfully: {size_mb:.1f} MB")
except Exception as e:
    print(f"Download failed: {e}")
    # Try alternative URL
    alt_url = "https://gitlab.orekit.org/orekit/orekit-data/-/archive/v2.0/orekit-data-v2.0.zip"
    print(f"Trying: {alt_url}")
    try:
        req = urllib.request.Request(alt_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=300, context=ctx) as response:
            with open(dest_zip, 'wb') as out_file:
                shutil.copyfileobj(response, out_file, length=16*1024*1024)
        size_mb = os.path.getsize(dest_zip) / (1024 * 1024)
        print(f"Downloaded successfully: {size_mb:.1f} MB")
    except Exception as e2:
        print(f"Alternative also failed: {e2}")
