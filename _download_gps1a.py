#!/usr/bin/env python3
"""Download GPS1A data and inspect contents"""
import urllib.request, ssl, tarfile, os

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

url = "https://isdc-data.gfz.de/grace-fo/Level-1A/JPL/INSTRUMENT/RL04/2024/gracefo_1A_2024-04-29_RL04.ascii.noLRI.tgz"
tgz_path = "data/gracefo_1A_2024-04-29_RL04.ascii.noLRI.tgz"

if not os.path.exists(tgz_path):
    print(f"Downloading {url} ...")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    r = urllib.request.urlopen(req, timeout=600, context=ctx)
    total = int(r.headers.get("Content-Length", 0))
    print(f"  Size: {total/1024/1024:.0f} MB")
    data = r.read()
    os.makedirs("data", exist_ok=True)
    with open(tgz_path, "wb") as f:
        f.write(data)
    print(f"  Saved to {tgz_path}")
else:
    print(f"Already downloaded: {tgz_path}")

# List contents
print("\nArchive contents:")
t = tarfile.open(tgz_path)
for name in sorted(t.getnames()):
    info = t.getmember(name)
    print(f"  {name} ({info.size/1024/1024:.1f} MB)")
t.close()
