#!/usr/bin/env python3
"""Check available GPS1A files for target dates"""
import urllib.request, ssl, re

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

url = 'https://isdc-data.gfz.de/grace-fo/Level-1A/JPL/INSTRUMENT/RL04/2024/'
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
r = urllib.request.urlopen(req, timeout=30, context=ctx)
html = r.read().decode()

targets = ['2024-04-29', '2024-04-30', '2024-05-01', '2024-05-02']
for line in html.split('\n'):
    for t in targets:
        if t in line:
            m = re.search(r'href="([^"]+)"', line)
            if m:
                print(m.group(1))
