#!/usr/bin/env python3
"""Check GPS1A ASCII YAML header and raw values"""
import tarfile

t = tarfile.open("data/gracefo_1A_2024-04-29_RL04.ascii.noLRI.tgz")
f = t.extractfile("GPS1A_2024-04-29_C_04.txt")
lines = f.read().decode("ascii", errors="replace")
t.close()

# Print YAML header
idx = lines.find("# End of YAML header")
print(lines[:idx])
print("=== End of YAML header ===")

# Print first 5 data lines
rest = lines[idx + 50:]
count = 0
for line in rest.split("\n"):
    if line.strip() and not line.startswith("#"):
        print(line)
        count += 1
        if count >= 5:
            break
