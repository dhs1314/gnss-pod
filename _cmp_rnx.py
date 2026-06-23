#!/usr/bin/env python3
"""Compare GPS1A vs GPS1B RINEX data to check if 1A has raw integer ambiguity"""
import tarfile

# Extract RINEX headers and first epoch from both
for level, tgz_name, member in [
    ("1A", "data/gracefo_1A_2024-04-29_RL04.ascii.noLRI.tgz", "GPS1A_2024-04-29_C_04.rnx"),
    ("1B", "data/gracefo_1B_2024-04-29_RL04.ascii.noLRI.tgz", "GPS1B_2024-04-29_C_04.rnx"),
]:
    t = tarfile.open(tgz_name)
    f = t.extractfile(member)
    lines = f.read().decode("ascii", errors="replace")
    t.close()

    # Print header
    header_end = lines.find("END OF HEADER")
    print(f"=" * 70)
    print(f"=== GPS{level} RINEX Header: {member}")
    print(f"=" * 70)
    print(lines[:header_end + 50])
    print()

    # Print first 2 epochs of data
    data = lines[header_end + 50:]
    epoch_count = 0
    for line in data.split("\n"):
        if line.strip():
            print(line)
            epoch_count += 1
            if epoch_count >= 6:
                break
    print()
    print()
