#!/usr/bin/env python3
"""Study GPS1A RINEX data structure before building loader"""
import tarfile

t = tarfile.open("data/gracefo_1A_2024-04-29_RL04.ascii.noLRI.tgz")
f = t.extractfile("GPS1A_2024-04-29_C_04.rnx")
data = f.read().decode("ascii", errors="replace")
t.close()

# Find END OF HEADER
hdr_end = data.find("END OF HEADER")
body = data[hdr_end + 50:]

# Parse a few epochs, understanding the 2-line-per-SV layout
lines = [l.strip() for l in body.split("\n") if l.strip()]
print("First 30 non-empty lines:")
for i, line in enumerate(lines[:30]):
    print(f"  [{i:3d}] {line[:120]}")

print("\n--- Parsing epochs ---")
epoch_idx = 0
i = 0
while i < len(lines) and epoch_idx < 3:
    line = lines[i]
    parts = line.split()
    # Epoch header: starts with 2-digit year
    if len(parts) >= 8 and len(parts[0]) == 2 and parts[0].isdigit():
        yr, mo, dy, hr, mn = [int(x) for x in parts[:5]]
        sec = float(parts[5])
        epoch_flag = int(parts[6])
        n_sv = int(parts[7])
        sv_list = [int(x) for x in parts[8:8+n_sv]]

        print(f"\nEpoch {epoch_idx}: {yr:02d}-{mo:02d}-{dy:02d} {hr:02d}:{mn:02d}:{sec:06.3f} flag={epoch_flag} n_sv={n_sv}")
        print(f"  SV list: {sv_list}")

        i += 1
        # Now parse n_sv SVs, each taking 2 lines (9 obs ÷ 5 per line)
        for sv_idx, sv in enumerate(sv_list):
            if i >= len(lines):
                break

            # Line 1: L1, L2, C1, P1, P2 (first 5 of 9 obs types)
            line1 = lines[i].split()
            # Line 2: LA, SA, S1, S2 (next 4 of 9 obs types)
            line2 = lines[i+1].split() if i+1 < len(lines) else []

            L1_cyc = float(line1[0]) if len(line1) > 0 else 0
            L2_cyc = float(line1[1]) if len(line1) > 1 else 0
            C1_m   = float(line1[2]) if len(line1) > 2 else 0
            P1_m   = float(line1[3]) if len(line1) > 3 else 0
            P2_m   = float(line1[4]) if len(line1) > 4 else 0

            # Check LLI (Loss of Lock Indicator) — last digit of phase value in RINEX
            # Format: F14.3 with last digit = LLI flag
            # Actually, RINEX uses the last digit after decimal as LLI

            print(f"  SV G{sv:02d}: L1={L1_cyc:14.5f}cyc  L2={L2_cyc:13.5f}cyc  "
                  f"C1={C1_m:13.3f}m  P1={P1_m:13.3f}m  P2={P2_m:13.3f}m  "
                  f"line2_parts={len(line2)}")

            i += 2  # Each SV takes 2 lines

        epoch_idx += 1
    else:
        i += 1
