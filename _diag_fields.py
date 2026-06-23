#!/usr/bin/env python3
"""Check ALL available fields in the BATCH pickle"""
import pickle
import numpy as np

gps1b = pickle.load(open("data/gracefo/2024/2024-04-29/GPS1B_2024-04-29_C_04.pkl", "rb"))

# Get all field names across all SVs and epochs
all_fields = set()
field_counts = {}
for gps_sod, sv_data in gps1b.items():
    for sv, rec in sv_data.items():
        for k in rec.keys():
            all_fields.add(k)
            field_counts[k] = field_counts.get(k, 0) + 1

print("All fields in BATCH pickle:")
for f in sorted(all_fields):
    print(f"  {f}: present in {field_counts[f]} records")

# Now check a few SVs at the first epoch for ALL fields
gps_sod = sorted(gps1b.keys())[0]
print(f"\nFirst epoch: gps_sod={gps_sod}")
for sv, rec in sorted(gps1b[gps_sod].items()):
    print(f"\n{sv} (all fields):")
    for k, v in sorted(rec.items()):
        if isinstance(v, float):
            print(f"  {k}: {v:.6f}")
        else:
            print(f"  {k}: {v}")
    if sv >= 'G03':
        break

# Check: is there a field with phase in CYCLES?
# L1_phase should be in cycles if it's raw phase, or in meters if smoothed
# Let's check L1_phase range for one SV across all epochs
print(f"\n--- L1_phase vs L1_range values ---")
for sv in ['G05', 'G07', 'G16']:
    l1_phase_vals = []
    l1_range_vals = []
    ca_phase_vals = []
    for gps_sod in sorted(gps1b.keys()):
        if sv in gps1b[gps_sod]:
            rec = gps1b[gps_sod][sv]
            if 'L1_phase' in rec:
                l1_phase_vals.append(float(rec['L1_phase']))
            if 'L1_range' in rec:
                l1_range_vals.append(float(rec['L1_range']))
            if 'CA_phase' in rec:
                ca_phase_vals.append(float(rec['CA_phase']))

    if l1_phase_vals:
        arr = np.array(l1_phase_vals)
        print(f"{sv} L1_phase: mean={np.mean(arr):.3f} std={np.std(arr):.3f} "
              f"range=[{np.min(arr):.3f}, {np.max(arr):.3f}] ({len(arr)} epochs)")
    if l1_range_vals:
        arr = np.array(l1_range_vals)
        print(f"{sv} L1_range: mean={np.mean(arr):.3f} std={np.std(arr):.3f} "
              f"range=[{np.min(arr):.3f}, {np.max(arr):.3f}] ({len(arr)} epochs)")
    if ca_phase_vals:
        arr = np.array(ca_phase_vals)
        print(f"{sv} CA_phase: mean={np.mean(arr):.3f} std={np.std(arr):.3f} "
              f"range=[{np.min(arr):.3f}, {np.max(arr):.3f}] ({len(arr)} epochs)")

# Check: is L1_phase actually in cycles but stored as float?
# If L1_phase is in cycles: values would be ~100 million
# If L1_phase is in meters: values would be ~20 million
# Let's check the magnitude
C = 299792458.0
F1 = 1575.42e6
LAM1 = C / F1  # ~0.1903 m

for sv in ['G05']:
    if sv in gps1b[gps_sod]:
        rec = gps1b[gps_sod][sv]
        if 'L1_phase' in rec:
            val = float(rec['L1_phase'])
            print(f"\n{sv} L1_phase={val:.3f}")
            print(f"  As meters: {val:.3f} m")
            print(f"  As cycles * lambda1: {val * LAM1:.3f} m (if val were cycles)")
            print(f"  Check: GPS range ~20,180 km. Does {val/1000:.0f} km make sense as range? {'YES (meter-level range)' if 15000 < val/1000 < 30000 else 'NO'}")
            print(f"  Check: Is {val:.0f} a reasonable cycle count? {'YES (~100M cycles)' if 50000000 < val < 200000000 else 'NO'}")

# The key question: does L1_phase differ from L1_range?
print(f"\n--- L1_phase vs L1_range comparison ---")
for sv in ['G05', 'G07']:
    if sv in gps1b[gps_sod]:
        rec = gps1b[gps_sod][sv]
        lp = rec.get('L1_phase')
        lr = rec.get('L1_range')
        cp = rec.get('CA_phase')
        cr = rec.get('CA_range')
        print(f"{sv}: L1_phase={lp}, L1_range={lr}, diff={float(lp)-float(lr) if (lp and lr) else 'N/A':.4f}")
        print(f"     CA_phase={cp}, CA_range={cr}")
