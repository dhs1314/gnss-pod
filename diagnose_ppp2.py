#!/usr/bin/env python3
"""Check if satellite clock is being applied"""
import pickle, sys, math, datetime as dt, numpy as np
sys.path.insert(0, '/workspace/gnss_pod/src')
from batch_v12 import parse_gnv1b
from run_ppp import GPS_SV_PLAN, satpos_from_sv
from ppp import PPPProcessor

gps1b = pickle.load(open('/workspace/gnss_pod/data/gracefo/2024/2024-04-29/gps1b_C.pkl', 'rb'))
gnv = parse_gnv1b('/workspace/gnss_pod/data/gracefo/2024/2024-04-29/GNV1B_2024-04-29_C_04.txt')
grace_pos = np.array(gnv[sorted(gnv.keys())[0]])
sv_plan = {int(r[0]): r for r in GPS_SV_PLAN}

# Check first obs
t = sorted(gps1b.keys())[0]
obs = gps1b[t]
sv_id = list(obs.keys())[0]
utc = dt.datetime(1980,1,6) + dt.timedelta(seconds=t)

# Get GPS position AND clock
sv_rec = sv_plan.get(int(sv_id[1:]))
print(f"SV plan record for {sv_id}: {sv_rec}")
print(f"  (prn, M0, Omega0, omega, ecc, inc, u0, sqrtA)")

# Check what satpos_from_sv returns
result = satpos_from_sv(sv_rec, utc)
if result:
    pos, vel = result
    print(f"Sat pos: {np.array(pos)}")
    print(f"Sat |pos|: {np.linalg.norm(np.array(pos)):.0f} m")
else:
    print("No position!")

# Key: GPS L1 = 21868738.004 m = 21868.738 km
# rho = |GPS_pos - GRACE_pos| = 32336793 m = 32336.8 km
# L1 - rho = -10468056 m = -10468.1 km

# This huge residual (-10468 km) MUST be clock offset!
# GPS satellite clock offset = (L1 - rho) / C
# = -10468056 / 299792458 = -34.93 ms

L1 = obs[sv_id]['L1']
L2 = obs[sv_id]['L2']
F1, F2 = 1575.42e6, 1227.60e6
ALPHA = F1*F1/(F1*F1-F2*F2)
BETA = -F2*F2/(F1*F1-F2*F2)
L_if = ALPHA*L1 + BETA*L2
rho = np.linalg.norm(np.array(pos) - grace_pos)
clock_bias = (L_if - rho) / 299792458 * 1000  # in ms
print(f"\nL_if = {L_if:.3f} m ({L_if/1000:.0f} km)")
print(f"rho  = {rho:.0f} m ({rho/1000:.0f} km)")
print(f"L_if - rho = {(L_if-rho)/1000:.0f} km")
print(f"Clock bias = {clock_bias:.3f} ms")
print(f"  (typical GPS satellite clock offset: ±1 ms, max ±2 ms)")
print(f"  This {abs(clock_bias):.1f} ms offset is {abs(clock_bias)/1:.1f}x larger than normal!")
print(f"  → Possible unit mismatch or data error")

# What if L1 is NOT in meters?
# If L1 is in km: L1 = 21868738 km → L_if = 21868737 km
# rho = 32337 km
# L_if - rho = 21868737 - 32337 = 21836400 km
# → Even worse

# If L1 is in microseconds:
# L_if_us = 21868737 μs
# L_if_m = 21868737 × 299792458 / 1e6 = 6.556 × 10^9 m = 6556083 km
# rho = 32337 km
# L_if - rho = 6523746 km → totally wrong

# So L1 IS in meters, and the clock bias is real
# But 35 ms is way too large for GPS satellite clock
# Something else must be wrong

# Let me check the GPS orbital radius vs the GPS L1 value
print(f"\n=== GPS orbital check ===")
print(f"L1 (GPS1B) = {L1:.3f} m = {L1/1000:.0f} km")
print(f"|GPS_pos| = {np.linalg.norm(np.array(pos)):.0f} m = {np.linalg.norm(np.array(pos))/1000:.0f} km")
print(f"Difference: {(L1 - np.linalg.norm(np.array(pos)))/1000:.0f} km")
print(f"  L1 is {L1/np.linalg.norm(np.array(pos))*100:.1f}% of orbital radius")
print(f"  Expected: 20265 km for GPS-to-ground range")
print(f"  GRACE orbital radius: 6857 km")
print(f"\n→ L1 ≈ 21869 km ≈ GPS orbital radius")
print(f"  This is satellite-to-Earth-center distance (not satellite-to-GRACE)")
print(f"  The L1 value in GPS1B = |GPS satellite position from Earth's center|?")
print(f"  No, that doesn't make physical sense either...")

# Wait - GPS L1 values are typically ~20000-25000 km
# And orbital radius is ~26560 km
# 21869 km is within the GPS orbital range
# But is L1 the range from GPS to GRACE or GPS to Earth's center?

# If L1 = range(GPS, GRACE) = 32367 km
# Then L1 should be ~32367 km
# But actual L1 = 21869 km ← wrong

# If L1 = |GPS_pos| = 26560 km
# Then L1 should be ~26560 km
# But actual L1 = 21869 km ← wrong by 4688 km

# If L1 is GPS pseudo-range in meters
# Then it's raw counts... and we already checked the unit conversion

# Let me check if maybe the GPS position is wrong
# GPS is at ~26560 km orbital radius
# L1 = 21869 km = 82.3% of expected
# This could mean the GPS position calculation is wrong
# OR the GPS orbital radius is stored in wrong units in GPS_SV_PLAN

# Check: what is the GPS orbital radius in the SV plan?
print(f"\n=== GPS SV plan orbital radius ===")
for r in GPS_SV_PLAN[:5]:
    prn, M0, Omega0, omega, ecc, inc, u0, sqrtA = r
    a_m = sqrtA**2  # assuming meters
    a_km = (sqrtA/1000)**2  # assuming km
    print(f"  PRN {prn}: sqrtA={sqrtA}")
    print(f"    if sqrtA in meters: a={a_m/1e6:.2f} km")
    print(f"    if sqrtA in km: a={a_km:.2f} km = {a_km*1000:.0f} m")