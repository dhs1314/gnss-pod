#!/usr/bin/env python3
"""Check exactly what batch_v12_fixed.py computes for one epoch"""
import pickle, numpy as np
from datetime import datetime, timedelta
from pathlib import Path
import sys
sys.path.insert(0, str(Path('.')))
from src.gps1b_loader import gps_sod_to_utc
from src.sp3_loader import get_gps_pos_from_sp3

C=299792458.0; F1=1575.42e6; F2=1227.60e6
F1_SQ,F2_SQ=F1*F1,F2*F2

gps_obs=pickle.load(open("data/gracefo/2024/2024-04-29/GPS1B_2024-04-29_C_04.pkl","rb"))
sp3=pickle.load(open("data/2024/120/igs_sp3_FIN.pkl","rb"))

GPS0=datetime(2000,1,1,12,0,0)
ref_orbit={}
for l in open("data/gracefo/2024/2024-04-29/GNV1B_2024-04-29_C_04.txt"):
    p=l.split()
    if len(p)<6: continue
    try:
        tg=float(p[0]); flag=p[2]
        if flag in ("C","E"): ref_orbit[GPS0+timedelta(seconds=tg)]=np.array([float(p[3]),float(p[4]),float(p[5])])
    except: continue

target=min(gps_obs.keys())+3600
for sod in sorted(gps_obs.keys()):
    if sod>=target: break

utc_batch=gps_sod_to_utc(sod)  # subtracts 18s
utc_direct=datetime(2000,1,1,12,0,0)+timedelta(seconds=sod)

print(f"sod={sod}")
print(f"batch utc={utc_batch}")
print(f"direct utc={utc_direct}")
print(f"SP3 epochs: {sp3['ts'][0]} to {sp3['ts'][-1]}")

# Reference position at utc_batch
orbit_ts=sorted(ref_orbit.keys())
t0=t1=None
for j,ti in enumerate(orbit_ts):
    if ti>=utc_batch: t1=ti; t0=orbit_ts[j-1] if j>0 else None; break
    t0=ti
if t1 is None: t0=t1=orbit_ts[-1]
if t0 is None: t0=orbit_ts[0]
dt0=(utc_batch-t0).total_seconds(); dtt=(t1-t0).total_seconds()
gr_batch=ref_orbit[t0]*(1-dt0/dtt)+ref_orbit[t1]*(dt0/dtt) if dtt else ref_orbit[t0]

# Reference at utc_direct
t0=t1=None
for j,ti in enumerate(orbit_ts):
    if ti>=utc_direct: t1=ti; t0=orbit_ts[j-1] if j>0 else None; break
    t0=ti
if t1 is None: t0=t1=orbit_ts[-1]
if t0 is None: t0=orbit_ts[0]
dt0=(utc_direct-t0).total_seconds(); dtt=(t1-t0).total_seconds()
gr_direct=ref_orbit[t0]*(1-dt0/dtt)+ref_orbit[t1]*(dt0/dtt) if dtt else ref_orbit[t0]

print(f"GRACE pos differs by: {np.linalg.norm(gr_batch-gr_direct):.1f} m")
print(f"(because reference orbit is interpolated at different UTC times)")

# LAM_IF from batch_v12
F_IF=(F1_SQ-F2_SQ)/(F1+F2)
LAM_IF=C/F_IF
print(f"\nF_IF={F_IF:.2f} Hz, LAM_IF={LAM_IF:.4f} m")
print(f"F1-F2={F1-F2:.2f} Hz, C/(F1-F2)={C/(F1-F2):.4f} m")

print(f"\n{'SV':5s} {'L_if_raw(km)':14s} {'L_if*LAM(km)':14s} {'rng(km)':10s} {'clk(km)':10s} {'O-C(m)':10s}")
for sv,rec in sorted(gps_obs[sod].items()):
    sp,clk,sv_v=get_gps_pos_from_sp3(sp3,sv,utc_batch)
    if sp is None: continue
    d=gr_batch-sp; rng=float(np.linalg.norm(d))
    if not(2e7<rng<5e7): continue
    L_if=float(rec["L_if"])
    L_if_scaled=L_if*LAM_IF
    mod_l=L_if_scaled-(rng+clk)  # KF model
    mod_p=float(rec["P_if"])-(rng+clk)
    print(f"{sv:5s} {L_if/1000:14.3f} {L_if_scaled/1000:14.3f} {rng/1000:10.1f} {clk/1000:10.1f} {mod_l:10.1f}")

# Now check what KF would do:
# The KF phase weight = 1/(0.003^2) ~ 111111
# The KF code weight = 1/(0.3^2) ~ 11
# Phase dominates. The phase says L_if_scaled = rho + clk_r
# Since L_if_scaled ~ 18,860 km and rho ~ 23,800 km
# The KF clk_r must be ~ -5,000 km = -0.017 seconds
print("\n--- Diagnosis ---")
print("The KF in batch_v12 sees phase observations at L_if*LAM_IF scale,")
print("which is ~0.862x the true range. The clock state absorbs the difference.")
print("Since phase weight >> code weight, the position is from phase only.")
print("But the geometry is compressed by 0.862x, giving wrong position.")
