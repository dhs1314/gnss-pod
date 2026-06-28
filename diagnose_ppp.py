#!/usr/bin/env python3
"""Diagnose PPP convergence failure for Apr29"""
import pickle, sys, math, datetime as dt, numpy as np
sys.path.insert(0, '/workspace/gnss_pod/src')
from batch_v12 import parse_gnv1b
from run_ppp import GPS_SV_PLAN, satpos_from_sv
from ppp import PPPProcessor

gps1b = pickle.load(open('/workspace/gnss_pod/data/gracefo/2024/2024-04-29/gps1b_C.pkl', 'rb'))
gnv = parse_gnv1b('/workspace/gnss_pod/data/gracefo/2024/2024-04-29/GNV1B_2024-04-29_C_04.txt')
grace_pos = np.array(gnv[sorted(gnv.keys())[0]])
sv_plan = {int(r[0]): r for r in GPS_SV_PLAN}

# Geometric range check
t = sorted(gps1b.keys())[0]
obs = gps1b[t]
sv_id = list(obs.keys())[0]
utc = dt.datetime(1980,1,6) + dt.timedelta(seconds=t)
sp = satpos_from_sv(sv_plan.get(int(sv_id[1:]), None), utc)
sat_pos = np.array(sp[0])
rho = np.linalg.norm(sat_pos - grace_pos)
L1 = obs[sv_id]['L1']
L2 = obs[sv_id]['L2']

# Ionospheric-free combination
F1, F2 = 1575.42e6, 1227.60e6
ALPHA = F1*F1/(F1*F1-F2*F2)
BETA = -F2*F2/(F1*F1-F2*F2)
L_if = ALPHA*L1 + BETA*L2

print(f"GRACE pos: {grace_pos} (|GRACE| = {np.linalg.norm(grace_pos):.0f} m)")
print(f"GPS {sv_id} pos: {sat_pos} (|GPS| = {np.linalg.norm(sat_pos):.0f} m)")
print(f"Geometric range: {rho:.0f} m = {rho/1000:.1f} km")
print(f"L1={L1:.3f} = {L1/1000:.1f} km")
print(f"L2={L2:.3f} = {L2/1000:.1f} km")
print(f"L_if={L_if:.3f} = {L_if/1000:.1f} km")
print(f"L_if - rho = {L_if-rho:.0f} m = {(L_if-rho)/1000:.1f} km")
print(f"\nInitial residuals (should be absorbed by clock):")
print(f"  (L_if - rho) = {L_if-rho:.0f} m ({abs(L_if-rho)/1000:.1f} km)")
print(f"  This will be estimated as receiver clock offset")

# Try PPPProcessor with verbose
kf = PPPProcessor(grace_pos.tolist())
records = []
for t_key in sorted(gps1b.keys())[:10]:
    obs_dict = gps1b[t_key]
    utc_dt = dt.datetime(1980,1,6) + dt.timedelta(seconds=t_key)
    for sv, rec in obs_dict.items():
        if rec.get('L1') is None: continue
        sp = satpos_from_sv(sv_plan.get(int(sv[1:]), None), utc_dt)
        if sp is None or sp[0] is None: continue
        sat_pos, sv_v = sp
        records.append({
            'time': utc_dt, 'sv': sv,
            'L1': rec['L1'], 'L2': rec['L2'],
            'P1': rec['P1'], 'P2': rec['P2'],
            'sat_pos': np.array(sat_pos),
            'sat_clock': 0.0,
            'el': float(rec.get('el', 45.0)),
            'az': float(rec.get('az', 0.0)),
            'sat_vel': np.zeros(3),
            'trop_dry': 0.0, 'trop_wet': 0.0,
        })

results = kf.process(records, ref_pos=grace_pos)
d3s = [r['d3'] for r in results if r.get('d3') is not None and not np.isnan(r['d3'])]
print(f"\nPPP results: {len(d3s)}/{len(results)} valid")
if d3s:
    print(f"RMS 3D: {math.sqrt(sum(x**2 for x in d3s)/len(d3s))*100:.1f} cm")
    print(f"Range: {min(abs(x) for x in d3s):.3f} ~ {max(abs(x) for x in d3s):.3f} m")
    for r in results[:3]:
        print(f"  {r['time'].strftime('%H:%M')}: dE={r.get('dE',0):.4f} dN={r.get('dN',0):.4f} dU={r.get('dU',0):.4f} d3={r.get('d3',0):.4f}")
else:
    print("ALL NaN!")
    # Check solve_epoch details
    print("\nChecking KFPPP internal state...")
    print(f"KFPPP.pos0 = {kf.pos0}")
    print(f"KFPPP.apply_rel = {kf.apply_rel}")
    print(f"KFPPP.apply_apc = {kf.apply_apc}")
    print(f"KFPPP.apply_trop = {kf.apply_trop}")