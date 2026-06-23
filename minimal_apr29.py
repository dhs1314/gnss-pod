#!/usr/bin/env python3
"""Minimal APR29 PPP test"""
import pickle, sys, math, datetime as dt, numpy as np
sys.path.insert(0, '/workspace/gnss_pod/src')
from batch_v12 import parse_gnv1b
from run_ppp import GPS_SV_PLAN, satpos_from_sv
from ppp import PPPProcessor

# Load data
gps1b = pickle.load(open('/workspace/gnss_pod/data/gracefo/2024/2024-04-29/gps1b_C.pkl', 'rb'))
gnv = parse_gnv1b('/workspace/gnss_pod/data/gracefo/2024/2024-04-29/GNV1B_2024-04-29_C_04.txt')

# GRACE position as numpy array
ts = sorted(gnv.keys())
grace_pos = np.array(gnv[ts[0]])  # shape (3,)
print(f"GPS1B={len(gps1b)}, GNV1B={len(gnv)}")
print(f"GRACE pos (numpy): {grace_pos}")
print(f"GRACE pos type: {grace_pos.dtype}, shape: {grace_pos.shape}")

# Build SV lookup
sv_plan = {int(r[0]): r for r in GPS_SV_PLAN}

def get_sat_pos(sv_id, utc):
    prn = int(sv_id[1:])
    if prn not in sv_plan: return None, None
    return satpos_from_sv(sv_plan[prn], utc)

# Prepare records
keys = sorted(gps1b.keys())[:50]
records = []
for t in keys:
    obs = gps1b[t]
    utc = dt.datetime(1980,1,6) + dt.timedelta(seconds=t)
    for sv_id, rec in obs.items():
        if rec.get('L1') is None: continue
        sp = get_sat_pos(sv_id, utc)
        if sp[0] is None: continue
        sat_pos, sv_v = sp
        records.append({
            'time': utc, 'sv': sv_id,
            'L1': rec['L1'], 'L2': rec['L2'],
            'P1': rec['P1'], 'P2': rec['P2'],
            'sat_pos': np.array(sat_pos),
            'sat_clock': 0.0,
            'el': float(rec.get('el', 45.0)),
            'az': float(rec.get('az', 0.0)),
            'sat_vel': np.zeros(3),
            'trop_dry': 0.0, 'trop_wet': 0.0,
        })

print(f"Records: {len(records)}")

# Run PPP
kf = PPPProcessor(grace_pos.tolist())
results = kf.process(records, ref_pos=grace_pos)

print(f"\nPPP epochs: {len(results)}")
if results:
    d3s = [r['d3'] for r in results if r.get('d3') is not None and not np.isnan(r['d3'])]
    print(f"Valid: {len(d3s)}/{len(results)}")
    if d3s:
        rms = math.sqrt(sum(x**2 for x in d3s)/len(d3s))*100
        dEs = [r.get('dE',0) for r in results]
        dNs = [r.get('dN',0) for r in results]
        dUs = [r.get('dU',0) for r in results]
        rms_e = math.sqrt(sum(x**2 for x in dEs)/len(dEs))*100
        rms_n = math.sqrt(sum(x**2 for x in dNs)/len(dNs))*100
        rms_u = math.sqrt(sum(x**2 for x in dUs)/len(dUs))*100
        print(f"RMS E/N/U/3D: {rms_e:.1f} / {rms_n:.1f} / {rms_u:.1f} / {rms:.1f} cm")
        print(f"d3 range: {min(abs(x) for x in d3s):.3f} ~ {max(abs(x) for x in d3s):.3f} m")
        print("First 5:")
        for r in results[:5]:
            t = r['time'].strftime('%H:%M:%S')
            print(f"  {t}: dE={r.get('dE',0):.4f} dN={r.get('dN',0):.4f} dU={r.get('dU',0):.4f} d3={r.get('d3',0):.4f}")
    else:
        print("ALL NaN d3!")
else:
    print("No results!")