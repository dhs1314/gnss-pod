#!/usr/bin/env python3
"""APR29 PPP Direct Test"""
from batch_v12 import parse_gnv1b
import pickle, sys, math, datetime as dt, numpy as np
sys.path.insert(0, '/workspace/gnss_pod/src')
from ppp import PPPProcessor
from run_ppp import satpos_from_sv

gnv = parse_gnv1b('/workspace/gnss_pod/data/gracefo/2024/2024-04-29/GNV1B_2024-04-29_C_04.txt')
gps1b = pickle.load(open('/workspace/gnss_pod/data/gracefo/2024/2024-04-29/gps1b_C.pkl', 'rb'))

keys_gnv = sorted(gnv.keys())
pos0 = gnv[keys_gnv[0]]
grace_pos = list(pos0)
print(f"GPS1B={len(gps1b)} epochs, GNV1B={len(gnv)} epochs")
print(f"GRACE pos: {grace_pos}")

keys = sorted(gps1b.keys())[:30]
records = []
for t in keys:
    obs_dict = gps1b[t]
    utc_dt = dt.datetime(1980,1,6) + dt.timedelta(seconds=t)
    for sv, rec in obs_dict.items():
        if rec.get('L1') is None: continue
        sp = satpos_from_sv(sv, utc_dt)
        if sp is None: continue
        sat_pos, sv_v = sp
        records.append({
            'time': utc_dt, 'sv': sv,
            'L1': rec['L1'], 'L2': rec['L2'],
            'P1': rec['P1'], 'P2': rec['P2'],
            'sat_pos': np.array(sat_pos),
            'sat_clock': 0.0,
            'el': rec.get('el', 45.0),
            'az': rec.get('az', 0.0),
            'sat_vel': np.array(sv_v) if sv_v else np.zeros(3),
            'trop_dry': 0.0, 'trop_wet': 0.0,
        })

print(f"PPP records: {len(records)}")

processor = PPPProcessor(grace_pos)
proc_results = processor.process(records, ref_pos=np.array(grace_pos))

print(f"\n结果: {len(proc_results)} epochs")
if proc_results:
    d3s = [r['d3'] for r in proc_results if r.get('d3') is not None]
    n = len(d3s)
    print(f"有效: {n}/{len(proc_results)}")
    if n > 0:
        rms_d3 = np.sqrt(np.nanmean([x**2 for x in d3s]))*100
        dEs = [r.get('dE', 0) for r in proc_results]
        dNs = [r.get('dN', 0) for r in proc_results]
        dUs = [r.get('dU', 0) for r in proc_results]
        rms_de = np.sqrt(np.nanmean([x**2 for x in dEs]))*100
        rms_dn = np.sqrt(np.nanmean([x**2 for x in dNs]))*100
        rms_du = np.sqrt(np.nanmean([x**2 for x in dUs]))*100
        print(f"RMS E/N/U/3D: {rms_de:.1f} / {rms_dn:.1f} / {rms_du:.1f} / {rms_d3:.1f} cm")
        print(f"d3范围: {min(abs(x) for x in d3s):.3f} ~ {max(abs(x) for x in d3s):.3f} m")
        print("前5个epoch:")
        for r in proc_results[:5]:
            t = r['time'].strftime('%H:%M:%S')
            print(f"  {t}: dE={r.get('dE',0):.4f} dN={r.get('dN',0):.4f} dU={r.get('dU',0):.4f} d3={r.get('d3',0):.4f} n={r.get('n_sat','?')}")
    else:
        print("ALL NaN!")
else:
    print("无结果!")