#!/usr/bin/env python3
"""APR29 PPP 测试 - 直接用PPPProcessor"""
import pickle, sys, math, datetime as dt
sys.path.insert(0, '/workspace/gnss_pod/src')

from ppp import PPPProcessor

# 加载数据
gps1b = pickle.load(open('/workspace/gnss_pod/data/gracefo/2024/2024-04-29/gps1b_C.pkl', 'rb'))
gnv_txt = open('/workspace/gnss_pod/data/gracefo/2024/2024-04-29/GNV1B_2024-04-29_C_04.txt').readlines()
sp3_pkl = '/workspace/gnss_pod/data/gracefo/2024/2024-04-29/igs_sp3.pkl'
sp3 = pickle.load(open(sp3_pkl, 'rb')) if __import__('os').path.exists(sp3_pkl) else None

print(f"GPS1B: {len(gps1b)} epochs")
print(f"GNV1B: {len(gnv_txt)} lines")
print(f"SP3: {len(sp3) if sp3 else 'None'} epochs")

# 解析GRACE位置
grace_pos = None
for line in gnv_txt:
    if line.startswith('#') or not line.strip(): continue
    parts = line.split()
    if len(parts) >= 8 and parts[0] not in ('#', 'HEADER', 'COMMENT'):
        try:
            gpX = float(parts[1]); gpY = float(parts[2]); gpZ = float(parts[3])
            if abs(gpX) > 1e6 and abs(gpY) > 1e6:
                grace_pos = [gpX, gpY, gpZ]
                print(f"GRACE C位置: X={gpX:.1f} Y={gpY:.1f} Z={gpZ:.1f}")
                break
        except: pass

if grace_pos is None:
    print("无法读取GRACE位置!")
    exit(1)

# 准备PPP记录 - 用broadcast星历
from run_ppp import satpos_from_sv, compute_relativityCorr, sagnac_correction

keys = sorted(gps1b.keys())[:30]
records = []

for t in keys:
    obs_dict = gps1b[t]
    utc = dt.datetime(1980,1,6) + dt.timedelta(seconds=t)
    
    for sv, rec in obs_dict.items():
        if rec.get('L1') is None: continue
        sat_pos = satpos_from_sv(sv, utc)
        if sat_pos is None: continue
        sp, sv_v = sat_pos
        
        # 天线PCV改正
        el_deg = rec.get('el', 45.0)
        el_rad = math.radians(el_deg)
        
        records.append({
            'time': utc,
            'sv': sv,
            'L1': rec['L1'], 'L2': rec['L2'],
            'P1': rec['P1'], 'P2': rec['P2'],
            'sat_pos': sp,
            'sat_clock': 0.0,
            'el': el_deg,
            'az': rec.get('az', 0.0),
            'sat_vel': sv_v,
            'trop_dry': 0.0,
            'trop_wet': 0.0,
        })

print(f"PPP records: {len(records)} obs")

# 运行PPP
processor = PPPProcessor(grace_pos)
proc_results = processor.process(records, ref_pos=grace_pos)

print(f"\nPPP结果: {len(proc_results)} epochs")
if proc_results:
    d3s = [r['d3'] for r in proc_results if r.get('d3') is not None]
    dEs = [r['dE'] for r in proc_results if r.get('dE') is not None]
    dNs = [r['dN'] for r in proc_results if r.get('dN') is not None]
    dUs = [r['dU'] for r in proc_results if r.get('dU') is not None]
    n = len(d3s)
    print(f"  有效: {n}/{len(proc_results)}")
    if n > 0:
        rms_d3 = math.sqrt(sum(x**2 for x in d3s)/n)*100
        rms_de = math.sqrt(sum(x**2 for x in dEs)/n)*100
        rms_dn = math.sqrt(sum(x**2 for x in dNs)/n)*100
        rms_du = math.sqrt(sum(x**2 for x in dUs)/n)*100
        print(f"  RMS E/N/U/3D: {rms_de:.1f} / {rms_dn:.1f} / {rms_du:.1f} / {rms_d3:.1f} cm")
        print(f"  d3范围: {min(abs(x) for x in d3s):.3f} - {max(abs(x) for x in d3s):.3f} m")
        print(f"  前5个epoch:")
        for r in proc_results[:5]:
            t = r['time'].strftime('%H:%M:%S')
            print(f"    {t}: dE={r.get('dE',0):.4f} dN={r.get('dN',0):.4f} dU={r.get('dU',0):.4f} d3={r.get('d3',0):.4f} n={r.get('n_sat','?')}")
else:
    print("  ALL NaN!")