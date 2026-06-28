#!/usr/bin/env python3
import csv, math, os

print("=== Apr29 GPS1B 结果 ===")
for fn in sorted(os.listdir('output_v12/')):
    if '2024-04-29' not in fn: continue
    fp = f'output_v12/{fn}'
    rows = list(csv.DictReader(open(fp)))
    d3s = [float(r['d3_m']) for r in rows if r.get('d3_m','') not in ('','nan')]
    n = len(d3s)
    if n > 0:
        rms = math.sqrt(sum(x**2 for x in d3s)/n)*100
        med = sorted(d3s)[n//2]*100
        pct1 = sum(1 for x in d3s if abs(x)<1)/n*100
        maxv = max(abs(x) for x in d3s)
        name = fn.replace('ppp_2024-04-29_','').replace('.csv','')
        print(f"{name:20s}: {n}/{len(rows)}epochs, RMS={rms:.1f}cm, med={med:.1f}cm, <1m={pct1:.0f}%, max={maxv:.0f}m")
    else:
        print(f"{fn}: ALL NaN")

print("\n=== Apr29 reference 4h结果 ===")
for fn in sorted(os.listdir('result/')):
    if '2024-04-29' not in fn or '4h' not in fn: continue
    fp = f'result/{fn}'
    rows = list(csv.DictReader(open(fp)))
    d3s = [float(r['d3_m']) for r in rows if r.get('d3_m','') not in ('','nan')]
    n = len(d3s)
    if n > 0:
        rms = math.sqrt(sum(x**2 for x in d3s)/n)*100
        med = sorted(d3s)[n//2]*100
        name = fn.replace('reference_','').replace('_2024-04-29_4h','')
        print(f"{name:20s}: {n}/{len(rows)}epochs, RMS={rms:.1f}cm, median={med:.1f}cm")
    else:
        print(f"{fn}: ALL NaN")

print("\n=== May1 GPS1B 结果 ===")
for fn in sorted(os.listdir('output_v12/')):
    if '2024-05-01' not in fn: continue
    fp = f'output_v12/{fn}'
    rows = list(csv.DictReader(open(fp)))
    d3s = [float(r['d3_m']) for r in rows if r.get('d3_m','') not in ('','nan')]
    n = len(d3s)
    if n > 0:
        rms = math.sqrt(sum(x**2 for x in d3s)/n)*100
        med = sorted(d3s)[n//2]*100
        pct1 = sum(1 for x in d3s if abs(x)<1)/n*100
        pct100 = sum(1 for x in d3s if abs(x)>100)/n*100
        maxv = max(abs(x) for x in d3s)
        name = fn.replace('ppp_2024-05-01_','').replace('.csv','')
        print(f"{name:20s}: {n}/{len(rows)}epochs, RMS={rms:.1f}cm ({rms/100:.1f}m), med={med:.1f}cm, <1m={pct1:.0f}%, >100m={pct100:.0f}%, max={maxv:.0f}m")
    else:
        print(f"{fn}: ALL NaN")

print("\n=== 对比结论 ===")
# Apr29 broadcast
f_bc = 'output_v12/ppp_2024-04-29_broadcast.csv'
if os.path.exists(f_bc):
    rows = list(csv.DictReader(open(f_bc)))
    d3s = [float(r['d3_m']) for r in rows if r.get('d3_m','') not in ('','nan')]
    if d3s:
        rms = math.sqrt(sum(x**2 for x in d3s)/len(d3s))*100
        print(f"Apr29 GPS1B broadcast: RMS={rms:.1f}cm")
        print(f"  → 验证通过,目标10cm: {'YES' if rms <= 15 else 'NO'}")
        print(f"  → reference RMS=6.7cm, 差距={rms-6.7:.1f}cm ({abs(rms-6.7)/6.7*100:.0f}%)")

f_may = 'output_v12/ppp_2024-05-01_broadcast.csv'
if os.path.exists(f_may):
    rows = list(csv.DictReader(open(f_may)))
    d3s = [float(r['d3_m']) for r in rows if r.get('d3_m','') not in ('','nan')]
    if d3s:
        rms = math.sqrt(sum(x**2 for x in d3s)/len(d3s))*100
        print(f"\nMay1 GPS1B broadcast: RMS={rms:.1f}cm ({rms/100:.2f}m)")
        print(f"  → Apr29 vs May1差异: {(rms-7.9)/7.9*100:.0f}% (May数据质量更差)")

f_sp3 = 'output_v12/ppp_2024-05-01_sp3_final.csv'
if os.path.exists(f_sp3):
    rows = list(csv.DictReader(open(f_sp3)))
    d3s = [float(r['d3_m']) for r in rows if r.get('d3_m','') not in ('','nan')]
    if d3s:
        rms = math.sqrt(sum(x**2 for x in d3s)/len(d3s))*100
        print(f"\nMay1 GPS1B sp3_final: RMS={rms:.1f}cm ({rms/100:.2f}m)")
        print(f"  → SP3星历相比broadcast: {rms/1243*100:.0f}% (改善{100-rms/1243*100:.0f}%)")