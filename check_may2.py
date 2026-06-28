#!/usr/bin/env python3
import csv, math, os, sys

print('=== May2 PPP结果 ===', flush=True)
for fn in sorted(os.listdir('output_v12/')):
    if '2024-05-02' not in fn: continue
    fp = f'output_v12/{fn}'
    rows = list(csv.DictReader(open(fp)))
    d3s = [float(r['d3_m']) for r in rows if r.get('d3_m','') not in ('','nan')]
    n = len(d3s)
    if n > 0:
        rms = math.sqrt(sum(x**2 for x in d3s)/n)*100
        d3s_s = sorted(d3s)
        med = d3s_s[n//2]*100
        pct1 = sum(1 for x in d3s if abs(x)<1)/n*100
        pct10 = sum(1 for x in d3s if abs(x)>10)/n*100
        pct100 = sum(1 for x in d3s if abs(x)>100)/n*100
        maxv = max(abs(x) for x in d3s)
        name = fn.replace('ppp_2024-05-02_','').replace('.csv','')
        print(f'{name}: {n}/{len(rows)}, RMS={rms:.1f}cm ({rms/100:.2f}m)')
        print(f'  median={med:.1f}cm, <1m={pct1:.0f}%, >10m={pct10:.0f}%, >100m={pct100:.0f}%, max={maxv:.0f}m')
    else:
        print(f'{fn}: ALL NaN ({len(rows)}行)')

# 对比Apr29和May1/2
print('\n=== Apr29 vs May1 vs May2 对比 ===', flush=True)
for date in ['2024-04-29', '2024-05-01', '2024-05-02']:
    for strat in ['broadcast', 'sp3_final']:
        fn = f'ppp_{date}_{strat}.csv'
        fp = f'output_v12/{fn}'
        if not os.path.exists(fp): continue
        rows = list(csv.DictReader(open(fp)))
        d3s = [float(r['d3_m']) for r in rows if r.get('d3_m','') not in ('','nan')]
        n = len(d3s)
        if n > 0:
            rms = math.sqrt(sum(x**2 for x in d3s)/n)*100
            d3s_s = sorted(d3s)
            med = d3s_s[n//2]*100
            pct1 = sum(1 for x in d3s if abs(x)<1)/n*100
            pct100 = sum(1 for x in d3s if abs(x)>100)/n*100
            maxv = max(abs(x) for x in d3s)
            print(f'{date} {strat}: {n}/{len(rows)}, RMS={rms:.1f}cm, med={med:.1f}cm, <1m={pct1:.0f}%, >100m={pct100:.0f}%, max={maxv:.0f}m')
        else:
            print(f'{date} {strat}: ALL NaN')