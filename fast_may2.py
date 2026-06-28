#!/usr/bin/env python3
"""快速PPP处理May2: 30s采样(原始60s)"""
import pickle, sys, math, time, os, csv
sys.path.insert(0, '/workspace/gnss_pod/src')

ODIR = '/workspace/gnss_pod/output_v12'
os.makedirs(ODIR, exist_ok=True)

def save_csv(results, out_path):
    if not results: return
    with open(out_path, 'w', newline='') as f:
        keys = list(results[0].keys())
        csv.DictWriter(f, fieldnames=keys).writeheader()
        csv.DictWriter(f, fieldnames=keys).writerows(results)

def main():
    may2_gps = 'data/gracefo/2024/2024-05-02/gps1b_C.pkl'
    may2_gnv = 'data/gracefo/2024/2024-05-02/gnv1b_C.pkl'
    may2_sp3 = 'data/gracefo/2024/2024-05-02/igs_sp3.pkl'
    
    print("加载数据...")
    gps1b = pickle.load(open(may2_gps, 'rb'))
    gnv = pickle.load(open(may2_gnv, 'rb'))
    sp3 = pickle.load(open(may2_sp3, 'rb')) if os.path.exists(may2_sp3) else None
    print(f"GPS1B: {len(gps1b)} epochs, GNV1B: {len(gnv)} epochs, SP3: {len(sp3) if sp3 else 0}")
    
    from run_ppp import run_ppp
    
    strategies = [('broadcast', None), ('sp3_final', sp3)]
    
    for strat, sp3_orbit in strategies:
        print(f"\n--- {strat} ---")
        t0 = time.time()
        results = run_ppp('2024-05-02', gps1b, gnv, strategy=strat, ref_orbit=sp3_orbit)
        dt = time.time() - t0
        
        n = len(results)
        d3s = [r['d3_m'] for r in results if r.get('d3_m') is not None and not math.isnan(r['d3_m'])]
        nv = len(d3s)
        
        out_path = f'{ODIR}/ppp_2024-05-02_{strat}.csv'
        save_csv(results, out_path)
        print(f"  {nv}/{n} epochs, {dt:.1f}s ({dt/n*1000:.0f}ms/epoch)")
        
        if nv > 0:
            rms = math.sqrt(sum(x**2 for x in d3s)/nv)*100
            d3s_s = sorted(d3s)
            med = d3s_s[nv//2]*100
            pct1 = sum(1 for x in d3s if abs(x)<1)/nv*100
            pct10 = sum(1 for x in d3s if abs(x)>10)/nv*100
            pct100 = sum(1 for x in d3s if abs(x)>100)/nv*100
            maxv = max(abs(x) for x in d3s)
            print(f"  RMS={rms:.1f}cm ({rms/100:.2f}m), median={med:.1f}cm")
            print(f"  <1m={pct1:.0f}%, >10m={pct10:.0f}%, >100m={pct100:.0f}%, max={maxv:.0f}m")
        else:
            print(f"  ALL NaN")

if __name__ == '__main__':
    main()