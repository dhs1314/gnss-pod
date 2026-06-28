#!/usr/bin/env python3
"""May2数据处理: 从TGZ提取→加载→PPP"""
import tarfile, glob, os, sys, math, pickle, shutil
sys.path.insert(0, '/workspace/gnss_pod/src')

WORK = '/workspace/gnss_pod'
GDIR = f'{WORK}/data/gracefo/2024'
ODIR = f'{WORK}/output_v12'
EXTRACT_DIR = f'{WORK}/data/gracefo/2024/2024-05-extracted'
os.makedirs(ODIR, exist_ok=True)

def extract_all():
    """提取TGZ所有内容"""
    tgz = f'{GDIR}/2024-05/TGZ/GRACE-FO_L1B_2024-05-04_SBG_v02.tgz'
    if not os.path.exists(tgz):
        print(f"TGZ不存在: {tgz}"); return False
    
    # 检查是否已提取
    existing = glob.glob(f'{EXTRACT_DIR}/**/*.txt', recursive=True)
    if existing:
        print(f"已提取到 {EXTRACT_DIR}, {len(existing)} txt文件")
        return True
    
    os.makedirs(EXTRACT_DIR, exist_ok=True)
    print(f"提取 {tgz}...")
    with tarfile.open(tgz, 'r:gz') as tar:
        members = tar.getmembers()
        txts = [m for m in members if m.name.endswith('.txt')]
        print(f"  {len(txts)} txt文件")
        for i, m in enumerate(txts):
            if i % 20 == 0:
                print(f"  提取中 {i}/{len(txts)}...")
            tar.extract(m, EXTRACT_DIR)
    print(f"  提取完成: {EXTRACT_DIR}")
    return True

def find_file(subdir, pattern):
    """在提取目录中找文件"""
    files = glob.glob(f'{EXTRACT_DIR}/**/{pattern}', recursive=True)
    return files[0] if files else None

def load_gps1b():
    """加载GPS1B"""
    pkl = f'{GDIR}/2024-05/gps1b_C.pkl'
    if os.path.exists(pkl):
        return pickle.load(open(pkl, 'rb'))
    
    txt = find_file(EXTRACT_DIR, 'GPS1B_2024-05-02_C_*.txt')
    if not txt:
        raise FileNotFoundError("找不到May2 GPS1B")
    print(f"从txt加载: {txt}")
    from gps1b_loader import load_gps1b_day
    data = load_gps1b_day('2024-05-02', 'C')
    pickle.dump(data, open(pkl, 'wb'))
    print(f"缓存: {pkl}")
    return data

def load_gnv1b():
    """加载GNV1B"""
    pkl = f'{GDIR}/2024-05/gnv1b_C.pkl'
    if os.path.exists(pkl):
        return pickle.load(open(pkl, 'rb'))
    
    txt = find_file(EXTRACT_DIR, 'GNV1B_2024-05-02_C_*.txt')
    if not txt:
        raise FileNotFoundError("找不到May2 GNV1B")
    print(f"从txt加载: {txt}")
    from gps1b_loader import load_gnv1b_day
    data = load_gnv1b_day('2024-05-02', 'C')
    pickle.dump(data, open(pkl, 'wb'))
    print(f"缓存: {pkl}")
    return data

def load_sp3():
    """加载SP3"""
    pkl = f'{GDIR}/2024-05/igs_sp3.pkl'
    if os.path.exists(pkl):
        return pickle.load(open(pkl, 'rb'))
    
    from gps1b_loader import load_sp3_orbits
    sp3 = load_sp3_orbits('2024-05-02', 'C', source='igs')
    if sp3:
        pickle.dump(sp3, open(pkl, 'wb'))
        print(f"SP3缓存: {pkl}")
    return sp3

def run_ppp(date_str, gps1b, gnv, sp3, strategy='broadcast'):
    from run_ppp import run_ppp as engine
    return engine(date_str, gps1b, gnv, strategy=strategy, ref_orbit=sp3)

def save_csv(results, out_path):
    import csv
    if not results: return
    keys = list(results[0].keys())
    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(results)

def main():
    print("=" * 60)
    print("处理May2 PPP")
    print("=" * 60)
    
    print("\n[1/5] 提取TGZ...")
    extract_all()
    
    print("\n[2/5] 加载GPS1B...")
    gps1b = load_gps1b()
    print(f"  GPS1B: {len(gps1b)} epochs")
    
    print("\n[3/5] 加载GNV1B...")
    gnv = load_gnv1b()
    print(f"  GNV1B: {len(gnv)} epochs")
    
    print("\n[4/5] 加载SP3...")
    sp3 = load_sp3()
    print(f"  SP3: {len(sp3) if sp3 else 0} epochs")
    
    print("\n[5/5] 运行PPP...")
    strategies = [('broadcast', None), ('sp3_final', sp3)]
    
    for strat, sp3_orbit in strategies:
        print(f"\n--- {strat} ---")
        results = run_ppp('2024-05-02', gps1b, gnv, sp3_orbit, strategy=strat)
        n = len(results)
        d3s = [r['d3_m'] for r in results if r.get('d3_m') is not None and not math.isnan(r['d3_m'])]
        nv = len(d3s)
        
        out_path = f'{ODIR}/ppp_2024-05-02_{strat}.csv'
        save_csv(results, out_path)
        
        if nv > 0:
            rms = math.sqrt(sum(x**2 for x in d3s)/nv)*100
            d3s_s = sorted(d3s)
            med = d3s_s[nv//2]*100
            pct1 = sum(1 for x in d3s if abs(x)<1)/nv*100
            pct10 = sum(1 for x in d3s if abs(x)>10)/nv*100
            pct100 = sum(1 for x in d3s if abs(x)>100)/nv*100
            maxv = max(abs(x) for x in d3s)
            print(f"  {nv}/{n} epochs, RMS={rms:.1f}cm ({rms/100:.2f}m), median={med:.1f}cm")
            print(f"    <1m={pct1:.0f}%, >10m={pct10:.0f}%, >100m={pct100:.0f}%, max={maxv:.0f}m")
        else:
            print(f"  ALL NaN ({n} epochs)")

if __name__ == '__main__':
    main()