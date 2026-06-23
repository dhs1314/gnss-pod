#!/usr/bin/env python3
"""
直接处理May2数据,绕过batch_v12.py的交互式input()
"""
import sys, os, math, pickle, tarfile, glob
sys.path.insert(0, '/workspace/gnss_pod/src')

GDIR = '/workspace/gnss_pod/data/gracefo/2024'

def extract_tgz(date_str, grace_id='C'):
    """手动提取tgz文件"""
    month = date_str[:7]  # 2024-05
    tgz_pattern = f'{GDIR}/{month}/*_TGZ*/*_TGZ*.tgz'
    tgz_files = glob.glob(tgz_pattern)
    if not tgz_files:
        # 搜索其他tgz
        tgz_files = glob.glob(f'{GDIR}/{month}/*.tgz')
    
    if not tgz_files:
        print(f"警告: 找不到 {month} 的tgz文件 (搜索: {tgz_pattern})")
        return False
    
    tgz = tgz_files[0]
    outdir = f'{GDIR}/{month}'
    print(f"提取 {tgz} → {outdir}...")
    
    try:
        with tarfile.open(tgz, 'r:gz') as tar:
            # 列出内容
            members = tar.getmembers()
            txt_files = [m for m in members if m.name.endswith('.txt')]
            print(f"  包含 {len(txt_files)} 个txt文件")
            
            # 检查目标文件是否已存在
            existing = sum(1 for m in txt_files if os.path.exists(os.path.join(outdir, os.path.basename(m.name))))
            if existing > 0:
                print(f"  已存在 {existing} 个文件,跳过提取")
                return True
            
            # 提取
            for member in txt_files:
                basename = os.path.basename(member.name)
                out_path = os.path.join(outdir, basename)
                if not os.path.exists(out_path):
                    member.name = basename  # 去掉路径前缀
                    tar.extract(member, outdir)
            print(f"  提取完成")
            return True
    except Exception as e:
        print(f"  提取失败: {e}")
        return False

def load_gps1b(date_str, grace_id='C'):
    """加载GPS1B数据"""
    pkl = f'{GDIR}/{date_str[:7]}/{date_str}/gps1b_{grace_id}.pkl'
    if os.path.exists(pkl):
        print(f"  从缓存加载: {pkl}")
        return pickle.load(open(pkl, 'rb'))
    
    # 直接从txt加载
    txt_pattern = f'{GDIR}/{date_str[:7]}/GPS1B_{date_str[:10]}_{grace_id}_*.txt'
    txt_files = glob.glob(txt_pattern)
    if not txt_files:
        raise FileNotFoundError(f"找不到GPS1B文件: {txt_pattern}")
    
    print(f"  从txt加载: {txt_files[0]}")
    from gps1b_loader import load_gps1b_day
    return load_gps1b_day(date_str, grace_id)

def load_gnv1b(date_str, grace_id='C'):
    """加载GNV1B数据"""
    pkl = f'{GDIR}/{date_str[:7]}/{date_str}/gnv1b_{grace_id}.pkl'
    if os.path.exists(pkl):
        print(f"  从缓存加载: {pkl}")
        return pickle.load(open(pkl, 'rb'))
    
    txt_pattern = f'{GDIR}/{date_str[:7]}/GNV1B_{date_str[:10]}_{grace_id}_*.txt'
    txt_files = glob.glob(txt_pattern)
    if not txt_files:
        raise FileNotFoundError(f"找不到GNV1B文件: {txt_pattern}")
    
    print(f"  从txt加载: {txt_files[0]}")
    from gps1b_loader import load_gnv1b_day
    return load_gnv1b_day(date_str, grace_id)

def load_sp3(date_str, grace_id='C'):
    """加载SP3精密星历"""
    pkl = f'{GDIR}/{date_str[:7]}/{date_str}/igs_sp3.pkl'
    if os.path.exists(pkl):
        print(f"  从缓存加载SP3: {pkl}")
        return pickle.load(open(pkl, 'rb'))
    
    from gps1b_loader import load_sp3_orbits
    sp3 = load_sp3_orbits(date_str, grace_id, source='igs')
    if sp3:
        # 保存缓存
        os.makedirs(os.path.dirname(pkl), exist_ok=True)
        pickle.dump(sp3, open(pkl, 'wb'))
        print(f"  SP3已缓存: {pkl}")
    return sp3

def run_ppp_direct(date_str, gps1b, gnv, sp3, strategy='broadcast'):
    """直接运行PPP,返回结果列表"""
    from run_ppp import run_ppp as ppp_engine
    return ppp_engine(date_str, gps1b, gnv, strategy=strategy, ref_orbit=sp3)

def save_csv(results, out_path):
    """保存结果到CSV"""
    import csv
    if not results: return
    keys = list(results[0].keys())
    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(results)
    print(f"  已保存: {out_path}")

def main():
    import glob
    
    date = '2024-05-02'
    grace = 'C'
    outdir = '/workspace/gnss_pod/output_v12'
    os.makedirs(outdir, exist_ok=True)
    
    print(f"=== 处理 {date} ===")
    
    # 1. 提取数据 (如果需要)
    print("\n[1/4] 提取数据...")
    extract_tgz(date, grace)
    
    # 2. 加载GPS1B
    print("\n[2/4] 加载GPS1B...")
    gps1b = load_gps1b(date, grace)
    print(f"  GPS1B: {len(gps1b)} epochs")
    
    # 3. 加载GNV1B
    print("\n[3/4] 加载GNV1B...")
    gnv = load_gnv1b(date, grace)
    print(f"  GNV1B: {len(gnv)} epochs")
    
    # 4. 处理
    print("\n[4/4] 运行PPP...")
    
    strategies = [
        ('broadcast', None),
        ('sp3_final', 'igs'),
    ]
    
    all_results = {}
    for strat_name, sp3_source in strategies:
        if sp3_source:
            print(f"\n  --- {strat_name} ---")
            sp3 = load_sp3(date, grace)
        else:
            sp3 = None
        
        results = run_ppp_direct(date, gps1b, gnv, sp3, strategy=strat_name)
        n = len(results)
        d3s = [r['d3_m'] for r in results if r.get('d3_m') is not None and not math.isnan(r['d3_m'])]
        nv = len(d3s)
        
        out_fn = f'ppp_{date[:10]}_{strat_name}.csv'
        out_path = f'{outdir}/{out_fn}'
        save_csv(results, out_path)
        
        if nv > 0:
            rms = math.sqrt(sum(x**2 for x in d3s)/nv)*100
            d3s_s = sorted(d3s)
            med = d3s_s[nv//2]*100
            pct1 = sum(1 for x in d3s if abs(x)<1)/nv*100
            pct10 = sum(1 for x in d3s if abs(x)>10)/nv*100
            pct100 = sum(1 for x in d3s if abs(x)>100)/nv*100
            maxv = max(abs(x) for x in d3s)
            print(f"  {strat_name}: {nv}/{n} epochs, RMS={rms:.1f}cm ({rms/100:.2f}m), median={med:.1f}cm")
            print(f"    <1m={pct1:.0f}%, >10m={pct10:.0f}%, >100m={pct100:.0f}%, max={maxv:.0f}m")
            all_results[strat_name] = {'n': nv, 'total': n, 'rms': rms, 'median': med, 'pct1': pct1, 'pct10': pct10, 'pct100': pct100, 'max': maxv}
        else:
            print(f"  {strat_name}: ALL NaN ({n} epochs)")
            all_results[strat_name] = {'n': 0, 'total': n, 'rms': float('nan')}
    
    print("\n=== 汇总 ===")
    for strat, r in all_results.items():
        if r['n'] > 0:
            print(f"{strat}: RMS={r['rms']:.1f}cm, median={r['median']:.1f}cm, <1m={r['pct1']:.0f}%, >100m={r['pct100']:.0f}%, max={r['max']:.0f}m")
        else:
            print(f"{strat}: ALL NaN")

if __name__ == '__main__':
    main()