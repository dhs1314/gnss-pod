#!/usr/bin/env python3
"""超简单测试: 验证May2 PPP是否正常工作"""
import pickle, sys, math, time
sys.path.insert(0, '/workspace/gnss_pod/src')

gps1b = pickle.load(open('data/gracefo/2024/2024-05-02/gps1b_C.pkl', 'rb'))
gnv = pickle.load(open('data/gracefo/2024/2024-05-02/gnv1b_C.pkl', 'rb'))
print(f"数据加载完成: GPS1B={len(gps1b)} epochs, GNV1B={len(gnv)} epochs")

# 检查第一个epoch的观测数量
keys = sorted(gps1b.keys())
first_key = keys[0]
obs_count = len(gps1b[first_key])
print(f"第一个epoch: t={first_key}, 观测数={obs_count}")

# 对比Apr29
apr = pickle.load(open('data/gracefo/2024/2024-04-29/gps1b_C.pkl', 'rb'))
apr_keys = sorted(apr.keys())
apr_obs = len(apr[apr_keys[0]])
print(f"Apr29第一个epoch: t={apr_keys[0]}, 观测数={apr_obs}")

# 只运行1个epoch,加超时
from run_ppp import run_ppp
one_key = {keys[0]: gps1b[keys[0]]}

print(f"\n运行1个epoch May2 PPP...")
t0 = time.time()
results = run_ppp('2024-05-02', one_key, gnv, strategy='broadcast', ref_orbit=None)
dt = time.time() - t0

d3s = [r['d3_m'] for r in results if r.get('d3_m') is not None and not math.isnan(r['d3_m'])]
print(f"1 epoch: {dt:.1f}s, {len(d3s)}/{len(results)} valid")

if d3s:
    rms = math.sqrt(sum(x**2 for x in d3s)/len(d3s))*100
    print(f"RMS={rms:.1f}cm ({rms/100:.2f}m)")
    print(f"→ 全天预计: {dt*len(gps1b)/60:.0f}分钟")
else:
    print("ALL NaN!")
    # 检查为什么
    first_obs = gps1b[keys[0]]
    sv0 = list(first_obs.keys())[0]
    d0 = first_obs[sv0]
    print(f"  第一个观测: SV={sv0}")
    print(f"  L1={d0.get('L1','N/A')}, L_if={d0.get('L_if','N/A')}, P1={d0.get('P1','N/A')}")
    print(f"  el={d0.get('el','N/A')}, az={d0.get('az','N/A')}")
    # 检查GNV1B GRACE位置
    gnv_keys = sorted(gnv.keys())
    gnv_first = gnv[gnv_keys[0]]
    if gnv_first:
        gp = gnv_first[0]
        print(f"  GRACE: X={gp.get('X',0):.1f}, Y={gp.get('Y',0):.1f}, Z={gp.get('Z',0):.1f}")
    print(f"  → PPP失败可能原因: 数据质量问题")