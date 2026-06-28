#!/usr/bin/env python3
import pickle
f1 = 'data/gracefo/2024/2024-05-01/GPS1B_2024-05-01_C_04.pkl'
f2 = 'data/gracefo/2024/2024-05-01/gps1b_C.pkl'
d1 = pickle.load(open(f1, 'rb'))
d2 = pickle.load(open(f2, 'rb'))
print(f'File 1: {len(d1)} epochs')
print(f'File 2: {len(d2)} epochs')
k1 = sorted(d1.keys())
k2 = sorted(d2.keys())
print(f'File1 range: {k1[0]} to {k1[-1]}')
print(f'File2 range: {k2[0]} to {k2[-1]}')

ce = sorted(set(k1) & set(k2))[0]
e1 = d1[ce]; e2 = d2[ce]
svs1 = sorted(e1.keys()); svs2 = sorted(e2.keys())
print(f'Epoch {ce}:')
print(f'  File1 SVs: {svs1}')
print(f'  File2 SVs: {svs2}')

# Compare fields for common SV
common_svs = sorted(set(svs1) & set(svs2))
if common_svs:
    sv = common_svs[0]
    r1 = e1[sv]; r2 = e2[sv]
    print(f'\n  {sv} fields comparison:')
    all_keys = sorted(set(list(r1.keys()) + list(r2.keys())))
    for k in all_keys:
        v1 = r1.get(k, 'N/A')
        v2 = r2.get(k, 'N/A')
        if isinstance(v1, float) and isinstance(v2, float):
            diff = v1 - v2
            flag = ' ***' if abs(diff) > 0.01 else ''
        else:
            diff = 'N/A'
            flag = ''
        print(f'    {k}: File1={v1}  File2={v2}  diff={diff}{flag}')
