#!/usr/bin/env python3
import pickle
data = pickle.load(open('data/gracefo/2024/2024-04-29/GPS1B_2024-04-29_C_04.pkl', 'rb'))
keys = sorted(data.keys())
print(f'N epochs: {len(keys)}')
print(f'First 5 keys: {keys[:5]}')
print(f'Last 5 keys: {keys[-5:]}')
print(f'Key type: {type(keys[0])}')
print()
k0 = keys[0]
print(f'First epoch keys: {sorted(data[k0].keys())}')
for sv in sorted(data[k0].keys())[:2]:
    rec = data[k0][sv]
    for k, v in rec.items():
        if isinstance(v, float):
            print(f'  {sv}.{k} = {v:.4f}')
        else:
            print(f'  {sv}.{k} = {v}')
