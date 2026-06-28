#!/usr/bin/env python3
import pickle
sp3 = pickle.load(open('data/2024/122/igs_sp3_FIN.pkl', 'rb'))
print('Keys:', list(sp3.keys()))
ts = sp3['ts']
print(f'ts type: {type(ts[0])}, first 3: {ts[:3]}')

t0 = ts[0]
# Check access patterns
if 'data' in sp3:
    print('Has "data" key')
    dk = list(sp3['data'].keys())[:3]
    print(f'First data keys: {dk}')
elif 'sat' in sp3:
    print('Has "sat" key')
else:
    # Try datetime key
    if t0 in sp3:
        print(f't0 is direct key, value type: {type(sp3[t0])}')
        if isinstance(sp3[t0], dict):
            k = list(sp3[t0].keys())
            print(f'SVs at t0: {k[:5]}...')
            sv0 = k[0]
            val = sp3[t0][sv0]
            print(f'{sv0}: {val}')
    else:
        print(f't0={t0} NOT a key')
        # Try all possible keys
        other_keys = [k for k in sp3.keys() if k != 'ts']
        print(f'Other keys: {other_keys[:5]}')
        if other_keys:
            k0 = other_keys[0]
            print(f'  key type: {type(k0)}')
            v0 = sp3[k0]
            print(f'  value type: {type(v0)}')
            if isinstance(v0, dict):
                print(f'  inner keys: {list(v0.keys())[:5]}')

# Check SP3 clock values for a few SVs
print("\n--- Testing get_gps_pos_from_sp3 ---")
from src.sp3_loader import get_gps_pos_from_sp3
from datetime import datetime
utc = datetime(2024, 5, 1, 0, 0, 0)
for sv in ['G03', 'G04', 'G06']:
    pos, clk, vel = get_gps_pos_from_sp3(sp3, sv, utc)
    if pos is not None:
        print(f'{sv}: pos=[{pos[0]/1000:.1f},{pos[1]/1000:.1f},{pos[2]/1000:.1f}]km clk={clk:.1f}m')
    else:
        print(f'{sv}: NO DATA')
