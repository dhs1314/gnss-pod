#!/usr/bin/env python3
"""直接读取GPS1B文件,找到数据行格式"""
fname = '/workspace/gnss_pod/data/GPS1B_2024-04-29_C_04.txt'
with open(fname, 'r') as f:
    lines = f.readlines()

# 打印第130-150行
print("=== 行130-150 ===")
for i in range(129, 150):
    print(f"[{i}] {lines[i].rstrip()[:180]}")

# 找prod_flag字段
print("\n=== 找数据行格式 ===")
for i, line in enumerate(lines):
    if i < 130: continue
    parts = line.strip().split()
    if len(parts) >= 8:
        print(f"行{i}: nfields={len(parts)}, 前8字段: {parts[:8]}")
        if i > 140: break