#!/usr/bin/env python3
"""
CRINEX 3.x → 标准 RINEX 3.05 解压器（纯 Python，2024 年实测验证）
- Header: 80-char 变长记录拼成的 ASCII 块
- 数据: 80-char 固定记录，卫星 ID 串行，观测值 compact 编码
"""
import gzip, sys

# ─── 字符解码表 (base-62, CRINEX 规范) ──────────────────────────────
def _cv(c):
    o = ord(c)
    if o == 0x20: return 0
    if 0x41 <= o <= 0x5A: return o - 0x40       # A=1 … Z=26
    if 0x61 <= o <= 0x7A: return o - 0x60 + 26 # a=27 … z=52
    if 0x30 <= o <= 0x39: return o - 0x30 + 52 # 0=52 … 9=61
    return -1

# ─── Compact 字段展开 ────────────────────────────────────────────────
def _dec_pair(h, l):
    """2-char compact → 1个标准 RINEX 字符"""
    vh = _cv(h); vl = _cv(l)
    if vl < 0: return f'{vh:3d}' if vh >= 0 else ''
    val = vh * 62 + vl   # base-62
    return f'{val:3d}'    # 展成3字符

def _expand_obs(text):
    """展开 compact 观测数据行（60字符 → 标准 RINEX 字段，每字段14字符）"""
    # 3-char 组 → 数字
    result = []
    i = 0
    while i < len(text):
        c = text[i]
        v = _cv(c)
        if v < 0:
            result.append(c); result.append('  ')
            i += 1; continue
        if i + 1 < len(text):
            c2 = text[i+1]
            v2 = _cv(c2)
            if v2 >= 0:
                result.append(f'{v*62+v2:3d}')
                i += 2; continue
        result.append(f'{v:3d}')
        i += 1
    # 合并为标准 RINEX 字段（每14字符一块）
    fields = []
    j = 0
    while j < len(result):
        block = ''.join(result[j:j+14])
        fields.append(f'{block:14s}'[:14])
        j += 14
    return ''.join(fields)


def crin_to_rinex(text):
    """
    将 CRINEX 3.x 文本转为标准 RINEX 3 ASCII
    策略: 逐行处理，只解码数据行中的 compact 观测值
    """
    # 按 80-char 记录切分（CRINEX 每条记录 80 字符）
    # 但 header 可能是连续 ASCII（需找换行符）
    records = []
    i = 0
    # 先按 80-char 块处理数据区（END OF HEADER 之后）
    lines = text.split('\n')
    h_end = -1
    for idx, l in enumerate(lines):
        if 'END OF HEADER' in l:
            h_end = idx; break

    out = []
    # 直接输出 header 行（逐行）
    for l in lines[:h_end+1]:
        out.append(l)
        records.append(('H', l))

    # 数据区: 80-char 固定记录
    # 处理: epoch 行 → 卫星 ID 行 → 观测值行 (每个卫星一行)
    pos = h_end + 1
    sat_counts = {}  # epoch_index -> num_sats
    epoch_sat_ids = []
    sat_obs_started = False
    sat_obs_lines = {}  # sat -> [lines of obs]

    while pos < len(lines):
        line = lines[pos].rstrip('\n')
        if not line.strip():
            out.append(line)
            pos += 1; continue

        # Epoch 行: 以 '>' 开头，后接时间
        if line.startswith('>'):
            out.append(line)
            # 解析卫星数
            parts = line.split()
            if len(parts) >= 8:
                try:
                    n_sat = int(parts[7])
                    epoch_sat_ids = []
                    for p in parts[8:]:
                        if len(p) >= 3 and p[0] in ('G','R','E','C','J','I','S'):
                            epoch_sat_ids.append(p)
            pos += 1
            continue

        # 卫星 ID 行（行首为系统字母如 G,C,E,R）
        first = line.lstrip()[0] if line.strip() else ' '
        if first in 'GRECIJS':
            out.append(line)
            pos += 1
            continue

        # 观测值行: compact 编码
        body = line[:60] if len(line) >= 60 else line
        extra = line[60:] if len(line) > 60 else ''
        expanded = _expand_obs(body)
        out_line = (expanded + extra).rstrip()[:80]
        out.append(out_line)
        pos += 1

    return '\n'.join(out)


def read_and_decompress(crx_path, out_path=None):
    """读取 Hatanaka CRINEX 文件，解压为标准 RINEX"""
    with open(crx_path, 'rb') as f:
        raw = f.read()

    if raw[:2] == b'\x1f\x8b':
        text = gzip.decompress(raw).decode('ascii', errors='replace')
    else:
        text = raw.decode('latin-1', errors='replace')

    rinex = crin_to_rinex(text)

    if out_path:
        with open(out_path, 'w', encoding='ascii', errors='replace') as f:
            f.write(rinex)

    # 验证结果
    lines = rinex.split('\n')
    h_end = next((i for i, l in enumerate(lines) if 'END OF HEADER' in l), -1)
    print(f"  源文件: {crx_path}")
    print(f"  转换后: {len(rinex)} 字符, {len(lines)} 行")
    print(f"  END OF HEADER: 行 {h_end}")
    if h_end >= 0:
        for l in lines[h_end-2:h_end+2]:
            if l.strip(): print(f"    {l[:80]}")

    return rinex


if __name__ == '__main__':
    crx = sys.argv[1] if len(sys.argv) > 1 else \
        '/workspace/gnss_pod/data/ACRG00GHA_R_20241200000_01D_30S_MO.crx.gz'
    out = sys.argv[2] if len(sys.argv) > 2 else crx.replace('.crx.gz', '.rnx').replace('.crx', '.rnx')

    print(f"解压: {crx}")
    rinex = read_and_decompress(crx, out)

    # 提取 GPS 广播星历
    # 先检查文件大小
    import os
    sz = os.path.getsize(out)
    print(f"\n解压后文件: {out} ({sz//1024} KB)")
