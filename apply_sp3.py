#!/usr/bin/env python3
"""Apply SP3 integration to GRACE-FO PPP section in main.py"""
import sys

SRC = '/workspace/gnss_pod/main.py'
src = open(SRC).read()
lines = src.split('\n')

# Find module-level import section (after 'from pathlib import Path')
import_insert_line = None
for i, l in enumerate(lines):
    if l.strip().startswith('from pathlib import Path'):
        import_insert_line = i
        break

if import_insert_line is None:
    print("ERROR: can't find import line"); sys.exit(1)

import_code = """\
from src.sp3_loader import load_whu_sp3_ultra_rapid, get_gps_pos_from_sp3
"""

# Find nav_records = generate_broadcast_nav in main() GRACE-FO section
nav_rec_line = None
for i, l in enumerate(lines):
    if 'nav_records = generate_broadcast_nav(start_date)' in l and '    nav_records' in l:
        nav_rec_line = i; break

if nav_rec_line is None:
    print(f"ERROR: nav_records line not found"); sys.exit(1)

# Find the GPS position block (satpos_from_navrec usage)
gps_block_start = None
for i in range(nav_rec_line, min(nav_rec_line+20, len(lines)):
    if '用广播星历计算 GPS 卫星位置' in lines[i] or 'satpos_from_navrec' in lines[i]:
        gps_block_start = i; break

if gps_block_start is None:
    print("ERROR: GPS block not found"); sys.exit(1)

print(f"Import at line {import_insert_line+1}")
print(f"nav_records at line {nav_rec_line+1}")
print(f"GPS block at line {gps_block_start+1}")

# Find end of GPS block (the r['sat_clock'] = 0.0 line)
gps_block_end = None
for i in range(gps_block_start, min(gps_block_start+15, len(lines)):
    if "r['sat_clock'] = 0.0" in lines[i] or "r['sat_pos'] = gps_pos" in lines[i]:
        gps_block_end = i; break

if gps_block_end is None:
    print("ERROR: GPS block end not found"); sys.exit(1)

print(f"GPS block: lines {gps_block_start+1}-{gps_block_end+1}")

# Build new lines
new_lines = []
new_lines += lines[:import_insert_line+1]  # up to and including import line
new_lines.append(import_code)
new_lines += lines[import_insert_line+1:nav_rec_line+1]  # between import and nav_records
new_lines.append("        sp3_data = load_whu_sp3_ultra_rapid(year, doy, str(data_dir))")
new_lines.append("")

# New GPS block
new_gps_block = """\
            sv = r['sv']
            # GPS 卫星位置（优先 WHU SP3 ultra-rapid，否则广播星历）
            if sp3_data is not None:
                gps_pos, gps_clk = get_gps_pos_from_sp3(sp3_data, sv, t)
                if gps_pos is None: continue
            else:
                nav_rec = next((n for n in nav_records
                               if n['sv'] == sv and
                               abs((n['time'] - t).total_seconds()) < 7200), None)
                if nav_rec is None: continue
                gps_pos = np.array(satpos_from_navrec(nav_rec, t)); gps_clk = 0.0
            r['sat_pos'] = gps_pos
            r['sat_clock'] = gps_clk"""

new_lines.append(new_gps_block)
new_lines.append("")
new_lines += lines[gps_block_end+1:]  # rest of file

new_src = '\n'.join(new_lines)

# Verify syntax
import ast
try:
    ast.parse(new_src)
    print("SYNTAX OK!")
except SyntaxError as e:
    print(f"SYNTAX ERROR at line {e.lineno}: {e.msg}")
    new_lines2 = new_src.split('\n')
    for i in range(max(0,e.lineno-3), min(len(new_lines2), e.lineno+3)):
        print(f"{i+1}: {repr(new_lines2[i][:100]})")
    sys.exit(1)

open(SRC, 'w').write(new_src)
print("Written!")
