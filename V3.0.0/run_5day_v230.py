"""V2.3.0 5-day validation — EKF orbit + Batch solver phase RMS.
Prints 3 columns per arc: EKF 3D RMS, EKF Phase RMS, Batch Phase RMS.
"""
import sys, os, pickle, subprocess, gzip, shutil
from pathlib import Path
from datetime import datetime
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).parent
DATA = ROOT / 'data'
CODE_DIR = DATA / 'CODE' / '2024'
CODE_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR = ROOT / 'results' / '5day_v230'
OUT_DIR.mkdir(parents=True, exist_ok=True)

DATES = ['2024-04-29', '2024-04-30', '2024-05-01', '2024-05-02', '2024-05-03']
GRACE = 'C'; INTERVAL = 30
ARC_HOURS = [0.17, 0.50]
ANTEX_FILE = DATA / "igs14.atx"
IERS_FILE = DATA / "IERS" / "eopc04_IAU2000.txt"

print("=" * 70)
print(f"V2.3.0 5-Day — EKF Orbit + Batch Phase RMS")
print("=" * 70)

# Step 1: Download CODE products
for d in DATES:
    dt = datetime.strptime(d, "%Y-%m-%d")
    doy = dt.strftime("%j"); y = dt.year; m = dt.month
    pref = f"COD0OPSFIN_{y}{doy}0000_01D"
    for fn_gz, fn_out in [
        (f"{pref}_05M_ORB.SP3.gz", f"{pref}_05M_ORB.SP3"),
        (f"{pref}_30S_CLK.CLK.gz", f"{pref}_30S_CLK.CLK"),
    ]:
        out_path = CODE_DIR / fn_out; gz_path = CODE_DIR / fn_gz
        if out_path.exists() or gz_path.exists(): continue
        try:
            url = f"http://ftp.aiub.unibe.ch/CODE/{y}/{fn_gz}"
            import urllib.request
            urllib.request.urlretrieve(url, gz_path)
        except: pass
    for dcb_fn in [f"P1P2{y%100:02d}{m:02d}.DCB", f"P1C1{y%100:02d}{m:02d}.DCB"]:
        dcb_p = CODE_DIR / dcb_fn
        if dcb_p.exists(): continue
        try:
            import urllib.request
            urllib.request.urlretrieve(f"http://ftp.aiub.unibe.ch/CODE/{y}/{dcb_fn}", dcb_p)
        except: pass

# Decompress
for f in CODE_DIR.glob("*.gz"):
    out = CODE_DIR / f.stem
    if not out.exists():
        with gzip.open(f) as gf, open(out, 'wb') as of:
            shutil.copyfileobj(gf, of)

print("Downloads: done")
print()

# Step 2: Run EKF + Batch on all dates
all_results = []

for date_str in DATES:
    dt = datetime.strptime(date_str, "%Y-%m-%d")
    doy = dt.strftime("%j"); y = dt.year; m = dt.month
    pref = f"COD0OPSFIN_{y}{doy}0000_01D"
    sp3_tgt = CODE_DIR / f"{pref}_05M_ORB.SP3"
    clk_tgt = CODE_DIR / f"{pref}_30S_CLK.CLK"
    dcbm = CODE_DIR / f"P1P2{y%100:02d}{m:02d}.DCB"
    dcbc = CODE_DIR / f"P1C1{y%100:02d}{m:02d}.DCB"
    if not dcbm.exists(): dcbm = CODE_DIR / "P1P22404.DCB"
    if not dcbc.exists(): dcbc = CODE_DIR / "P1C12404.DCB"

    gps1b = DATA / f"GPS1B_{date_str}_{GRACE}_04.pkl"
    if not gps1b.exists():
        gps1b = DATA / f"gracefo/{y}/{date_str}/GPS1B_{date_str}_{GRACE}_04.pkl"
    if not sp3_tgt.exists() or not clk_tgt.exists() or not gps1b.exists():
        continue

    for arc_h in ARC_HOURS:
        chi2 = '25' if arc_h < 0.3 else '100'
        cmd = [
            sys.executable, str(ROOT / 'run_sequential_pod.py'),
            '--date', date_str, '--hours', str(arc_h),
            '--interval', str(INTERVAL), '--grace-id', GRACE,
            '--dynamics-mode', 'simplified',
            '--sp3-file', str(sp3_tgt),
            '--clk-file', str(clk_tgt),
            '--dcb-file', str(dcbm),
            '--dcb-p1c1-file', str(dcbc),
            '--antex-file', str(ANTEX_FILE),
            '--iers-c04', str(IERS_FILE),
            '--enable-phase-windup', '--enable-relativity',
            '--ar-min-epochs', '6', '--gravity-nmax', '90',
            '--chi2-threshold', chi2,
            '--batch-lsq-v2',  # Framework v3 Batch Solver
        ]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=600,
                env={**os.environ,
                     'JAVA_HOME': r'C:\Program Files\JetBrains\PyCharm Community Edition 2024.3.5\jbr',
                     'OREKIT_DATA_PATH': str(DATA / 'orekit')})

            rms_3d = phase_ekf = phase_batch = code_rms = rej = nsv = 0
            for line in r.stdout.split('\n'):
                if '3D RMS vs GNV1B' in line:
                    rms_3d = float(line.split()[4])
                if 'Phase:' in line and 'RMS=' in line and 'accepted' in line:
                    phase_ekf = float(line.split('RMS=')[1].split('m')[0])
                if 'Code:' in line and 'accepted' in line:
                    code_rms = float(line.split('RMS=')[1].strip().split('m')[0])
                if 'Rejected:' in line and 'Final' not in line:
                    try: rej = int(line.split()[-1])
                    except: pass
                if 'Final SV count:' in line:
                    nsv = int(line.split()[-1])
                if 'Phase RMS:  EKF=' in line:
                    # "Phase RMS:  EKF=0.276m  Batch=0.160m"
                    parts = line.replace('m','').split()
                    for j, p in enumerate(parts):
                        if p.startswith('EKF='): phase_ekf = float(p[4:])
                        if p.startswith('Batch='): phase_batch = float(p[6:])

            print(f"  {date_str} {arc_h}h | RMS={rms_3d:.3f}m | "
                  f"Ph_EKF={phase_ekf:.3f}m Ph_B={phase_batch:.3f}m | "
                  f"SVs={nsv}")
            all_results.append({
                'date': date_str, 'arc_h': arc_h,
                'rms_3d': rms_3d,
                'phase_ekf': phase_ekf,
                'phase_batch': phase_batch,
                'code': code_rms, 'rej': rej, 'n_sv': nsv,
            })
        except Exception as e:
            print(f"    ERROR: {e}")

if not all_results:
    print("No results."); sys.exit(1)

# Step 3: Summary
print(f"\n{'='*70}")
print(f"V2.3.0 5-Day Summary ({len(all_results)} arcs)")
print(f"{'='*70}")

for arc_h in ARC_HOURS:
    subset = [r for r in all_results if r['arc_h'] == arc_h]
    if not subset: continue
    rms = [r['rms_3d'] for r in subset]
    pe = [r['phase_ekf'] for r in subset]
    pb = [r['phase_batch'] for r in subset]
    imp = [(e - b)/e*100 for e, b in zip(pe, pb) if e > 0 and b > 0]
    print(f"\n{arc_h}h ({len(subset)} dates):")
    print(f"  3D RMS:   mean={np.mean(rms):.3f}m  median={np.median(rms):.3f}m  "
          f"best={np.min(rms):.3f}m  worst={np.max(rms):.3f}m")
    print(f"  EKF Phase: mean={np.mean(pe):.3f}m  Batch Phase: mean={np.mean(pb):.3f}m")
    if imp:
        print(f"  Batch improvement: mean={np.mean(imp):+.1f}%  best={np.max(imp):+.1f}%  worst={np.min(imp):+.1f}%")

# Table
print(f"\n{'Date':<12} {'Arc':>5} {'3D_RMS':>7} {'Ph_EKF':>7} {'Ph_Bch':>7} {'B_Impr':>6} {'SVs':>4}")
for r in all_results:
    imp = (r['phase_ekf'] - r['phase_batch'])/r['phase_ekf']*100 if r['phase_ekf'] > 0 and r['phase_batch'] > 0 else 0
    print(f"{r['date']:<12} {r['arc_h']:>5.2f}h {r['rms_3d']:>7.3f} {r['phase_ekf']:>7.3f} {r['phase_batch']:>7.3f} {imp:>+6.1f}% {r['n_sv']:>4}")

# Plot: 3 bars per arc date
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8), gridspec_kw={'height_ratios': [2, 1]})

w = 0.25
for i, arc_h in enumerate(ARC_HOURS):
    subset = [r for r in all_results if r['arc_h'] == arc_h]
    x = np.arange(len(subset)) + i * (len(subset) + 0.5)
    rms_v = [r['rms_3d'] for r in subset]
    pb_v = [r['phase_batch'] for r in subset]
    pe_v = [r['phase_ekf'] for r in subset]
    labels = [r['date'][5:] for r in subset]

    ax1.bar(x - w, rms_v, w, color='#607D8B', alpha=0.8, label=f'{arc_h}h 3D RMS' if i == 0 else '')
    ax1.bar(x, pe_v, w, color='#FF9800', alpha=0.8, label=f'{arc_h}h EKF Phase' if i == 0 else '')
    ax1.bar(x + w, pb_v, w, color='#2196F3', alpha=0.8, label=f'{arc_h}h Batch Phase' if i == 0 else '')

    for xi, rv, pv in zip(x, rms_v, pb_v):
        imp = (r['phase_ekf'] - r['phase_batch'])/r['phase_ekf']*100 if r['phase_ekf'] > 0 and r['phase_batch'] > 0 else 0
        ax2.bar(xi + i*0.05, [imp], 0.25, color='#4CAF50' if imp > 0 else '#F44336', alpha=0.7)

ax1.set_ylabel('RMS [m]', fontsize=12)
ax1.set_title(f'V2.3.0 5-Day — EKF Orbit + Batch Solver Phase RMS', fontsize=13)
ax1.legend(fontsize=8, loc='upper left', ncol=3)
ax1.grid(True, alpha=0.3)
all_x = sorted(set(r['date'][5:] for r in all_results))
ax1.set_xticks(range(len(all_x)))
ax1.set_xticklabels(all_x, fontsize=9)

ax2.set_ylabel('Phase Improvement [%]', fontsize=11)
ax2.set_xlabel('Date', fontsize=11)
ax2.axhline(y=0, color='gray', linewidth=0.5)
ax2.grid(True, alpha=0.3)
ax2.set_xticks(range(len(all_x)))
ax2.set_xticklabels(all_x, fontsize=9)

plt.tight_layout()
png_path = OUT_DIR / '5day_v230.png'
plt.savefig(str(png_path), dpi=120)
pickle.dump(all_results, open(str(OUT_DIR / '5day_v230.pkl'), 'wb'))
print(f"\n  Report: {png_path}")
print(f"  Data:   {OUT_DIR / '5day_v230.pkl'}")
