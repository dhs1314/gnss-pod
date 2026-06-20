"""Multi-date validation — loops run_sequential_pod.py for each date.

Reuses the existing fully-tested pipeline: no new code for data loading.

Usage:
    py -3.12 run_multiday_v2.py --dates 2024-04-29 2024-04-30 2024-05-01 \
      --hours 0.17,0.5 --interval 30 --grace-id C
"""
import sys, os, pickle, subprocess
from pathlib import Path
from datetime import datetime
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).parent
DATA = ROOT / 'data'

# ── Parse ──
import argparse
p = argparse.ArgumentParser()
p.add_argument('--dates', nargs='+', required=True, help='YYYY-MM-DD ...')
p.add_argument('--hours', default='0.17,0.5', help='comma-separated arc lengths')
p.add_argument('--interval', type=float, default=30.0)
p.add_argument('--grace-id', default='C')
p.add_argument('--mission', default='GRACE-FO')
p.add_argument('--sat-id', default='C')
p.add_argument('--gravity-nmax', type=int, default=90)
p.add_argument('--chi2-017', type=float, default=25.0)
p.add_argument('--chi2-05', type=float, default=100.0)
args = p.parse_args()

arc_hours = [float(h) for h in args.hours.split(',')]

print("="*70)
print(f"Multi-Date Validation: {len(args.dates)} days, arcs={arc_hours}h")
print("="*70)

all_rms = {}
all_phase_rms = {}
all_code_rms = {}
all_rej = {}

for date_str in args.dates:
    y, m, d = [int(x) for x in date_str.split('-')]
    doy = (datetime(y, m, d) - datetime(y, 1, 1)).days + 1

    # Check that GPS1B exists (either .pkl or .rnx)
    rinex = DATA / f"GPS1B_{date_str}_{args.grace_id}_04.pkl"
    rnx   = DATA / f"GPS1B_{date_str}_{args.grace_id}_04.rnx"
    if not rinex.exists() and not rnx.exists():
        # Check gracefo subdir
        rinex2 = DATA / 'gracefo' / str(y) / date_str / f"GPS1B_{date_str}_{args.grace_id}_04.pkl"
        if rinex2.exists():
            rinex = rinex2
        else:
            print(f"  {date_str}: NO GPS1B DATA, skip")
            continue

    # Check GNV1B
    gnv = DATA / 'gracefo' / str(y) / date_str / f"GNV1B_{date_str}_{args.grace_id}_04.txt"
    gnv2 = DATA / f"GNV1B_{date_str}_{args.grace_id}_04.txt"
    if not gnv.exists() and not gnv2.exists():
        print(f"  {date_str}: NO GNV1B, skip")
        continue

    # Check SP3 for this date
    sp3_pkl = DATA / f"CODE/{y}/COD0OPSFIN_{y}{doy:03d}0000_01D_05M_ORB.pkl"
    sp3_txt  = DATA / f"CODE/{y}/COD0OPSFIN_{y}{doy:03d}0000_01D_05M_ORB.SP3"
    if not sp3_pkl.exists() and not sp3_txt.exists():
        print(f"  {date_str}: NO SP3 for DOY={doy:03d}, skip")
        continue

    clk_path = DATA / f"CODE/{y}/COD0OPSFIN_{y}{doy:03d}0000_01D_30S_CLK.CLK"
    if not clk_path.exists():
        print(f"  {date_str}: NO CLK for DOY={doy:03d}, skip")
        continue

    dcb_p1p2 = DATA / f"CODE/{y}/P1P2{doy:03d}.DCB"
    dcb_p1c1 = DATA / f"CODE/{y}/P1C1{doy:03d}.DCB"

    date_rms = {}
    date_phase = {}
    date_code = {}
    date_rej = {}

    for arc_h in arc_hours:
        chi2 = args.chi2_017 if arc_h < 0.3 else args.chi2_05

        # Build command
        cmd = [
            sys.executable, str(ROOT / 'run_sequential_pod.py'),
            '--date', date_str, '--hours', str(arc_h),
            '--interval', str(args.interval),
            '--grace-id', args.grace_id,
            '--dynamics-mode', 'simplified',
            '--sp3-file', str(sp3_txt if sp3_txt.exists() else sp3_pkl),
            '--clk-file', str(clk_path),
            '--dcb-file', str(dcb_p1p2) if dcb_p1p2.exists() else '',
            '--antex-file', str(DATA / 'igs14.atx'),
            '--iers-c04', str(DATA / 'IERS' / 'eopc04_IAU2000.txt'),
            '--enable-phase-windup', '--enable-relativity',
            '--ar-min-epochs', '6',
            '--gravity-nmax', str(args.gravity_nmax),
            '--chi2-threshold', str(chi2),
        ]
        # Filter out empty args
        cmd = [c for c in cmd if c]

        print(f"  {date_str} {arc_h}h (chi2={chi2})...", end=' ', flush=True)
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600,
                                     env={**os.environ,
                                          'JAVA_HOME': r'C:\Program Files\JetBrains\PyCharm Community Edition 2024.3.5\jbr',
                                          'OREKIT_DATA_PATH': str(DATA / 'orekit')})
            # Parse 3D RMS from output
            for line in result.stdout.split('\n'):
                if '3D RMS vs GNV1B' in line:
                    # "  3D RMS vs GNV1B: 0.293 m  (mean=0.266, max=0.567)"
                    parts = line.split()
                    rms = float(parts[4])
                    date_rms[arc_h] = rms
                    print(f"{rms:.3f}m")
                if 'Phase:' in line and 'RMS=' in line:
                    parts = line.split('RMS=')
                    if len(parts) > 1:
                        date_phase[arc_h] = float(parts[1].split('m')[0])
                if 'Code:' in line and 'RMS=' in line:
                    parts = line.split('RMS=')
                    if len(parts) > 1:
                        date_code[arc_h] = float(parts[1].strip().split('m')[0])
                if 'Rejected:' in line:
                    date_rej[arc_h] = int(line.split()[-1])

            if arc_h not in date_rms:
                print("FAIL (no RMS line)")
                if result.stderr:
                    print(f"    ERR: {result.stderr[:200]}")
        except subprocess.TimeoutExpired:
            print("TIMEOUT")
        except Exception as e:
            print(f"ERROR: {e}")

    if date_rms:
        all_rms[date_str] = date_rms
        all_phase_rms[date_str] = date_phase
        all_code_rms[date_str] = date_code
        all_rej[date_str] = date_rej

# ── Summary ──
print(f"\n{'='*70}")
print("Summary")
print(f"{'='*70}")
for arc_h in arc_hours:
    vals = [all_rms[d][arc_h] for d in all_rms if arc_h in all_rms[d]]
    if vals:
        v = np.array(vals)
        print(f"  {arc_h}h: n={len(v)}  mean={np.mean(v):.3f}m  "
              f"median={np.median(v):.3f}m  std={np.std(v):.3f}m  "
              f"best={np.min(v):.3f}m  worst={np.max(v):.3f}m")

# ── Save ──
out_dir = ROOT / 'results' / 'multiday'
out_dir.mkdir(parents=True, exist_ok=True)
pickle.dump({
    'rms': dict(all_rms), 'phase': dict(all_phase_rms),
    'code': dict(all_code_rms), 'rej': dict(all_rej),
    'arc_hours': arc_hours, 'args': vars(args),
}, open(str(out_dir / 'multiday_results.pkl'), 'wb'))
print(f"\nSaved: {out_dir / 'multiday_results.pkl'}")

# ── Plot ──
if len(all_rms) >= 2:
    fig, axes = plt.subplots(len(arc_hours), 1, figsize=(12, 3*len(arc_hours)), squeeze=False)
    for idx, arc_h in enumerate(arc_hours):
        ax = axes[idx][0]
        dates_sorted = sorted(all_rms.keys())
        rms_vals = [all_rms[d].get(arc_h, np.nan) for d in dates_sorted]
        phase_vals = [all_phase_rms.get(d, {}).get(arc_h, np.nan) for d in dates_sorted]
        code_vals = [all_code_rms.get(d, {}).get(arc_h, np.nan) for d in dates_sorted]

        x = range(len(dates_sorted))
        ax.plot(x, rms_vals, 'o-', color='#2196F3', markersize=7, linewidth=1.5, label='3D RMS')
        ax.plot(x, phase_vals, 's--', color='#4CAF50', markersize=4, alpha=0.7, label='Phase RMS')
        ax.plot(x, code_vals, '^--', color='#FF9800', markersize=4, alpha=0.7, label='Code RMS')
        ax.axhline(y=np.nanmean(rms_vals), color='red', linestyle=':', alpha=0.5,
                   label=f'Mean {np.nanmean(rms_vals):.3f}m')
        ax.set_xticks(list(x))
        ax.set_xticklabels([d[5:] for d in dates_sorted], rotation=45, fontsize=8)
        ax.set_ylabel(f'{arc_h}h [m]')
        ax.set_title(f'GRACE-FO {args.grace_id} — {arc_h}h Arc POD')
        ax.legend(fontsize=7, loc='upper left')
        ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(str(out_dir / 'multiday_accuracy.png'), dpi=150, bbox_inches='tight')
    print(f"Saved plot: {out_dir / 'multiday_accuracy.png'}")
