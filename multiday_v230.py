"""V2.3.0 Multi-date validation (Phase 13.0).

Processes all available dates with GPS1B+GNV1B data.
Reuses CODE products from 2024-04-29 for all dates.
Runs EKF (0.17h/0.5h) + Batch solver on each date.
"""
import sys, os, pickle, math
from pathlib import Path
from datetime import datetime, timedelta
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict

ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'src'))

DATA = ROOT / 'data'
OUT_DIR = ROOT / 'results' / 'multiday_v230'
OUT_DIR.mkdir(parents=True, exist_ok=True)

GRACE_ID = 'C'; INTERVAL = 30; GRAV_NMAX = 90
GPS_SOD_START = 767620800; MJD_J2000 = 51544.5

# ── Find all dates with GPS1B + GNV1B ──
print("Scanning available dates...")
available_dates = []
for f in sorted(os.listdir(str(DATA))):
    if f.startswith('GPS1B_') and f.endswith(f'_{GRACE_ID}_04.pkl'):
        d = f[6:16]
        y, m, day = [int(x) for x in d.split('-')]
        gnv = DATA / 'gracefo' / str(y) / d / f'GNV1B_{d}_{GRACE_ID}_04.txt'
        gnv2 = DATA / f'GNV1B_{d}_{GRACE_ID}_04.txt'
        if gnv.exists() or gnv2.exists():
            available_dates.append(d)

# Limit to ~10 dates spread across the available range for speed
if len(available_dates) > 12:
    step = max(1, len(available_dates) // 10)
    available_dates = available_dates[::step][:12]
    available_dates.append('2024-04-29')  # always include reference date
    available_dates = sorted(set(available_dates))

ARC_HOURS = [0.17, 0.5]
print(f"Dates to process: {len(available_dates)}")
for d in available_dates[:5]:
    print(f"  {d}")
if len(available_dates) > 5:
    print(f"  ... and {len(available_dates)-5} more")

# ── Load shared products (CODE from 2024-04-29 reused for all) ──
from gravity_model import read_icgem_gfc
gfc_path = DATA / 'gravity' / 'GGM05C.gfc'
Cnm, Snm, Nmax, GM_grav, R_grav = read_icgem_gfc(str(gfc_path))
Nmax = min(Nmax, GRAV_NMAX)

sp3 = pickle.load(open(str(DATA / "CODE/2024/COD0OPSFIN_20241200000_01D_05M_ORB.pkl"), 'rb'))
from precision_products import (read_rinex_clk, read_antex, load_code_dcb_pair,
                                 compute_dcb_if_correction, setup_iers_from_c04)
clk_data = read_rinex_clk(str(DATA / "CODE/2024/COD0OPSFIN_20241200000_01D_30S_CLK.CLK"))
antex = read_antex(str(DATA / "igs14.atx"))
dcb_pair = load_code_dcb_pair(str(DATA / "CODE/2024/P1C12404.DCB"),
                               str(DATA / "CODE/2024/P1P22404.DCB"))
setup_iers_from_c04(str(DATA / "IERS/eopc04_IAU2000.txt"))
print(f"Products loaded: SP3={len(sp3['ts'])}ep CLK={len(clk_data)}SV grav Nmax={Nmax}")

# ── EKF config ──
ekf_config_template = {
    'Cd': 2.2, 'CR': 1.3, 'area_drag': 0.68, 'area_srp': 3.4, 'mass': 580.0,
    'bodies': ['Sun', 'Moon'], 'Cnm': Cnm, 'Snm': Snm,
    'GM_grav': GM_grav, 'R_grav': R_grav, 'gravity_nmax': Nmax,
    'sigma_acc_process': 1e-3, 'tau_emp': 600.0, 'sigma_emp_ss': 1e-8,
    'sigma_zwd_rw': 1e-9, 'sigma_phase': 0.20, 'sigma_code': 0.30,
    'el_min': 0.087, 'use_phase_windup': True, 'use_relativity': True,
    'use_cycle_slip': False, 'ar_min_epochs': 6,
    'antex_data': antex, 'dcb_data': dcb_pair,
    'elev_exp_phase': 1.0, 'mw_max_epochs': 200,
}

all_results = []  # {date, arc_h, rms_3d, phase_ekf, phase_batch}


def process_one_arc(date_str, arc_hours):
    """Run EKF for one arc, return metrics dict or None."""
    from run_sequential_pod import load_gnv1b, compute_epoch_geometry, interpolate_ref
    from sequential_filter import SequentialEKF
    from coordinates import ecef_to_eci, eci_to_ecef
    from troposphere import saastamoinen_zhd
    from batch_solver import BatchLinearSolver

    y, m, d_ = [int(x) for x in date_str.split('-')]

    # Load GPS1B
    gps1b_path = DATA / f'GPS1B_{date_str}_{GRACE_ID}_04.pkl'
    if not gps1b_path.exists():
        gps1b_path = DATA / f'GPS1B_{date_str}_{GRACE_ID}_04.rnx'
        if not gps1b_path.exists():
            return None
    if str(gps1b_path).endswith('.pkl'):
        gps1b = pickle.load(open(str(gps1b_path), 'rb'))
    else:
        from gps1b_rnx_loader import load_gps1b_rnx
        gps1b = load_gps1b_rnx(str(gps1b_path))

    gps_start = min(gps1b.keys())
    gps_end = gps_start + int(arc_hours * 3600)
    epochs = sorted(set(g for g in gps1b.keys()
                        if gps_start <= g <= gps_end
                        and abs((g - gps_start) % INTERVAL) <= 2.0))
    if len(epochs) < 6:
        return None

    # Load GNV1B
    gnv_path = DATA / 'gracefo' / str(y) / date_str / f'GNV1B_{date_str}_{GRACE_ID}_04.txt'
    gnv2 = DATA / f'GNV1B_{date_str}_{GRACE_ID}_04.txt'
    ref_orbit, ref_vel = load_gnv1b(str(gnv_path if gnv_path.exists() else gnv2))

    # Per-SV code bias
    N_BIAS = min(60, len(epochs))
    sv_p_res = {}
    J2000 = datetime(2000, 1, 1, 12, 0, 0)
    for gps_sod in epochs[:N_BIAS]:
        ref_pos = interpolate_ref(ref_orbit, gps_sod)
        if ref_pos is None: continue
        utc_dt = J2000 + timedelta(seconds=gps_sod)
        recs = gps1b.get(int(gps_sod), gps1b.get(gps_sod, {}))
        from run_sequential_pod import get_sat_geometry
        for sv_id, rec in recs.items():
            if 'P_if' not in rec: continue
            sat_pos, sat_clk, rho_corr = get_sat_geometry(sp3, sv_id, utc_dt, ref_pos, clk_data)
            if sat_pos is None: continue
            dcb_corr = compute_dcb_if_correction(dcb_pair, sv_id)
            P_r = float(rec['P_if']) + dcb_corr + sat_clk - rho_corr
            sv_p_res.setdefault(sv_id, []).append(P_r)
    sv_bias = {}
    for sv, vals in sv_p_res.items():
        if len(vals) >= 3: sv_bias[sv] = float(np.median(vals))
    sv_bias_ref = float(np.mean(list(sv_bias.values()))) if sv_bias else 0.0
    for sv in sv_bias: sv_bias[sv] -= sv_bias_ref

    # Init EKF
    r0_ecef = interpolate_ref(ref_orbit, epochs[0])
    v0_ecef = interpolate_ref(ref_vel, epochs[0])
    if r0_ecef is None or v0_ecef is None: return None
    mjd_start = MJD_J2000 + epochs[0] / 86400.0
    r0_eci, v0_eci = ecef_to_eci(r0_ecef, v0_ecef, mjd_start)

    chi2 = 25 if arc_hours < 0.3 else 100
    ekf_cfg = dict(ekf_config_template)
    ekf_cfg['chi2_threshold'] = chi2
    ekf_cfg['elev_exp_code'] = 1.0 if arc_hours < 0.3 else 0.70
    ekf_cfg['clock_rw'] = 0.0004 if arc_hours < 0.3 else 0.001

    ekf = SequentialEKF(ekf_cfg)
    state = ekf.initialize(r0_eci, v0_eci, mjd_start, epochs[0])

    dr_vals = []; pass1_geometry = []
    SEC_PER_DAY = 86400.0; C_L = 299792458.0; OE = 7.2921151467e-5

    for i_ep, gps_sod in enumerate(epochs):
        mjd_utc = MJD_J2000 + gps_sod / 86400.0
        if i_ep > 0:
            mjd_prev = MJD_J2000 + epochs[i_ep-1] / 86400.0
            state = ekf.predict(state, gps_sod, mjd_prev, mjd_prev + 69.184/86400.0)
        rcv_ecef, _ = eci_to_ecef(state.r_eci, state.v_eci, mjd_utc)
        ep_data = compute_epoch_geometry(gps_sod, gps1b, sp3, rcv_ecef, clk_data)
        if not ep_data: continue
        doy = (datetime(y, m, d_) - datetime(y, 1, 1)).days + 1
        state, stats = ekf.process_epoch(state, ep_data, sp3, sv_bias, sv_bias_ref,
                                          mjd_utc, mjd_utc + 69.184/86400.0, doy)
        r_ecef, _ = eci_to_ecef(state.r_eci, state.v_eci, mjd_utc)
        r_gnv = interpolate_ref(ref_orbit, gps_sod)
        if r_gnv is not None:
            dr_vals.append(np.linalg.norm(r_ecef - r_gnv))

        # Store geometry for batch solver
        lat_rad = np.arcsin(rcv_ecef[2] / np.linalg.norm(rcv_ecef))
        h_m = np.linalg.norm(rcv_ecef) - 6378137.0; zhd = saastamoinen_zhd(lat_rad, h_m)
        for d in ep_data:
            sv = d['sv']; se = np.asarray(d['sat_pos'], dtype=float)
            sc = float(d.get('sat_clk', 0)); el = float(d.get('el', 0.5))
            rho = np.linalg.norm(se - rcv_ecef)
            sag = (OE/C_L)*(se[0]*rcv_ecef[1] - se[1]*rcv_ecef[0])
            mf = 1.001/np.sqrt(0.002001+np.sin(el)**2)
            dcb_c = compute_dcb_if_correction(dcb_pair, sv)
            d['_geo_full'] = rho + sag - sc + zhd * mf
            d['_obs_code'] = float(d.get('P_if_raw', 0)) + dcb_c - sv_bias.get(sv, 0.0)
            d['_obs_phase'] = float(d.get('L_if_raw', 0)) - sv_bias.get(sv, 0.0)
        pass1_geometry.append(ep_data)

    rms_3d = np.sqrt(np.mean([d**2 for d in dr_vals])) if dr_vals else 0
    ekf_phase = stats['rms_phase'] if rms_3d > 0 else 0

    # Batch solver
    bls = BatchLinearSolver(pass1_geometry, sigma_phase=0.20, sigma_code=0.30)
    bls_sol = bls.solve()
    batch_phase = bls_sol['rms_phase']
    batch_code = bls_sol['rms_code']

    return {
        'rms_3d': rms_3d, 'phase_ekf': ekf_phase,
        'phase_batch': batch_phase, 'code_batch': batch_code,
        'n_sv': state.n_sv, 'n_epochs': len(dr_vals),
    }


# ── Process all dates ──
for date_str in available_dates:
    print(f"\n{date_str}:")
    for arc_h in ARC_HOURS:
        metrics = process_one_arc(date_str, arc_h)
        if metrics:
            imp = (metrics['phase_ekf'] - metrics['phase_batch'])/metrics['phase_ekf']*100 if metrics['phase_ekf'] > 0 else 0
            all_results.append({
                'date': date_str, 'arc_h': arc_h,
                'rms_3d': metrics['rms_3d'],
                'phase_ekf': metrics['phase_ekf'],
                'phase_batch': metrics['phase_batch'],
                'improvement': imp,
                'n_sv': metrics['n_sv'],
                'n_epochs': metrics['n_epochs'],
            })
            print(f"  {arc_h}h: EKF={metrics['rms_3d']:.3f}m  "
                  f"phase: {metrics['phase_ekf']:.3f}→{metrics['phase_batch']:.3f}m ({imp:+.0f}%)  "
                  f"SVs={metrics['n_sv']}")

# ── Statistics ──
print(f"\n{'='*70}")
print(f"Multi-Date Summary ({len(all_results)} arcs across {len(available_dates)} dates)")
print(f"{'='*70}")

for arc_h in ARC_HOURS:
    subset = [r for r in all_results if r['arc_h'] == arc_h]
    if not subset: continue
    rms_v = np.array([r['rms_3d'] for r in subset])
    imp_v = np.array([r['improvement'] for r in subset])
    print(f"\n{arc_h}h ({len(subset)} dates):")
    print(f"  EKF 3D RMS:   mean={np.mean(rms_v):.3f}m  median={np.median(rms_v):.3f}m  "
          f"best={np.min(rms_v):.3f}m  worst={np.max(rms_v):.3f}m")
    print(f"  Batch Phase improvement: mean={np.mean(imp_v):+.1f}%  "
          f"best={np.max(imp_v):+.1f}%  worst={np.min(imp_v):+.1f}%")
    n_pos = sum(1 for v in imp_v if v > 0)
    print(f"  Batch improves EKF phase: {n_pos}/{len(subset)} dates")

# ── Save ──
pickle.dump(all_results, open(str(OUT_DIR / 'multiday_v230.pkl'), 'wb'))

# ── Plot ──
fig, axes = plt.subplots(len(ARC_HOURS), 1, figsize=(14, 5*len(ARC_HOURS)), squeeze=False)
dates_plt = sorted(set(r['date'] for r in all_results))

for idx, arc_h in enumerate(ARC_HOURS):
    ax = axes[idx][0]
    subset = [r for r in all_results if r['arc_h'] == arc_h]
    if not subset: continue

    date_to_rms = {r['date']: r['rms_3d'] for r in subset}
    rms_vals = [date_to_rms.get(d, np.nan) for d in dates_plt]

    colors = ['#2196F3' if not np.isnan(v) and v < 1.0 else '#FF9800' for v in rms_vals]
    ax.bar(range(len(dates_plt)), rms_vals, color=colors, alpha=0.7, edgecolor='#333')
    ax.axhline(y=np.nanmean(rms_vals), color='red', linestyle='--', alpha=0.5,
               label=f'Mean: {np.nanmean(rms_vals):.3f}m')
    ax.axhline(y=0.5, color='green', linestyle=':', alpha=0.4)
    ax.set_ylabel(f'{arc_h}h 3D RMS [m]', fontsize=12)
    ax.set_title(f'V2.3.0 EKF — {arc_h}h Arc POD ({len(subset)} dates)', fontsize=13)
    ax.set_xticks(range(len(dates_plt)))
    ax.set_xticklabels([d[5:] for d in dates_plt], rotation=45, fontsize=9)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    # Annotate values
    for i, v in enumerate(rms_vals):
        if not np.isnan(v):
            ax.text(i, v+0.05, f'{v:.2f}', ha='center', fontsize=7, rotation=90)

plt.tight_layout()
png_path = OUT_DIR / 'multiday_v230.png'
plt.savefig(str(png_path), dpi=150, bbox_inches='tight')
print(f"\n  Plot: {png_path}")
