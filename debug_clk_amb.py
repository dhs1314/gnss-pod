"""Test clock+ambiguity estimation with fixed orbit (GNV1B reference)."""
import sys, pickle, numpy as np
from pathlib import Path
from datetime import datetime, timedelta
sys.path.insert(0, 'src')
from gps1b_rnx_loader import load_gps1b_rnx
from sp3_loader import get_gps_pos_from_sp3 as _sp3_get
from orbit_dynamics import total_acc_eci, GM, OMEGA_E
from coordinates import ecef_to_eci, eci_to_ecef
from gravity_model import read_icgem_gfc
from orbit_integrator import integrate_orbit_eci_with_stm

C_LIGHT = 299792458.0
SEC_PER_DAY = 86400.0
MJD_J2000 = 51544.5
J2000 = datetime(2000, 1, 1, 12, 0, 0)
F1, F2 = 1575.42e6, 1227.60e6
F1_SQ, F2_SQ = F1 * F1, F2 * F2
ALPHA = F1_SQ / (F1_SQ - F2_SQ)
BETA = -F2_SQ / (F1_SQ - F2_SQ)

def load_gnv1b(filepath):
    orbit = {}
    with open(filepath, encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            parts = line.split()
            if len(parts) < 6: continue
            try:
                t, flag = float(parts[0]), parts[2]
                if flag not in ('C', 'E'): continue
                X, Y, Z = float(parts[3]), float(parts[4]), float(parts[5])
                if abs(X) < 1e3: continue
                orbit[t] = np.array([X, Y, Z])
            except: continue
    return orbit

def interpolate_ref(ref_orbit, gps_sod):
    ts = sorted(ref_orbit.keys())
    t0 = t1 = None
    for i, ti in enumerate(ts):
        if ti >= gps_sod: t1 = ti; t0 = ts[i-1] if i > 0 else None; break
        t0 = ti
    if t1 is None: t0 = t1 = ts[-1]
    if t0 is None: t0 = ts[0]
    if t0 == t1: return ref_orbit[t0]
    a = (gps_sod - t0) / (t1 - t0)
    return ref_orbit[t0] * (1 - a) + ref_orbit[t1] * a

def get_sat_geometry(sp3, sv, utc_dt, rcv_pos):
    pos, clk, vel = _sp3_get(sp3, sv, utc_dt)
    if pos is None or abs(clk) > 0.1 * C_LIGHT:
        return None, None, None
    rho = float(np.linalg.norm(pos - rcv_pos))
    for _ in range(5):
        travel_time = rho / C_LIGHT
        tx_dt = utc_dt - timedelta(seconds=travel_time)
        pos_tx, clk_tx, vel_tx = _sp3_get(sp3, sv, tx_dt)
        if pos_tx is None: break
        rho_new = float(np.linalg.norm(pos_tx - rcv_pos))
        if abs(rho_new - rho) < 1e-8:
            pos, clk = pos_tx, clk_tx; rho = rho_new; break
        pos, clk = pos_tx, clk_tx; rho = rho_new
    if not (1.8e7 < rho < 2.8e7): return None, None, None
    sag = (OMEGA_E / C_LIGHT) * (pos[0] * rcv_pos[1] - pos[1] * rcv_pos[0])
    return pos, clk, rho + sag

# Load data
dp = Path('data')
date_str = '2024-04-29'; y, m, d = 2024, 4, 29; grace_id = 'C'
doy = (datetime(y, m, d) - datetime(y, 1, 1)).days + 1

gnv_path = dp / 'gracefo' / str(y) / date_str / f'GNV1B_{date_str}_{grace_id}_04.txt'
ref_orbit = load_gnv1b(str(gnv_path))
print(f'GNV1B: {len(ref_orbit)} epochs')

rnx_path = dp / f'GPS1B_{date_str}_{grace_id}_04.rnx'
gps1b_raw = load_gps1b_rnx(str(rnx_path))
print(f'GPS1B: {len(gps1b_raw)} epochs')

sp3_pkl = dp / str(y) / f'{doy:03d}' / 'igs_sp3_FIN.pkl'
if not sp3_pkl.exists(): sp3_pkl = dp / str(y) / f'{doy:03d}' / 'igs_sp3.pkl'
sp3 = pickle.load(open(str(sp3_pkl), 'rb'))
print(f'SP3: {len(sp3["ts"])} epochs')

# Select epochs (5 min at 30s)
hours = 0.0833
interval = 30.0
gps_sod_start = min(gps1b_raw.keys())
gps_sod_end = gps_sod_start + hours * 3600
epochs = []
for gps_sod in sorted(gps1b_raw.keys()):
    if not (gps_sod_start <= gps_sod <= gps_sod_end): continue
    dt_ep = gps_sod - gps_sod_start
    nearest = round(dt_ep / interval) * interval
    if abs(dt_ep - nearest) > max(2.0, interval * 0.1): continue
    epochs.append(gps_sod)
print(f'Epochs: {len(epochs)}')

# Compute sv_bias from first epochs
N_BIAS = min(60, len(epochs))
sv_p_residuals = {}
epoch_geo = {}
for gps_sod in epochs:
    utc_dt = J2000 + timedelta(seconds=gps_sod)
    ref_pos = interpolate_ref(ref_orbit, gps_sod)
    ep_data = []
    recs = gps1b_raw.get(gps_sod, {})
    for sv_id, rec in recs.items():
        if 'L_if' not in rec: continue
        sat_pos, sat_clk, rho_corr = get_sat_geometry(sp3, sv_id, utc_dt, ref_pos)
        if sat_pos is None: continue
        el = np.arcsin(abs(sat_pos[2] - ref_pos[2]) / rho_corr)
        if el < 0.087: continue
        L_if_raw = float(rec['L_if']); P_if_raw = float(rec['P_if'])
        L_r = L_if_raw + sat_clk - rho_corr
        P_r = P_if_raw + sat_clk - rho_corr
        ep_data.append({'sv': sv_id, 'sat_pos': sat_pos, 'sat_clk': sat_clk,
                        'rho_corr': rho_corr, 'el': el,
                        'L_r': L_r, 'P_r': P_r,
                        'L_if_raw': L_if_raw, 'P_if_raw': P_if_raw})
        if len(epoch_geo) < N_BIAS:
            sv_p_residuals.setdefault(sv_id, []).append(P_r)
    if ep_data: epoch_geo[gps_sod] = ep_data

sv_bias = {}
for sv, p_vals in sv_p_residuals.items():
    if len(p_vals) >= 3:
        sv_bias[sv] = float(np.median(p_vals))
print(f'SV code biases: {len(sv_bias)} SVs')

all_svs = set()
for gps_sod in epochs:
    for d in epoch_geo.get(gps_sod, []):
        if d['sv'] in sv_bias: all_svs.add(d['sv'])
sv_list = sorted(all_svs)
N_SV = len(sv_list); N_EPOCH = len(epochs)
sv_to_idx = {sv: i for i, sv in enumerate(sv_list)}
print(f'SVs: {N_SV}, Epochs: {N_EPOCH}')

W_PHASE = 1.0 / 0.01**2
W_CODE = 1.0 / 0.30**2

# Solve for clocks + ambiguities only (fixed orbit = GNV1B)
N_CLK_AMB = N_EPOCH + N_SV
N_full = np.zeros((N_CLK_AMB, N_CLK_AMB))
b_full = np.zeros(N_CLK_AMB)

n_phase = 0; n_code = 0; total_res_sq = 0.0
for i_ep, gps_sod in enumerate(epochs):
    ed_list = epoch_geo.get(gps_sod, [])
    ref_pos = interpolate_ref(ref_orbit, gps_sod)

    for d in ed_list:
        sv = d['sv']
        if sv not in sv_to_idx or sv not in sv_bias: continue
        i_sv = sv_to_idx[sv]
        sat_pos = d['sat_pos']; sat_clk = d['sat_clk']

        # Modeled range (from GNV1B — fixed)
        rho_model = float(np.linalg.norm(sat_pos - ref_pos))
        sag = (OMEGA_E / C_LIGHT) * (sat_pos[0] * ref_pos[1] - sat_pos[1] * ref_pos[0])
        rho_corr_model = rho_model + sag

        # Phase
        obs = d['L_if_raw'] + sat_clk - rho_corr_model - sv_bias[sv]
        h = np.zeros(N_CLK_AMB)
        h[i_ep] = 1.0          # clock
        h[N_EPOCH + i_sv] = 1.0  # ambiguity
        N_full += W_PHASE * np.outer(h, h)
        b_full += W_PHASE * obs * h
        n_phase += 1
        total_res_sq += W_PHASE * obs**2

        # Code
        obs = d['P_if_raw'] + sat_clk - rho_corr_model - sv_bias[sv]
        h = np.zeros(N_CLK_AMB)
        h[i_ep] = 1.0
        # no ambiguity for code
        N_full += W_CODE * np.outer(h, h)
        b_full += W_CODE * obs * h
        n_code += 1
        total_res_sq += W_CODE * obs**2

# Add weak prior on first clock epoch
N_full[0, 0] += 1.0 / 10000.0

# Solve
try:
    dx = np.linalg.solve(N_full, b_full)
except np.linalg.LinAlgError:
    print('Singular!')
    dx = np.linalg.lstsq(N_full, b_full, rcond=None)[0]

clk = dx[:N_EPOCH]
B_sv = dx[N_EPOCH:]

# Recompute residuals
residuals = []
for i_ep, gps_sod in enumerate(epochs):
    ed_list = epoch_geo.get(gps_sod, [])
    ref_pos = interpolate_ref(ref_orbit, gps_sod)
    for d in ed_list:
        sv = d['sv']
        if sv not in sv_to_idx or sv not in sv_bias: continue
        i_sv = sv_to_idx[sv]
        sat_pos = d['sat_pos']; sat_clk = d['sat_clk']
        rho_model = float(np.linalg.norm(sat_pos - ref_pos))
        sag = (OMEGA_E / C_LIGHT) * (sat_pos[0] * ref_pos[1] - sat_pos[1] * ref_pos[0])
        rho_corr_model = rho_model + sag

        # Phase
        obs = d['L_if_raw'] + sat_clk - rho_corr_model - sv_bias[sv]
        model = clk[i_ep] + B_sv[i_sv]
        residuals.append(obs - model)

post_rms = np.sqrt(np.mean(np.array(residuals)**2))
print(f'\nClocks + Ambiguities only (orbit fixed to GNV1B):')
print(f'  n_phase={n_phase}, n_code={n_code}')
print(f'  Clock range: [{np.min(clk):.1f}, {np.max(clk):.1f}] m')
print(f'  Amb range: [{np.min(B_sv):.1f}, {np.max(B_sv):.1f}] m')
print(f'  Post-fit residual RMS: {post_rms:.4f} m')
print(f'  Weighted RMS: {np.sqrt(total_res_sq/(n_phase+n_code)):.4f}')
