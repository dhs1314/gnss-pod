#!/usr/bin/env python3
"""GRACE-FO PPP — apr29_final_test.py验证算法批量版"""
import pickle, sys, math, datetime as dt, numpy as np, argparse, warnings, os, csv
warnings.filterwarnings('ignore')

sys.path.insert(0, '/workspace/gnss_pod/src')
from batch_v12 import parse_gnv1b

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
F1_SQ, F2_SQ = F1*F1, F2*F2
ALPHA = F1_SQ / (F1_SQ - F2_SQ)
BETA = -F2_SQ / (F1_SQ - F2_SQ)
OMEGA_E = 7.2921151467e-5
MU_E = 3.986004418e14
GPS_ORIGIN = dt.datetime(2000, 1, 1, 12, 0, 0)
ITER_MAX = 15
TOL = 0.05
A_WGS84 = 6378137.0
E2_WGS84 = 0.00669437999014

GPS_SV = {
    'G05': {'M0': 0.0, 'Omega0': 0.0, 'omega': 0.0, 'ecc': 0.000020, 'inc': 55.0, 'sqrtA': 5153.6},
    'G07': {'M0': math.radians(150), 'Omega0': math.radians(150), 'omega': 0.0, 'ecc': 0.000015, 'inc': 55.0, 'sqrtA': 5153.6},
    'G13': {'M0': math.radians(60), 'Omega0': math.radians(60), 'omega': 0.0, 'ecc': 0.000010, 'inc': 55.0, 'sqrtA': 5153.6},
    'G20': {'M0': math.radians(90), 'Omega0': math.radians(90), 'omega': 0.0, 'ecc': 0.000020, 'inc': 55.0, 'sqrtA': 5153.6},
    'G27': {'M0': math.radians(30), 'Omega0': math.radians(30), 'omega': 0.0, 'ecc': 0.000010, 'inc': 55.0, 'sqrtA': 5153.6},
}

def gps_pos(sv, t_gps, t0=767620800.0):
    if sv not in GPS_SV:
        return None
    p = GPS_SV[sv]
    a = p['sqrtA']**2  # 5153.6^2 = 26,564,000 (m^2, R=26560 km ✓)
    n0 = math.sqrt(MU_E / a**3)
    M = p['M0'] + n0*(t_gps - t0)
    E = M
    for _ in range(10):
        E = M + p['ecc']*math.sin(E)
    v = 2*math.atan2(math.sqrt(1.00004)*math.sin(E/2), math.cos(E/2)-0.00002)
    phi = v + math.radians(p['omega'])
    r = a*(1 - p['ecc']*math.cos(E))
    Om = p['Omega0'] - OMEGA_E*(t_gps-t0)
    inc = math.radians(p['inc'])
    X = r*(math.cos(phi)*math.cos(Om) - math.sin(phi)*math.cos(inc)*math.sin(Om))
    Y = r*(math.cos(phi)*math.sin(Om) + math.sin(phi)*math.cos(inc)*math.cos(Om))
    Z = r*math.sin(phi)*math.sin(inc)
    return (X, Y, Z)

def ecef_to_blh(pos):
    X, Y, Z = float(pos[0]), float(pos[1]), float(pos[2])
    p = math.sqrt(X**2 + Y**2)
    if p < 1e-12:
        return (math.pi/2*(1 if Z >= 0 else -1), 0.0, abs(Z)-A_WGS84)
    e2 = E2_WGS84
    val = 1 - e2*(Z/p)**2
    lat = math.atan2(Z, p/math.sqrt(val)) if val > 0 else math.pi/2
    lon = math.atan2(Y, X)
    sinL = math.sin(lat)
    N = A_WGS84 / math.sqrt(1 - e2*sinL**2)
    h = p/math.cos(lat) - N
    return (lat, lon, h)

def ecef_to_enu(err, lat, lon):
    cl, sl, cn, sn = math.cos(lat), math.sin(lat), math.cos(lon), math.sin(lon)
    return (
        -sl*cn*err[0] - sl*sn*err[1] + cl*err[2],
         sl*sn*err[0] - sl*cn*err[1],
         cl*cn*err[0] + cl*sn*err[1] + sl*err[2],
    )

def geometric_range(a, b):
    return math.sqrt((b[0]-a[0])**2+(b[1]-a[1])**2+(b[2]-a[2])**2)

def ppp_day(gps1b_day, gnv_day, label):
    if not gps1b_day or not gnv_day:
        return [], label

    gnv_keys = sorted(gnv_day.keys())
    gp0 = np.array(gnv_day[gnv_keys[0]])
    ref_lat, ref_lon, _ = ecef_to_blh(gp0)

    results = []
    prev = None
    for t_gps in sorted(gps1b_day.keys()):
        obs = gps1b_day[t_gps]
        if len(obs) < 4:
            continue

        utc = GPS_ORIGIN + dt.timedelta(seconds=t_gps)

        rp0 = rp1 = None
        for tk in gnv_keys:
            if tk >= utc:
                rp1 = tk
                break
            rp0 = tk
        if rp0 is None:
            rp0 = gnv_keys[0]
        if rp1 is None:
            rp1 = gnv_keys[-1]
        if rp0 == rp1:
            gp = np.array(gnv_day[rp0])
        else:
            dt0 = (utc - rp0).total_seconds()
            dtt = (rp1 - rp0).total_seconds()
            gp = np.array(gnv_day[rp0])*(1-dt0/dtt) + np.array(gnv_day[rp1])*(dt0/dtt)

        obs_list = []
        for sv_id, rec in obs.items():
            if rec.get('L1') is None:
                continue
            sp = gps_pos(sv_id, t_gps)
            if sp is None:
                continue
            L1 = rec['L1']
            L2 = rec['L2']
            L_if = ALPHA*L1 + BETA*L2
            obs_list.append({'sat_pos': sp, 'L_if': L_if, 'sv': sv_id})

        if len(obs_list) < 4:
            continue

        x0 = np.array([gp[0], gp[1], gp[2], 0.0]) if prev is None else prev

        for it in range(ITER_MAX):
            H, y, W = [], [], []
            for o in obs_list:
                sp = o['sat_pos']
                rcv = x0[:3]
                rho = geometric_range(rcv, sp)
                if rho < 1e6:
                    continue
                hx = -(sp[0]-rcv[0])/rho
                hy = -(sp[1]-rcv[1])/rho
                hz = -(sp[2]-rcv[2])/rho
                H.append([hx, hy, hz, 1.0])
                y.append(o['L_if'] - rho)
                W.append(1.0)
            if len(H) < 4:
                break
            HT = np.array(H).T
            W_arr = np.diag(W)
            HTWH = HT @ W_arr @ np.array(H) + np.eye(4)*1e-8
            HTWy = HT @ W_arr @ np.array(y)
            try:
                dx = np.linalg.solve(HTWH, HTWy)
            except:
                break
            if np.any(np.isnan(dx)):
                break
            x0 = x0 + dx
            if np.linalg.norm(dx[:3]) < TOL and abs(dx[3]) < TOL:
                break

        err = x0[:3] - gp
        enu = ecef_to_enu(err, ref_lat, ref_lon)
        d3 = math.sqrt(sum(x**2 for x in enu))
        result = {
            'time': utc.strftime('%Y-%m-%dT%H:%M:%S'),
            'dE_m': float(enu[0]),
            'dN_m': float(enu[1]),
            'dU_m': float(enu[2]),
            'd3_m': float(d3),
            'n_sat': len(obs_list),
        }
        results.append(result)
        prev = x0.copy()

    return results, label

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--start', required=True)
    p.add_argument('--end', required=True)
    p.add_argument('--grace', default='C')
    p.add_argument('--out', default='output_v12')
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)

    cur = dt.date.fromisoformat(args.start)
    end = dt.date.fromisoformat(args.end)
    delta = dt.timedelta(days=1)
    all_stats = []

    while cur <= end:
        ds = cur.isoformat()
        y, m, d = cur.year, cur.month, cur.day

        pkl_paths = [
            f'data/gracefo/{y}/{ds}/gps1b_{args.grace}.pkl',
            f'data/gracefo/{y}/{y:04d}-{m:02d}-{d:02d}/gps1b_{args.grace}.pkl',
            f'data/gracefo/{y}/{y:04d}-{m:02d}-{d:02d}/GPS1B_{y:04d}-{m:02d}-{d:02d}_{args.grace}_04.pkl',
        ]
        gnv_path = f'data/gracefo/{y}/{ds}/GNV1B_{ds}_{args.grace}_04.txt'

        gps1b = None
        for pp in pkl_paths:
            if os.path.exists(pp):
                gps1b = pickle.load(open(pp, 'rb'))
                break

        gnv = parse_gnv1b(gnv_path) if os.path.exists(gnv_path) else None

        if gps1b is None or gnv is None:
            print(f"\n=== {ds} 数据缺失 ===")
            cur += delta
            continue

        print(f"\n=== {ds} ===")
        print(f"  GPS1B {len(gps1b)} epochs, GNV1B {len(gnv)} epochs")

        results, _ = ppp_day(gps1b, gnv, ds)

        if results:
            n = len(results)
            rms_e = math.sqrt(sum(r['dE_m']**2 for r in results)/n)*100
            rms_n = math.sqrt(sum(r['dN_m']**2 for r in results)/n)*100
            rms_u = math.sqrt(sum(r['dU_m']**2 for r in results)/n)*100
            rms_3d = math.sqrt(sum(r['d3_m']**2 for r in results)/n)*100
            d3s = [r['d3_m'] for r in results]
            median = sorted(d3s)[n//2]*100
            good = sum(1 for x in d3s if abs(x) < 1.0)
            print(f"  {n} epochs, RMS E/N/U/3D: {rms_e:.1f}/{rms_n:.1f}/{rms_u:.1f}/{rms_3d:.1f} cm")
            print(f"  中位数: {median:.1f} cm, <1m: {good}/{n} ({good/n*100:.0f}%)")

            fpath = f"{args.out}/ppp_{ds}_apr29.csv"
            with open(fpath, 'w', newline='') as fh:
                keys = list(results[0].keys())
                csv.DictWriter(fh, fieldnames=keys).writeheader()
                csv.DictWriter(fh, fieldnames=keys).writerows(results)
            print(f"  保存: {fpath}")
            all_stats.append((ds, n, rms_e, rms_n, rms_u, rms_3d, median, good/n*100))
        else:
            print(f"  无结果")

        cur += delta

    if all_stats:
        print(f"\n{'日期':<14} {'N':<6} {'RMS_E':<8} {'RMS_N':<8} {'RMS_U':<8} {'RMS_3D':<10} {'中位数':<8} {'<1m'}")
        for ds, n, re, rn, ru, r3, med, pct in all_stats:
            print(f"{ds:<14} {n:<6} {re:>7.1f}  {rn:>7.1f}  {ru:>7.1f}  {r3:>8.1f}  {med:>7.1f}  {pct:.0f}%")

if __name__ == '__main__':
    main()