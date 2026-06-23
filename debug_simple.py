#!/usr/bin/env python3
"""Debug PPP convergence for Apr29"""
import pickle, sys, datetime as dt, math, numpy as np
sys.path.insert(0, '/workspace/gnss_pod/src')

C = 299792458.0; F1, F2 = 1575.42e6, 1227.60e6
F1_SQ, F2_SQ = F1*F1, F2*F2
ALPHA = F1_SQ/(F1_SQ-F2_SQ); BETA = -F2_SQ/(F1_SQ-F2_SQ)
MU_E = 3.986004418e14; OMEGA_E = 7.2921151467e-5
GPS_ORIGIN = dt.datetime(1980, 1, 6)
A_WGS84 = 6378137.0; E2_WGS84 = 0.00669437999014

GPS_SV = [(1,0.0,0.0,0.0,2e-5,55.0,0.0,5153.6),(3,30.0,30.0,60.0,1.5e-5,54.8,30.0,5153.6),(5,90.0,90.0,0.0,2e-5,55.0,0.0,5153.6),(7,150.0,150.0,0.0,1.5e-5,54.8,30.0,5153.6),(8,180.0,180.0,0.0,1e-5,55.0,0.0,5153.6),(9,210.0,210.0,0.0,2e-5,55.0,0.0,5153.6),(10,240.0,240.0,0.0,1.5e-5,54.8,0.0,5153.6),(11,270.0,270.0,0.0,1e-5,55.0,0.0,5153.6),(13,300.0,300.0,0.0,1.5e-5,54.8,0.0,5153.6),(14,330.0,330.0,0.0,1e-5,55.0,0.0,5153.6),(15,0.0,0.0,0.0,2e-5,55.0,0.0,5153.6),(17,30.0,30.0,60.0,1e-5,54.8,30.0,5153.6),(18,90.0,90.0,0.0,2e-5,55.0,0.0,5153.6),(19,120.0,120.0,0.0,1.5e-5,54.8,0.0,5153.6),(20,150.0,150.0,0.0,1e-5,55.0,0.0,5153.6),(21,180.0,180.0,0.0,2e-5,55.0,0.0,5153.6),(22,210.0,210.0,0.0,1.5e-5,54.8,0.0,5153.6),(23,240.0,240.0,0.0,1e-5,55.0,0.0,5153.6),(24,270.0,270.0,0.0,2e-5,54.8,0.0,5153.6),(25,300.0,300.0,0.0,1.5e-5,55.0,0.0,5153.6),(26,330.0,330.0,0.0,1e-5,54.8,0.0,5153.6),(27,0.0,0.0,0.0,2e-5,55.0,0.0,5153.6),(28,30.0,30.0,60.0,1.5e-5,54.8,30.0,5153.6),(29,90.0,90.0,0.0,1e-5,55.0,0.0,5153.6),(30,120.0,120.0,0.0,2e-5,54.8,0.0,5153.6),(31,150.0,150.0,0.0,1e-5,55.0,0.0,5153.6),(32,180.0,180.0,0.0,1.5e-5,54.8,0.0,5153.6)]
SV_PLAN = {r[0]: r for r in GPS_SV}

def gps_pos(sv_id, t_gps, t0=767620800.0):
    prn = int(sv_id[1:]) if isinstance(sv_id, str) else sv_id
    if prn not in SV_PLAN: return None
    r = SV_PLAN[prn]
    _, M0_d, Omega0_d, omega_d, ecc, inc_d, u0_d, sqrtA = r
    a = sqrtA**2; n0 = math.sqrt(MU_E / a**3)
    M = math.radians(M0_d) + n0*(t_gps - t0)
    E = M
    for _ in range(10): E = M + ecc*math.sin(E)
    v = 2*math.atan2(math.sqrt(1.00004)*math.sin(E/2), math.cos(E/2)-0.00002)
    phi = v + math.radians(omega_d); r_m = a*(1 - ecc*math.cos(E))
    Om = math.radians(Omega0_d) - OMEGA_E*(t_gps - t0)
    inc = math.radians(inc_d)
    X = r_m*(math.cos(phi)*math.cos(Om) - math.sin(phi)*math.cos(inc)*math.sin(Om))
    Y = r_m*(math.cos(phi)*math.sin(Om) + math.sin(phi)*math.cos(inc)*math.cos(Om))
    Z = r_m*math.sin(phi)*math.sin(inc)
    return (X, Y, Z)

def geo_range(a, b):
    return math.sqrt((b[0]-a[0])**2+(b[1]-a[1])**2+(b[2]-a[2])**2)

def ecef_to_blh(pos):
    X, Y, Z = float(pos[0]), float(pos[1]), float(pos[2])
    p = math.sqrt(X**2+Y**2)
    if p < 1e-12: return (math.pi/2*(1 if Z>=0 else -1), 0.0, abs(Z)-A_WGS84)
    e2 = E2_WGS84
    val = 1 - e2*(Z/p)**2
    lat = math.atan2(Z, p/math.sqrt(val)) if val > 0 else math.pi/2*(1 if Z>=0 else -1)
    lon = math.atan2(Y, X); sinL = math.sin(lat)
    return (lat, lon, 0.0)

def ecef_to_enu(err, lat, lon):
    cl, sl, cn, sn = math.cos(lat), math.sin(lat), math.cos(lon), math.sin(lon)
    return (-sl*cn*err[0]-sl*sn*err[1]+cl*err[2], sl*sn*err[0]-sl*cn*err[1], cl*cn*err[0]+cl*sn*err[1]+sl*err[2])

# Load data
gps1b = pickle.load(open('data/gracefo/2024/2024-04-29/gps1b_C.pkl', 'rb'))
from batch_v12 import parse_gnv1b
gnv = parse_gnv1b('data/gracefo/2024/2024-04-29/GNV1B_2024-04-29_C_04.txt')

gnv_keys = sorted(gnv.keys())
grace_ref = np.array(gnv[gnv_keys[0]])
ref_lat, ref_lon, _ = ecef_to_blh(grace_ref)
print(f"GRACE: {grace_ref}")
print(f"B/L: lat={math.degrees(ref_lat):.4f}, lon={math.degrees(ref_lon):.4f}")

# First epoch
t_gps = sorted(gps1b.keys())[0]
utc = GPS_ORIGIN + dt.timedelta(seconds=t_gps)
obs_dict = gps1b[t_gps]
print(f"\nFirst epoch: t_gps={t_gps}, UTC={utc}, {len(obs_dict)} sats")

# GRACE ref
rp0 = rp1 = None
for tk in gnv_keys:
    if tk >= utc: rp1 = tk; break
    rp0 = tk
if rp0 is None: rp0 = gnv_keys[0]
if rp1 is None: rp1 = gnv_keys[-1]
if rp0 == rp1: gp = np.array(gnv[rp0])
else:
    dt0 = (utc - rp0).total_seconds()
    dtt = (rp1 - rp0).total_seconds()
    gp = np.array(gnv[rp0])*(1-dt0/dtt) + np.array(gnv[rp1])*(dt0/dtt)
print(f"GRACE ref: {gp}")

# Build obs list
obs_list = []
for sv_id, rec in obs_dict.items():
    if rec.get('L1') is None or rec.get('el') is None: continue
    sp = gps_pos(sv_id, t_gps)
    if sp is None: continue
    L1 = rec['L1']; L2 = rec['L2']
    L_if = ALPHA*L1 + BETA*L2
    el = math.radians(max(rec['el'], 5.0))
    rho = geo_range(gp, sp)
    obs_list.append({'sv': sv_id, 'sat_pos': sp, 'L_if': L_if, 'el': el, 'rho0': rho})

print(f"Obs list: {len(obs_list)} sats")
for o in obs_list[:3]:
    print(f"  {o['sv']}: |GPS|={np.linalg.norm(o['sat_pos'])/1000:.0f}km, rho0={o['rho0']/1000:.0f}km, L_if={o['L_if']/1000:.0f}km, L_if-rho0={(o['L_if']-o['rho0'])/1000:.0f}km")

# Solve
x0 = np.array([gp[0], gp[1], gp[2], 0.0])
print(f"\nInitial x0: {x0[:3]}")
for it in range(15):
    H, y, W = [], [], []
    for o in obs_list:
        sp = o['sat_pos']; rcv = x0[:3]
        rho = geo_range(rcv, sp)
        if rho < 1e6:
            o['rho_iter'] = rho
            continue
        sigma = 0.003/max(math.sin(o['el']), 0.1)
        w = 1.0/sigma**2
        hx = -(sp[0]-rcv[0])/rho; hy = -(sp[1]-rcv[1])/rho; hz = -(sp[2]-rcv[2])/rho
        H.append([hx, hy, hz, 1.0])
        y.append(o['L_if'] - rho)
        W.append(w)
        o['rho_iter'] = rho
    if len(H) < 4: print(f"  iter {it}: <4 obs"); break
    HT = np.array(H).T; W_arr = np.diag(W)
    HTWH = HT @ W_arr @ np.array(H) + np.eye(4)*1e-8
    HTWy = HT @ W_arr @ np.array(y)
    try: dx = np.linalg.solve(HTWH, HTWy)
    except: print(f"  iter {it}: solve failed"); break
    if np.any(np.isnan(dx)) or np.any(np.isinf(dx)): print(f"  iter {it}: NaN"); break
    x0 = x0 + dx
    dr = np.linalg.norm(dx[:3])
    print(f"  iter {it}: |dx|={dr:.4f}m, clock={dx[3]:.2f}m, pos=({x0[0]:.0f},{x0[1]:.0f},{x0[2]:.0f})")
    if dr < 0.05 and abs(dx[3]) < 0.05: print(f"  → Converged!"); break

pos_est = x0[:3]
err = pos_est - gp
enu = ecef_to_enu(err, ref_lat, ref_lon)
d3 = math.sqrt(sum(x**2 for x in enu))
print(f"\nFinal: dE={enu[0]*100:.2f}cm, dN={enu[1]*100:.2f}cm, dU={enu[2]*100:.2f}cm, d3={d3*100:.2f}cm")