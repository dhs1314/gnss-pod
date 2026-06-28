#!/usr/bin/env python3
"""Debug single PPP epoch to find root cause of large errors."""
import pickle, numpy as np
from datetime import datetime, timedelta

C=299792458.0; F1=1575.42e6; F2=1227.60e6
F1_SQ=F1*F1; F2_SQ=F2*F2
OMEGA_E=7.2921151467e-5; J2000=datetime(2000,1,1,12,0,0)
A=6378137.0; E2=0.00669437999014

sp3=pickle.load(open('data/2024/120/igs_sp3_FIN.pkl','rb'))
gps1b=pickle.load(open('data/gracefo/2024/2024-04-29/gps1b_C.pkl','rb'))
ref={}
for l in open('data/gracefo/2024/2024-04-29/GNV1B_2024-04-29_C_04.txt'):
    if l.startswith('#') or not l.strip(): continue
    p=l.split()
    if len(p)<6: continue
    try:
        t=float(p[0]); flag=p[2]
        if flag in ('C','E'): ref[t]=np.array([float(p[3]),float(p[4]),float(p[5])])
    except: pass

# Find stable epoch (1 hour in)
target_gps = min(gps1b.keys()) + 3600
for gps_sod in sorted(gps1b.keys()):
    if gps_sod >= target_gps: break

utc = J2000 + timedelta(seconds=gps_sod - 18)
print(f'GPS J2000={gps_sod}, UTC={utc}')

# GRACE position
ts_ref=sorted(ref.keys()); t0=t1=None
for i,ti in enumerate(ts_ref):
    if ti>=gps_sod: t1=ti; t0=ts_ref[i-1] if i>0 else None; break
    t0=ti
if t1 is None: t0=t1=ts_ref[-1]
if t0 is None: t0=ts_ref[0]
if t0==t1: grace=ref[t0]
else: a=(gps_sod-t0)/(t1-t0); grace=ref[t0]*(1-a)+ref[t1]*a
print(f'GRACE: r={np.linalg.norm(grace)/1000:.1f}km')

def blh(p):
    X,Y,Z=float(p[0]),float(p[1]),float(p[2])
    r=np.sqrt(X**2+Y**2)
    if r<1e-6: lat=np.pi/2 if Z>=0 else -np.pi/2; lon=0.0
    else:
        lat=np.arctan2(Z,r)
        for _ in range(10):
            sl=np.sin(lat); N=A/np.sqrt(1-E2*sl**2)
            new_lat=np.arctan2(Z+E2*N*sl,r)
            if abs(new_lat-lat)<1e-15: lat=new_lat; break
            lat=new_lat
        lon=np.arctan2(Y,X)
    return np.array([lat,lon,0.0])

lat,lon,_=blh(grace)
R_enu=np.array([[-np.sin(lon),np.cos(lon),0],
                [-np.sin(lat)*np.cos(lon),-np.sin(lat)*np.sin(lon),np.cos(lat)],
                [np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)]])
print(f'lat={np.degrees(lat):.2f} lon={np.degrees(lon):.2f}')

# SP3 positions for each SV
print(f'\n{"SV":5s} {"rng(km)":>9s} {"L_if(km)":>10s} {"P_if(km)":>10s} {"clk(km)":>9s} {"O-C(km)":>9s} {"el":>8s}')
for sv,rec in sorted(gps1b[gps_sod].items()):
    sts=sp3['ts']; seps=sp3['epochs']; st0=st1=None
    for i,ti in enumerate(sts):
        if ti>=utc: st1=ti; st0=sts[i-1] if i>0 else None; break
        st0=ti
    if st1 is None: st0=st1=sts[-1]
    if st0 is None: st0=sts[0]
    sp0=seps[st0].get(sv); sp1=seps[st1].get(sv)
    if sp0 is None and sp1 is None: continue
    if sp0 is None: sat_pos=np.array(sp1[:3]); sat_clk=sp1[3]
    elif sp1 is None: sat_pos=np.array(sp0[:3]); sat_clk=sp0[3]
    else:
        aa=(utc-st0).total_seconds()/(st1-st0).total_seconds()
        sat_pos=np.array(sp0[:3])*(1-aa)+np.array(sp1[:3])*aa
        sat_clk=sp0[3]*(1-aa)+sp1[3]*aa

    if abs(sat_clk)>0.1*C: continue

    los=sat_pos-grace; rng=np.linalg.norm(los)
    e_enu=R_enu@(los/rng); el=np.arcsin(np.clip(e_enu[2],-1,1))
    if el<np.radians(7): continue

    L_if=rec['L_if']; P_if=rec['P_if']
    oc=(L_if - rng + sat_clk)/1000
    print(f'{sv:5s} {rng/1000:9.1f} {L_if/1000:10.3f} {P_if/1000:10.3f} {sat_clk/1000:9.1f} {oc:9.1f} {np.degrees(el):8.1f}')
