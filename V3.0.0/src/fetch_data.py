"""
GNSS PPP 数据获取模块（纯标准库版，无外部依赖）
===============================================
仅依赖: numpy, urllib (stdlib)
"""

import os, sys, zipfile, shutil, warnings
from pathlib import Path
from datetime import datetime, timedelta
from typing import Optional, Tuple, List, Dict
import urllib.request
import numpy as np

C       = 299792458.0
MU_E    = 3.9860050e14
OMEGA_E = 7.2921151467e-5

IGS_STATIONS = {
    "wuhn": {"name": "Wuhan",    "lat": 30.5317, "lon": 114.3572},
    "bjfs": {"name": "Beijing",  "lat": 39.6086, "lon": 115.8925},
}


def ymd_to_doy(year, month, day):
    return (datetime(year, month, day) - datetime(year, 1, 1)).days + 1

def doy_to_ymd(year, doy):
    d = datetime(year, 1, 1) + timedelta(days=doy - 1)
    return d.year, d.month, d.day

def datetime_to_mjd(dt):
    return (dt - datetime(1858, 11, 17)).total_seconds() / 86400.0

def blh_to_ecef(lat, lon, h=0.0):
    a, f = 6378137.0, 1/298.257223563
    N = a / np.sqrt(1 - 2*f + f*f * np.sin(lat)**2)
    return np.array([
        (N+h)*np.cos(lat)*np.cos(lon),
        (N+h)*np.cos(lat)*np.sin(lon),
        (N*(1-f*f) + h)*np.sin(lat)
    ])

def ecef_to_blh(pos):
    X, Y, Z = pos
    a, f = 6378137.0, 1/298.257223563
    e2 = 2*f - f*f
    p = np.sqrt(X*X + Y*Y)
    lon = np.arctan2(Y, X)
    lat = np.arctan2(Z, p/np.sqrt(1-e2*(Z/p)**2))
    for _ in range(5):
        N = a/np.sqrt(1-e2*np.sin(lat)**2)
        H = p/np.cos(lat) - N
        lat = np.arctan2(Z, p*(1-e2*N/(N+H)))
    return lat, lon, H

def compute_satpos_nav(brdc, ts):
    """从广播星历计算卫星 ECEF 位置（ICD-200 算法）"""
    toe   = brdc.get('toe',  0.0)
    sqrtA = brdc.get('sqrtA', 0.0)
    ecc   = brdc.get('e',     0.0)
    M0    = brdc.get('M0',    0.0)
    i0    = brdc.get('i0',    0.0)
    omega = brdc.get('omega', 0.0)
    Omega0= brdc.get('Omega0',0.0)
    dN    = brdc.get('dN',   0.0)
    iDot  = brdc.get('iDot',  0.0)
    OmegaDot= brdc.get('OmegaDot',0.0)
    cuc, cus = brdc.get('cuc',0.0), brdc.get('cus',0.0)
    crc, crs = brdc.get('crc',0.0), brdc.get('crs',0.0)
    cic, cis = brdc.get('cic',0.0), brdc.get('cis',0.0)
    af0, af1 = brdc.get('af0',0.0), brdc.get('af1',0.0)

    a  = sqrtA**2
    n0 = np.sqrt(MU_E/a**3)
    n  = n0 + dN

    gps_origin = datetime(1980, 1, 6)
    gps_secs   = (ts - gps_origin).total_seconds()
    t_k = (gps_secs % 604800) - toe
    if t_k >  302400: t_k -= 604800
    if t_k < -302400: t_k += 604800

    M = M0 + n*t_k
    E = M
    for _ in range(10): E = M + ecc*np.sin(E)
    sinE, cosE = np.sin(E), np.cos(E)
    v = np.arctan2(np.sqrt(1-ecc**2)*sinE, cosE-ecc)
    phi = v + omega

    du = cuc*np.cos(2*phi) + cus*np.sin(2*phi)
    dr = crc*np.cos(2*phi) + crs*np.sin(2*phi)
    di = cic*np.cos(2*phi) + cis*np.sin(2*phi)

    u = phi + du
    r = a*(1 - ecc*cosE) + dr
    i = i0 + iDot*t_k + di
    xp = r*np.cos(u)
    yp = r*np.sin(u)
    Om = Omega0 + (OmegaDot - OMEGA_E)*t_k - OMEGA_E*toe

    pos = np.array([
        xp*np.cos(Om) - yp*np.cos(i)*np.sin(Om),
        xp*np.sin(Om) + yp*np.cos(i)*np.cos(Om),
        yp*np.sin(i)
    ])
    clk = af0 + af1*t_k
    return pos, clk


class RinexReader:
    """读取 RINEX 2.11 观测文件"""
    def __init__(self, filepath):
        self.filepath = filepath
        self.version  = 2.11
        self.marker   = ""
        self.obs_types= []
        self.interval= 30.0
        self.epochs   = []
        self._parse()

    def _parse(self):
        with open(self.filepath, encoding='utf-8', errors='replace') as f:
            lines = f.readlines()
        h_end = next((i for i,l in enumerate(lines) if 'END OF HEADER' in l[60:]), len(lines))
        self._parse_header(lines[:h_end+1])
        self._parse_body(lines[h_end+1:])

    def _parse_header(self, lines):
        for line in lines:
            tag = line[60:].strip()
            if tag == 'MARKER NAME': self.marker = line[:20].strip()
            elif tag == '# / TYPES OF OBS' or tag == 'OBS TYPES':
                n = int(line[:6].strip())
                toks = line[6:].split()
                self.obs_types = toks[:n]
            elif tag == 'INTERVAL':
                try: self.interval = float(line[:10].strip())
                except: pass

    def _parse_body(self, lines):
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if not line or line.startswith('COMMENT'): i+=1; continue
            parts = line.split()
            if len(parts) < 7: i+=1; continue
            try:
                yr = int(parts[0])
                if yr<80: yr+=2000
                elif yr<100: yr+=1900
                mo,dy,hr,mn = int(parts[1]),int(parts[2]),int(parts[3]),int(parts[4])
                sec = float(parts[5])
                epoch = datetime(yr,mo,dy,hr,mn,int(sec),int((sec%1)*1e6))
                i+=1; sat_obs={}
                for _ in range(100):
                    if i>=len(lines): break
                    raw = lines[i]
                    if raw.strip() and raw[0].isdigit(): break
                    pos=0; k=0
                    obs_vals={}
                    while pos<len(raw) and k<len(self.obs_types):
                        tok=raw[pos:pos+14].strip()
                        try: obs_vals[self.obs_types[k]]=float(tok) if tok not in ('','0','0.0') else float('nan')
                        except: obs_vals[self.obs_types[k]]=float('nan')
                        pos+=14; k+=1
                    sv_id=raw[:3].strip()
                    if sv_id: sat_obs[sv_id]=obs_vals
                    i+=1
                self.epochs.append({'time':epoch,'obs':sat_obs})
            except: i+=1; continue

    def __iter__(self): return iter(self.epochs)
    def __len__(self): return len(self.epochs)


class NavReader:
    """读取 RINEX 广播星历文件"""
    def __init__(self, filepath):
        self.filepath=filepath
        self.records=[]
        self._parse()

    def _parse(self):
        with open(self.filepath, encoding='utf-8', errors='replace') as f:
            lines=f.readlines()
        i=0
        while i<len(lines) and 'END OF HEADER' not in lines[i][60:]: i+=1
        i+=1
        while i<len(lines):
            line=lines[i].strip()
            if not line or line.startswith('COMMENT'): i+=1; continue
            parts=line.split()
            if len(parts)<8: i+=1; continue
            try:
                yr=int(parts[0]); mo=int(parts[1]); dy=int(parts[2])
                hr=int(parts[3]); mn=int(parts[4]); sc=float(parts[5])
                sv=parts[6]
                if yr<80: yr+=2000
                elif yr<100: yr+=1900
                epoch=datetime(yr,mo,dy,hr,mn,int(sc))
                i+=1; p1=lines[i].strip().split(); i+=1
                p2=lines[i].strip().split(); i+=1
                p3=lines[i].strip().split(); i+=1
                def f(lst,idx,d=0.0):
                    try: return float(lst[idx])
                    except: return d
                self.records.append({
                    'sv':sv,'time':epoch,
                    'M0':f(p1,0),'e':f(p1,1),'sqrtA':f(p1,2),
                    'i0':f(p1,3),'Omega0':f(p1,4),'omega':f(p1,5),
                    'OmegaDot':f(p1,6),'iDot':f(p1,7),'dN':f(p1,8),
                    'cuc':f(p1,9),'cus':f(p1,10),'crc':f(p1,11),'crs':f(p1,12),
                    'cic':f(p1,13),'cis':f(p1,14),'toe':f(p1,15),
                    'tgd':f(p1,16),'af0':f(p1,17),'af1':f(p1,18),'af2':f(p1,19),
                })
            except: i+=1; continue

    def get_nav(self, sv, time, max_age=7200):
        best,best_dt=None,float('inf')
        for r in self.records:
            if r['sv'][:3]!=sv[:3]: continue
            dt=abs((r['time']-time).total_seconds())
            if dt<best_dt and dt<max_age:
                best_dt=dt; best=r
        return best


class Sp3Reader:
    """读取 IGS SP3 精密星历文件"""
    def __init__(self, filepath):
        self.filepath=filepath
        self.epochs={}
        self.svs=[]
        self._parse()

    def _parse(self):
        with open(self.filepath, encoding='utf-8', errors='replace') as f:
            lines=f.readlines()
        cur_t=None; svs=set()
        for line in lines:
            if line[0]=='*':
                yr=int(line[3:7]);mo=int(line[8:10]);d=int(line[11:13])
                h=int(line[14:16]);mn=int(line[17:19]);sc=float(line[20:31])
                cur_t=datetime(yr,mo,d,h,mn,int(sc))
                self.epochs[cur_t]=[]
            elif line[0] in 'PV':
                try:
                    sv=line[2:5].strip()
                    x,y,z,clk=float(line[6:22]),float(line[22:38]),float(line[38:54]),float(line[54:70])
                    svs.add(sv)
                    self.epochs[cur_t].append({'sv':sv,'pos':np.array([x,y,z]),'clock':clk*1e-6*C})
                except: pass
        self.svs=sorted(svs)

    def get_sat_pos(self, sv, time):
        times=sorted(self.epochs.keys())
        if not times: return np.zeros(3),0.0
        t0=t1=None
        for t in times:
            if t<=time: t0=t
            if t>=time and t1 is None: t1=t
        if t0 is None: t0=t1
        if t1 is None: t1=t0
        def get(t):
            for s in self.epochs.get(t,[]):
                if s['sv']==sv: return s['pos'].copy(),s['clock']
            return None,None
        p0,c0=get(t0); p1,c1=get(t1)
        if p0 is None: return np.zeros(3),0.0
        if p1 is None or np.allclose(p0,np.zeros(3)): return p0,c0 or 0.0
        w=(time-t0).total_seconds()/max((t1-t0).total_seconds(),1)
        return (1-w)*p0+w*p1, (1-w)*(c0 or 0)+w*(c1 or 0)


class IGSDataDownloader:
    """IGS 数据下载器（urllib 版）"""
    def __init__(self, data_dir="./data", timeout=60):
        self.data_dir=Path(data_dir); self.data_dir.mkdir(exist_ok=True)
        self.timeout=timeout

    def _download(self, url, local):
        try:
            urllib.request.urlretrieve(url, local)
            return Path(local).stat().st_size>1000
        except Exception as e:
            warnings.warn(f"下载失败 {url}: {e}"); return False

    def _decompress(self, gz_path, out_path):
        try:
            import gzip
            with gzip.open(gz_path,'rb') as fi:
                with open(out_path,'wb') as fo: shutil.copyfileobj(fi,fo)
            return True
        except: return False

    def download_rinex(self, station, year, doy):
        yy=str(year)[2:]; yyyy=str(year)
        ddir=self.data_dir/yyyy/f"{doy:03d}"; ddir.mkdir(exist_ok=True,parents=True)
        base="https://cddis.nasa.gov/archive/gnss/data/daily"
        for ext in ['.rnx.gz','.rnx','.crx.gz']:
            suf=ext.replace('.gz',''); fname=f"{station}{doy:03d}0.{yy}{ext}"
            url=f"{base}/{yyyy}/{doy:03d}/{yy}o/{fname}"
            gz=ddir/fname; local=ddir/fname.replace('.gz','')
            if local.exists() and local.stat().st_size>1000: return str(local)
            if self._download(url,gz):
                if ext.endswith('.gz') and self._decompress(gz,local): return str(local)
                elif not ext.endswith('.gz'): return str(gz)
        return None

    def download_nav(self, year, doy):
        yy=str(year)[2:]; yyyy=str(year)
        ddir=self.data_dir/yyyy/f"{doy:03d}"; ddir.mkdir(exist_ok=True,parents=True)
        base="https://cddis.nasa.gov/archive/gnss/data/daily"
        fname=f"brdc{doy:03d}0.{yy}n.gz"
        url=f"{base}/{yyyy}/{doy:03d}/{yy}n/{fname}"
        gz=ddir/fname; local=ddir/fname.replace('.gz','')
        if self._download(url,gz):
            self._decompress(gz,local); return str(local)
        return None


if __name__=='__main__':
    print("数据模块 OK (stdlib+urllib)")
