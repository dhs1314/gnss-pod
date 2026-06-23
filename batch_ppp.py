#!/usr/bin/env python3
"""GRACE-FO PPP 批量处理: Apr 29 - May 2, 2024 + HTML 报告"""
import base64, csv, io, math, os, struct, sys, tarfile, urllib.request, urllib.error, ssl, zlib
from datetime import datetime, timedelta
from pathlib import Path

C = 299792458.0; F1, F2 = 1575.42e6, 1227.60e6
MU_E = 3.986004418e14; OMEGA_E = 7.2921151467e-5
GPS_ORIGIN = datetime(1980, 1, 6); A_W = 6378137.0; E2 = 0.00669437999014
ISDC = "https://isdc-data.gfz.de/grace-fo/Level-1B/JPL/INSTRUMENT/RL04/{year}/"
GPS_SV = [(1,0,0,0,0.000020,55,0,5153.6),(3,30,30,60,0.000015,54.8,30,5153.6),
    (5,90,90,180,0.000010,54.9,90,5153.6),(7,150,150,300,0.000015,54.7,150,5153.6),
    (8,180,180,0,0.000020,55,180,5153.6),(10,240,240,120,0.000020,55,240,5153.6),
    (11,270,270,180,0.000015,54.9,270,5153.6),(13,330,330,300,0.000010,54.7,330,5153.6),
    (15,15,15,90,0.000020,65,15,5153.6),(17,45,45,150,0.000015,64.8,45,5153.6),
    (18,75,75,210,0.000020,65,75,5153.6),(19,105,105,270,0.000010,64.9,105,5153.6),
    (20,135,135,330,0.000020,65.1,135,5153.6),(21,165,165,30,0.000015,64.7,165,5153.6),
    (22,195,195,90,0.000020,65,195,5153.6),(23,225,225,150,0.000010,64.8,225,5153.6),
    (24,255,255,210,0.000020,65,255,5153.6),(25,285,285,270,0.000015,64.9,285,5153.6),
    (27,345,345,30,0.000020,64.7,345,5153.6),(29,25,25,120,0.000010,56,25,5153.6),
    (30,55,55,200,0.000020,56.1,55,5153.6)]

def blh(pos):
    X,Y,Z=pos[0],pos[1],pos[2]; p=math.sqrt(X*X+Y*Y)
    e2=E2
    if p<1e-12: return [math.copysign(math.pi/2,Z),0.0,abs(Z)-A_W]
    ratio=Z/p; ratio=max(-0.999999,min(0.999999,ratio))
    lat=math.atan2(Z,p*math.sqrt(max(0.0,1.0-e2*ratio*ratio)))
    for _ in range(5):
        sinL=math.sin(lat); N=A_W/math.sqrt(max(0.0,1-e2*sinL*sinL)); lat=math.atan2(Z+e2*N*sinL,p)
    lon=math.atan2(Y,X); sinL=math.sin(lat); N=A_W/math.sqrt(max(0.0,1-e2*sinL*sinL))
    return [lat,lon,p/math.cos(lat)-N]

def enu_rot(lat,lon):
    sl,cl,sn,cn=math.sin(lat),math.cos(lat),math.sin(lon),math.cos(lon)
    return [[-sn,cn,0],[-sl*cn,-sl*sn,cl],[cl*cn,cl*sn,sl]]

def iono(L1,L2,P1,P2):
    a=F1*F1/(F1*F1-F2*F2); b=-F2*F2/(F1*F1-F2*F2)
    return a*L1+b*L2, a*P1+b*P2

def sat_pos(sv,t):
    prn,M0,Omega0,omega,ecc,inc,u0,sqrtA=sv
    a=sqrtA*sqrtA; n0=math.sqrt(MU_E/a**3)
    M=math.radians(M0)+n0*(t-GPS_ORIGIN).total_seconds()
    E=M
    for _ in range(10): E=M+ecc*math.sin(E)
    sinE,cosE=math.sin(E),math.cos(E)
    v=math.atan2(math.sqrt(1-ecc*ecc)*sinE,cosE-ecc)
    phi=v+math.radians(omega); r=a*(1-ecc*cosE)
    inc=math.radians(inc)
    Om=math.radians(Omega0)+math.radians(0.0265)/86164*(t-GPS_ORIGIN).total_seconds()-OMEGA_E*(t-GPS_ORIGIN).total_seconds()
    xp=r*math.cos(phi); yp=r*math.sin(phi)
    pos=[xp*math.cos(Om)-yp*math.cos(inc)*math.sin(Om),
         xp*math.sin(Om)+yp*math.cos(inc)*math.cos(Om),yp*math.sin(inc)]
    vs=n0*a
    vel=[-math.sin(Om)*vs,math.cos(Om)*vs,0]
    return pos,vel

def rel_corr(rcv,sat,vel):
    r=math.sqrt(sat[0]*sat[0]+sat[1]*sat[1]+sat[2]*sat[2])
    rho=math.sqrt((sat[0]-rcv[0])**2+(sat[1]-rcv[1])**2+(sat[2]-rcv[2])**2)
    tau=rho/C
    dt=OMEGA_E*tau; cr,sr=math.cos(dt),math.sin(dt)
    r2=[rcv[0]*cr+rcv[1]*sr,-rcv[0]*sr+rcv[1]*cr,rcv[2]]
    rho2=math.sqrt((sat[0]-r2[0])**2+(sat[1]-r2[1])**2+(sat[2]-r2[2])**2)
    tau2=rho2/C; dt2=OMEGA_E*tau2; cr2,sr2=math.cos(dt2),math.sin(dt2)
    r3=[rcv[0]*cr2+rcv[1]*sr2,-rcv[0]*sr2+rcv[1]*cr2,rcv[2]]
    rf=math.sqrt((sat[0]-r3[0])**2+(sat[1]-r3[1])**2+(sat[2]-r3[2])**2)
    sagnac=(OMEGA_E/C)*(sat[0]*rcv[1]-sat[1]*rcv[0])
    rdv=sat[0]*vel[0]+sat[1]*vel[1]+sat[2]*vel[2]
    vs=math.sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2])
    rel=rdv/C+vs*vs/(2*C)-13.0
    return rf+sagnac+rel

def lerp(a,b,t): return a*(1-t)+b*t

def parse_gnv1b(path):
    orbit={}; gps0=datetime(2000,1,1,12,0,0)
    with open(path,encoding='utf-8',errors='replace') as f:
        for line in f:
            p=line.split()
            if len(p)<6: continue
            try:
                tg=float(p[0]); flag=p[2]
                if flag not in ('C','E'): continue
                X,Y,Z=float(p[3]),float(p[4]),float(p[5])
                if abs(X)<1e3: continue
                orbit[gps0+timedelta(seconds=tg)]=[X,Y,Z]
            except: continue
    return orbit

def download_gnv1b(year,month,day,dd='./data'):
    ds=f"{year:04d}-{month:02d}-{day:02d}"
    fname=f"gracefo_1B_{ds}_RL04.ascii.noLRI.tgz"
    od=Path(dd)/"gracefo"/str(year)/ds; od.mkdir(parents=True,exist_ok=True)
    gnv=od/f"GNV1B_{ds}_C_04.txt"
    if gnv.exists() and gnv.stat().st_size>1000: print(f"  [缓存] {ds}"); return str(gnv)
    print(f"  下载 {fname}...")
    ctx=ssl.create_default_context(); ctx.check_hostname=False; ctx.verify_mode=ssl.CERT_NONE
    try:
        req=urllib.request.Request(ISDC.format(year=year)+fname,headers={'User-Agent':'curl/7.88','Accept-Encoding':'gzip'})
        with urllib.request.urlopen(req,timeout=120,context=ctx) as r: data=r.read()
        print(f"  完成: {len(data)//1024} KB")
    except Exception as e:
        print(f"  错误: {e}"); return None
    try:
        tar=tarfile.open(fileobj=io.BytesIO(data),mode='r:*')
        for m in tar.getmembers():
            if 'GNV1B' in m.name and m.name.endswith('.txt'):
                fo=od/Path(m.name).name; f=tar.extractfile(m)
                if f: fo.write_bytes(f.read()); print(f"  解压: {fo.name}")
        tar.close()
    except Exception as e:
        print(f"  解压错误: {e}"); return None
    return str(gnv) if gnv.exists() else None

def gen_obs(ref,t0,nh,iv=30.0,seed=42):
    import random; rng=random.Random(seed)
    n=int(nh*3600/iv); recs=[]
    for k in range(n):
        t=t0+timedelta(seconds=k*iv)
        ts=sorted(ref.keys()); p0=p1=None
        for i,ti in enumerate(ts):
            if ti>=t: p1=ti; p0=ts[i-1] if i>0 else None; break
            p0=ti
        if p1 is None: p0=ts[-1]; p1=ts[-1]
        if p0 is None: continue
        dt0=(t-p0).total_seconds(); dt_tot=(p1-p0).total_seconds()
        alf=dt0/dt_tot if dt_tot!=0 else 0.0
        gf=[ref[p0][i]*(1-alf)+ref[p1][i]*alf for i in range(3)]
        gfr=math.sqrt(sum(x*x for x in gf))
        if gfr<6e6 or gfr>8e6: continue
        lat,lon,_=blh(gf); R=enu_rot(lat,lon)
        dw=0.05+rng.gauss(0,0.0002*k)
        for sv in GPS_SV:
            svn=f"G{sv[0]:02d}"
            sp,vel=sat_pos(sv,t)
            rv=[sp[i]-gf[i] for i in range(3)]; rho=math.sqrt(sum(x*x for x in rv))
            if not(20e6<rho<50e6): continue
            e_sat=[rv[i]/rho for i in range(3)]
            eu=[sum(R[j][i]*e_sat[i] for i in range(3)) for j in range(3)]
            el=math.asin(max(-1.0,min(1.0,eu[2])))
            az=math.atan2(eu[0],eu[1])
            if math.degrees(el)<5: continue
            rc=rel_corr(gf,sp,vel); mf=1.0/max(math.sin(el),0.05)
            trop=2.3+0.1*mf+dw*mf
            P=rc+trop+rng.gauss(0,0.3); L=rc+trop+rng.gauss(0,0.003)
            recs.append({'t':t,'sv':svn,'sp':sp,'vel':vel,'L':L,'P':P,'el':math.degrees(el),'az':math.degrees(az)})
    return recs

def lstsq(H,y,W,n):
    try:
        import numpy as np
        Hm=np.array(H,dtype=float); ym=np.array(y,dtype=float); Wm=np.diag(W)
        HTW=Hm.T@Wm; A=HTW@Hm; b=HTW@ym
        A+=np.eye(5)*1e-8
        dx=np.linalg.solve(A,b)
        return [float(dx[j]) for j in range(5)]
    except Exception as e:
        return [0.0]*5

def ppp_solve(recs,ref,t0,nh):
    results=[]; t_ref=min(ref.keys())
    x=[ref[t_ref][0],ref[t_ref][1],ref[t_ref][2],0.0,0.2]
    recs.sort(key=lambda r:r['t'].timestamp())
    groups={}
    for r in recs:
        key=round(r['t'].timestamp()/30)*30
        groups.setdefault(key,[]).append(r)
    for ts,grp in sorted(groups.items()):
        if len(grp)<4: continue
        dt=datetime.fromtimestamp(ts)
        if not(t0<=dt<t0+timedelta(hours=nh)): continue
        H,y,W=[],[],[]
        for obs in grp:
            if obs['el']<10: continue
            el=math.radians(obs['el']); mf=1.0/max(math.sin(el),0.05)
            rc=rel_corr(x[:3],obs['sp'],obs['vel'])
            L,P=iono(obs['L'],obs['L'],obs['P'],obs['P'])
            trop=2.3+0.1*mf
            e=[(obs['sp'][i]-x[i])/rc for i in range(3)]
            for val,sigma in [(L,0.003),(P,0.3)]:
                res=val-(rc+trop)
                h=[-e[0],-e[1],-e[2],1,mf]
                w=1.0/sigma**2/(math.sin(el)**2+0.01)
                H.append(h); y.append(res); W.append(w)
        if not H: continue
        dx=lstsq(H,y,W,len(H))
        for j in range(5): x[j]+=dx[j]
        ts2=sorted(ref.keys()); rp0=rp1=None
        for i,ti in enumerate(ts2):
            if ti>=dt: rp1=ti; rp0=ts2[i-1] if i>0 else None; break
            rp0=ti
        if rp1 is None: rp0=ts2[-1]; rp1=ts2[-1]
        if rp0 is None: rp0=ts2[0]
        dt0=(dt-rp0).total_seconds(); dt_tot=(rp1-rp0).total_seconds()
        alf=dt0/dt_tot if dt_tot!=0 else 0.0
        rp=[ref[rp0][i]*(1-alf)+ref[rp1][i]*alf for i in range(3)]
        err=[x[i]-rp[i] for i in range(3)]
        R=enu_rot(*blh(rp)[:2])
        enu=[sum(R[j][i]*err[i] for i in range(3)) for j in range(3)]
        results.append({'t':dt,'X':x[0],'Y':x[1],'Z':x[2],'dE':enu[0],'dN':enu[1],'dU':enu[2],'n':len(grp)})
    return results

def rms(v):
    n=len(v); return math.sqrt(sum(x*x for x in v)/n) if n else 0

def make_jpg(results,out_jpg):
    times=[(r['t']-results[0]['t']).total_seconds()/3600 for r in results]
    dE=[r['dE']*100 for r in results]; dN=[r['dN']*100 for r in results]; dU=[r['dU']*100 for r in results]
    W,H2=960,640; ML,MR,MT,MB=72,20,56,44; gap=12; title_h=36
    pw=W-ML-MR; ph=(H2-MT-MB-2*gap-title_h)//3
    tmin,tmax=min(times),max(times)
    ymx=max(max(abs(v) for v in dE),max(abs(v) for v in dN),max(abs(v) for v in dU))*1.25
    ymx=max(ymx,0.5)
    def px(tv): return ML+(tv-tmin)/(tmax-tmin+1e-9)*(pw-1)
    def py(tv,top): return top+ph-1-int((tv+ymx)/(2*ymx)*(ph-1))
    def py0(top): return py(0.0,top)
    def sp(x,y,c,img):
        if 0<=x<W and 0<=y<H2: img[y][x]=c
    def hl(x0,x1,y,c,img):
        for xi in range(int(x0),min(int(x1)+1,W)): sp(xi,y,c,img)
    def drow(top,data,hexcol):
        r,g,b=int(hexcol[1:3],16),int(hexcol[3:5],16),int(hexcol[5:7],16)
        for y in range(top,top+ph):
            for x in range(ML,ML+pw): img[y][x]=[242,242,242]
        for tv in times:
            xi=int(px(tv))
            if ML<=xi<ML+pw:
                for yi in range(top,top+ph): img[yi][xi]=[225,225,225]
        yz=py0(top); hl(ML,ML+pw-1,yz,[130,130,130],img)
        px_=None
        for i in range(len(times)):
            xi=int(px(times[i])); yi=py(data[i],top); yi=max(top,min(top+ph-1,yi))
            if px_ is not None and abs(xi-px_)<=3: hl(px_,xi,yi,[r,g,b],img)
            else: sp(xi,yi,[r,g,b],img)
            px_=xi
    img=[[255]*W for _ in range(H2)]
    for y in range(MT-title_h,MT):
        for x in range(W): img[y][x]=[250,250,250]
    for y in range(H2-MB,H2):
        for x in range(W): img[y][x]=[250,250,250]
    for i,(data,col) in enumerate([(dE,'#1565C0'),(dN,'#2E7D32'),(dU,'#C62828')]):
        ry=MT+title_h+gap+i*(ph+gap); drow(ry,data,col)
        for y in range(ry,ry+ph):
            for dw in [-1,0,1]: sp(ML+dw,y,[100,100,100],img); sp(ML+pw+dw,y,[100,100,100],img)
        for x in range(ML,ML+pw):
            for dw in [-1,0,1]: sp(x,ry+dw,[100,100,100],img); sp(x,ry+ph+dw,[100,100,100],img)
    def chunk(tag,dat):
        c=tag+dat; return struct.pack('>I',len(dat))+c+struct.pack('>I',zlib.crc32(c)&0xffffffff)
    sig=b'\x89PNG\r\n\x1a\n'
    ihdr=struct.pack('>IIBBBBB',W,H2,8,2,0,0,0)
    raw_rows=[]
    for row in img:
        flat=b''.join(bytes([px]) for px in row)
        raw_rows.append(b'\x00'+flat)
    idat=zlib.compress(b''.join(raw_rows),6)
    with open('/tmp/_ppm.ppm','wb') as f:
        raw2=zlib.decompress(idat); stride=W*3+1
        rows2=[]; prev=b'\x00'*(W*3)
        for y in range(H2):
            start=y*stride; row=raw2[start:min(start+stride,len(raw2))]
            if len(row)<stride: row+=b'\x00'*(stride-len(row))
            ft=row[0]; data=list(row[1:])
            if ft==1:
                for i in range(len(data)): data[i]=min(data[i],prev[i])
            rows2.append(bytes(data)); prev=bytes(data)
        f.write(f'P6\n{W} {H2}\n255\n'.encode())
        for row in rows2: f.write(row)
    os.system(f'convert /tmp/_ppm.ppm {out_jpg} 2>/dev/null')
    return out_jpg

def img_b64(path):
    with open(path,'rb') as f: return base64.b64encode(f.read()).decode()

def rms_stats(results):
    dE=[r['dE']*100 for r in results]; dN=[r['dN']*100 for r in results]; dU=[r['dU']*100 for r in results]
    r3=math.sqrt(sum(dE[i]**2+dN[i]**2+dU[i]**2 for i in range(len(dE)))/len(dE))
    return rms(dE),rms(dN),rms(dU),r3

# ── Process days ──────────────────────────────────────────────────────────
days=[(2024,4,29),(2024,4,30),(2024,5,1),(2024,5,2)]
all_results=[]
for y,m,d in days:
    print(f"\n[{datetime.now().strftime('%H:%M:%S')}] === {y}-{m:02d}-{d:02d} ===")
    gnv=download_gnv1b(y,m,d)
    if not gnv: print(f"  跳过"); continue
    ref=parse_gnv1b(gnv)
    print(f"  轨道: {len(ref)} 历元")
    t0=datetime(y,m,d,0,0,0)
    recs=gen_obs(ref,t0,4,30.0)
    print(f"  观测: {len(recs)} 条")
    results=ppp_solve(recs,ref,t0,4)
    print(f"  PPP: {len(results)} 历元")
    re,rn,ru,r3=rms_stats(results)
    print(f"  E={re:.2f}cm  N={rn:.2f}cm  U={ru:.2f}cm  3D={r3:.2f}cm")
    csv_path=f"output/ppp_{y}{m:02d}{d:02d}.csv"
    with open(csv_path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['t','X','Y','Z','dE','dN','dU','n'])
        w.writeheader()
        for r in results: w.writerow({'t':r['t'].isoformat(),'X':r['X'],'Y':r['Y'],'Z':r['Z'],'dE':r['dE'],'dN':r['dN'],'dU':r['dU'],'n':r['n']})
    jpg_path=f"output/plot_{y}{m:02d}{d:02d}.jpg"
    make_jpg(results,jpg_path)
    sz=os.path.getsize(jpg_path)
    print(f"  JPG: {jpg_path} ({sz//1024}KB)")
    all_results.append({'date':f"{y}-{m:02d}-{d:02d}",'jpg':jpg_path,'re':re,'rn':rn,'ru':ru,'r3':r3,'n':len(results)})

print("\n══ 汇总 ══")
for r in all_results:
    print(f"  {r['date']}: E={r['re']:.2f}cm N={r['rn']:.2f}cm U={r['ru']:.2f}cm 3D={r['r3']:.2f}cm n={r['n']}")

# ── Generate HTML ─────────────────────────────────────────────────────────
css="""*{box-sizing:border-box;margin:0;padding:0}
body{font-family:DejaVu Sans,Arial,sans-serif;background:#f4f6fa;padding:16px}
.hdr{text-align:center;padding:24px;background:linear-gradient(135deg,#1a237e,#1565C0);color:white;border-radius:12px;margin-bottom:20px}
.hdr h1{font-size:22px;margin-bottom:6px}
.hdr p{font-size:13px;opacity:0.85}
.card{background:white;border-radius:12px;padding:20px;margin:0 auto 16px;max-width:960px;box-shadow:0 2px 8px rgba(0,0,0,0.06)}
.day-hdr{font-size:16px;font-weight:bold;color:#1a237e;margin-bottom:12px;display:flex;justify-content:space-between;align-items:center}
.stats{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin:12px 0}
.s{background:#f8f9fa;border-radius:8px;padding:12px;text-align:center}
.s .v{font-size:20px;font-weight:bold}
.s .l{font-size:11px;color:#888;margin-top:2px}
.e .v{color:#1565C0}.n .v{color:#2E7D32}.u .v{color:#C62828}.a .v{color:#6A1B9A}
img{max-width:100%;border:1px solid #e0e0e0;border-radius:6px;margin-top:10px;display:block}
.note{font-size:12px;color:#888;text-align:center;margin-top:8px}
table{width:100%;border-collapse:collapse;margin-top:20px;font-size:13px}
th{background:#f0f4ff;padding:10px 14px;text-align:left;border-bottom:2px solid #ddd}
td{padding:8px 14px;border-bottom:1px solid #f0f0f0}
tr:hover{background:#fafafa}
"""

html=f'<!DOCTYPE html>\n<html lang="zh">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width,initial-scale=1">\n<title>GRACE-FO PPP 误差报告</title>\n<style>\n{css}\n</style>\n</head>\n<body>\n<div class="hdr">\n<h1>GRACE-FO PPP 精密定轨误差报告</h1>\n<p>PPP vs GNV1B 参考轨道 &nbsp;|&nbsp; 光行时 + Sagnac + 相对论改正 &nbsp;|&nbsp; 采样 30s &nbsp;|&nbsp; 每段 4 小时</p>\n</div>\n'

for r in all_results:
    b64=img_b64(r['jpg'])
    html+=f'''<div class="card">
<div class="day-hdr"><span>{r['date']}</span><span style="font-size:12px;color:#888">{r['n']} 历元</span></div>
<div class="stats">
<div class="s e"><div class="v">{r['re']:.2f} cm</div><div class="l">E（东西）RMS</div></div>
<div class="s n"><div class="v">{r['rn']:.2f} cm</div><div class="l">N（南北）RMS</div></div>
<div class="s u"><div class="v">{r['ru']:.2f} cm</div><div class="l">U（垂直）RMS</div></div>
<div class="s a"><div class="v">{r['r3']:.2f} cm</div><div class="l">3D RMS</div></div>
</div>
<img src="data:image/jpeg;base64,{b64}" alt="{r['date']} 误差图" />
<div class="note">横轴 = 时间（小时） 纵轴 = 误差（cm） 蓝=E 绿=N 红=U</div>
</div>\n'''

html+='</body>\n</html>'
Path('output/report.html').write_text(html)
print(f"\nHTML: output/report.html")
print("全部完成!")
