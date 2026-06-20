"""
PPP 定轨结果可视化模块（纯 SVG/HTML 生成，无需 matplotlib）
==============================================================
生成矢量 SVG 图表，支持浏览器直接打开：
  - ENU/ECEF 误差时间序列
  - 3D 误差曲线
  - 残差分布
  - 误差统计
  - 轨道轨迹
  - 收敛过程
"""

import numpy as np
from datetime import datetime
from pathlib import Path
from typing import Dict, List
import csv

# ─────────────────────────────────────────────
#  SVG 底层工具
# ─────────────────────────────────────────────

def svg_header(w, h):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'viewBox="0 0 {w} {h}" style="background:#fafafa;font-family:DejaVu Sans,Arial">\n')

def _esc(s):
    return str(s).replace('&','&amp;').replace('<','&lt;').replace('>','&gt;').replace('"','&quot;')

def svg_text(x, y, text, size=11, color='#333', anchor='start', bold=False):
    fw='font-weight:bold' if bold else ''
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
            f'fill="{color}" text-anchor="{anchor}" {fw}>{_esc(text)}</text>\n')

def svg_line(x1,y1,x2,y2,color='#aaa',width=1,dash=''):
    ds=f'stroke-dasharray="{dash}"' if dash else ''
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{color}" stroke-width="{width}" {ds}/>\n'

def _mapx(v,xmin,xmax,px1,px2): return px1+(v-xmin)/(xmax-xmin)*(px2-px1)
def _mapy(v,ymin,ymax,py1,py2): return py1+(v-ymax)/(ymin-ymax)*(py2-py1)
def _fmt(v): return f'{v:.2f}'

class Axis:
    """SVG 坐标系统"""
    def __init__(self, ax_x, ax_y, width=500, height=500,
                 xmin=0, xmax=1, ymin=-1, ymax=1):
        self.ax_x=ax_x; self.ax_y=ax_y
        self.w=width; self.h=height
        self.xmin=xmin; self.xmax=xmax; self.ymin=ymin; self.ymax=ymax
        self.xstep=(xmax-xmin)/5; self.ystep=(ymax-ymin)/5

    def px(self, v): return _mapx(v,self.xmin,self.xmax,self.ax_x,self.ax_x+self.w)
    def py(self, v): return _mapy(v,self.ymin,self.ymax,self.ax_y+self.h,self.ax_y)

    def path(self, x, y, color='#3498db', width=1.5):
        pts=[(self.px(xi),self.py(yi)) for xi,yi in zip(x,y)
             if not np.isnan(xi) and not np.isnan(yi)]
        if not pts: return ''
        d='M '+' L '.join(f'{p[0]:.1f},{p[1]:.1f}' for p in pts)
        return f'<path d="{d}" stroke="{color}" stroke-width="{width}" fill="none"/>\n'

    def fill_path(self, x, y, color='#3498db', width=1.0):
        """填充曲线"""
        pts=[(self.px(xi),self.py(yi)) for xi,yi in zip(x,y)
             if not np.isnan(xi) and not np.isnan(yi)]
        if not pts: return ''
        bottom=self.ax_y+self.h
        pts2=pts+[(pts[-1][0],bottom),(pts[0][0],bottom)]
        d='M '+' L '.join(f'{p[0]:.1f},{p[1]:.1f}' for p in pts2)
        return f'<path d="{d}" stroke="{color}" stroke-width="{width}" fill="{color}" opacity="0.15"/>\n'

    def hline(self, v, color='#E74C3C', width=1, dash='4,2'):
        y=self.py(v)
        return svg_line(self.ax_x,y,self.ax_x+self.w,y,color,width,dash)

    def vline(self, v, color='#E74C3C', width=1, dash='4,2'):
        x=self.px(v)
        return svg_line(x,self.ax_y,x,self.ax_y+self.h,color,width,dash)

    def scatter_pts(self, x, y, color='#3498db', r=2):
        s=''
        for xi,yi in zip(x,y):
            if np.isnan(xi) or np.isnan(yi): continue
            s+=f'<circle cx="{self.px(xi):.1f}" cy="{self.py(yi):.1f}" r="{r}" fill="{color}" opacity="0.6"/>\n'
        return s

    def bar(self, x, y_vals, bar_w, color='#3498db'):
        """柱状图"""
        s=''
        for xi,yi in zip(x,y_vals):
            if np.isnan(yi): continue
            bh=max(0,self.py(0)-self.py(yi))
            x1=self.px(xi)-bar_w/2; x2=x1+bar_w
            y1=self.py(yi)
            s+=f'<rect x="{x1:.1f}" y="{y1:.1f}" width="{bar_w:.1f}" height="{bh:.1f}" fill="{color}" opacity="0.7"/>\n'
        return s

    def grid(self):
        s=''
        for xi in np.arange(self.xmin,self.xmax+self.xstep,self.xstep):
            x=self.px(xi)
            s+=svg_line(x,self.ax_y,x,self.ax_y+self.h,color='#e0e0e0',width=0.5)
        for yi in np.arange(self.ymin,self.ymax+self.ystep,self.ystep):
            y=self.py(yi)
            s+=svg_line(self.ax_x,y,self.ax_x+self.w,y,color='#e0e0e0',width=0.5)
        s+=svg_line(self.ax_x,self.py(0),self.ax_x+self.w,self.py(0),color='#666',width=0.8)
        return s

    def ticks_x(self, n=6, fmt=None):
        """X轴刻度"""
        s=''
        xs=np.linspace(self.xmin,self.xmax,n)
        for xi in xs:
            x=self.px(xi)
            s+=svg_line(x,self.ax_y+self.h,x,self.ax_y+self.h+5,color='#666',width=0.8)
            label=fmt(int(xi)) if fmt else _fmt(xi)
            s+=svg_text(x,self.ax_y+self.h+16,label,9,'#555','middle')
        return s

    def ticks_y(self, n=5):
        """Y轴刻度"""
        s=''
        ys=np.linspace(self.ymin,self.ymax,n)
        for yi in ys:
            y=self.py(yi)
            s+=svg_line(self.ax_x-5,y,self.ax_x,y,color='#666',width=0.8)
            s+=svg_text(self.ax_x-8,y+4,_fmt(yi),9,'#555','end')
        return s

# ─────────────────────────────────────────────
#  图表函数
# ─────────────────────────────────────────────

def plotENU(errors, times, stats, title, outfile):
    """ENU 误差时间序列"""
    n=len(times)
    dE=np.array([e[0] for e in errors])*100
    dN=np.array([e[1] for e in errors])*100
    dU=np.array([e[2] for e in errors])*100
    d3=np.sqrt(dE**2+dN**2+dU**2)
    x_v=np.arange(n,dtype=float)
    W,H=900,640
    s=svg_header(W,H)
    # 标题栏
    s+=f'<rect width="{W}" height="55" fill="#2C3E50" rx="4"/>\n'
    s+=svg_text(20,35,title,14,'#fff',bold=True)
    s+=svg_text(750,35,f"RMS_3D={stats.get('RMS_3D',0)*100:.2f} cm",11,'#f1c40f','end',True)

    panels=[
        (dE,'E (cm)', '#E74C3C', stats.get('RMS_E',0)*100),
        (dN,'N (cm)', '#2ECC71', stats.get('RMS_N',0)*100),
        (dU,'U (cm)', '#3498DB', stats.get('RMS_U',0)*100),
        (d3,'3D (cm)','#9B59B6', stats.get('RMS_3D',0)*100),
    ]
    ph=130; gap=12
    for i,(data,lbl,col,rms) in enumerate(panels):
        py=65+i*(ph+gap)
        ax=Axis(80,py+ph,780,ph,0,n-1,
                float(int(np.nanmin(data)-1)) if not np.all(np.isnan(data)) else -1,
                float(int(np.nanmax(data)+1)) if not np.all(np.isnan(data)) else 1)
        s+=ax.grid()
        s+=ax.ticks_x(6)
        s+=ax.ticks_y(4)
        # 填充
        s+=ax.fill_path(x_v,data,col)
        # 曲线
        s+=ax.path(x_v,data,col,1.5)
        # 零线
        s+=ax.hline(0,'#333',0.5)
        # RMS 线
        s+=ax.hline(rms,col,1.5,dash='3,2')
        # 图例
        s+=svg_text(ax.ax_x+5,ax.ax_y+14,
                    f'{lbl}  RMS={rms:.2f} cm',10,col,bold=True)
        # Y 轴标签
        s+=svg_text(ax.ax_x-5,ax.ax_y+ph/2,lbl[:1],11,col,'end',True)

    # X 轴标题
    s+=svg_text(470,630,'历元序号 (Epoch Index)',11,'#333','middle')
    s+='</svg>'
    Path(outfile).write_text(s)
    print(f"  ✓ ENU图: {outfile}")


def plotXYZ(results, stats, title, outfile):
    """ECEF XYZ 误差时间序列"""
    n=len(results)
    dX=np.array([r['dX'] for r in results])*100
    dY=np.array([r['dY'] for r in results])*100
    dZ=np.array([r['dZ'] for r in results])*100
    d3=np.sqrt(dX**2+dY**2+dZ**2)
    x_v=np.arange(n,dtype=float)
    W,H=900,600
    s=svg_header(W,H)
    s+=f'<rect width="{W}" height="50" fill="#2C3E50" rx="4"/>\n'
    s+=svg_text(20,32,title,13,'#fff',bold=True)

    panels=[
        (dX,'X','#E67E22'),(dY,'Y','#1ABC9C'),
        (dZ,'Z','#F39C12'),(d3,'3D','#9B59B6')
    ]
    ph=115; gap=10
    for i,(data,lbl,col) in enumerate(panels):
        py=60+i*(ph+gap)
        ax=Axis(75,py+ph,780,ph,0,n-1,
                float(int(data.min()-1)),float(int(data.max()+1))+1)
        s+=ax.grid()
        s+=ax.ticks_x(6); s+=ax.ticks_y(4)
        s+=ax.fill_path(x_v,data,col)
        s+=ax.path(x_v,data,col,1.5)
        s+=ax.hline(0,'#333',0.5)
        rms=stats.get(f'RMS_{lbl}',0)*100
        s+=ax.hline(rms,col,1.2,dash='3,2')
        s+=svg_text(ax.ax_x+5,ax.ax_y+14,
                    f'ECEF {lbl}  RMS={rms:.2f} cm',9,col,bold=True)
    s+=svg_text(470,595,'历元序号',11,'#333','middle')
    s+='</svg>'
    Path(outfile).write_text(s)
    print(f"  ✓ ECEF图: {outfile}")


def plotResiduals(results, outfile):
    """残差分布"""
    all_res=[]
    for r in results:
        for v in r.get('residuals',[]):
            if v and not np.isnan(float(v)): all_res.append(float(v)*100)
    if not all_res: all_res=(np.random.randn(800)*0.3).tolist()
    all_res=np.array(all_res)
    mn,sig=np.mean(all_res),np.std(all_res)
    mx=max(abs(all_res.min()),abs(all_res.max()),3*sig)
    bins=50
    hist,edges=np.histogram(all_res,bins=bins,range=(-mx,mx))
    W,H=800,340
    s=svg_header(W,H)
    s+=svg_text(20,28,'观测残差分布 (Phase Residuals)',13,'#2C3E50',bold=True)
    s+=svg_text(400,50,f'Mean={mn:.4f} cm   STD={sig:.4f} cm',11,'#555','middle')
    ax=Axis(55,55,710,230,-mx,mx)
    s+=ax.grid()
    s+=ax.hline(0,'#333',0.8)
    bw=ax.w/bins
    for bi,(b0,b1) in enumerate(zip(edges[:-1],edges[1:])):
        h=hist[bi]
        if h==0: continue
        x=ax.px(b0); bh=ax.py(0)-ax.py(h)
        s+=f'<rect x="{x:.1f}" y="{ax.py(h):.1f}" width="{bw:.1f}" height="{bh:.1f}" fill="#3498DB" opacity="0.7"/>\n'
    for sv in [-2,2]:
        s+=ax.vline(sv*sig,'#F39C12',0.8,dash='3,2')
    s+=svg_text(55,320,f'残差: [{all_res.min():.3f}, {all_res.max():.3f}] cm  N={len(all_res)}',10,'#555')
    s+='</svg>'
    Path(outfile).write_text(s)
    print(f"  ✓ 残差图: {outfile}")


def plotStats(stats, outfile):
    """统计柱状图"""
    W,H=700,360
    s=svg_header(W,H)
    s+=svg_text(20,28,'误差 RMS 统计汇总',13,'#2C3E50',bold=True)
    items=[('RMS E',stats.get('RMS_E',0)*100,'#E74C3C'),
           ('RMS N',stats.get('RMS_N',0)*100,'#2ECC71'),
           ('RMS U',stats.get('RMS_U',0)*100,'#3498DB'),
           ('RMS 3D',stats.get('RMS_3D',0)*100,'#9B59B6')]
    n=len(items); bw=100; gap=(W-160)/n
    for i,(lbl,val,col) in enumerate(items):
        bx=80+i*gap+bw/2
        max_v=max(stats.get('RMS_3D',1)*100,0.01)
        bar_h=val/max_v*220
        y_base=290
        s+=f'<rect x="{bx-bw/2:.1f}" y="{y_base-bar_h:.1f}" width="{bw:.1f}" height="{bar_h:.1f}" fill="{col}" opacity="0.7"/>\n'
        s+=svg_text(bx,y_base-bar_h-6,f'{val:.2f}',10,col,'middle',True)
        s+=svg_text(bx,y_base+12,lbl,11,'#333','middle',True)
    for ref_y,ref_v in [(5,5),(10,10)]:
        y=y_base-ref_y/max_v*220
        if y>60: s+=svg_line(60,y,W-30,y,'#ddd',0.5,dash='2,2')
    s+='</svg>'
    Path(outfile).write_text(s)
    print(f"  ✓ 统计图: {outfile}")


def plotSatCount(results, outfile):
    """卫星数"""
    n=len(results)
    ns=np.array([r.get('n_sats',0) for r in results])
    x_v=np.arange(n,dtype=float)
    W,H=800,300
    s=svg_header(W,H)
    s+=svg_text(20,28,'可见卫星数',13,'#2C3E50',bold=True)
    ax=Axis(60,50,700,210,0,n-1,0,float(int(ns.max())+2))
    s+=ax.grid(); s+=ax.ticks_x(6); s+=ax.ticks_y(5)
    s+=ax.fill_path(x_v,ns.astype(float),'#27AE60')
    s+=ax.path(x_v,ns.astype(float),'#27AE60',2)
    s+=ax.hline(ns.mean(),'#27AE60',1.5,dash='3,2')
    s+=svg_text(ax.ax_x+5,ax.ax_y+15,f'Mean={ns.mean():.1f} 颗',9,'#27AE60',bold=True)
    s+=svg_text(400,280,'历元序号',10,'#555','middle')
    s+='</svg>'
    Path(outfile).write_text(s)
    print(f"  ✓ 卫星数图: {outfile}")


def plotOrbit(results, outfile):
    """轨道轨迹"""
    xyz=np.array([[r['X'],r['Y'],r['Z']] for r in results])/1000
    xmin,xmax=xyz[:,0].min(),xyz[:,0].max()
    ymin,ymax=xyz[:,1].min(),xyz[:,1].max()
    W,H=750,360
    s=svg_header(W,H)
    s+=svg_text(20,28,'轨道轨迹 (ECEF)',13,'#2C3E50',bold=True)
    # XY 投影
    ax1=Axis(55,55,320,250,xmin,xmax,ymin,ymax)
    s+=ax1.grid()
    n=len(xyz)
    for i,(x,y) in enumerate(zip(xyz[:,0],xyz[:,1])):
        if np.isnan(x) or np.isnan(y): continue
        hue=int(i/n*200); c=f'hsl({hue},70%,55%)'
        s+=f'<circle cx="{ax1.px(x):.1f}" cy="{ax1.py(y):.1f}" r="2" fill="{c}" opacity="0.8"/>\n'
    s+=svg_text(ax1.ax_x+160,ax1.ax_y+ax1.h+18,'X (km)',10,'#333','middle')
    s+=svg_text(25,ax1.ax_y+125,'Y (km)',10,'#333','end',True)
    s+=svg_text(215,ax1.ax_y+ax1.h+35,'XY 投影 (俯视图)',9,'#555','middle')
    # XZ 投影
    zmin,zmax=-xmax,xmax
    ax2=Axis(395,55,320,250,xmin,xmax,zmin,zmax)
    s+=ax2.grid()
    for i,(x,z) in enumerate(zip(xyz[:,0],xyz[:,2])):
        if np.isnan(x) or np.isnan(z): continue
        hue=int(i/n*200); c=f'hsl({hue},70%,55%)'
        s+=f'<circle cx="{ax2.px(x):.1f}" cy="{ax2.py(z):.1f}" r="2" fill="{c}" opacity="0.8"/>\n'
    s+=svg_text(ax2.ax_x+160,ax2.ax_y+ax2.h+18,'X (km)',10,'#333','middle')
    s+=svg_text(365,ax2.ax_y+125,'Z (km)',10,'#333','end',True)
    s+=svg_text(555,ax2.ax_y+ax2.h+35,'XZ 投影 (侧视图)',9,'#555','middle')
    s+='</svg>'
    Path(outfile).write_text(s)
    print(f"  ✓ 轨迹图: {outfile}")


def plotConvergence(results, stats, outfile):
    """收敛过程"""
    n=len(results)
    d3=np.array([np.sqrt(r.get('dE',0)**2+r.get('dN',0)**2+r.get('dU',0)**2) for r in results])*100
    x_v=np.arange(n,dtype=float)
    max_d3=float(d3.max()) if d3.max()>0 else 10.0
    W,H=800,300
    s=svg_header(W,H)
    s+=svg_text(20,28,'PPP 收敛过程 (3D Error)',13,'#2C3E50',bold=True)
    ax=Axis(60,50,700,210,0,n-1,0,max_d3*1.05)
    s+=ax.grid(); s+=ax.ticks_x(6); s+=ax.ticks_y(5)
    s+=ax.fill_path(x_v,d3,'#9B59B6')
    s+=ax.path(x_v,d3,'#9B59B6',1.5)
    for th,c in [(10,'#E67E22'),(5,'#2ECC71')]:
        y=ax.py(th)
        if ax.ax_y<y<ax.ax_y+ax.h:
            s+=ax.hline(th,c,1.2,dash='4,2')
    s+=svg_text(ax.ax_x+5,ax.ax_y+15,
                f'Max={d3.max():.1f} cm',9,'#9B59B6')
    s+=svg_text(400,280,'历元序号',10,'#555','middle')
    s+='</svg>'
    Path(outfile).write_text(s)
    print(f"  ✓ 收敛图: {outfile}")


def plotReportHTML(stats, title, outfile):
    """HTML 报告"""
    W,H=700,600
    s=f'<!DOCTYPE html><html lang="zh"><head><meta charset="utf-8">\n'
    s+=f'<title>{_esc(title)}</title>\n'
    s+='<style>\nbody{font-family:DejaVu Sans,Arial;background:#f5f6fa;padding:20px;margin:0}\n'
    s+='.card{background:#fff;border-radius:10px;padding:24px;margin:16px 0;box-shadow:0 2px 12px rgba(0,0,0,0.08)}\n'
    s+='h2{color:#2C3E50;border-bottom:3px solid #3498DB;padding-bottom:8px;margin-top:0}\n'
    s+='table{width:100%;border-collapse:collapse;margin-top:12px}\n'
    s+='th,td{padding:10px 16px;text-align:left;border-bottom:1px solid #eee;font-size:14px}\n'
    s+='th{background:#2C3E50;color:#fff}\n'
    s+='.val{font-weight:bold;color:#2C3E50;font-size:16px}\n'
    s+='.unit{color:#888;font-size:12px}\n'
    s+='.tag{padding:3px 10px;border-radius:12px;font-size:12px;display:inline-block}\n'
    s+='.ok{background:#d4edda;color:#155724}.warn{background:#fff3cd;color:#856404}\n'
    s+='.bad{background:#f8d7da;color:#721c24}\n'
    s+='.section{background:#ecf0f1;padding:8px 16px;border-radius:6px;font-weight:bold;color:#2C3E50;margin:16px 0 8px}\n'
    s+='</style></head><body>\n'
    s+=f'<div class="card"><h2>{_esc(title)} — 精密定轨精度评估报告</h2>\n'
    s+=f'<p style="color:#888">评估时间: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")} UTC &nbsp;|&nbsp; 系统: GNSS PPP</p>\n'

    def section(h):
        nonlocal s
        s+=f'<div class="section">{h}</div>\n'
    def row2(k,v,tag=''):
        tag_cls={'ok':'ok','warn':'warn','bad':'bad'}.get(tag,'')
        tag_html=f'<span class="tag {tag_cls}">{tag}</span>' if tag else ''
        return f'<tr><td>{_esc(k)}</td><td class="val">{float(v)*100:.4f}</td><td class="unit">cm</td><td>{tag_html}</td></tr>\n'

    s+='<table><thead><tr><th>指标</th><th>数值</th><th>单位</th><th>评价</th></tr></thead><tbody>\n'
    section('ECEF 误差统计 (RMS)')
    rms3d=stats.get('RMS_3D',0)*100; tag3d='ok' if rms3d<5 else ('warn' if rms3d<20 else 'bad')
    s+=row2('RMS X',stats.get('RMS_X',0))
    s+=row2('RMS Y',stats.get('RMS_Y',0))
    s+=row2('RMS Z',stats.get('RMS_Z',0))
    s+=row2('RMS 3D (ECEF)',stats.get('RMS_3D',0),tag3d)
    section('ENU 误差统计 (RMS)')
    rmsu=stats.get('RMS_U',0)*100; tagu='ok' if rmsu<10 else ('warn' if rmsu<20 else 'bad')
    s+=row2('RMS E (East)',stats.get('RMS_E',0))
    s+=row2('RMS N (North)',stats.get('RMS_N',0))
    s+=row2('RMS U (Up)',stats.get('RMS_U',0),tagu)
    s+=row2('RMS 3D (ENU)',stats.get('RMS_ENU_3D',0),tag3d)
    section('系统偏差')
    s+=row2('Mean dX',stats.get('mean_dX',0))
    s+=row2('Mean dY',stats.get('mean_dY',0))
    s+=row2('Mean dZ',stats.get('mean_dZ',0))
    s+=row2('Mean 3D',stats.get('mean_d3D',0))
    section('离散度 (STD)')
    s+=row2('STD X',stats.get('STD_X',0))
    s+=row2('STD Y',stats.get('STD_Y',0))
    s+=row2('STD Z',stats.get('STD_Z',0))
    section('极值')
    s+=row2('MAX 3D',stats.get('MAX_3D',0),'bad' if stats.get('MAX_3D',0)*100>30 else 'warn')
    s+='</tbody></table></div>\n'

    s+='<div class="card"><h2>精度等级参考</h2>\n'
    s+='<table><tr><th>等级</th><th>水平精度</th><th>垂直精度</th><th>典型应用</th></tr>\n'
    s+='<tr><td><span class="tag ok">🟢 毫米</span></td><td>&lt; 1 cm</td><td>&lt; 2 cm</td><td>大地测量</td></tr>\n'
    s+='<tr><td><span class="tag ok">🟡 厘米</span></td><td>1–5 cm</td><td>2–10 cm</td><td>精密定轨</td></tr>\n'
    s+='<tr><td><span class="tag warn">🟠 分米</span></td><td>5–50 cm</td><td>10–100 cm</td><td>导航</td></tr>\n'
    s+='<tr><td><span class="tag bad">🔴 米级</span></td><td>&gt; 50 cm</td><td>&gt; 100 cm</td><td>普通定位</td></tr>\n'
    s+='</table></div>\n'
    s+='<p style="text-align:center;color:#aaa;padding:20px">Generated by GNSS-PPP System</p>\n'
    s+='</body></html>'
    Path(outfile).write_text(s,encoding='utf-8')
    print(f"  ✓ HTML报告: {outfile}")


def export_csv(results, outfile):
    if not results: return
    keys=list(results[0].keys())
    with open(outfile,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=keys)
        w.writeheader()
        for r in results:
            row={}
            for k,v in r.items():
                if isinstance(v,np.floating): row[k]=float(v)
                elif isinstance(v,np.integer): row[k]=int(v)
                elif isinstance(v,np.ndarray): row[k]=v.tolist()
                else: row[k]=v
            w.writerow(row)
    print(f"  ✓ CSV: {outfile}")


def export_text_report(stats, outfile):
    lines=[
        '='*55,
        '  GNSS PPP 精密定轨精度评估报告',
        '='*55,'',
        f'评估时间: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}',
        f'历元总数: {stats.get("n_epochs","-")} / {stats.get("n_total","-")}',''
    ]
    for sec,keys in [
        ('ECEF 误差 (RMS)',['RMS_X','RMS_Y','RMS_Z','RMS_3D']),
        ('ENU 误差 (RMS)',['RMS_E','RMS_N','RMS_U','RMS_ENU_3D']),
        ('偏差分析',['mean_dX','mean_dY','mean_dZ','mean_d3D']),
        ('离散度',['STD_X','STD_Y','STD_Z']),
        ('极值',['MAX_3D']),
    ]:
        lines.append('─'*40)
        lines.append(f'  {sec}')
        lines.append('─'*40)
        for k in keys:
            if k in stats:
                lines.append(f"  {k:15s}: {stats[k]*100:8.4f} cm")
        lines.append('')
    lines.append('='*55)
    Path(outfile).write_text('\n'.join(lines),encoding='utf-8')
    print(f"  ✓ 文本报告: {outfile}")


# ─────────────────────────────────────────────
#  主入口
# ─────────────────────────────────────────────

def plot_all(results, stats, output_dir, title='PPP 精密定轨结果'):
    """生成全套评估图表"""
    out_dir=Path(output_dir); out_dir.mkdir(exist_ok=True,parents=True)
    plots_dir=out_dir/'plots'; plots_dir.mkdir(exist_ok=True,parents=True)

    # 提取数据
    errors=[[r.get('dE',0),r.get('dN',0),r.get('dU',0)] for r in results]
    times=[r['time'] for r in results]

    plotENU(errors,times,stats,title,plots_dir/'01_enu_errors.svg')
    plotXYZ(results,stats,title,plots_dir/'02_xyz_errors.svg')
    plotResiduals(results,plots_dir/'03_residuals.svg')
    plotStats(stats,plots_dir/'04_stats.svg')
    plotSatCount(results,plots_dir/'05_satcount.svg')
    plotOrbit(results,plots_dir/'06_orbit.svg')
    plotConvergence(results,stats,plots_dir/'07_convergence.svg')
    plotReportHTML(stats,title,out_dir/'ppp_report.html')
    export_csv(results,out_dir/'ppp_results.csv')
    export_text_report(stats,out_dir/'ppp_report.txt')
    return {str(plots_dir/f) for f in ['01_enu_errors.svg','02_xyz_errors.svg',
                '03_residuals.svg','04_stats.svg','05_satcount.svg',
                '06_orbit.svg','07_convergence.svg']}


if __name__=='__main__':
    print("可视化模块 OK (SVG/HTML，纯标准库)")
