# V2.2.2 软件架构文档

**版本**: V2.2.2 — Sequential EKF + WL-AR PPP
**日期**: 2026-06-10
**目标**: 低轨卫星 (GRACE-FO) 精密定轨，精度优于 1m (0.5h 弧段)

---

## 1. 软件架构

```
┌──────────────────────────────────────────────────────────┐
│                  run_sequential_pod.py                    │
│                  (主控入口 / 数据管线)                      │
├──────────────────────────────────────────────────────────┤
│  数据加载          │  精密产品          │  动力学模型        │
│  gps1b_rnx_loader │  sp3_loader        │  gravity_model    │
│  gps1b_loader     │  precision_products│  orbit_dynamics   │
│  (GNV1B 参考轨道)  │  (DCB/ANTEX/C04)   │  third_body       │
│                   │                    │  srp / empirical  │
├──────────────────────────────────────────────────────────┤
│                     orbit_integrator                      │
│                (RK4/DP8 + STM + 灵敏度矩阵)                │
├──────────────────────────────────────────────────────────┤
│  coordinates  │  troposphere  │  measurement_corrections  │
│  (ECEF↔ECI)   │  (ZHD+GMF)    │  (PCO/Wind-up/Relativity) │
├──────────────────────────────────────────────────────────┤
│                  sequential_filter.py                     │
│  ┌─────────────┐  ┌──────────────┐  ┌──────────────────┐ │
│  │  EKFState   │  │ SequentialEKF│  │    MWBuffer      │ │
│  │  状态容器    │  │  predict()   │  │  WL 模糊度固定    │ │
│  │  remove_sv()│  │ process_epoch│  │  自校准 b_r_wl    │ │
│  └─────────────┘  └──────────────┘  └──────────────────┘ │
└──────────────────────────────────────────────────────────┘
```

### 1.1 模块依赖关系

```
run_sequential_pod.py
 ├── gps1b_rnx_loader.py    → GPS1B RINEX 观测数据加载
 ├── sp3_loader.py           → CODE SP3 精密星历 + CLK 钟差
 ├── precision_products.py   → DCB, ANTEX PCO, IERS C04 ERP
 ├── gravity_model.py        → ICGEM GFC 重力场模型读取 + 球谐综合
 ├── orbit_dynamics.py
 │    ├── third_body.py      → 日/月/行星第三体摄动 (JPL DE)
 │    ├── srp.py             → 太阳光压 (Cannonball)
 │    ├── empirical.py       → 经验加速度 (RTN)
 │    ├── solid_tides.py     → 固体潮修正
 │    └── relativity_orbit.py→ 相对论轨道修正
 ├── orbit_integrator.py     → 数值积分 (RK4/DP8) + STM
 ├── coordinates.py          → ECEF ↔ ECI 坐标变换
 ├── troposphere.py          → Saastamoinen ZHD + GMF 映射
 ├── measurement_corrections.py → PCO / Wind-up / Relativity
 ├── cycle_slip.py           → TurboEdit MW+GF 周跳探测
 └── sequential_filter.py
      ├── EKFState           → 状态向量容器 + SV 增删
      ├── SequentialEKF      → EKF 预测/更新主逻辑
      └── MWBuffer           → WL 模糊度累积 + 自校准固定
```

---

## 2. 状态向量与滤波参数

### 2.1 状态向量 (ECI 坐标系)

```
x = [r_x, r_y, r_z,  v_x, v_y, v_z,  aR, aT, aN,  zwd,  clk,  amb_1..amb_N]
     ├── 位置(3) ──┤ ├── 速度(3) ──┤ ├─ 经验RTN(3) ─┤  ZWD  钟差   IF模糊度(N个)
     I_R = 0:3      I_V = 3:6       I_EMP = 6:9      I_ZWD=9 I_CLK=10 I_AMB=11+
```

**总维数**: N_BASE(11) + N_sv

### 2.2 动力学模型

| 项目 | 模型 | 参数 |
|------|------|------|
| 重力场 | GGM05C 静态 | Nmax=90 |
| 第三体 | 日/月 质点引力 | JPL DE 近似 |
| 太阳光压 | Cannonball | CR=1.3, A/m=0.00586 |
| 大气阻力 | Cannonball | Cd=2.2, A/m=0.00117 |
| 经验加速度 | RTN Gauss-Markov | τ=600s, σ_ss=1e-8 m/s² |

### 2.3 过程噪声

| 状态分量 | 模型 | 噪声参数 |
|---------|------|---------|
| 位置/速度 | 未建模加速度映射 | σ_acc = 1e-3 m/s² |
| 经验 RTN | Gauss-Markov | q = σ_ss² × (1 - exp(-2dt/τ)) |
| ZWD | 随机游走 | q = 1e-9 × dt |
| 钟差 | 随机游走 (TCXO) | q = 0.001 × dt |
| 模糊度 | 常数 | q = 0 |

### 2.4 测量噪声

| 观测量 | σ | R |
|-----------|---|-----|
| 载波相位 (IF) | 0.01 m | 1e-4 m² |
| 伪距 (IF) | 0.30 m | 0.09 m² |

### 2.5 初始协方差

| 状态 | P0 | σ |
|------|-----|-----|
| 位置 | 100 m² | 10 m |
| 速度 | 1.0 (m/s)² | 1 m/s |
| 经验 RTN | 1e-12 (m/s²)² | 1e-6 m/s² |
| ZWD | 0.25 m² | 0.5 m |
| 钟差 | 1e10 m² | 1e5 m |
| 模糊度 | 100 m² | 10 m |

---

## 3. 相位模糊度固定 (WL-AR)

### 3.1 MW 组合

```
MW = (L1_cyc - L2_cyc) - (f1·P1 + f2·P2) / ((f1+f2) · λ_wl)   [cycles]

其中:
  λ_wl = c / (f1 - f2) ≈ 0.862 m  (宽巷波长)
  f1 = 1575.42 MHz,  f2 = 1227.60 MHz
```

MW 组合消除了几何距离、钟差、对流层和一阶电离层，仅保留:

```
MW = N_wl + b_r_wl - b_s_wl + noise
```

### 3.2 自校准接收机偏差

```
b_r_wl = median{ fractional_part(mean_MW_sv) }  对所有 SV (≥3颗, ≥min_epochs历元)
```

### 3.3 WL 固定判定

```
N_w_float = mean_MW_sv - b_r_wl
N_w_fixed = round(N_w_float)
判定: |N_w_float - N_w_fixed| ≤ 0.35 cycles  AND  std(MW) ≤ 0.35 cycles
```

### 3.4 IF 模糊度约束

```
B_if = λ_nl · N1 + coeff_w · N_w

其中:
  λ_nl = c / (f1 + f2) ≈ 0.107 m  (窄巷波长)
  coeff_w = c·f2 / (f1² - f2²) ≈ 0.3776 m/cycle
```

WL 固定后:
- 模糊度状态修正: `B_if ← B_if + coeff_w · (N_w_fixed - N_w_float)`
- 协方差收紧: `P_amb ← min(P_amb, 0.10 m²)`

### 3.5 SV 裁剪

超过 `prune_timeout` (默认 1800s) 未被观测的 SV 从状态中移除，避免死星模糊度污染协方差矩阵。

---

## 4. 输入输出

### 4.1 命令行参数

```
run_sequential_pod.py
  ──date YYYY-MM-DD           # 处理日期 (必需)
  ──hours FLOAT               # 弧段长度 [h] (default: 0.5)
  ──interval FLOAT            # 采样间隔 [s] (default: 30)
  ──grace-id C|D              # GRACE-FO 卫星 ID (default: C)
  ──data-dir PATH             # 数据根目录 (default: ./data)
  ──gravity-nmax INT          # 重力场最大阶数 (default: 90)
  ──sigma-acc FLOAT           # 未建模加速度 σ [m/s²] (default: 1e-3)
  ──tau-emp FLOAT             # 经验加速度时间常数 [s] (default: 600)
  ──sigma-phase FLOAT         # 相位 σ [m] (default: 0.01)
  ──sigma-code FLOAT          # 伪距 σ [m] (default: 0.30)
  ──sp3-file PATH             # CODE SP3 精密星历
  ──clk-file PATH             # CODE CLK 精密钟差
  ──antex-file PATH           # IGS ANTEX 天线文件
  ──dcb-file PATH             # CODE P1P2 DCB 文件
  ──iers-c04 PATH             # IERS C04 ERP 文件
  ──ar-min-epochs INT         # WL 固定最小历元数 (default: 10)
  ──enable-phase-windup       # 启用相位缠绕修正
  ──enable-relativity         # 启用相对论 Shapiro 修正
  ──enable-cycle-slip         # 启用 TurboEdit 周跳探测
  ──enable-all-corrections    # 启用所有 Phase 2.3 修正
```

### 4.2 输入数据

| 数据 | 文件 | 来源 | 说明 |
|------|------|------|------|
| GNV1B | `GNV1B_YYYY-MM-DD_C_04.txt` | ISDC/GFZ | 参考轨道 (1s) |
| GPS1B RINEX | `GPS1B_YYYY-MM-DD_C_04.rnx` | ISDC/GFZ | GPS 观测 (10s) |
| SP3 | `COD0OPSFIN_20241200000_01D_05M_ORB.SP3` | CODE | 精密星历 (5min) |
| CLK | `COD0OPSFIN_20241200000_01D_30S_CLK.CLK` | CODE | 精密钟差 (30s) |
| GGM05C | `GGM05C.gfc` | ICGEM | 重力场模型 |
| ANTEX | `igs14.atx` | IGS | 卫星天线 PCO/PCV |
| DCB P1C1 | `P1C1YYYY.DCB` | CODE | 码偏差 P1-C1 |
| DCB P1P2 | `P1P2YYYY.DCB` | CODE | 码偏差 P1-P2 |
| IERS C04 | `eopc04_IAU2000.txt` | IERS | 地球定向参数 |

### 4.3 输出

| 输出 | 文件 | 内容 |
|------|------|------|
| 轨道结果 | `results/sequential_ekf/seq_YYYY-MM-DD_C_X.Xh.pkl` | 历元时间、ECEF 位置/速度、经验加速度、ZWD、钟差、滤波统计 |
| 控制台输出 | stdout | 每历元状态报告 + MW 固定日志 + 最终精度统计 |

### 4.4 输出 pickle 结构

```python
{
    'epochs': [gps_sod, ...],         # GPS 秒 (N 历元)
    'r_ecef': [(N,3) ndarray],        # ECEF 位置 [m]
    'v_ecef': [(N,3) ndarray],        # ECEF 速度 [m/s]
    'a_rtn':  [(N,3) ndarray],        # RTN 经验加速度 [m/s²]
    'zwd':    [float, ...],           # 天顶湿延迟 [m]
    'clk':    [float, ...],           # 接收机钟差 [m]
    'stats':  [dict, ...],            # 每历元统计
    'n_sv':   [int, ...],             # 每历元 SV 数
    'r_gnv':  [(N,3) ndarray],        # GNV1B 参考位置 [m]
}
```

### 4.5 精度评估

程序输出 3D RMS 误差 (vs GNV1B 参考轨道):

```
-- Results --
  3D RMS vs GNV1B: 0.803 m  (mean=0.651, max=1.664)
  Mean aR=2.114e-11 aT=-3.883e-11 aN=-8.202e-13 m/s^2
  Mean ZWD=-0.0462 m
  Phase: 216 accepted, RMS=0.2270m
  Code:  216 accepted, RMS=0.5853m
  Rejected: 0
  Final SV count: 12
```

---

## 5. 精度验证结果

**测试条件**: 2024-04-29, GRACE-FO C, GGM05C Nmax=90, sigma_acc=1e-3
**测量修正**: PCO + Wind-up + Relativity + DCB + IERS C04 (无周跳探测)
**对比基准**: GNV1B 官方定轨产品 (GPS-based reduced-dynamic)

### 5.1 浮点 PPP vs WL-AR PPP

| 弧长 | 历元数 | 浮点 PPP (Phase 2.3) | WL-AR (Phase 3.0) | 精度提升 |
|------|--------|---------------------|-------------------|---------|
| 0.17h (10min) | 21 | 0.936 m | **0.803 m** | 14.2% |
| 0.5h (30min) | 61 | 1.816 m | **1.005 m** | 44.7% |
| 2h | 241 | 6.63 m | **6.06 m** | 8.6% |

### 5.2 WL 固定统计

| 弧长 | 总 SV 数 | WL 固定 SV | 固定率 | b_r_wl [cyc] | MW std (中位数) |
|------|---------|-----------|--------|-------------|----------------|
| 0.17h | 12 | 11 | 91.7% | +0.225 | 0.036 cyc |
| 0.5h | 16 | 15 | 93.8% | +0.219 | 0.054 cyc |

### 5.3 典型滤波统计

| 指标 | 0.17h | 0.5h |
|------|-------|------|
| 相位残差 RMS | 0.227 m | 0.454 m |
| 伪距残差 RMS | 0.585 m | 0.754 m |
| 更新拒绝率 | 0% | 0.2% |
| 经验 aT 均值 | -3.88e-11 m/s² | -1.41e-11 m/s² |
| ZWD 均值 | -0.046 m | +0.095 m |

### 5.4 误差退化分析 (2h+ 弧段)

```
Epoch     Time    |dr|      相位通过率    n_sv
────────────────────────────────────────────────
  96      0.8h    3.4m      9/11 (82%)     18
 144      1.2h    2.2m      9/9  (100%)    18
 168      1.4h    4.3m      7/7  (100%)    17
 192      1.6h    7.9m     10/11 (91%)     17
 216      1.8h   12.4m      6/9  (67%)     17
 240      2.0h   12.2m      4/10 (40%)     18
```

**退化原因**: GGM05C 静态重力场缺少固体潮、海潮、高阶重力场修正，未被建模的加速度 (~0.01 m/s²) 在 1h+ 弧段上累积导致位置漂移。

---

## 6. 已知限制

1. **弧长限制**: GGM05C-only 动力学的有效弧长约 1 小时；2h+ 需要完整动力学 (固体潮/海潮/150+阶)
2. **NL 模糊度**: WL 固定仅约束 1 个自由度，N1 仍为浮点，每个 SV 贡献 ~0.3m 位置误差
3. **周跳探测**: TurboEdit GF 阈值 (0.02m) 对 GRACE 数据过于敏感，需要调优
4. **24h 全日弧**: SV 裁剪已实现，但动力学模型不足导致 ~2h 后发散
5. **卫星 WL 偏差**: G20 等少数卫星残差 >0.35cyc，需外部 OSB 产品修正

---

## 7. 典型运行命令

```bash
# 0.5h 弧段 (推荐)
py -3.12 run_sequential_pod.py \
  --date 2024-04-29 --hours 0.5 --interval 30 --grace-id C \
  --sp3-file data/CODE/2024/COD0OPSFIN_20241200000_01D_05M_ORB.SP3 \
  --clk-file data/CODE/2024/COD0OPSFIN_20241200000_01D_30S_CLK.CLK \
  --antex-file data/igs14.atx \
  --dcb-file data/CODE/2024/P1P22404.DCB \
  --iers-c04 data/IERS/eopc04_IAU2000.txt \
  --ar-min-epochs 6 --enable-phase-windup --enable-relativity

# 0.17h 短弧 (快速测试)
py -3.12 run_sequential_pod.py \
  --date 2024-04-29 --hours 0.17 --interval 30 --grace-id C \
  --sp3-file data/CODE/2024/COD0OPSFIN_20241200000_01D_05M_ORB.SP3 \
  --clk-file data/CODE/2024/COD0OPSFIN_20241200000_01D_30S_CLK.CLK \
  --antex-file data/igs14.atx --dcb-file data/CODE/2024/P1P22404.DCB \
  --iers-c04 data/IERS/eopc04_IAU2000.txt \
  --ar-min-epochs 5 --enable-phase-windup --enable-relativity
```
