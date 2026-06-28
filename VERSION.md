# GRACE-FO PPP-AR POD — Version 2.2.4

**Release Date**: 2026-06-16

---

## 1. 软件标识

| 项目 | 值 |
|------|-----|
| 版本号 | **V2.2.4** |
| 发布日 | 2026-06-16 |
| 适用卫星 | GRACE-FO C / D |
| 观测系统 | GPS L1/L2 (双频 IF) |
| 滤波方法 | 序贯扩展卡尔曼滤波 (Sequential EKF) |
| 模糊度固定 | WL PPP-AR (自标定) + NL SD bootstrapping |
| 动力学 | GGM05C N=90 (简化) / Orekit v13.1.5 (可选) |
| 许可 | Internal research code |

---

## 2. 软件依赖

| 组件 | 版本 | 备注 |
|------|------|------|
| Python | 3.12.6 | |
| numpy | 2.2.6 | |
| scipy | 1.17.1 | 仅 `linalg.solve` |
| matplotlib | 3.10.9 | 仅全天线评估 |
| astropy | (any) | ECI↔ECEF 坐标变换 |
| orekit-jpype | 13.1.5 | 可选 (`--dynamics-mode orekit`) |
| Java (JBR) | OpenJDK 21.0.6 | 可选, PyCharm 自带 |

---

## 3. 外部数据产品

| 产品 | 目录 | 来源 | 延迟 |
|------|------|------|------|
| GPS1B RINEX | `data/GPS1B_*.rnx` | JPL/GRACE-FO L1B | 1-3天 |
| GNV1B 参考轨道 | `data/GNV1B_*.txt` | JPL/GRACE-FO L1B | 1-3天 |
| CODE Final SP3 | `data/CODE/2024/COD*ORB.SP3` | CODE/CDDIS | 12-18天 |
| CODE 30s CLK | `data/CODE/2024/COD*CLK.CLK` | CODE/CDDIS | 12-18天 |
| CODE P1P2 DCB | `data/CODE/2024/P1P2*.DCB` | CODE | 月 |
| CODE P1C1 DCB | `data/CODE/2024/P1C1*.DCB` | CODE | 月 |
| IGS14.atx | `data/igs14.atx` | IGS | 年 |
| IERS C04 ERP | `data/IERS/eopc04_IAU2000.txt` | IERS | 30天 |
| GGM05C.gfc | `data/gravity/GGM05C.gfc` | ICGEM | 静态 |

---

## 4. 软件架构

```
run_sequential_pod.py          ← 主入口: CLI + 数据加载 + EKF 控制循环
├── src/sequential_filter.py   ← EKF 核心 (SequentialEKF + EKFState)
│   ├── src/orbit_dynamics.py  ← 力模型: GGM05C + 第三体 + SRP + Drag
│   ├── src/gravity_model.py   ← ICGEM .gfc 重力场读取
│   ├── src/third_body.py      ← 日月解析星历
│   ├── src/srp.py            ← Cannonball 光压
│   ├── src/empirical.py      ← RTN 经验力 (ECI 变换)
│   ├── src/troposphere.py    ← Saastamoinen ZHD + GMF 映射函数
│   ├── src/coordinates.py    ← ECI↔ECEF (astropy)
│   ├── src/ambiguity.py      ← MW 组合 + WL/NL 常量
│   └── src/cycle_slip.py     ← TurboEdit MW+GF 周跳检测
├── src/orbit_integrator.py   ← RK4 积分 + STM 传播
├── src/batch_lsq.py          ← BatchAmbiguityResolver (全弧段)
├── src/batch_solver.py       ← BatchLinearSolver (clock+ZWD+amb)
├── src/batch_orbit.py        ← BatchOrbitLSQ (实验性)
├── src/orekit_bridge.py      ← Orekit NumericalPropagator 桥接
├── src/precision_products.py ← CODE 产品读取器 (CLK/ANTEX/DCB/IERS)
├── src/measurement_corrections.py ← PCO + 相位缠绕 + 相对论
├── src/sp3_loader.py         ← SP3 星历解析
├── src/gps1b_rnx_loader.py   ← RINEX 读取
├── eval_day2.py              ← 全天线评估 (6时段 × 0.5h → 18窗口)
└── fullday_assessment.py     ← 全天线评估 (实验性)
```

### 状态向量 (ECI, 维度 = 11 + N_sv)

```
索引  | 符号      | 维度 | 说明
------|----------|------|------
0:2   | r        | 3    | 卫星 ECI 位置 [m]
3:5   | v        | 3    | 卫星 ECI 速度 [m/s]
6:8   | aR,aT,aN | 3    | 经验 RTN 加速度 (Gauss-Markov, τ=600s)
9     | zwd      | 1    | 天顶湿延迟 [m] (随机游走, σ=1e-9 m/√s)
10    | clk      | 1    | 接收机钟差 [m] (随机游走, σ=0.02-0.03 m/√s)
11..N | amb_k    | N_sv | IF 相位模糊度 per SV [m] (常数)
```

### EKF 处理流程 (每历元)

```
1. predict():  轨道积分 (DOPRI8, 10s sub-steps, 30s EKF step)
               + 状态协方差传播 (STM Φ + 过程噪声 Q)
   Q 参数:     σ_acc=1e-3 m/s², τ_emp=600s, Q_amb=0

2. 每颗 SV 按序处理:
   a. 卫星天线 PCO 修正 + ECI 转换
   b. 相对论 Shapiro 修正 (L1/L2)
   c. 相位缠绕修正 (Wu et al. 1993)
   d. 码观测更新 → 锚定 clock (无 amb 项)
   e. 相位观测更新 → 锚定 ambiguity (有 amb 项)
   f. MW 积累 → WL 固定 (≥6 epochs, σ<0.35 cyc)
   g. WL 渐进收紧: P_amb 100→0.10→0.04→0.01 m²
   h. NL SD 整数 bootstrapping (Phase 5/6)

3. 接收机端 WL 偏差自标定 (b_r_wl):
   中位数 (小数部分 of Σ MW per SV)

4. 死卫星剪枝: timeout=60s
```

### V2.2.4 关键参数 (自适应)

```python
# 仰角加权
elev_exp_phase = 1.0                              # 相位: 完整仰角依赖
elev_exp_code  = 1.0 if hours < 0.3 else 0.70     # 码: 短弧紧, 长弧松

# Clock 过程噪声 (自适应)
clock_rw = 0.0004 if hours < 0.3 else 0.001       # σ_rw = √clock_rw

# 离群值检测 (自适应)
chi2_threshold = 25 if hours < 0.3 else 100        # 短弧紧, 长弧松

# WL 模糊度固定
ar_min_epochs = 6                                  # 30s 数据最优值
```

### 文件清单 (`src/`, 共 25 个 .py 文件)

| 文件 | ~行数 | 描述 |
|------|------|------|
| `sequential_filter.py` | 1140 | EKF 核心 |
| `orbit_dynamics.py` | 300 | 力模型 |
| `orbit_integrator.py` | 560 | 轨道积分 + STM |
| `gravity_model.py` | 200 | .gfc 读取 |
| `coordinates.py` | 80 | ECI↔ECEF |
| `ambiguity.py` | 80 | MW 组合 |
| `cycle_slip.py` | 200 | 周跳检测 |
| `batch_lsq.py` | 380 | BatchAmbiguityResolver |
| `batch_solver.py` | 220 | BatchLinearSolver |
| `batch_orbit.py` | 300 | BatchOrbitLSQ |
| `orekit_bridge.py` | 500 | Orekit 桥接 |
| `precision_products.py` | 400 | CODE 产品读取 |
| `measurement_corrections.py` | 200 | PCO/相位缠绕/相对论 |
| `empirical.py` | 100 | RTN 变换 |
| `troposphere.py` | 100 | 对流层 |
| `third_body.py` | 80 | 第三体 |
| `srp.py` | 80 | 光压 |
| `solid_tides.py` | 120 | 固体潮 |
| `relativity_orbit.py` | 60 | 相对论 |
| `sp3_loader.py` | 200 | SP3 |
| `gps1b_rnx_loader.py` | 200 | RINEX |
| `gps1b_loader.py` | 60 | 旧版加载器 |
| `ppp.py` | — | 单点定位 |
| `fetch_data.py` | — | 数据下载工具 |
| `plotting.py` | — | 绘图工具 |

---

## 5. 已验证精度

**测试条件**: 2024-04-29, GRACE-FO C, GPS only, GGM05C Nmax=90, 30s EKF step  
**参照**: JPL GNV1B 参考轨道 (~2cm 精度)

| 弧段 | Float PPP | V2.2.4 | **Arc AR (Phase 12)** | 累计提升 |
|------|-----------|--------|------|----------|
| 0.17h (10min, 20 epochs) | 0.936 m | 0.293 m | **0.246 m** | **-74%** |
| 0.5h (30min, 60 epochs) | 1.817 m | 0.986 m | **0.986 m** | **-46%** |

### Arc AR 运行命令 (Phase 12.0, 0.246m at 0.17h)

```powershell
py -3.12 run_sequential_pod.py `
  --date 2024-04-29 --hours 0.17 --interval 30 --grace-id C `
  ... (same data flags) --arc-ar
```
  --date 2024-04-29 --hours 0.17 --interval 30 --grace-id C `
  --dynamics-mode simplified `
  --sp3-file data/CODE/2024/COD0OPSFIN_20241200000_01D_05M_ORB.SP3 `
  --clk-file data/CODE/2024/COD0OPSFIN_20241200000_01D_30S_CLK.CLK `
  --dcb-file data/CODE/2024/P1P22404.DCB `
  --antex-file data/igs14.atx `
  --iers-c04 data/IERS/eopc04_IAU2000.txt `
  --enable-phase-windup --enable-relativity `
  --ar-min-epochs 6 --gravity-nmax 90
```

**0.5h (最佳精度)**:
```powershell
py -3.12 run_sequential_pod.py `
  --date 2024-04-29 --hours 0.5 --interval 30 --grace-id C `
  --dynamics-mode simplified `
  --sp3-file data/CODE/2024/COD0OPSFIN_20241200000_01D_05M_ORB.SP3 `
  --clk-file data/CODE/2024/COD0OPSFIN_20241200000_01D_30S_CLK.CLK `
  --dcb-file data/CODE/2024/P1P22404.DCB `
  --antex-file data/igs14.atx `
  --iers-c04 data/IERS/eopc04_IAU2000.txt `
  --enable-phase-windup --enable-relativity `
  --ar-min-epochs 6 --chi2-threshold 100 --gravity-nmax 90
```

**全天线评估**:
```powershell
py -3.12 eval_day2.py
# 输出: results/accuracy_2024-04-29.png
```

**Orekit 动力学** (追加):
```powershell
$env:JAVA_HOME = "C:\Program Files\JetBrains\PyCharm Community Edition 2024.3.5\jbr"
$env:OREKIT_DATA_PATH = "data/orekit"
# ... 同上, 追加 --dynamics-mode orekit
```

### CLI 参数完整列表

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--date` | (必需) | YYYY-MM-DD |
| `--hours` | 0.5 | 弧段长度 [h] |
| `--interval` | 30.0 | EKF 步长 [s] |
| `--data-dir` | ./data | 数据根目录 |
| `--grace-id` | C | 卫星 ID (C 或 D) |
| `--gravity-nmax` | 90 | GGM05C 截断阶数 |
| `--sigma-acc` | 1e-3 | 未建模加速度 [m/s²] |
| `--sigma-phase` | 0.20 | 相位 σ 天顶 [m] |
| `--sigma-code` | 0.30 | 码 σ 天顶 [m] |
| `--chi2-threshold` | 25 | χ² 检验阈值 |
| `--tau-emp` | 600.0 | 经验力相关时间 [s] |
| `--ar-min-epochs` | 10 | WL 固定前最少历元数 |
| `--enable-phase-windup` | False | 相位缠绕修正 |
| `--enable-relativity` | False | 相对论 Shapiro |
| `--enable-cycle-slip` | False | MW+GF 周跳检测 |
| `--dynamics-mode` | simplified | simplified / orekit |
| `--sp3-file` | auto | SP3 文件路径 |
| `--clk-file` | auto | CLK 文件路径 |
| `--dcb-file` | auto | DCB 文件路径 |
| `--antex-file` | auto | ANTEX 文件路径 |
| `--iers-c04` | auto | IERS C04 文件路径 |
| `--batch-ar` | False | 两步批处理 AR |
| `--batch-lsq` | False | BatchLSQ (实验性) |

---

## 6. 版本历史

| 版本 | 日期 | 关键改进 | 0.17h RMS |
|------|------|----------|-----------|
| V2.2.1 | 2026-05 | 基线 Float PPP, GGM05C N=90 | 0.936 m |
| V2.2.2 | 2026-06 | + CODE 精密产品 + WL-AR + NL SD | 0.815 m |
| V2.2.3 | 2026-06-16 | + 仰角加权 + 自适应 code/chi2 | 0.323 m |
| **V2.2.4** | **2026-06-16** | **+ 自适应 clock_rw + 码平滑缓冲区** | **0.293 m** |

### V2.2.4 变更明细

| 文件 | 行 | 变更 |
|------|-----|------|
| `src/sequential_filter.py` | +6 | `+_clock_rw` 配置 (cfg.get) |
| | +5 | `+_cp_buf` 码平滑缓冲区 (per-SV deque) |
| | ±1 | `self.sigma_zwd_rw` 默认保持 1e-9 (1e-10 回退) |
| | ±5 | `_build_process_noise`: `Q[10,10] = self._clock_rw * dt` |
| | +5 | `_cp_buf.pop(sv)` 在周跳/剪枝时清理 |
| | +8 | `amb_init`: 中位数替代单历元码相 |
| `run_sequential_pod.py` | +1 | `clock_rw: 0.0004 if hours<0.3 else 0.001` |
| `references/POD_Improvement_Roadmap.md` | +80 | Phase 7-10 方案 |
| `memory/phase7_1_results.md` | — | 新增记忆文件 |
| `VERSION.md` | — | 本文档 (V2.2.4 改写) |

### EKF 架构已知限制

- 新 SV 加入需 ~10 epochs (~5 min) 收敛模糊度
- 长弧段 (>0.5h) clock drift ~0.3 m/h → 绝对 N1 参考不可靠
- 无整数钟产品 → 无法用 LAMBDA 全模糊度搜索
- EKF 方法极限: ~0.25 m (0.17h), ~0.80 m (0.5h)

### 数据可用性

- 2024-04-29 完整 GPS1B 数据 (8640 epochs, 24h, 10s 采样)
- GNV1B 参考轨道 (86400 epochs, 1s 采样)
- CODE Final SP3+CLK+DCB (已下载)
- 其他日期 (2024-04-30 ~ 2024-05-08) L1B 数据已归档在 `data/gracefo/2024/`

---

## 7. 数据目录结构

```
d:\prj\gnss_pod\data\
├── GPS1B_2024-04-29_C_04.pkl             # 24h GPS L1/L2/P1/P2 观测
├── GPS1B_2024-04-29_C_04.rnx             # RINEX 源文件 (备份)
├── GNV1B_2024-04-29_C_04.txt             # 参考轨道 (ASCII, ~23MB)
├── igs14.atx                              # IGS 天线文件
├── gravity/
│   └── GGM05C.gfc                         # 重力场 (Nmax=360, 使用 ≤90)
├── CODE/2024/
│   ├── COD0OPSFIN_20241200000_01D_05M_ORB.SP3   # GPS 精密轨道 (5min)
│   ├── COD0OPSFIN_20241200000_01D_05M_ORB.pkl   # SP3 pickle 缓存
│   ├── COD0OPSFIN_20241200000_01D_30S_CLK.CLK   # GPS 精密钟差 (30s)
│   ├── P1P22404.DCB                              # P1-P2 DCB
│   └── P1C12404.DCB                              # P1-C1 DCB
├── IERS/
│   └── eopc04_IAU2000.txt                 # 地球定向参数
├── orekit/                                # Orekit 数据 (可选)
│   ├── UTC-TAI.history
│   ├── eopc04_IAU2000.txt
│   └── GGM05C.gfc
├── gracefo/2024/2024-04-29/               # GRACE-FO L1B 原始文档
│   ├── GNV1B_2024-04-29_C_04.txt
│   ├── GPS1B_2024-04-29_C_04.txt
│   ├── CLK1B_2024-04-29_C_04.txt
│   └── ... (ACC1B, SCA1B, KBR1B, etc.)
└── 2024/                                   # (空目录, 未用)
```
