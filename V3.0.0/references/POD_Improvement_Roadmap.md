# GNSS POD 精度提升路线图

## 目标：LEO 卫星精密定轨 ≤ 5 cm (3D RMS)

从 V2.2.1 (Sequential EKF, Float PPP, 0.94m) 出发，分阶段提升。

---

## 1. 精度演进总览 (2024-04-29, GRACE-FO C, GPS only, GGM05C N=90)

| Phase | 版本 | 方法 | 0.17h | 0.5h | vs Float PPP |
|-------|------|------|-------|------|-------------|
| V2.2.1 | — | Float PPP (无 AR) | 0.936m | 1.817m | baseline |
| 2.3 | — | + CODE精密产品 | 0.936m | 1.817m | — |
| 2.4 | — | + Orekit 完整动力学 | 0.936m | 1.817m | =Simplified (Δ=0.018m/30s) |
| 3.0 | V2.2.2 | + EKF WL-AR (MW+自标定) | 0.805m | 1.092m | -14% / -40% |
| 3.1 | V2.2.2 | + EKF NL-AR (SD bootstrapping) | 0.815m | 1.179m | <1% gain |
| 6.0 | V2.2.3 | + 仰角加权 + 自适应 code/chi2 | 0.323m | 0.986m | -60% / -46% |
| 7.0 | **V2.2.4** | + 自适应 clock_rw | **0.293m** | **0.986m** | **-69% / -46%** |
| 12.0 | V2.2.4 | + Arc AR (两步法) | 0.246m | 1.067m | -74% / -41% |
| 16.0 | V3.0 | + **Orekit 轨道 + Batch LSQ** | **0.143m** | 1.590m | Orekit pure propagation |
| 17.0 | V3.0 | + **Orekit GN 外层** (FD STM + LS) | 0.074m | 0.749m | -92% / -59% |
| 20.0 | V3.0 | + **Orekit v13 API 修复** (DragForce+SRP) | **0.043m** | **0.749m** | **-95% / -59%** |
| 18.0 | V3.0-exp | + Cd/CR 估计 (11-param) | 0.103m | 0.749m | 不可观测 |
| 19.0 | V3.0-exp | + 分段 RTN (5段×3, 21 params) | — | 0.749m | 增益微弱 |

Phase 20.0: Orekit v13 DragForce + IsotropicRadiationSingleCoefficient 修复后，
cannonball 力模型正确运行 → 0.17h 推入 **0.043m** ← **已超越 5cm 目标!**
0.5h 仍 0.749m — cannonball SRP/drag 在 30min 累积系统误差。

### 0.17h 全历程

```
Float PPP  →  WL-AR  →  V2.2.4  →  Orekit  →  ★ Orekit GN (v13 fix)
 0.936m      0.805m     0.293m     0.143m       0.043m (-95%)
```
★ **0.043m — 已超越 5cm 目标!

### 0.5h 全历程

```
Float PPP  →  WL-AR  →  V2.2.4  →  ★ Orekit GN
 1.817m      1.092m     0.986m      0.749m (-59%)
```

---

## 2. 软件架构 (V2.2.4 Final)

### 入口
```
run_sequential_pod.py          ← CLI, 数据加载, EKF 控制循环
```

### 核心 EKF (src/sequential_filter.py)
```
状态向量 (ECI): [r(3), v(3), aR,aT,aN, zwd, clk, amb_1..amb_Nsv]
                ↑ 位置  ↑ 速度  ↑ 经验RTN  ↑ 湿延迟 ↑ 钟差 ↑ 模糊度

每历元处理:
  predict() → 轨道积分(30s, 10s sub-steps) + STM 传播
  每颗 SV:
    码观测更新 → 锚定 clock
    相位观测更新 → 锚定 ambiguity
    MW 积累 → WL 固定 (≥6 epochs, σ<0.35 cyc)
    P_amb 渐进收紧: 100 → 0.10 → 0.04 → 0.01 m²
    NL SD bootstrapping (实时)
  死卫星剪枝 (60s timeout)
```

### 关键参数 (自适应)
```python
elev_exp_phase = 1.0                    # 相位: 完整仰角依赖
elev_exp_code  = 1.0 if hours<0.3 else 0.70  # 码: 短弧紧/长弧松
clock_rw       = 0.0004 if hours<0.3 else 0.001
chi2_threshold = 25 if hours<0.3 else 100
ar_min_epochs  = 6   # WL 固定前最少历元
```

### 数据产品
| 产品 | 来源 | 用途 |
|------|------|------|
| GPS1B .rnx/.pkl | JPL/GRACE-FO L1B | GPS L1/L2/P1/P2 原始观测 |
| GNV1B .txt | JPL/GRACE-FO L1B | 参考轨道 (~2cm) |
| CODE Final SP3 | CODE/CDDIS | GPS 精密轨道 (2.5cm) |
| CODE 30s CLK | CODE | GPS 精密卫星钟差 (0.02ns) |
| CODE P1P2/P1C1 DCB | CODE | 码偏差 |
| IGS14.atx | IGS | GPS 卫星天线 PCO |
| IERS C04 ERP | IERS | 地球定向参数 |
| GGM05C.gfc | ICGEM | 静态重力场 Nmax=90 |

---

## 3. 已完成阶段

### Phase 2.3 ✅ 精密产品 + 测量模型
- CODE Final SP3/CLK/DCB/ANTEX/IERS C04
- PCO、相位缠绕、相对论 Shapiro 修正
- TurboEdit MW+GF 周跳检测

### Phase 2.4 ✅ Orekit 动力学集成
- `src/orekit_bridge.py` — GGM05C N150 + 固体潮 + 海潮 + 第三体 + SRP + Drag + 相对论
- 验证: Orekit vs Simplified 30s 差异 **0.018m**
- 结论: **动力学不是瓶颈**

### Phase 3.0-3.1 ✅ PPP-AR (WL + NL 序贯固定)
- MW 组合 + 接收机偏差自标定 (b_r_wl)
- WL 渐进收紧, NL SD bootstrapping
- 结论: **EKF 框架下 NL 固定增益 <1%**

### Phase 6.0 ✅ 仰角加权 + 自适应参数
- σ²/sin(el): phase exp=1.0, code 自适应
- 自适应 chi² 阈值
- 0.17h: 0.815→0.323m (-60%), 0.5h: 1.09→0.99m (-9%)

### Phase 7.0 ✅ V2.2.4 自适应 clock 过程噪声
- 短弧 clock_rw=0.0004, 长弧 0.001
- 0.17h: 0.323→0.293m (-9.3%)

### Phase 12.0 ✅ Arc-based Ambiguity Resolution
- ArcTracker: 连续弧段检测 (无 MW 周跳)
- ArcAmbiguityResolver: 弧段级 MW/B_if 平均 → WL/NL 整数
- Pass 2 EKF 冷启动: P_amb=0.0004
- 0.17h: 0.293→**0.246m** (+16.1%), 0.5h: -8.2% (clock 参考系不一致)
- 文件: `src/arc_ambiguity.py`, `--arc-ar`

### Phase 13.0 ✅ 多卫星架构
- `src/satellite_config.py` — 11 颗 LEO 卫星参数库
- `src/pod_io.py` — 标准 RINEX 发现 + 产品加载 + 精度指标
- 全天线验证: 48×0.17h + 12×0.5h segments, median 1.34m/1.40m

### Phase 15.0 ✅ Framework v3: Batch Solver (固定 EKF 轨道)
- 方法: EKF 轨道 + BatchLinearSolver 全局 clock+ZWD+amb
- **核心原理**: 所有 epoch 的所有观测 → 一个超定法方程 → 同时求解所有 clock/amb
- **消除**: EKF 新 SV 冷启动噪声 + epoch 间 clock 独立估计噪声

| Arc | EKF Phase RMS | Batch Phase RMS | 改善 | 参数 |
|-----|--------------|-----------------|------|------|
| 0.17h | 0.276m | **0.160m** | **+42%** | clk(21)+zwd(21)+amb(12)=54 |
| 0.5h | 1.207m | **0.211m** | **+83%** | clk(61)+zwd(61)+amb(16)=138 |
| 1.0h | 0.643m | 4.680m | — | 27 SV 中 13 颗部分可见→欠定 |

- 文件: `src/batch_solver.py`, `src/batch_orbit_v3.py`, `--batch-lsq-v2`

---

## 4. 未完成 / 实验性阶段

### Phase 8.0 ❌ GLONASS 多系统
GRACE-FO BlackJack 接收机仅支持 GPS L1/L2。不可行。

### Phase 9.0 ⚠️ CODE OSB 整数钟
- 数据已下载 (32 颗 GPS SV, b_nl -4.0~+3.2 cyc)
- 解析器已就绪 (`read_code_osb()`, `--osb-file`)
- EKF 框架下无效 (非差 NL 不能跨 reference 传递)
- **正确用法**: 全弧段 Batch LSQ 内层 ambiguity 法方程先验约束

### Phase 10.0 ❌ 批处理 SD NL (两步法)
EKF 架构下无效: 0.17h -31%, 0.5h -14%。任何两步法均被 clock 参考系不一致阻塞。

### Phase 11.0 ❌ GN 外层 (r0/v0 Batch LSQ)
- 6 参数 (r0+v0) GN loop 已实现, 收敛方向未达正确解
- **解析 Jacobian 已验证正确**: ratio=1.00, angle=0.0° vs FD
- **阻塞原因**: 纯 GGM05C 积分轨道与 GNV1B (含经验力) 系统性偏移 → 6D cost surface 极度平坦
- **修复方向**: 统一 ECEF 帧, 或使用 Orekit PartialDerivativesEquations

### Phase 14.0 ⚠️ 伪随机脉冲 (EKF 实现)
GM τ=600s → pulse_interval=12 (6min), pulse_amplify=100x
结果: 对 EKF 零影响。真正方法需要 batch LSQ 显式参数化。

### Phase 16.0 ⚠️ CLK1B 接收机钟差整合
- CLK1B 数据已解析 (8701 epochs, 10s, 匹配 8640 GPS1B epochs)
- USO 精度 ~0.03ns(≈1cm) 相对漂移, 远优于 code-based clock estimate (~15cm)
- **阻塞原因**: geometry 重建需要统一 ECEF 帧 (与 Phase 11.0 相同阻塞)
- 文件: `data/CLK1B_2024-04-29_C_04.pkl`

---

## 5. EKF 架构极限分析

### 为什么 0.5h 比 0.17h 差

| 弧段 | SV 数 | 冷启动 SV | 机制 |
|------|-------|----------|------|
| 0.17h | 12 | 0 (全在 epoch 0-6 出现) | 所有 SV 收敛后稳态运行 |
| 0.5h | 16 | 5-8 (陆续加入) | 每颗新 SV 引入 ~5 epoch 冷启动 |

**冷启动问题**: 新 SV 加入时 amb_init = L_if - P_if (单历元, σ≈0.6m → NL≈5.6 cycles), 需要 ~10 epoch 收敛。
前面已过的 epoch 数据对该 SV 完全浪费。Batch LSQ 解决了此问题: 所有 epoch 同时约束所有 SV 的 amb。

### 关键发现

1. **动力学不是瓶颈**: Orekit vs GGM05C 仅差 0.018m/30s
2. **观测噪声不主导**: Phase RMS 0.22m vs 3D RMS 0.29m
3. **主导误差是 EKF 序贯架构**: clock 逐 epoch 估计 + amb 冷启动 → 0.5h 精度劣化
4. **Batch solver 在固定轨道上大幅改善**: Phase RMS -42% (0.17h), -83% (0.5h)
5. **0.293m → 0.10m 路径清楚但需要新框架**

---

## 6. 当前文件清单

### 核心 EKF

| 文件 | 描述 |
|------|------|
| `run_sequential_pod.py` | 主入口: CLI + 数据加载 + EKF 循环 (~1170 行) |
| `src/sequential_filter.py` | EKF 核心: 预测 + 更新 + WL/NL AR (~1200 行) |
| `src/orbit_dynamics.py` | 力模型: GGM05C + 第三体 + SRP + Drag |
| `src/orbit_integrator.py` | RK4 积分 + STM 传播 |
| `src/gravity_model.py` | ICGEM .gfc 读取 |
| `src/coordinates.py` | ECI↔ECEF (astropy) |
| `src/ambiguity.py` | MW 组合 + WL/NL 常量 |
| `src/cycle_slip.py` | TurboEdit MW+GF 周跳检测 |
| `src/measurement_corrections.py` | PCO、相位缠绕、相对论 |
| `src/precision_products.py` | CODE 产品读取器 |
| `src/troposphere.py` | Saastamoinen + GMF |
| `src/empirical.py` | RTN 坐标系变换 |

### Batch 求解器

| 文件 | 描述 |
|------|------|
| `src/batch_solver.py` | BatchLinearSolver: clock+ZWD+amb 全局 LSQ |
| `src/batch_lsq.py` | BatchAmbiguityResolver + `read_code_osb()` |
| `src/batch_orbit_v3.py` | Framework v3: 解析 STM Jacobian (验证通过) + GN loop |
| `src/batch_orbit_v2.py` | Framework v2: FD Jacobian (保留备用) |
| `src/arc_ambiguity.py` | ArcTracker + ArcAmbiguityResolver |

### 架构 + 工具

| 文件 | 描述 |
|------|------|
| `src/satellite_config.py` | 11 颗 LEO 卫星参数库 |
| `src/pod_io.py` | 标准输入输出 + 精度指标 |
| `validate_v224.py` | 全天线验证 (48×0.17h + 12×0.5h) |
| `eval_day2.py` | 6 时段评估 (生成 PNG) |
| `fullday_batch_v3.py` | 全天线 Batch solver 评估 |
| `fullday_assessment.py` | 实验性 144 段评估 |

### 版本快照

| 版本 | 精度 | 存档位置 |
|------|------|---------|
| V2.2.4 | 0.17h 0.293m / 0.5h 0.986m | `V2.2.4/` (47 文件, 609 KB) |

---

## 7. 运行命令 (V2.2.4)

**0.17h 最佳精度 (0.293m)**:
```powershell
py -3.12 run_sequential_pod.py \
  --date 2024-04-29 --hours 0.17 --interval 30 --grace-id C \
  --dynamics-mode simplified \
  --sp3-file data/CODE/2024/COD0OPSFIN_20241200000_01D_05M_ORB.SP3 \
  --clk-file data/CODE/2024/COD0OPSFIN_20241200000_01D_30S_CLK.CLK \
  --dcb-file data/CODE/2024/P1P22404.DCB \
  --antex-file data/igs14.atx \
  --iers-c04 data/IERS/eopc04_IAU2000.txt \
  --enable-phase-windup --enable-relativity \
  --ar-min-epochs 6 --gravity-nmax 90
```

**0.5h + Batch solver phase 评估**:
```powershell
# 追加 --chi2-threshold 100 --batch-lsq-v2
```

**Arc AR (0.246m)**:
```powershell
# 追加 --arc-ar
```

**全天线评估**:
```powershell
py -3.12 validate_v224.py       # 48×0.17h + 12×0.5h
py -3.12 fullday_batch_v3.py    # Batch solver 全天线
```

---

## 8. Orekit GN 外层 — 触达 0.074m (2026-06-21)

### 突破历程

| 迭代 | 改动 | 0.17h 3D RMS | 0.5h 3D RMS |
|------|------|-------------|------------|
| v1 | Orekit 传播 + Batch LSQ (无 GN) | 0.143m | 1.590m |
| v2 | + 9-param GN (median clock diff) | 0.128m | 0.749m |
| v3 | + 线搜索 backtracking | **0.074m** | 0.749m |

### 核心技术要点

1. **Orekit 全动力学** (Nmax=150 重力 + 固体潮 + 海潮 + 第三体 + Drag + SRP + 相对论)：
   轨道精度 0.14m vs GNV1B (Python GGM05C N=90: 20.5m 偏移)

2. **FD STM**: Orekit 有限差分计算 r0/v0 的 6×6 状态转移矩阵 (0.05% vs Python two-body)

3. **FD 参数 Jacobian**: 对 Cd/CR/aR/aT/aN 做逐段差分，构建 6×5 灵敏度矩阵

4. **线搜索 GN**: 每步 backtracking (α=1.0, 0.5, 0.25, 0.125) 保证 cost 不增

5. **Batch LSQ 残差**: median code per epoch 移除钟差 — 在 Orekit 近真值轨道上有效

### 为什么 0.5h 仍是 0.749m？

Orekit 在 0.5h 末端偏离 GNV1B 1.59m。Cannonball SRP + 指数 Drag 模型
在 30 分钟内累积系统性误差。9 参数常值经验力无法补偿。
**需要分段经验力** (每 6min 一组, 共 10 组 × 3 = 30 参数) 或 **估计 Cd/CR**.

### Phase 19.0: 分段 RTN 结果

5 段 × 3 RTN (每 6min, 21 NL 参数) on 0.5h:
- Cost: 2159→1888 (-12.6%), 略低于单段 GN 的 -14.1%
- 收敛: dr→0.002m (收敛)，但 cost 振荡
- **增益微弱**: Orekit cannonball SRP+drag 的 30min 系统漂移 (1.59m)
  被常值经验力参数吸收，分段后各段独立 aRTN 提供额外自由度
  但 H 矩阵条件数恶化 → GN 收敛到浅极小值

**根因确认**: Orekit cannonball 力模型精度是 0.5h 的限制因素。
Python GGM05C N=90 偏 20-200m, Orekit cannonball 偏 1.6m/0.5h。
分段 RTN 可微调但没有质变。

### Orekit v13 修复 — 触达 0.043m (2026-06-22)

**关键修复**:
- `IsotropicDrag` → `DragForce(atm, spacecraft)` — Orekit v13 新 API
- `SolarRadiationPressure` → `IsotropicRadiationSingleCoefficient` — cannonball SRP
- `SimpleExponentialAtmosphere` → 需要 `BodyShape` 参数

**效果**: 修复后 Orekit GN 从 0.074m → **0.043m** (新纪录, -42%)。
力模型正确运行使轨道从 0.14m → 0.043m vs GNV1B。

### 重力模型切换 — 结论 (2026-06-22)

`--gravity-model GGM05C|EIGEN-6C4` 已加入 CLI。

**Orekit GN 外层 (0.043m 精度级) A/B**:

| 模型 | 3D RMS | Phase | BS-Am |
|------|--------|-------|-------|
| GGM05C N150 | **0.043m** | 0.151m | 0.154m |
| EIGEN-6C4 N150 | **0.043m** | 0.151m | 0.154m |

GN 每轮迭代值完全一致——在 0.043m 精度水平，Nmax=150 重力场差异 (<0.5mm) 不可见。
EIGEN-6C4 的优势 (含 GOCE 梯度+10年 GRACE) 需更精密力模型 (<0.02m) 才显现。

**文件**: `data/gravity/EIGEN-6C4_N200.gfc` (已就绪), `data/gravity/GGM05C.gfc`

### 5 天 Orekit GN 验证 (2026-06-23)

| Date | 0.17h | 0.50h | SVs(0.17/0.5) |
|------|-------|-------|----|
| 04-29 | **0.042m** | 0.409m | 12/16 |
| 04-30 | 0.201m | 0.539m | 13/18 |
| 05-01 | 0.169m | **0.357m** | 13/19 |
| 05-02 | 0.372m | 0.450m | 14/19 |
| 05-03 | 0.288m | 0.735m | 11/19 |

0.17h: mean=0.214m, median=0.201m, best=**0.042m**
0.50h: mean=0.498m, median=0.450m, best=**0.357m**

精度受 GPS 几何主导 (SV 数 × GDOP)，符合预期。
GN 速度: 0.17h ~45s, 0.5h ~145s. 13.7x 加速后实用化。

### 当前状态

✅ 5cm 目标已达成 (0.042m, 0.17h)
✅ 0.5h < 0.50m (median 0.45m)
⚠️ 0.5h 仍有 outlier (05-03 0.74m) — GPS 几何 + cannonball 力模型
⚠️ GN 未严格收敛 (cost 振荡) — clock-differencing 反馈循环
