# GRACE-FO PPP-AR POD — Version 2.3.0

**Release Date**: 2026-06-20  |  **Frozen**: 2026-06-20

---

## 1. 软件标识

| 项目 | 值 |
|------|-----|
| 版本号 | **V2.3.0** |
| 发布日 | 2026-06-20 |
| 基于 | V2.2.4 (2026-06-16) |
| 适用卫星 | GRACE-FO C / D + 10 种 LEO 卫星 (配置库) |
| 观测系统 | GPS L1/L2 (双频 IF) |
| 滤波方法 | 序贯扩展卡尔曼滤波 (Sequential EKF) |
| 模糊度固定 | WL PPP-AR (自标定) + Arc AR (弧段级) + Batch Solver (全弧段) |
| 动力学 | GGM05C N=90 (简化) / Orekit v13.1.5 (可选) |
| 许可 | Internal research code |

### V2.3.0 新增功能 (vs V2.2.4)

| 功能 | 文件 | 描述 |
|------|------|------|
| **Framework v3 Batch Solver** | `src/batch_solver.py`, `src/batch_orbit_v3.py` | 固定 EKF 轨道, 全局 clock+ZWD+amb LSQ |
| **Arc-based AR** | `src/arc_ambiguity.py` | 弧段级 MW/B_if 平均 → 整数固定 |
| **Satellite Config DB** | `src/satellite_config.py` | 11 颗 LEO 卫星物理参数库 |
| **POD I/O** | `src/pod_io.py` | 标准 RINEX 发现 + 产品加载 + 精度指标 |
| **CLK1B 钟差解析** | `run_sequential_pod.py --clk1b` | GRACE-FO USO 钟差 ~0.03ns 精度 |
| **CODE OSB 解析** | `src/batch_lsq.py::read_code_osb()` | SINEX BIAS 格式, 32 颗 GPS SV |
| **Multi-date validation** | `run_multiday_v2.py` | 多日期批量运行 |
| **全天线 Batch 评估** | `fullday_batch_v3.py` | 24h × 12 Batch solver 弧段 |
| **解析 STM Jacobian** | `src/batch_orbit_v3.py` | r0/v0 6 维 Jacobian 已验证 (ratio=1.00, angle=0.0° vs FD) |

---

## 2. 软件依赖

| 组件 | 版本 | 备注 |
|------|------|------|
| Python | 3.12.6 | |
| numpy | 2.2.6 | |
| scipy | 1.17.1 | 仅 `linalg.solve` |
| matplotlib | 3.10.9 | 仅评估脚本 |
| astropy | (any) | ECI↔ECEF 坐标变换 |
| orekit-jpype | 13.1.5 | 可选 (`--dynamics-mode orekit`) |
| Java (JBR) | OpenJDK 21.0.6 | 可选, PyCharm 自带 |

---

## 3. 外部数据产品

| 产品 | 目录 | 来源 | 延迟 |
|------|------|------|------|
| GPS1B RINEX | `data/GPS1B_*.rnx` | JPL/GRACE-FO L1B | 1-3天 |
| GNV1B 参考轨道 | `data/GNV1B_*.txt` | JPL/GRACE-FO L1B | 1-3天 |
| CLK1B 接收机钟差 | `data/CLK1B_*.txt` | JPL/GRACE-FO L1B | 1-3天 |
| CODE Final SP3 | `data/CODE/2024/COD*ORB.SP3` | CODE/CDDIS | 12-18天 |
| CODE 30s CLK | `data/CODE/2024/COD*CLK.CLK` | CODE/CDDIS | 12-18天 |
| CODE P1P2 DCB | `data/CODE/2024/P1P2*.DCB` | CODE | 月 |
| CODE P1C1 DCB | `data/CODE/2024/P1C1*.DCB` | CODE | 月 |
| CODE OSB | `data/CODE/2024/COD*OSB.BIA` | CODE | 12-18天 |
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
│   ├── src/empirical.py       ← RTN 经验力 (ECI 变换)
│   ├── src/troposphere.py     ← Saastamoinen ZHD + GMF 映射函数
│   ├── src/coordinates.py     ← ECI↔ECEF (astropy)
│   ├── src/ambiguity.py       ← MW 组合 + WL/NL 常量
│   ├── src/cycle_slip.py      ← TurboEdit MW+GF 周跳检测
│   └── src/measurement_corrections.py ← PCO + 相位缠绕 + 相对论
├── src/orbit_integrator.py    ← RK4 积分 + STM 传播
├── src/batch_solver.py        ← BatchLinearSolver (全弧段 clock+ZWD+amb)
├── src/batch_lsq.py           ← BatchAmbiguityResolver + CODE OSB 解析
├── src/batch_orbit_v3.py      ← Framework v3: 解析 STM Jacobian + GN loop
├── src/batch_orbit_v2.py      ← Framework v2: FD Jacobian (备用)
├── src/arc_ambiguity.py       ← ArcTracker + ArcAmbiguityResolver (Phase 12)
├── src/orekit_bridge.py       ← Orekit NumericalPropagator 桥接
├── src/precision_products.py  ← CODE 产品读取器 (CLK/ANTEX/DCB/IERS)
├── src/satellite_config.py    ← 11颗卫星物理参数库 (Phase 13)
├── src/pod_io.py              ← 标准 I/O 抽象
├── validate_v224.py           ← 全天线验证 (48×0.17h + 12×0.5h)
├── fullday_batch_v3.py        ← 全天线 Batch solver 评估
├── eval_day2.py               ← 6时段评估脚本
├── run_multiday_v2.py         ← 多日期批量运行
└── run_v230.ps1               ← 一键运行 PowerShell 脚本
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

2. 每颗 SV 按序处理:
   a. 卫星天线 PCO 修正 + ECI 转换
   b. 相对论 Shapiro 修正 (L1/L2)
   c. 相位缠绕修正 (Wu et al. 1993)
   d. 码观测更新 → 锚定 clock (无 amb 项)
   e. 相位观测更新 → 锚定 ambiguity (有 amb 项)
   f. MW 积累 → WL 固定 (≥6 epochs, σ<0.35 cyc)
   g. WL 渐进收紧: P_amb 100→0.10→0.04→0.01 m²
   h. NL SD 整数 bootstrapping (实时)

3. 接收机端 WL 偏差自标定 (b_r_wl):
   中位数 (小数部分 of Σ MW per SV)

4. 死卫星剪枝: timeout=60s
```

---

## 5. 已验证精度

**测试条件**: 2024-04-29, GRACE-FO C, GPS only, GGM05C Nmax=90, 30s EKF step  
**参照**: JPL GNV1B 参考轨道 (~2cm 精度)

| Arc | Float PPP | V2.2.4 EKF | Arc AR (Phase 12) | Batch Solver Phase | 说明 |
|-----|-----------|-----------|-------------------|-------------------|------|
| 0.17h (10min) | 0.936 m | **0.293 m** | 0.246 m | 0.160 m | 轨道 0.293m; Phase 0.16m (Batch) |
| 0.5h (30min) | 1.817 m | **0.986 m** | 1.067 m | 0.211 m | 轨道 0.986m; Phase 0.21m (Batch) |

**累计提升 (0.17h)**:
```
Float PPP  →  EKF WL-AR  →  V2.2.3  →  V2.2.4 EKF  →  Arc AR  →  Batch Phase
 0.936m       0.805m         0.323m     0.293m         0.246m      0.160m (phase)
               (-14%)        (-60%)     (-69%)         (-74%)
```

### 运行命令 (V2.3.0)

**EKF 最佳 (0.293m at 0.17h)**:
```powershell
py -3.12 run_sequential_pod.py `
  --date 2024-04-29 --hours 0.17 --interval 30 --grace-id C `
  --dynamics-mode simplified `
  --sp3-file ../data/CODE/2024/COD0OPSFIN_20241200000_01D_05M_ORB.SP3 `
  --clk-file ../data/CODE/2024/COD0OPSFIN_20241200000_01D_30S_CLK.CLK `
  --dcb-file ../data/CODE/2024/P1P22404.DCB `
  --antex-file ../data/igs14.atx `
  --iers-c04 ../data/IERS/eopc04_IAU2000.txt `
  --enable-phase-windup --enable-relativity `
  --ar-min-epochs 6 --gravity-nmax 90
```

**EKF + Batch Solver (Phase RMS 评估)**:
```powershell
# 追加 --batch-lsq-v2
```

**Arc AR (0.246m)**:
```powershell
# 追加 --arc-ar
```

**CODE OSB (实验性)**:
```powershell
# 追加 --osb-file ../data/CODE/2024/COD0OPSFIN_20241200000_01D_01D_OSB.BIA
```

**CLK1B 接收机钟差 (实验性)**:
```powershell
# 追加 --clk1b ../data/CLK1B_2024-04-29_C_04.pkl
```

**全天线评估**:
```powershell
py -3.12 validate_v224.py              # 48×0.17h + 12×0.5h → PNG
py -3.12 fullday_batch_v3.py           # 12×0.17h Batch solver → PNG
```

**一键运行**:
```powershell
powershell -File run_v230.ps1 -Mode 017       # 仅 0.17h
powershell -File run_v230.ps1 -Mode batch     # 0.17h + Batch solver
powershell -File run_v230.ps1 -Mode all       # 0.17h + 0.5h
powershell -File run_v230.ps1 -Mode fullday   # 全天线验证
```

---

## 6. 版本历史

| 版本 | 日期 | 关键改进 | 0.17h RMS |
|------|------|----------|-----------|
| V2.2.1 | 2026-05 | 基线 Float PPP, GGM05C N=90 | 0.936 m |
| V2.2.2 | 2026-06 | + CODE 精密产品 + WL-AR + NL SD | 0.815 m |
| V2.2.3 | 2026-06-16 | + 仰角加权 + 自适应 code/chi2 | 0.323 m |
| V2.2.4 | 2026-06-16 | + 自适应 clock_rw + 码平滑缓冲区 | 0.293 m |
| **V2.3.0** | **2026-06-20** | **+ Arc AR + Batch Solver + 多卫星架构** | **0.246m** |

### V2.3.0 变更明细 (vs V2.2.4)

**新增文件 (7)**:
| 文件 | 描述 |
|------|------|
| `src/arc_ambiguity.py` | ArcTracker + ArcAmbiguityResolver (弧段级 AR) |
| `src/batch_orbit_v3.py` | Framework v3: 解析 STM Jacobian + GN loop |
| `src/batch_orbit_v2.py` | Framework v2: FD Jacobian (备用) |
| `src/satellite_config.py` | 11 颗 LEO 卫星参数库 |
| `src/pod_io.py` | 标准 I/O 抽象 |
| `fullday_batch_v3.py` | 全天线 Batch solver 评估 |
| `run_multiday_v2.py` | 多日期批量运行 |

**修改文件 (2)**:
| 文件 | 变更 |
|------|------|
| `src/sequential_filter.py` | +mw_max_epochs (弧段级 MW), +CLK1B/OSB 配置 |
| `run_sequential_pod.py` | +--arc-ar, --batch-lsq-v2, --clk1b, --osb-file, --batch-v3 |

### CLI 参数完整列表 (V2.3.0)

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
| `--chi2-threshold` | auto | χ² 检验阈值 (自适应) |
| `--tau-emp` | 600.0 | 经验力相关时间 [s] |
| `--ar-min-epochs` | 6 | WL 固定前最少历元数 |
| `--enable-phase-windup` | False | 相位缠绕修正 |
| `--enable-relativity` | False | 相对论 Shapiro |
| `--enable-cycle-slip` | False | MW+GF 周跳检测 |
| `--dynamics-mode` | simplified | simplified / orekit |
| `--sp3-file` | auto | SP3 文件路径 |
| `--clk-file` | auto | CLK 文件路径 |
| `--dcb-file` | auto | DCB 文件路径 |
| `--antex-file` | auto | ANTEX 文件路径 |
| `--iers-c04` | auto | IERS C04 文件路径 |
| `--arc-ar` | False | Arc-based AR (Phase 12) |
| `--batch-lsq-v2` | False | Batch Solver (Framework v3) |
| `--batch-lsq` | False | 两步批处理 AR (Phase 10) |
| `--osb-file` | None | CODE OSB BIA 文件 (Phase 9) |
| `--clk1b` | None | CLK1B 接收机钟差 pkl |

### EKF 架构已知限制

- 新 SV 加入需 ~10 epochs (~5 min) 收敛模糊度
- 长弧段 (>0.5h) clock drift ~0.3 m/h → 绝对 N1 参考不可靠
- 无整数钟产品 → 无法用 LAMBDA 全模糊度搜索
- EKF 方法极限: ~0.25 m (0.17h), ~0.80 m (0.5h)
- Batch solver 改善 phase RMS 42-83%, 但无法改善轨道 (需 GN 外层)

### 触达 0.10m 的剩余路径

三个组件各自可用但阻塞于 ECI↔ECEF 帧一致性:
1. BatchLinearSolver (phase RMS 0.16-0.21m) ✅
2. CODE OSB (32 颗 GPS SV) ✅
3. 解析 STM Jacobian (ratio=1.00 vs FD) ✅

修复方案: 全域统一帧 (ECEF 旋转矩阵 或 ECI 全帧) ← 详细见 `references/POD_Improvement_Roadmap.md`
