# GRACE-FO PPP-AR POD: 输入输出规范与软件架构 (V3.2)

## 1. 软件架构

```
                           +----------------------------+
                           |         CLI                |
                           |  eval_5day_orekit.py       |
                           |  run_sequential_pod.py     |
                           +------------+---------------+
                                        |
               +------------------------+------------------------+
               |                        |                        |
               v                        v                        v
     +------------------+   +--------------------+   +--------------------+
     | (1) 数据加载层    |   | (2) 算法核心层      |   | (3) 输出与验证层    |
     |                  |   |                    |   |                    |
     | GPS 观测加载      |   | EKF 序贯滤波       |   | 3D/3V RMS 统计     |
     | GPS 轨道/钟差     |   | Batch 线性求解     |   | 精度图表 (PNG)      |
     | 天线/DCB/EOP     |   | Orekit GN 外层     |   | 结果存档 (PKL)      |
     | 广播/精密可选     |   | 模糊度固定 (WL/NL)  |   | QC 质量报告        |
     | 重力场模型        |   | 数据质量管理       |   |                    |
     | 卫星宏模型参数    |   | 伪距自主初轨       |   |                    |
     +------------------+   +--------------------+   +--------------------+
```

### 1.1 模块依赖

```
eval_5day_orekit.py
+-- src/code_orbit.py            # 伪距自主初轨
+-- src/sequential_filter.py     # EKF 序贯滤波
|   +-- src/orbit_dynamics.py    # 力模型
|   +-- src/orbit_integrator.py  # RK4 积分 + STM
|   +-- src/gravity_model.py     # ICGEM 重力场
|   +-- src/coordinates.py       # ECI-ECEF 变换
|   +-- src/cycle_slip.py        # TurboEdit 周跳
|   +-- src/ambiguity.py         # MW/WL/NL
|   +-- src/measurement_corrections.py  # PCO/缠绕/相对论
|   +-- src/troposphere.py       # 对流层
+-- src/batch_solver.py          # 批量线性求解
+-- src/batch_orbit_v3.py        # GN 精密定轨外层
|   +-- src/orekit_bridge.py     # Orekit v13 (Java)
|   +-- src/empirical.py         # RTN 经验力
|   +-- src/srp.py               # 光压
|   +-- src/third_body.py        # 第三体
|   +-- src/solid_tides.py       # 固体潮
+-- src/data_quality.py          # 4 层 QC
+-- src/precision_products.py    # SP3/CLK/DCB/ANTEX/IERS
+-- src/sp3_loader.py            # SP3 解析
+-- src/fetch_data.py            # 广播星历
+-- src/satellite_config.py      # 卫星参数 DB
+-- src/gracefo_macro.py         # Box-wing 宏模型
+-- run_sequential_pod.py        # 单弧段入口

外部依赖: orekit-jpype 13.1.5, numpy, matplotlib, jpype1, astropy
```

---

## 2. 输入规范

### 2.1 必选输入

| 类别 | 参数 | 格式 | 说明 |
|------|------|------|------|
| LEO GPS 观测 | GPS1B pkl | `{gps_sod: {sv: {L1,L2,P1,P2,...}}}` | 双频载波相位+伪距 |
| 日期范围 | `--dates` | `YYYY-MM-DD,...` | 多天批量 |
| 弧段长度 | `--hours` | `0.17,0.50,...` | 小时 |
| 卫星标识 | `--grace-id` | `C`/`D` 等 | satellite_config.py 中定义 |

### 2.2 GPS 轨道/钟差 (二选一)

| 选择 | 参数 | 格式 | 精度 |
|------|------|------|------|
| **精密星历** (默认) | SP3 + CLK 文件 | SP3-c 5min + RINEX CLK 30s | ~2.5cm 轨道 |
| **广播星历** | `--broadcast` | IGS RINEX 2.11 nav -> pkl 缓存 | ~1m 轨道 |

切换: `--broadcast` 自动从 `data/BRDC/` 加载预计算广播轨道.

### 2.3 辅助产品

| 产品 | 路径 | 用途 |
|------|------|------|
| DCB 码偏差 | `data/CODE/{year}/P1P2{YYMM}.DCB` | 卫星硬件延迟 |
| ANTEX 天线 | `data/igs14.atx` | GPS PCO/PCV |
| IERS EOP | `data/IERS/eopc04_IAU2000.txt` | 地球定向参数 |
| OSB 偏差 (可选) | `data/CODE/{year}/*_OSB.BIA` | 非差 WL/NL 偏差 |

### 2.4 卫星动力学参数

通过 `src/satellite_config.py` 管理，数据库中含 11 颗 LEO 卫星:

| 参数 | GRACE-FO C | 说明 |
|------|-----------|------|
| mass | 580.0 kg | 卫星质量 |
| area_drag | 0.68 m^2 | 阻力截面积 |
| area_srp | 3.4 m^2 | 光压截面积 |
| CD | 2.2 | 阻力系数 |
| CR | 1.3 | 光压系数 |

添加新卫星只需在 SATELLITE_DB 中增加条目.

### 2.5 力模型参数

| 参数 | 默认值 | 可选 |
|------|--------|------|
| 重力场 | GGM05C Nmax=150 | EIGEN-6C4 Nmax=200 / 任意 .gfc |
| 固体潮 | IERS 2010 | -- |
| 海潮 | FES2004 N=50 | -- |
| 第三体 | Sun+Moon | -- |
| 光压 | Cannonball | Box-wing (GRACE-FO 8 面板) |
| 阻力 | SimpleExponential | NRLMSISE00, Harris-Priester |
| 相对论 | Schwarzschild | -- |

### 2.6 EKF + 批处理参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| sigma_phase | 0.20 m | 相位噪声 (天顶) |
| sigma_code | 0.30 m | 码噪声 (天顶) |
| sigma_acc_process | 1e-3 m/s^2 | 未建模加速度 |
| chi2_threshold | 25 (<0.3h) | 卡方检验 |
| el_min | 5 deg | 高度截止角 |
| clock_rw | 0.0004 (<0.3h) | 钟差随机游走 |
| ar_min_epochs | 6 | MW 宽巷最少历元 |
| GN max_iter | 6 | GN 外循环迭代 |
| GN damping | 0.5 | Levenberg-Marquardt |

---

## 3. 输出规范

### 3.1 精密轨道 (每历元)

| 字段 | 格式 | 单位 |
|------|------|------|
| gps_sod | float | 秒 |
| r_ecef | float[3] | 米 |
| v_ecef | float[3] | 米/秒 |
| r_eci | float[3] | 米 (可选) |
| v_eci | float[3] | 米/秒 (可选) |
| clk | float | 米 |

### 3.2 精度统计 (每弧段)

| 字段 | 说明 |
|------|------|
| rms_3d | 3D 位置 RMS (m) |
| rms_3v | 3D 速度 RMS (mm/s) |
| phase_ekf | EKF 相位残差 RMS (m) |
| phase_batch | Batch 相位残差 RMS (m) |
| phase_gn | GN 外层相位残差 RMS (m) |
| n_sv | 活跃卫星数 |
| avg_cov | 平均覆盖率 (%) |
| qc_score / qc_grade | 质量评分 (0-1, A-F) |

### 3.3 可视化

3 面板 PNG: (a) 3D+3V RMS 曲线, (b) Phase RMS 曲线, (c) 覆盖率-精度散点图.

文件名: `results/{N}day_{start}_{end}[_BRDC]_orekit/{...}.png`

---

## 4. 处理流程

```
Step 0: 运动学 WLS
  P_if -> r_ecef (~5-10m)

Step 1: 伪距 GN (--code-arc-hours)
  P_if (长弧) + 动力学 -> 轨道表 (~0.2m)
  异常检测: |resid| > 5m 剔除

Step 2: EKF 序贯滤波
  码(锚定钟差) -> 相位(锚定模糊度)
  MW -> WL 固定 -> (可选) OSB NL
  保护: TurboEdit > Chi2 > Gap > MW > 覆盖率

Step 3: Batch AR + Orekit GN
  3a: Batch 线性求解 (clk+zwd+amb)
  3b: GN 外循环 x6
    Orekit 传播 -> 重建几何 -> Batch -> Jacobian -> Newton + 线搜索
```

广播模式 (--broadcast): Step 3b 跳过，用 EKF 轨道评估.

### 数据质量保护链 (8 层)

| 层 | 检查 | 阈值 | 动作 |
|----|------|------|------|
| 1 | 观测范围 | 15,000-35,000 km | 跳过 |
| 2 | 高度角 | el < 5 deg | 跳过卫星 |
| 3 | Chi2 | > 25/100 | 拒绝观测 |
| 4 | TurboEdit | GF>0.08m, MW>4sigma | 重置 amb |
| 5 | Gap | >2 epoch | 分裂弧段 |
| 6 | MW | sigma>10cyc / sigma>0.30cyc | 拒绝 SV |
| 7 | 覆盖率 | <40% | 拒绝 SV |
| 8 | QC | <0.5 | 标记弧段 |

---

## 5. 数据格式 (可替换接口)

### 5.1 GPS1B 观测

```
{gps_sod_int: {sv: {
    'L1': float,   # L1 相位 [m]
    'L2': float,   # L2 相位 [m]
    'P1': float,   # P1 伪距 [m]
    'P2': float,   # P2 伪距 [m]
    'L_if': float, # IF 组合相位
    'P_if': float, # IF 组合伪距
    'L1_cyc': float, # L1 相位 [cycle]
    'L2_cyc': float, # L2 相位 [cycle]
    'L1_SNR': float, # L1 SNR [dB-Hz]
}}
```

### 5.2 GPS 轨道

```
{'ts': [datetime,...],
 'epochs': {datetime: {sv: [x_m, y_m, z_m, clk_m]}}}
```

clk_m 单位为米. 精密轨道和广播星历共用此格式.

### 5.3 输出轨道

当前为 Python pickle. 可扩展 SP3-c / CSV / RINEX 3.04.

---

## 6. 运行命令

```powershell
$env:OREKIT_DATA_PATH = 'd:\prj\gnss_pod\data\orekit'
$env:JAVA_HOME = '...'

# 标准精密定轨
python eval_5day_orekit.py --dates 2024-04-29,...,2024-05-08 --hours 0.17,0.50

# 广播星历
python eval_5day_orekit.py --broadcast --skip-code-orbit --dates ... --hours 0.17

# 自洽模式 (无 GNV1B)
python eval_5day_orekit.py --code-arc-hours 1.0 --dates ... --hours 0.17
```

---

## 7. 精度总览

| 版本 | 0.17h | 0.50h | 改进 |
|------|-------|-------|------|
| V2.2.4 | 0.293m | 0.986m | 自适应 clock_rw |
| V3.0.0 | 0.043m* | 0.409m | Orekit GN |
| V3.1.0 | 0.169m | 0.493m | QC + gap + IRLS |
| **V3.2.0** | **0.169m** | **0.493m** | 自洽初轨 + 速度 + BRDC |
| V3.2 BRDC | **1.355m** | -- | 真实 IGS 广播星历 |

* 单天最优

### BRDC vs CODE

| | CODE | BRDC | 退化 |
|---|------|------|------|
| Mean | 0.169m | 1.355m | 8x |
| Best | 0.047m | 0.428m | 9x |
