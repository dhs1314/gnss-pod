# GRACE-FO 精密定轨 (POD) 完整处理流程

## 目标：LEO 卫星轨道精度 ≤ 5 cm (3D RMS)

---

## 总览

```
┌──────────────────────────────────────────────────────────────────────┐
│                        POD Pipeline Overview                          │
├──────────────────────────────────────────────────────────────────────┤
│                                                                       │
│  Phase 1: 数据获取                                                    │
│    GPS1B → ACC1B → SCA1B → CLK1B → AOD1B  (GRACE-FO星上数据)         │
│    SP3 → CLK → ERP → SINEX → ANTEX → DCB   (IGS外部产品)             │
│                                                                       │
│  Phase 2: 数据预处理                                                  │
│    GPS1B二进制 → RINEX → 数据质量检查 → 粗差剔除                      │
│    时间对齐 → 坐标系转换 → 天线相位中心校正                            │
│                                                                       │
│  Phase 3: 测量模型                                                    │
│    电离层-free组合 → PCO/PCV → 相位缠绕 → 对流层改正 → 相对论        │
│                                                                       │
│  Phase 4: 力学模型                                                    │
│    重力场 → 固体潮/海潮 → 第三体 → 大气阻力 → 光压 → 经验加速度       │
│                                                                       │
│  Phase 5: 轨道确定                                                    │
│    初始轨道 → 浮点解 → 模糊度固定 → 固定解 → 残差编辑 → 最终解       │
│                                                                       │
│  Phase 6: 质量评估                                                    │
│    轨道重叠检验 → KBR残差 → SLR检核 → 精度预算                        │
│                                                                       │
└──────────────────────────────────────────────────────────────────────┘
```

---

## Phase 1: 数据获取

### 1.1 星上数据 (GRACE-FO, 卫星 C/D)

| 数据产品 | 作用 | 采样率 | 下载地址 |
|---------|------|--------|---------|
| **GPS1B** | GPS L1/L2 相位+伪距观测值 | 10 s | [GFZ ISDC](https://isdc-data.gfz.de/grace-fo) |
| **ACC1B** | 非引力加速度 (三轴) | 1 s | [GFZ ISDC](https://isdc-data.gfz.de/grace-fo) |
| **SCA1B** | 卫星姿态四元数 | 1 s | [GFZ ISDC](https://isdc-data.gfz.de/grace-fo) |
| **CLK1B** | 接收机时钟校正 | 10 s | [GFZ ISDC](https://isdc-data.gfz.de/grace-fo) |
| **AOD1B** | 大气/海洋去混频产品 | 3 h | [GFZ ISDC](https://isdc-data.gfz.de/grace-fo) |

```bash
# 下载示例 (GFZ ISDC)
wget -r -np -nH --cut-dirs=3 \
  https://isdc-data.gfz.de/grace-fo/level-1b/jpl/rl03/2026/01/GPS1B_2026-01-01_C_03.dat
wget -r -np -nH --cut-dirs=3 \
  https://isdc-data.gfz.de/grace-fo/level-1b/jpl/rl03/2026/01/ACC1B_2026-01-01_C_03.dat
wget -r -np -nH --cut-dirs=3 \
  https://isdc-data.gfz.de/grace-fo/level-1b/jpl/rl03/2026/01/SCA1B_2026-01-01_C_03.dat
```

### 1.2 外部精密产品 (IGS/CODE/IERS)

| 数据产品 | 作用 | 精度 | 延迟 | 下载地址 |
|---------|------|------|------|---------|
| **IGS Final SP3** | GPS卫星精密轨道 | ~2.5 cm | 12-18 天 | [CDDIS](https://cddis.nasa.gov) |
| **IGS Final CLK** | GPS卫星精密时钟 (30s) | ~0.02 ns | 12-18 天 | [CDDIS](https://cddis.nasa.gov) |
| **IERS C04 ERP** | 极移 xp,yp / UT1-UTC / LOD | ~0.05 mas | 30 天 | [IERS](https://www.iers.org) |
| **IGS SINEX** | 地面站坐标+速度 | ~2 mm | 每周 | [CDDIS](https://cddis.nasa.gov) |
| **IGS ANTEX** | GPS卫星+地面站天线PCO/PCV | ~1 mm | 不定期 | [IGS](https://igs.org) |
| **CODE DCB** | 差分码偏差 (P1-C1, P1-P2) | ~0.1 ns | 每月 | [CODE](http://ftp.aiub.unibe.ch) |

```bash
# IGS Final 产品文件名示例
# 轨道: igs21234.sp3.Z
# 时钟: igs21234.clk.Z
# ERP:  igs21237.erp.Z
# SINEX: igs21P2123.snx.Z
```

### 1.3 辅助数据

| 数据 | 说明 |
|------|------|
| **重力场模型** | EGM2008 / EIGEN-6C4 (ICGEM 格式) |
| **海潮模型** | FES2014b (球谐系数 50阶) |
| **LEO 天线 PCO/PCV** | GRACE-FO 在轨标定 PCV 文件 |
| **卫星宏模型** | 质量、面积、光学属性 (用于阻力和光压建模) |

---

## Phase 2: 数据预处理

### 2.1 GPS1B 二进制 → RINEX 转换

```
GPS1B 二进制 (.dat)
    ↓ Bin2AsciiLevel1 工具
ASCII 观测值 + 导航电文
    ↓ 格式重组
RINEX 3.04 观测文件 (.rnx) + 导航文件 (.nav)
```

**RINEX 3.04 观测值映射：**

| GPS1B 字段 | RINEX 3 观测码 | 说明 |
|-----------|---------------|------|
| L1 phase (cycles) | `L1C` | L1 载波相位 |
| L2 phase (cycles) | `L2W` | L2 载波相位 |
| L1 pseudorange (m) | `C1C` | L1 C/A 码伪距 |
| L2 pseudorange (m) | `C2W` | L2 P/Y 码伪距 |

### 2.2 数据时间对齐

```
GPS1B 时间标签 → CLK1B 校正 → GPS 时间
SCA1B 姿态数据 → 线性内插到 GPS 观测历元
ACC1B 加速度 → 降采样 + 内插到 GPS 观测历元
IGS SP3/CLK → 内插到 GPS 观测历元 (通常用 9阶 Lagrange)
```

### 2.3 粗差剔除 & 质量检查

```python
for epoch in arc:
    for sat in visible_sats:
        # 1. SNR 阈值
        if L1_SNR < 30 dB-Hz or L2_SNR < 25 dB-Hz: REJECT
        
        # 2. 伪距-相位一致性 (码减相位差异)
        if abs((C1C - L1C*lambda_L1) - median_mw) > 30 m: REJECT
        
        # 3. 电离层变化率 (GF组合时间差分)
        if abs(delta_GF) > 0.02 m/s: REJECT  (可能周跳)
        
        # 4. 仰角截止
        if elevation < 10°: REJECT
```

### 2.4 周跳检测 (TurboEdit 算法)

```
Melbourne-Wübbena (MW) 组合:
  MW = (L1 - L2) - (f1*P1 + f2*P2)/(f1 + f2) * (f1-f2)/c
  
  周跳判定: |MW - <MW>_sliding| > 4*sigma_MW

Geometry-Free (GF) 组合:
  GF = L1 - L2
  
  周跳判定: |delta_GF| > threshold  (时间差分跳变)
```

---

## Phase 3: 测量模型

### 3.1 观测方程

**无电离层组合 (Ionosphere-Free, IF)** — 消除一阶电离层：

$$
\Phi_{IF} = \frac{f_1^2}{f_1^2 - f_2^2}\Phi_1 - \frac{f_2^2}{f_1^2 - f_2^2}\Phi_2
$$

$$
P_{IF} = \frac{f_1^2}{f_1^2 - f_2^2}P_1 - \frac{f_2^2}{f_1^2 - f_2^2}P_2
$$

**完整观测方程 (IF 组合)：**

$$
\begin{aligned}
P_{IF} &= \rho_{LEO}^{GPS} + c(\delta t_r - \delta t^s) + T_{trop} + \varepsilon_P \\
\Phi_{IF} &= \rho_{LEO}^{GPS} + c(\delta t_r - \delta t^s) + T_{trop} + \lambda_{IF} N_{IF} + \varepsilon_\Phi
\end{aligned}
$$

### 3.2 测量改正项 (依次施加)

#### 改正 1: 天线相位中心偏移 (PCO)

```python
# GPS 发射天线 PCO (从 ANTEX 文件)
pco_gps_sat = read_antex("igs14.atx", prn, "L1", "L2")

# LEO 接收天线 PCO (GRACE-FO 在轨标定值)
pco_leo_L1 = [-0.0005, 0.0002, 0.4513]  # 示例 (m), 星固坐标系
pco_leo_L2 = [-0.0003, 0.0001, 0.4487]
```

#### 改正 2: 天线相位中心变化 (PCV)

```python
# 在轨标定 PCV (球谐展开)
def pcv_correction(elevation, azimuth, pcv_coeffs, max_degree=8):
    """在轨残差 PCV 模型 (Kang et al., 2021)"""
    corr = 0
    for n in range(max_degree + 1):
        for m in range(n + 1):
            P_nm = legendre(n, m, sin(elevation))
            corr += P_nm * (pcv_coeffs.C[n,m] * cos(m*azimuth) 
                          + pcv_coeffs.S[n,m] * sin(m*azimuth))
    return corr
```

#### 改正 3: 相位缠绕 (Phase Wind-Up)

只有载波相位需要此改正：

```python
def phase_wind_up(leo_pos, gps_pos, leo_attitude, gps_attitude):
    """
    Wu et al. (1993) 相位缠绕改正
    量级: 0.5-1 周 (L1) = ~10-19 cm → 必须改正
    """
    # 有效偶极子方向
    d_leo = leo_attitude_body_z  # LEO 天线方向 (朝向天顶)
    d_gps = gps_attitude_body_z  # GPS 天线方向 (朝向地心)
    
    # 信号传播方向
    k = unit_vector(gps_pos - leo_pos)
    
    # 有效偶极子
    D_leo = d_leo - k * dot(k, d_leo) + cross(k, leo_attitude_y)
    D_gps = d_gps - k * dot(k, d_gps) - cross(k, gps_attitude_y)
    
    # 相位缠绕角
    phi = atan2(dot(cross(k, D_leo), D_gps), dot(D_leo, D_gps))
    delta_phi = round((prev_phi - phi) / (2*pi))
    return (delta_phi + phi - prev_phi) / (2*pi)  # cycles
```

#### 改正 4: 对流层延迟

```python
def troposphere_correction(leo_pos, gps_pos, epoch):
    """
    对流层天顶延迟 (Saastamoinen 模型) + 映射函数 (Vienna 1)
    量级: 2-3 m (天顶) → 10-50 m (低仰角)
    """
    # 天顶延迟
    zhd = saastamoinen_zhd(pressure, latitude, height)
    zwd = ztd - zhd  # 湿延迟作为待估参数
    
    # 映射函数: Vienna 1 (VMF1) 或 GMF
    mf_h = vmf1_hydrostatic(epoch, latitude, height, elevation)
    mf_w = vmf1_wet(epoch, latitude, elevation)
    
    return zhd * mf_h + zwd * mf_w
```

#### 改正 5: 相对论效应

```python
def relativity_correction(gps_pos, gps_vel, leo_pos, leo_vel):
    """
    Shapiro 信号延迟 (广义相对论)
    量级: ~2 cm (GPS→LEO)
    """
    R_gps = norm(gps_pos)
    R_leo = norm(leo_pos)
    rho   = norm(gps_pos - leo_pos)
    
    # Shapiro 公式
    dt_rel = 2 * MU / c**3 * log((R_gps + R_leo + rho) 
                                  / (R_gps + R_leo - rho))
    return c * dt_rel  # 距离改正 (m)
```

### 3.3 GPS 卫星位置/时钟内插

```python
def interpolate_gps_orbit_clock(sp3_file, clk_file, epoch, prn):
    """9 阶 Lagrange 内插 GPS 精密轨道和 30s 时钟"""
    # SP3: 读取 10 个节点 (epoch 前后各 5, 间隔 900s)
    # CLK: 读取 10 个节点 (epoch 前后各 5, 间隔 30s)
    
    orbit_eci = lagrange_interp_9th(sp3_positions, times, epoch)
    clock_corr = lagrange_interp_9th(clk_values, clk_times, epoch)
    
    return orbit_eci, clock_corr
```

---

## Phase 4: 力学模型

### 4.1 轨道传播方程

卫星运动的 Newton 方程 (ECI 惯性系):

$$
\ddot{\mathbf{r}} = -\frac{GM}{r^3}\mathbf{r} + \mathbf{a}_\text{non-spherical} + \mathbf{a}_\text{tides} + \mathbf{a}_\text{3rd-body} + \mathbf{a}_\text{drag} + \mathbf{a}_\text{srp} + \mathbf{a}_\text{rel} + \mathbf{a}_\text{emp}
$$

### 4.2 各力模型配置

| 力模型 | 配置参数 | 备注 |
|--------|---------|------|
| 中心引力 | $GM = 398600.4415\ \text{km}^3\text{/s}^2$ | WGS84 |
| 非球形引力 | EGM2008, 120×120 | 球谐展开 |
| 固体潮 | IERS 2010, 频率相关 | k20=0.30190 |
| 海潮 (FES2014b) | 50×50 | 必须加载! |
| 极潮 | IERS 2010 | 固体极+海极 |
| 大气/海洋去混频 | AOD1B RL06, 100×100 | STOKES 系数改正 |
| 月球 | DE430 星历 | JPL 行星历 |
| 太阳 | DE430 星历 | 含光压 |
| 大气阻力 | DTM2000, Cd=1.0-1.5 | 用ACC数据约束 |
| 太阳辐射压 | Cr=0.7-1.0, Cannonball | 含地球遮挡 |
| 广义相对论 | Schwarzschild + 一阶 | |
| 经验加速度 | RTN 分段常值, 15 min | σ=5×10⁻⁹ m/s² |

### 4.3 ACC1B 加速度计数据处理

```
原始 ACC1B (1 Hz, 卫星固连系)
    ↓ 1. 偏差+尺度因子校正
    ↓ 2. 低通滤波 (cutoff ~0.1 Hz, 去除高频噪声)
    ↓ 3. 降采样到 GPS 观测历元 (10s)
    ↓ 4. SCA1B姿态 → 转换到 ECI
校准后非引力加速度 (ECI)
    ↓ 用于大气阻力/光压模型验证
    ↓ 或直接代入轨道积分 (替代阻力+光压模型)
```

### 4.4 RTN 经验加速度配置

```python
# 经验加速度参数化
n_intervals = int(24 * 3600 / 900)  # = 96 段/天

for interval in range(n_intervals):
    for direction in ['R', 'T', 'N']:
        # 每个方向每段一个常值加速度
        param = add_parameter(
            name=f"emp_{direction}_{interval}",
            prior_value=0.0,
            prior_sigma={  # 先验约束
                'R': 1e-9,   # 径向 (最弱)
                'T': 5e-9,   # 沿轨 (最强，需要较大自由度)
                'N': 1e-9,   # 法向
            }[direction],
            lower_bound=-1e-6,
            upper_bound=+1e-6
        )
```

---

## Phase 5: 轨道确定

### 5.1 两阶段估计策略

```
阶段 1: 浮点解 (Float Solution)
  ├── 伪距 + 载波相位 (iono-free 组合)
  ├── 待估参数: 轨道初值 (6) + 经验加速度 (288) + 钟差 (每历元)
  │               + 模糊度 (浮点, 每弧段) + 天顶湿延迟 (每2h)
  ├── Batch LS / EKF
  └── 输出: 浮点轨道 ~10 cm 精度

阶段 2: 模糊度固定 + 固定解 (Fixed Solution)
  ├── 从浮点解提取宽巷+窄巷模糊度
  ├── LAMBDA 固定整数模糊度
  ├── 重新运行 Batch LS (已固定的模糊度作为约束)
  └── 输出: 固定轨道 ~2-3 cm 精度 ✅
```

### 5.2 模糊度固定流程 (AR)

```
输入: 浮点解 PPP 模糊度 N_float + 协方差 Q_N

Step 1: 宽巷固定 (Melbourne-Wübbena)
  N_wl = (L1 - L2) - (f1*P1 + f2*P2)/(f1+f2) * (f1-f2)/c
  固定判定: |N_wl - round(N_wl)| < 0.25 周 AND σ_Nwl < 0.15 周

Step 2: 星间单差 (消除接收机端误差)
  SD_ambiguities = satellite_b - satellite_ref (对每个历元的模糊度差分)

Step 3: LAMBDA 搜索
  Z^T * Q_N * Z  去相关变换
  Integer Least Squares 搜索
  Ratio Test: R = ||N_2nd_best||^2 / ||N_best||^2 > 3.0

Step 4: 部分固定 (Partial Fixing)
  若 Ratio < 3.0: 按成功率排序, 固定高成功率模糊度子集

Step 5: 约束更新
  固定解 = 浮点解 + Q_N*inv(Q_N_fixed) * (N_fixed - N_float)
```

### 5.3 BatchLSEstimator 配置

```java
// Orekit 批最小二乘 POD 配置
BatchLSEstimator estimator = new BatchLSEstimator(
    optimizer,        // Levenberg-Marquardt
    propagatorBuilder // 含全部力学模型
);

// 添加测量值 (遍历所有历元、所有可见 GPS 卫星)
for (Observation obs : preprocessedObservations) {
    // L1/L2 伪距 → iono-free 组合
    Pseudorange prIF = new Pseudorange(
        obs.epoch, obs.pr_L1, obs.pr_L2,
        sigma_PR_IF, 1.0,
        CombinationType.IONO_FREE,
        satSystem, prn);
    estimator.addMeasurement(prIF);
    
    // L1/L2 载波相位 → iono-free 组合
    Phase phaseIF = new Phase(
        obs.epoch, obs.ph_L1, obs.ph_L2,
        sigma_PH_IF, 1.0,
        CombinationType.IONO_FREE,
        satSystem, prn);
    estimator.addMeasurement(phaseIF);
}

// 估计
estimator.estimate();
```

### 5.4 待估参数汇总

| 参数类型 | 数量 (24h弧段) | 先验 σ | 说明 |
|---------|:----------:|--------|------|
| 初始位置 (3D) | 3 | 1 m | ECI 坐标 |
| 初始速度 (3D) | 3 | 0.001 m/s | ECI 坐标 |
| 经验加速度 R | 96 | 1×10⁻⁹ m/s² | 径向 |
| 经验加速度 T | 96 | 5×10⁻⁹ m/s² | 沿轨 |
| 经验加速度 N | 96 | 1×10⁻⁹ m/s² | 法向 |
| 接收机钟差 | ~8640 | 100 m | 逐历元 (10s) |
| GPS 模糊度 | ~100-200 | 10 cycles | 每弧段每卫星 |
| Cd (可选) | 1 | 0.1 | 阻力系数 |
| Cr (可选) | 1 | 0.05 | 光压系数 |
| 对流层湿延迟 | 12 | 0.01 m | 每2小时 |
| **总计** | **~9000+** | | |

---

## Phase 6: 质量评估

### 6.1 轨道重叠检验 (每日弧段 6h 重叠)

```
弧段1: [T00:00, T30:00]    24h + 前后 3h padding
弧段2: [T24:00, T54:00]
             └── 重叠: [T24:00, T30:00], 6小时

轨道差异 RMS (3D):
  期望: < 2 cm (固定解), < 5 cm (浮点解)
```

### 6.2 KBR/星间测距残差 (仅 GRACE-FO)

```
测量: KBR 星间测距 (精度 ~1 μm/s 距离变化率)
与 POD 轨道预测的星间距离差异:
  期望残差 RMS < 1 mm/s
```

### 6.3 SLR (卫星激光测距) 外部检核

```
独立于 GPS 的外部检核手段:
  SLR 残差 RMS < 2-3 cm → 5cm 轨道精度成立
```

### 6.4 精度预算表

| 指标 | 目标值 | 检验方法 |
|------|--------|---------|
| 轨道重叠 3D RMS | < 2 cm (固定) / < 5 cm (浮点) | 6h重叠弧段 |
| KBR 残差 RMS | < 1 mm/s | 星间测距 |
| SLR 残差 RMS | < 2 cm | 独立外部检核 |
| 伪距残差 RMS | < 50 cm | 观测值拟合 |
| 载波相位残差 RMS | < 5 mm | 观测值拟合 |
| 模糊度固定率 | > 90% | LAMBDA Ratio > 3.0 |

---

## 完整 Pipeline 总结 (一句话版)

```
GPS1B+ACC1B+SCA1B → 预处理(RINEX+周跳+粗差) → 测量建模(IF组合+PCV+PCO+
缠绕+对流层+相对论) + IGS精密产品(SP3+CLK+ERP) → Reduced-Dynamic力学模型
(重力场+潮汐+阻力+光压+RTN经验加速度) → 浮点解 Batch LS → LAMBDA模糊度固定
→ 固定解 Batch LS → 轨道重叠+KBR+SLR检核 → 5cm精度轨道 ✅
```

---

## 参考资料

1. Jäggi et al. (2006). Reduced-dynamic orbit determination. *Advances in Space Research*, 38(11).
2. Kang et al. (2021). GRACE-FO antenna phase center modeling and precise orbit determination. *Remote Sensing*, 13(21), 4204.
3. Kroes et al. (2005). Precise GRACE baseline determination using GPS. *GPS Solutions*, 9.
4. Teunissen (1995). The least-squares ambiguity decorrelation adjustment. *Journal of Geodesy*, 70(1).
5. Wu et al. (1993). Effects of antenna orientation on GPS carrier phase. *Manuscripta Geodaetica*, 18.
6. GRACE-FO Level 1 Data Product User Handbook (JPL D-56922).
7. IERS Conventions 2010 (Petit & Luzum, eds.).
8. Orekit Documentation: https://www.orekit.org/
