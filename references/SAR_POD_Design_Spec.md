# 遥感SAR卫星精密定轨系统 — 技术方案与设计文档

> **文档版本**: v1.0  
> **日期**: 2026-06-07  
> **用途**: 作为 AI 编程智能体（Coding Agent）的输入方案，指导完整软件实现  
> **参考验证数据**: GRACE-FO L1B 导航数据

---

## 1. 项目概述

### 1.1 任务目标

设计并实现一个模块化、可配置的遥感SAR卫星精密定轨（POD）软件系统，支持：

- **三种处理时延模式**：实时（Real-Time）、近实时（Near-Real-Time, NRT）、延时（Delayed）
- **双导航系统**：GPS、北斗（BDS），支持单独或联合定轨
- **双动力学后端**：简化动力学模型（自研） / Orekit 库（可选），通过配置切换
- **姿态数据利用**：输入姿态四元数，用于精确天线 PCO/PCV 旋转校正
- **处理时效约束**：
  - 近实时模式：2小时数据，地面处理 ≤ 2分钟
  - 延时模式：24-48小时数据，地面处理 ≤ 5分钟

### 1.2 输入数据

| 数据类型 | 格式 | 来源 | 备注 |
|---------|------|------|------|
| 星载GPS观测 | L1/L2 伪距 + 载波相位 | 导航接收机 | 10s采样（GRACE-FO L1B格式参考）|
| 星载北斗观测 | B1/B3 伪距 + 载波相位 | 导航接收机 | 同GPS格式 |
| 姿态四元数 | q0, q1, q2, q3 | 卫星姿控系统 | 用于天线PCO旋转 |
| GPS广播星历 | RINEX 导航文件 | 接收机解码 / IGS | 实时模式必需 |
| 北斗广播星历 | RINEX 导航文件 | 接收机解码 / MGEX | 实时模式必需 |
| 精密轨道 (SP3) | SP3-c/d 格式 | IGS/MGEX | 近实时/延时模式 |
| 精密钟差 (CLK) | RINEX CLK 格式 | IGS/MGEX | 近实时/延时模式 |
| 天线模型 (ANTEX) | ATX 格式 | IGS/MGEX | PCO/PCV校正 |
| DCB偏差 | DCB/BSX 格式 | CODE/CAS | 频间偏差校正 |
| ERP地球定向参数 | ERP 格式 | IERS | EOP数据 |
| SINEX文件 | SINEX 格式 | IGS/MGEX | 测站坐标（可选）|
| 重力场模型 | GFC/SHA/EGM格式 | GFZ/ICGEM | 静态 + 时变重力场 |
| 海潮模型 | OTIS/netCDF网格 | AVISO/TPXO | 海潮负荷形变计算 |
| 行星历表 | BSP二进制 | JPL | DE440日月行星位置 |
| 大气密度模型 | 软件模型(非文件) | NRLMSISE-00 | 代码内置，无外部文件 |
| 卫星宏模型 | JSON/XML | 自建 | 卫星各面面积 + 光学特性 |

### 1.2.1 数据下载来源详解

> 以下逐一列出每种输入数据的下载地址、文件命名规范、更新频率和适用场景。
> FTP链接建议用专业下载工具（wget/lftp/FileZilla）；HTTPS可直接浏览器访问。

---

#### ① 星载导航接收机原始观测数据 + 姿态数据

**说明**：这是星上的原始数据，实际操作中由卫星地面站接收下传，不对外公开下载。对于算法验证，使用 **GRACE-FO L1B 公开数据** 作为替身。

| 项目 | 详情 |
|------|------|
| **数据类型** | GPS L1A/L1B 观测 + SCA1B 姿态四元数 |
| **数据层级** | Level-1B（已做时间校正、降采样到10s） |
| **下载地址** | **NASA PODAAC** — https://podaac.jpl.nasa.gov/dataset/GRACEFO_L1B_ASCII_GRAZ_AUX_RL04 |
| **备选入口** | https://podaac-tools.jpl.nasa.gov/drive/files/allData/gracefo/L1B/ |
| **文件命名示例** | `GRA1A_2024-001-00_A_GNV1B.txt` (轨道导航，含GPS L1B观测) |
|   | `GRA1A_2024-001-00_B_SCA1B.txt` (星敏姿态四元数) |
| **更新频率** | 静态存档，已有2024.05~2026全年数据 |
| **数据格式** | ASCII文本，固定宽度字段，每字段有头说明 |
| **关键文件** | GNV1B（GPS导航解 + 观测值）、SCA1B（姿态）、VGN1B（矢量导航解） |
| **使用脚本** | 需自写RINEX转换器（GNV1B → RINEX OBS），或使用GFZ提供的转换工具 |

**GRACE-FO L1B 数据字段关键说明**：

```
GNV1B文件（导航/观测数据）:
  - GPS_time: GPS连续秒（自1980-01-06 00:00:00）
  - PRN: 卫星PRN编号 (G01-G32)
  - L1_Phase, L2_Phase: 载波相位 (单位: cycles)
  - C1_Code, C2_Code: 伪距 (单位: m)
  - LLI_L1, LLI_L2: 失锁指示器 (Loss of Lock Indicator)

SCA1B文件（星敏姿态）:
  - GPS_time: GPS连续秒
  - q0, q1, q2, q3: 姿态四元数 (ECI→SAT body frame)
  - 精度: ~0.01° (milli-degree级)
```

---

#### ② GPS广播星历（BRDC / Broadcast Navigation）

**用途**：实时模式的卫星轨道计算；近实时/延时模式的初始轨道初值

| 项目 | 详情 |
|------|------|
| **数据格式** | RINEX 2.11 / 3.0x 导航文件（`*.n` 或 `*.rnx`） |
| **下载地址** | **CDDIS** — https://cddis.nasa.gov/archive/gnss/data/daily/YYYY/DDD/YYn/ |
| **备选地址** | **IGN** — ftp://igs.ign.fr/pub/igs/data/YYYY/DDD/ |
|   | **武汉大学 IGS数据中心** — ftp://igs.gnsswhu.cn/pub/gps/data/daily/ |
| **文件命名** | `brdcDDD0.YYn` (RINEX 2) — 其中 `DDD`=年积日, `YY`=年份后两位 |
|   | `BRD400GLR_S_YYYYDDD0000_01D_MN.rnx` (RINEX 3) |
| **示例** | `brdc0010.24n` = 2024年第001天的GPS广播星历 |
| **更新频率** | 每日一个文件，约在UTC 23:00发布当日完整文件 |
| **实时获取** | IGS实时数据流 (NTRIP caster: `ntrip.igs-ip.net:2101`) |
| **文件大小** | ~80-200KB/天 |

**下载命令示例**：
```bash
# CDDIS需要注册账号后使用wget
wget -r -np -nH --cut-dirs=3 -A '*n' \
  ftp://cddis.nasa.gov/gnss/data/daily/2024/001/24n/
```

---

#### ③ 北斗广播星历（BDS Broadcast Navigation）

**用途**：实时模式北斗卫星轨道计算

| 项目 | 详情 |
|------|------|
| **数据格式** | RINEX 3.0x 导航文件，包含BDS CNAV星历 |
| **下载地址** | **CDDIS MGEX** — https://cddis.nasa.gov/archive/gnss/data/campaign/mgex/daily/rinex3/YYYY/DDD/YYp/ |
| **文件命名** | `BRD400DLR_S_YYYYDDD0000_01D_MN.rnx`（Multi-GNSS广播星历） |
| **更新频率** | 每日一个文件 |
| **注意事项** | 所有MGEX站点的广播星历都包含GPS+BDS+GLONASS+GALILEO，无需单独下载北斗 |

---

#### ④ GPS精密轨道 (SP3)

**IGS分析中心**：CODE(AIUB)、GFZ、JPL、MIT、ESA、NRCan等

| 产品级别 | 延迟 | 精度(轨道) | 下载地址 | 文件命名示例 |
|---------|------|-----------|---------|-------------|
| **IGU (Ultra-Rapid)** | T+3-9h | ~5cm | `cddis.nasa.gov/gnss/products/` + `WWWW/` | `iguWWWWD_HH.sp3.Z` |
| **IGR (Rapid)** | T+17-41h | ~2.5cm | `cddis.nasa.gov/gnss/products/` + `WWWW/` | `igrWWWWD.sp3.Z` |
| **IGS (Final)** | T+12-18d | ~2.5cm | `cddis.nasa.gov/gnss/products/` + `WWWW/` | `igsWWWWD.sp3.Z` |
| **CODE Final** | T+14d | ~2cm | `ftp.aiub.unibe.ch/CODE/YYYY/` | `COD0OPSFIN_YYYYDDD0000_01D_05M_ORB.SP3.gz` |

**命名规则**：`WWWW`=GPS周, `D`=星期几(0-6), `HH`=UTC时

**CDDIS下载路径**：
```
https://cddis.nasa.gov/archive/gnss/products/
├── WWWW/                    # GPS周年积 → 每种子目录
│   ├── igsWWWWD.sp3.Z      # IGS Final (5min间隔)
│   ├── igrWWWWD.sp3.Z      # IGS Rapid
│   └── iguWWWWD_00.sp3.Z   # IGS Ultra-Rapid (00h/06h/12h/18h UTC)
```

**下载命令示例**：
```bash
# IGS Final SP3 for GPS Week 2296, Day 0
wget https://cddis.nasa.gov/archive/gnss/products/2296/igs22960.sp3.Z
gunzip igs22960.sp3.Z
```

---

#### ⑤ GPS精密钟差 (CLK)

**与SP3配套使用，必须同一分析中心的同一天产品**

| 产品级别 | 延迟 | 精度(钟差) | 下载地址 | 文件命名示例 |
|---------|------|-----------|---------|-------------|
| **IGU** | T+3-9h | ~3ns RMS | 同SP3路径 | `iguWWWWD_HH.clk.Z` |
| **IGR** | T+17-41h | ~0.05ns RMS | 同SP3路径 | `igrWWWWD.clk.Z` |
| **IGS Final** | T+12-18d | ~0.02ns RMS | 同SP3路径 | `igsWWWWD.clk.Z` |
| **CODE Final** | T+14d | ~0.02ns | `ftp.aiub.unibe.ch/CODE/YYYY/` | `COD0OPSFIN_YYYYDDD0000_01D_30S_CLK.CLK.gz` |

**注意**：RINEX CLK格式为30s间隔，包含每颗GPS卫星的精密钟差($10^{-14}$~$10^{-13}$s精度)。

---

#### ⑥ 北斗精密产品 (SP3 + CLK)

**来源**：MGEX多系统分析中心

| 分析中心 | 简称 | 系统 | 精度(轨道) | 延迟 | 下载地址 | 文件命名模式 |
|---------|------|------|-----------|------|---------|-------------|
| **武汉大学** | WUM | GPS+GLO+GAL+BDS | 3-5cm (BDS-3) | T+12-18d | `cddis.nasa.gov/gnss/products/mgex/` | `WUM0MGXFIN_*` |
| **GFZ波茨坦** | GBM | GPS+GLO+GAL+BDS | 2-5cm (BDS-3) | T+12-18d | `cddis.nasa.gov/gnss/products/mgex/` | `GBM0MGXFIN_*` |
| **CODE伯尔尼** | COD | GPS+GLO+GAL+BDS | 2-3cm (BDS-3 MEO) | T+12-14d | `ftp.aiub.unibe.ch/CODE_MGEX/` | `COM0OPSFIN_*` |

**下载示例**：
```bash
# GFZ MGEX Final (包含GPS+BDS)
wget https://cddis.nasa.gov/archive/gnss/products/mgex/YYYY/GBM0MGXFIN_YYYYDDD0000_01D_05M_ORB.SP3.gz
wget https://cddis.nasa.gov/archive/gnss/products/mgex/YYYY/GBM0MGXFIN_YYYYDDD0000_01D_30S_CLK.CLK.gz
```

**实时SSR改正流**（用于实时模式BDS精密轨道）：
```
NTRIP Caster: ntrip.gnsslab.cn:2101  (武汉大学)
Mount Point: SSR00CNE0  (BDS+GPS SSR改正)
```

---

#### ⑦ 天线模型 (ANTEX)

| 项目 | 详情 |
|------|------|
| **数据格式** | ANTEX 1.4 格式（`*.atx`） |
| **最新版本** | igs20.atx（IGS参考框架IGS20） |
| **旧版** | igs14.atx（IGS14框架） |
| **下载地址** | https://files.igs.org/pub/station/general/igs20.atx |
| **备用地址** | ftp://igs.org/pub/station/general/igs20.atx |
| **文件大小** | ~28MB，包含5000+天线校准数据 |
| **内容** | GPS/BDS/GAL/GLO卫星发射天线PCO/PCV + 地面站接收天线PCO/PCV |
| **更新频率** | 约每年一次或新版IGS框架发布时 |
| **BDS可用性** | igs20.atx包含C08-C46等BDS-2/3卫星天线，但BDS GEO卫星(如C01-C05)数据较不完整 |

**关键内容**：
```
TYPE / SERIAL NO 行:
  G01  → GPS PRN 01 卫星天线（Block IIR等）
  C19  → BDS PRN 19 卫星天线 (BDS-3 MEO)
  C08  → BDS PRN 08 卫星天线 (BDS-3 IGSO)

PCO值: 卫星本体坐标系下的相位中心偏移 (mm), 对L1/L2/B1/B3分别给出
PCV值: 天底角NADIR相关的变化 (mm), -90°到+90°范围
```

---

#### ⑧ DCB (差分码偏差)

**GPS DCB**：

| 项目 | 详情 |
|------|------|
| **数据格式** | IONEX DCB格式（`*.DCB`） |
| **来源** | CODE (伯尔尼大学天文研究所) |
| **下载地址** | `ftp://ftp.aiub.unibe.ch/CODE/YYYY/` → `P1P2YYMM.DCB` 或 `P1C1YYMM.DCB` |
| **文件命名** | `P1P22401.DCB`（2024年1月GPS P1-P2 DCB） |
| **更新频率** | 每月一个文件 |
| **精度** | ~0.02ns (σ) |
| **内容** | GPS卫星P1-C1, P1-P2 DCB (纳秒)，以及接收机DCB |

**北斗 DCB (MGEX DCB)**：

| 项目 | 详情 |
|------|------|
| **数据格式** | Bias-SINEX (BSX) 格式 |
| **来源** | 中科院CAS (Wang et al.) |
| **下载地址** | http://pub.ionosphere.cn/product/dcb/ |
| **备选地址** | CAS发布在CDDIS MGEX目录 |
| **更新频率** | 每月 |
| **精度** | ~0.05ns |
| **内容** | BDS-2/3卫星C2I-C6I, C2I-C7I DCB，以及IONO bias (TGD) |

**注意**：星载接收机不需要接收机端DCB（因为只有一个接收机），但必须应用发射卫星端DCB。

---

#### ⑨ ERP (地球定向参数)

| 项目 | 详情 |
|------|------|
| **数据格式** | IERS格式（`*.erp`）或`*.txt` |
| **来源** | IERS (国际地球自转服务) |
| **下载地址** | **IERS C04** — https://hpiers.obspm.fr/iers/eop/eopc04/eopc04.1962-now |
| **备用地址** | ftp://cddis.nasa.gov/products/iers/ |
| **文件命名** | `eopc04_IAU2000.62-now`（持续更新的长文件） |
| **更新频率** | 每日更新（但延时数据约滞后30天） |
| **关键字段** | x_pole (arcsec), y_pole (arcsec), UT1-UTC (s), LOD (s), dX, dY (天极偏移) |
| **精度** | UT1-UTC ~10μs, 极移 ~0.03mas（最终值） |
| **注意** | 快速ERP（Bulletin A）可从IERS获取，但精度较低（0.1-0.2mas） |

**数据格式示例**：
```
YYYY MM DD MJD     x_pole  y_pole  UT1-UTC   LOD     dX_sun  dY_sun
2024 01 01 60310  0.12345 -0.45678  0.012345  0.001234  0.000123 0.000456
```

---

#### ⑩ SINEX (参考框架 / 测站坐标)

| 项目 | 详情 |
|------|------|
| **数据格式** | SINEX 2.0x 格式（`*.snx`） |
| **用途** | 将精密轨道从IGS参考框架转换到WGS84/ITRF时需要 |
| **来源** | IGS Reference Frame Working Group |
| **下载地址** | https://cddis.nasa.gov/archive/gnss/products/ → `WWWW/`目录 |
|   | ftp://igs.org/pub/product/ → `WWWW/`目录 |
| **文件命名** | `igsWWWW.snx` (IGS Final) 或 `igrWWWW.snx` (IGS Rapid) |
| **更新频率** | 每周一个文件 |
| **是否需要** | **SAR卫星定轨场景下通常不需要**（参考框架只有cm级影响），但在与外部分析中心轨道比对时需要 |

---

#### ⑪ 重力场模型 (EIGEN-6S4 / GOCO06s)

| 项目 | 详情 |
|------|------|
| **数据格式** | ICGEM格式（`*.gfc`）或SHA/EGM格式 |
| **用途** | 地球静态重力场球谐系数 $C_{nm}$, $S_{nm}$ |
| **下载地址** | **ICGEM (GFZ波茨坦)** — http://icgem.gfz-potsdam.de/tom_longtime |
| **直达链接** | http://icgem.gfz-potsdam.de/getmodel/gfc/eigen-6s4.gfc |
| **文件大小** | ~25MB (Nmax=240) |
| **阶数** | 220×220（可截断使用，如80×80用于近实时） |
| **时变模型** | C20, C30, C40的时变系数已在ICGEM文件中包含（drift项） |

**备选重力场模型**：
| 模型 | 最高阶 | NRT推荐 | 延时推荐 | 下载 |
|------|--------|---------|---------|------|
| EIGEN-6S4 | 300 | 80×80 | 150×150 | `eigen-6s4.gfc` |
| GOCO06s | 300 | 80×80 | 150×150 | `goco06s.gfc` |
| EGM2008 | 2190 | 100×100 | 200×200 | `egm2008.gfc` |

**ICGEM下载命令**：
```bash
wget "http://icgem.gfz-potsdam.de/getmodel/gfc/eigen-6s4.gfc" -O eigen_6s4.gfc
```

---

#### ⑫ 海潮模型 (FES2014b)

| 项目 | 详情 |
|------|------|
| **数据格式** | NetCDF网格（全球海潮振幅+相位） |
| **用途** | 计算海潮负荷对低轨卫星轨道的影响 |
| **下载地址** | **AVISO** — https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes.html |
| **注册要求** | 需要注册AVISO账号（免费） |
| **文件大小** | ~2GB |
| **包含潮汐** | M2, S2, N2, K2, K1, O1, P1, Q1等11个主潮波 |
| **替代模型** | FES2004（公开无需注册，精度略低，可用于NRT模式） |
|   | TPXO9 (Oregon State) — https://www.tpxo.net/global |
| **注意** | 海潮模型是近实时/延时模式的可选项；忽略它可能引入~1-2cm误差 |

**简化方案**：如果下载FES2014b全网格太繁琐，可以在延时模式下启用（高精度），近实时模式下用FES2004或关闭。

---

#### ⑬ 行星历表 (DE440 / DE441)

| 项目 | 详情 |
|------|------|
| **数据格式** | JPL BSP (Binary SPK) 二进制格式 |
| **用途** | 计算月球、太阳及大行星精确位置（用于N体摄动计算） |
| **下载地址** | **JPL SSD** — https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/ |
| **直达链接** | https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de440.bsp |
| **文件大小** | ~100MB |
| **覆盖时间** | 1549-2650年 |
| **读取方式** | 使用`jplephem`(Python)或Orekit内置SPICE reader |
| **精度** | 月球位置 ~0.1mas, 太阳位置 ~0.01mas |
| **注意** | 对SAR卫星5cm定轨，日月的DE440历表是必需的（N体摄动贡献~1-2cm） |

---

#### ⑭ 卫星宏模型文件 (Satellite Macro Model)

##### 概念

卫星宏模型（Macro Model）是精密定轨非保守力建模的核心输入之一。它不是可下载的公开数据，而是**根据目标卫星的几何设计、表面材料和姿态特性自行构建的物理模型**，用于精确计算大气阻力和太阳辐射压（SRP）产生的非保守力加速度。

在LEO轨道上（400-1400 km），大气阻力和SRP是除地球重力外最重要的摄动力：
- **大气阻力**：对500km轨道SAR卫星，可产生 ~10⁻⁷ m/s² 量级的持续摄动，24h累积位置误差可达 **百米级**
- **太阳辐射压**：对800kg级卫星，可产生 ~10⁻⁸~10⁻⁷ m/s² 量级的加速度，需在cm级定轨中建模

##### Box-Wing 模型

SAR卫星受限于成本和运力，通常不具备加速度计（ACC），因此必须采用**Box-Wing几何光学模型**：将卫星简化为多个平板面的组合（"Box"为星本体六个面，"Wing"为太阳翼和SAR天线展开面），每个面用面积、法向量、光学系数描述。

```
物理原理：

1. 大气阻力加速度：
   a_drag = -0.5 * Cd * (ρ * |v_rel|² * A_eff) / m * v̂_rel

   其中 A_eff = Σᵢ Aᵢ * max(cos θᵢ, 0)
   - Cd：阻力系数（待估参数，典型值 2.0~2.5）
   - ρ：大气密度（来自 NRLMSISE-00 / DTM-2013 模型）
   - v_rel：卫星相对大气的速度矢量
   - θᵢ：面板 i 法向量与来流方向的夹角
   - A_eff：有效迎风面积（各面板投影面积之和）

2. 太阳辐射压加速度：
   a_srp = Cr * (S₀ / c) * (A_eff / m) * r̂_sun * ν_shadow

   各面板贡献分解：
   a_srp_panel = -(S₀/c) * A * cos(θ) * [(1 - ρ_s) * r̂_sun + 2*(ρ_s/3 + ρ_d)*n̂]

   其中：
   - Cr：光压系数（待估参数，典型值 1.0~1.5）
   - S₀ = 1361 W/m²（太阳常数）
   - c = 3e8 m/s（光速）
   - ν_shadow：地影因子（0=全影，1=全光照；Cannonball/Conical模型）
   - ρ_s：镜面反射系数 (specular / reflectivity)
   - ρ_d：漫反射系数 (diffuse)
   - ρ_a = 1 - ρ_s - ρ_d：吸收系数（absorbed energy 转化为热辐射，贡献各向同性）
   - n̂：面板法向量
   - r̂_sun：太阳方向单位矢量
   - θ：面板法向量与太阳方向的夹角

3. 能量守恒约束：ρ_s + ρ_d + ρ_a = 1.0
   → 对每个面板必须保证 reflectivity + diffuse ≤ 1.0
```

##### 面板属性含义

| 属性 | 符号 | 物理含义 | 对阻力的影响 | 对光压的影响 |
|------|------|---------|:-----------:|:-----------:|
| `area_m2` | A | 面板几何面积 | 决定迎风截面积 | 决定受光截面积 |
| `normal_vector` | n̂ | 面板法向量（星体系） | 决定来流方向投影 cosθ | 决定太阳方向投影 cosθ |
| `reflectivity_coeff` | ρ_s | 镜面反射系数 | — | 镜面反射分量（2倍冲量） |
| `diffuse_coeff` | ρ_d | 漫反射系数 | — | 漫反射分量（2/3倍冲量） |
| `absorption_coeff` | ρ_a | 吸收系数 (=1-ρ_s-ρ_d) | — | 吸收后热辐射（1倍冲量） |

##### 星体坐标系 (S/C Body Frame)

卫星宏模型中所有面板法向量定义在 **卫星本体坐标系** 下：

```
星体坐标系定义 (以对地定向SAR卫星为例)：

     +Z (天顶方向，背离地心)
      ↑
      |
      +----→ +X (沿飞行方向，沿轨道切线)
     /
   +Y (轨道面法线，右手定则)

面板命名约定：
  +X 面 → 沿飞行方向前侧      (RAM方向，大气阻力主导面)
  -X 面 → 沿飞行方向后侧      (Wake方向)
  +Y 面 → 轨道面正法向        (SAR通常侧视方向)
  -Y 面 → 轨道面负法向
  +Z 面 → 天顶方向            (SRP主导面，受照时间长)
  -Z 面 → 天底/向地方向       (SAR天线展开方向)
  +Solar Wing → 太阳翼 +Y面  (大面积，SRP主力)
  -Solar Wing → 太阳翼 -Y面
  SAR Antenna → 大型SAR天线  (大面积，对阻力和光压都敏感)
```

##### 姿态依赖的有效面积计算

由于SAR卫星在轨姿态变化（侧视成像、太阳翼跟踪），宏模型必须与 **姿态四元数** 联用，实时计算有效截面积：

```python
def compute_effective_area(panels, attitude_quaternion, direction_eci):
    """
    根据卫星当前姿态计算对某方向的有效截面积
    
    Args:
        panels: 面板列表 [{area, normal_vector_sbf, ...}, ...]
        attitude_quaternion: 星体→ECI 四元数 [q0, q1, q2, q3]
        direction_eci: 来流/光照方向单位矢量 (ECI系)
                       对阻力：v_rel方向；对光压：r_sun方向
    
    Returns:
        A_eff: 有效截面积 (m²)
    """
    R_sbf_to_eci = quaternion_to_rotation_matrix(attitude_quaternion)
    A_eff = 0.0
    for panel in panels:
        # 面板法向量从星体系转到ECI
        n_eci = R_sbf_to_eci @ panel.normal_vector_sbf
        # 投影面积 = 面积 × cos(入射角)，仅计算正投影
        cos_theta = max(np.dot(n_eci, direction_eci), 0.0)
        A_eff += panel.area * cos_theta
    return A_eff
```

**关键设计原则**：
- 阻力计算使用 `direction = v̂_rel`（相对大气的来流方向）
- 光压计算使用 `direction = -r̂_sun`（太阳光入射反向）
- SAR天线是大型展开结构（典型面积 10~30 m²），其投影面积**随侧视角变化显著**，不能用常值面积近似

##### SAR卫星的特殊考虑

相比普通LEO卫星，SAR卫星的宏模型构建面临特殊挑战：

| 挑战 | 原因 | 建模策略 |
|------|------|---------|
| **大SAR天线** | 展开天线面积可达15-30m²，与星本体相当 | 单独作为面板建模，随姿态精确旋转 |
| **太阳翼跟踪** | 太阳翼法向量相对星体变化 | 动态计算或按模式分段建模 |
| **侧视成像姿态机动** | 侧视角可达20-55°，改变阻力/SRP剖面 | 将SAR天线面设为独立面板，由姿态驱动 |
| **缺少ACC数据** | 无加速度计直接测量非保守力 | Box-Wing模型是唯一手段，参数须可估 |
| **多面板几何复杂** | 结构不规则（馈源、展开机构等） | 简化等效面板，面积守恒，光学系数等效 |

##### 完整JSON示例

```json
{
  "satellite_name": "SAR_SAT_01",
  "description": "C-band SAR satellite, sun-synchronous dawn-dusk orbit, ~500km altitude",
  "mass_kg": 1200.0,
  "mass_sigma_kg": 5.0,
  "center_of_mass_offset_sbf": [0.05, 0.0, -0.15],
  "panels": [
    {
      "name": "+X (RAM方向 / 沿飞行前侧)",
      "area_m2": 6.25,
      "reflectivity_coeff": 0.30,
      "diffuse_coeff": 0.40,
      "normal_vector_sbf": [1.0, 0.0, 0.0],
      "note": "阻力主导面，SRP贡献在晨昏轨道较小"
    },
    {
      "name": "-X (Wake方向 / 沿飞行后侧)",
      "area_m2": 6.25,
      "reflectivity_coeff": 0.25,
      "diffuse_coeff": 0.45,
      "normal_vector_sbf": [-1.0, 0.0, 0.0],
      "note": "尾流方向，阻力贡献为零"
    },
    {
      "name": "+Y (轨道面正法向 / 右侧视方向)",
      "area_m2": 8.00,
      "reflectivity_coeff": 0.20,
      "diffuse_coeff": 0.50,
      "normal_vector_sbf": [0.0, 1.0, 0.0],
      "note": "SAR侧视方向，成像时朝向地面"
    },
    {
      "name": "-Y (轨道面负法向 / 左侧视方向)",
      "area_m2": 8.00,
      "reflectivity_coeff": 0.20,
      "diffuse_coeff": 0.50,
      "normal_vector_sbf": [0.0, -1.0, 0.0]
    },
    {
      "name": "+Z (天顶方向)",
      "area_m2": 8.00,
      "reflectivity_coeff": 0.15,
      "diffuse_coeff": 0.55,
      "normal_vector_sbf": [0.0, 0.0, 1.0],
      "note": "SRP主导面，晨昏轨道中持续受照"
    },
    {
      "name": "-Z (天底/向地方向 / SAR天线安装面)",
      "area_m2": 8.00,
      "reflectivity_coeff": 0.20,
      "diffuse_coeff": 0.60,
      "normal_vector_sbf": [0.0, 0.0, -1.0],
      "note": "安装SAR天线的星本体底部"
    },
    {
      "name": "Solar Wing +Y (太阳翼受照面)",
      "area_m2": 12.00,
      "reflectivity_coeff": 0.05,
      "diffuse_coeff": 0.15,
      "normal_vector_sbf": [0.0, 1.0, 0.0],
      "note": "太阳能电池板，吸收率高(~0.8)；法向量随太阳跟踪变化"
    },
    {
      "name": "Solar Wing -Y (太阳翼背面)",
      "area_m2": 12.00,
      "reflectivity_coeff": 0.10,
      "diffuse_coeff": 0.30,
      "normal_vector_sbf": [0.0, -1.0, 0.0],
      "note": "背面热控涂层，光学特性不同"
    },
    {
      "name": "SAR Antenna (合成孔径雷达天线)",
      "area_m2": 25.00,
      "reflectivity_coeff": 0.08,
      "diffuse_coeff": 0.65,
      "normal_vector_sbf": [0.0, 0.0, -1.0],
      "note": "大型展开天线，面积最大，对阻力和SRP均敏感；法向量天底方向"
    }
  ],
  "summary": {
    "total_panel_area_m2": 93.50,
    "mass_area_ratio_kg_per_m2": 12.83,
    "reference_drag_area_m2": 15.0,
    "reference_srp_area_m2": 35.0
  }
}
```

##### 构建流程与验证

宏模型的构建需要卫星研制方的配合，典型流程：

```
卫星CAD模型 → 提取外表面几何 → 简化为平板面集合
    ↓
各面板表面材料 → 实验室测量光学系数 (ρ_s, ρ_d)
    ↓
组装JSON宏模型 → 与在轨遥测比对
    ↓
初始定轨验证 → 调整Cd/Cr先验值 → 固化宏模型
```

如果没有详细CAD模型，可用的简化方案：
1. **单面板Isotropic模型**：只用总质量和常值面积质量比（A/m），对阻力用球形等效，对SRP用各向同性——精度最低
2. **六面Box模型**：星本体用六个标准面等效，太阳翼和天线分别附加——NRT模式推荐
3. **完整Box-Wing模型**：8-12个面板，包含天线细节——延时模式精度最优

---

### 1.2.2 各处理模式数据需求汇总

| 数据 | 实时 | 近实时 (NRT) | 延时 (Delayed) |
|------|:--:|:--:|:--:|
| 星载观测 (RINEX OBS) | ✓ | ✓ | ✓ |
| 姿态四元数 | ✓ | ✓ | ✓ |
| GPS广播星历 | ✓ | ✓(初值) | ✓(初值) |
| BDS广播星历 | ✓(如果用) | ✓(初值) | ✓(初值) |
| **GPS SP3 精密轨道** | ✗ (或用SSR流) | ✓ (IGU/IGR) | ✓ (IGS Final) |
| **BDS SP3 精密轨道** | ✗ | ✓ (GBM/WUM) | ✓ (GBM/WUM/COD) |
| **GPS CLK 精密钟差** | ✗ | ✓ | ✓ |
| **BDS CLK 精密钟差** | ✗ | ✓ | ✓ |
| **ANTEX 天线模型** | ✗ | ✓ | ✓ |
| **DCB (GPS)** | ✗ | ✓ | ✓ |
| **DCB (BDS)** | ✗ | ✓ | ✓ |
| **ERP 地球定向参数** | ✓(预报) | ✓ | ✓ |
| **SINEX** | ✗ | ✗ | ✓(比对用) |
| **重力场模型 (EIGEN-6S4)** | ✗ | ✓ (80×80) | ✓ (150×150) |
| **海潮模型 (FES2014b)** | ✗ | ✓(可选) | ✓ |
| **DE440 行星历表** | ✗ | ✓ | ✓ |
| **卫星宏模型** | ✗ | ✓ | ✓ |

### 1.3 输出数据

| 输出 | 格式 | 说明 |
|------|------|------|
| 轨道状态向量 | SP3 格式 | 15min间隔，含精度标识 |
| 协方差矩阵 | 自定义二进制/JSON | 每个历元6×6协方差 |
| 残差文件 | CSV | 观测值-计算值残差，用于QA |
| 处理日志 | 文本 | 迭代次数、收敛状态、告警 |
| 精度评估报告 | JSON | 3D RMS、重叠弧段检验、Sigma0 |

---

## 2. 系统总体架构

### 2.1 三层架构

```
┌─────────────────────────────────────────────────────────┐
│                    配置与调度层                          │
│  config.yaml / config.json → 解析 → 调度引擎            │
└──────────────────────┬──────────────────────────────────┘
                         │
┌─────────────────────────────────────────────────────────┐
│                    核心定轨引擎层                          │
│  ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐│
│  │观测预处理  │ │力学模型   │ │数值积分   │ │参数估计   ││
│  │ObservPrep │ │Dynamics  │ │Integrator│ │Estimator ││
│  └──────────┘ └──────────┘ └──────────┘ └──────────┘│
│  ┌──────────┐ ┌──────────┐ ┌──────────┐              │
│  │模糊度固定 │ │质量评估   │ │Orekit适配 │              │
│  │AR_Module │ │QAModule  │ │OrekitAdapt│              │
│  └──────────┘ └──────────┘ └──────────┘              │
└──────────────────────┬──────────────────────────────────┘
                         │
┌─────────────────────────────────────────────────────────┐
│                    数据接入与输出层                        │
│  数据加载器 ← 原始文件    输出写入器 → SP3/CSV/JSON      │
└─────────────────────────────────────────────────────────┘
```

### 2.2 模块依赖关系

```
DataLoader → ObservPrep → Dynamics → Integrator → Estimator → OutputWriter
                ↓             ↓           ↓           ↓
             AR_Module ← QAModule (残差分析驱动AR决策)
             OrekitAdapt (可选后端，替代Dynamics+Integrator)
```

---

## 3. 三种处理模式详细设计

### 3.1 实时模式（Real-Time）

**场景**：星上实时导航或地面快速定位  
**数据窗口**：滑动窗口 5-10 分钟  
**处理时延**：< 1s（星上）/ < 10s（地面）  
**动力学模型**：简化动力学（计算量极小）

| 组件 | 配置 |
|------|------|
| 轨道模型 | 二体问题 + J2 摄动 |
| 经验力 | 无（或1个常数加速度参数） |
| 观测值 | 伪距 only（不用法相位） |
| 钟差处理 | 接收机钟差作为白噪声估计 |
| 积分器 | 解析解（二体）或简单RK4 |
| 估计方法 | 卡尔曼滤波（扩展EKF） |

**算法选择**：扩展卡尔曼滤波（EKF），单历元更新，无平滑。

### 3.2 近实时模式（Near-Real-Time）

**场景**：地面接收2小时数据，快速产出轨道  
**数据窗口**：2小时连续弧段  
**处理时延**：≤ 2分钟  
**动力学模型**：Reduced-Dynamic（中等复杂度）

| 组件 | 配置 |
|------|------|
| 轨道模型 | EIGEN-6S4 50×50（或简化80×80） |
| 经验力 | RTN 15min分段常值，3参数/段 |
| 观测值 | 伪距 + 相位（浮点模糊度） |
| 钟差处理 | 接收机钟差随机游走，GPS钟差用精密CLK |
| 积分器 | Gauss-Jackson 8阶（快速） |
| 估计方法 | 最小二乘批处理（LSQ），2-3次迭代 |

**关键优化**：
- 降采样到30s（原始10s → 30s，减少2/3观测）
- 截止高度角10°，减少低质量观测
- 并行：2小时弧段分4段（每段30min）并行处理，后拼接

### 3.3 延时模式（Delayed）

**场景**：获取完整数据后高精度定轨  
**数据窗口**：24-48小时弧段  
**处理时延**：≤ 5分钟  
**动力学模型**：完整 Reduced-Dynamic

| 组件 | 配置 |
|------|------|
| 轨道模型 | EIGEN-6S4 150×150 + 时变C20-C40 |
| 经验力 | RTN 6min分段常值，强过程噪声约束 |
| 观测值 | 伪距 + 相位 + PPP-AR模糊度固定 |
| 钟差处理 | 接收机钟差随机游走，GPS/BDS钟差精密CLK |
| 积分器 | Gauss-Jackson 8阶 或 DOPRI8（可选） |
| 估计方法 | LSQ批处理 + LAMBDA模糊度固定，3-5次迭代 |

**关键优化**：
- 完整观测（10s采样，不降采样）
- 截止高度角5°（高精度要求）
- PPP-AR：先浮点解，再LAMBDA固定，Ratio>3.0接受
- 24h弧段并行：分12段（每段2h）并行，后重叠弧段拼接

---

## 4. 核心算法模块设计

### 4.1 模块1：数据接入与对齐架构（DataInputArchitecture）

> **设计目标**：适配各类星载导航接收机数据格式，解析与输入接口完全解耦，明确数据对齐策略和时间基准。
> **核心原则**：三层解耦（Parser → 统一模型 → InputManager），以观测数据GPST时间为对齐基准。

---

#### 4.1.1 总体架构：三层解耦设计

```
┌─────────────────────────────────────────────────────────────────────┐
│  Layer 3：数据输入管理器（DataInputManager）                     │  ← POD算法唯一调用入口
│  职责：数据对齐、时间插值、统一查询接口                           │
│  关键方法：load_and_align() → AlignedDataBuffer                 │
├─────────────────────────────────────────────────────────────────────┤
│  Layer 2：统一数据模型（Unified Data Model）                      │  ← 所有数据的归一化结构
│  职责：定义 ObsEpoch, Sp3Orbit, ClkData, AttitudeQuat 等        │
│  关键点：所有时间标签统一为 GPST（GPS Time，连续秒）              │
├─────────────────────────────────────────────────────────────────────┤
│  Layer 1：格式解析器（FormatParser）——适配器模式                  │  ← 格式相关，可独立扩展
│  职责：解析特定格式原始文件 → 归一化数据对象（Layer 2 结构）       │
│  扩展方式：新增 Parser 类，实现 FormatParser 接口，注册到 Factory  │
│  已有实现：RinexObsParser, GraceFOObsParser, Sp3Parser...      │
└─────────────────────────────────────────────────────────────────────┘
          ↑ Raw Files (RINEX, GRACE-FO L1B, SP3, CLK, SCA1B...)
```

**解耦要点**：
- Layer 1 只负责"格式解析"，不关心POD算法需要什么
- Layer 2 只定义"数据长什么样"，不关心数据从哪来、怎么对齐
- Layer 3 只负责"数据对齐和查询"，不关心具体文件格式

---

#### 4.1.2 时间基准与对齐策略（核心设计）

##### 时间基准（Time Reference Benchmark）

**GPST（GPS Time）** 是系统唯一时间基准。

| 属性 | 说明 |
|------|------|
| **定义** | 1980-01-06 00:00:00 GPST 起算的连续秒数（不含闰秒） |
| **与UTC关系** | GPST = UTC + 闰秒偏移（GPS以来累计18闰秒，2026年=18s） |
| **系统内使用** | 所有Layer 2数据对象的时间标签、POD状态向量时间、积分时间均用GPST |

**各数据源时间标签 → GPST 转换规则**：

| 数据源 | 原始时间标签 | 转换方法 |
|---------|-------------|---------|
| RINEX OBS | GPS Week + GPS Seconds of Week | 直接：`t_gpst = week*604800 + sow` |
| RINEX NAV | GPST（Toe） | 直接 |
| GRACE-FO L1B | GPS Seconds（GPST） | 直接 |
| SP3 | GPST（文件头MJD→GPST） | `t_gpst = (mjd - 44244.0) * 86400` |
| CLK | GPST（年积日+秒-of-day） | 直接 |
| ERP | UTC 0h | `t_gpst = utc_to_gpst(utc)`（需闰秒表）|
| SCA1B姿态 | GPS Seconds（GPST） | 直接 |
| ANTEX | 无时间标签（天线校准常数） | N/A |

---

##### 主时间轴（Primary Time Grid）定义

**主时间轴 = 星载接收机观测历元集合**，以GPST表示。

**选择理由**：
1. 观测数据是POD的"驱动数据"——残差计算和设计矩阵评估都以观测历元为基准
2. 观测历元数量最多（每小时360个@10s采样），精度要求最高
3. 辅助数据（SP3/CLK/ERP等）都是"背景产品"，插值到观测历元是最自然的处理方式

**主时间轴生成流程**：
```
Step 1: 解析观测文件（Layer 1: RinexObsParser / GraceFOObsParser / ...）
Step 2: 提取所有观测历元时间标签 → 转为GPST
Step 3: 排序、去重 → 得到主时间轴 T_primary = [t₁, t₂, ..., t_N] (GPST)
Step 4: T_primary 作为后续所有对齐操作的基准
```

---

##### 数据对齐总流程

```
输入：各数据源原始文件（格式各异，时间标签各异）
  │
  ▼
[Layer 1: 格式解析器]
  每个Parser读取对应格式 → 输出Layer 2归一化对象（时间标签已转GPST）
  │
  ▼
[Layer 3: DataInputManager.load_and_align()]
  │
  ├─ Step A: 确定主时间轴 T_primary（来自观测数据ObsEpoch集合）
  │
  ├─ Step B: 对每个 t ∈ T_primary，插值所有辅助数据源：
  │    ├─ SP3精密轨道 → 插值到 t（Lagrange 7点）
  │    ├─ CLK精密钟差 → 插值到 t（线性插值）
  │    ├─ ERP地球定向 → 插值到 t（线性插值）
  │    ├─ 姿态四元数   → 插值到 t（Slerp球面线性插值）
  │    ├─ DCB差分偏差 → 查表/线性插值到 t
  │    └─ ANTEX天线    → 查表（不插值，按卫星PRN查）
  │
  ├─ Step C: 构建 AlignedDataBuffer（按T_primary索引的对齐数据缓冲区）
  │
  ▼
输出：AlignedDataBuffer（POD算法直接查询使用）
```

---

#### 4.1.3 Layer 1：格式解析器设计（适配器模式）

**核心接口**：

```python
# src/data/format_parser.py
from abc import ABC, abstractmethod
from typing import Iterable, Dict, Any, List, Type
from pathlib import Path

class FormatParser(ABC):
    """格式解析器抽象基类。"""
    @abstractmethod
    def parse(self, file_path: Path, **kwargs) -> Iterable[Any]:
        """解析文件，返回归一化数据对象迭代器。"""
        ...
    @property
    @abstractmethod
    def supported_formats(self) -> List[str]:
        """返回此Parser支持的文件格式标识列表"""
        ...

class ParserRegistry:
    """Parser注册表（工厂模式）。"""
    _registry: Dict[str, Type[FormatParser]] = {}
    @classmethod
    def register(cls, format_id: str, parser_cls: Type[FormatParser]):
        cls._registry[format_id] = parser_cls
    @classmethod
    def get_parser(cls, format_id: str) -> FormatParser:
        return cls._registry[format_id]()
```

**具体Parser实现（已规划）**：

| Parser类 | 支持的格式 | 输出Layer 2对象 | 说明 |
|---------|------------|----------------|------|
| `RinexObsParser` | RINEX v2.11 / v3.0x OBS | `ObsEpoch` | 通用格式，最常用 |
| `GraceFOObsParser` | GRACE-FO L1B ASCII (GNV1B) | `ObsEpoch` | NASA JPL格式，验证用 |
| `SentinelObsParser` | Sentinel-1 GNSS RINEX | `ObsEpoch` | ESA哥白尼计划 |
| `TerraSarObsParser` | TerraSAR-X GNSS DAS | `ObsEpoch` | DLR/ASL格式 |
| `FengYunObsParser` | 风云三号GNSS二进制 | `ObsEpoch` | 中国气象卫星 |
| `RinexNavParser` | RINEX v2/v3 NAV (GPS/BDS/GLO/GAL) | `NavData` | 广播星历 |
| `Sp3Parser` | SP3-c / SP3-d | `Sp3Orbit` | 精密轨道 |
| `ClkParser` | RINEX CLK (30s/5s) | `ClkData` | 精密钟差 |
| `AntexParser` | ANTEX 1.4 | `AntexData` | 天线相位中心 |
| `ErpParser` | IERS ERP (eopc04 / bulletin_a) | `ErpData` | 地球定向参数 |
| `AttitudeParser` | SCA1B / 通用四元数文本 | `AttitudeQuat` | 姿态四元数 |
| `DcbParser` | CODE .DCB / CAS .BSX | `DcbData` | 差分码偏差 |

---

#### 4.1.4 Layer 2：统一数据模型（归一化结构）

> 所有Layer 1 Parser的输出、Layer 3的输入，均使用以下归一化数据类。
> 时间标签统一为GPST（`t_gpst: float`，单位：秒）。

```python
# src/data/unified_model.py
from dataclasses import dataclass, field
from typing import Dict, List, Optional
import numpy as np

GPST_EPOCH_AS_MJD = 44244.0

@dataclass
class ObsEpoch:
    """一个历元的星载GNSS观测数据（归一化）。时间基准：GPST。"""
    t_gpst: float
    recv_time_internal: float
    satellites: Dict[str, 'SatObservation']
    recv_clock_offset: float = 0.0
    integration_time: float = 0.0

@dataclass
class SatObservation:
    sat_id: str
    sys: str
    L1: Optional[float] = None
    L2: Optional[float] = None
    C1: Optional[float] = None
    C2: Optional[float] = None
    D1: Optional[float] = None
    D2: Optional[float] = None
    B1: Optional[float] = None
    B2: Optional[float] = None
    B3: Optional[float] = None
    LLI: Dict[str, int] = field(default_factory=dict)
    sig_strength: Dict[str, float] = field(default_factory=dict)
    elevation_deg: Optional[float] = None
    azimuth_deg: Optional[float] = None

@dataclass
class NavData:
    """GPS/北斗广播星历（归一化）。时间基准：GPST。"""
    sat_id: str; sys: str
    t_oc: float; af0: float; af1: float; af2: float
    t_oe: float; sqrt_a: float; e: float; i0: float
    Omega0: float; omega: float; M0: float
    delta_n: float; idot: float; Omega_dot: float
    tgd1: Optional[float] = None; tgd2: Optional[float] = None
    health: int = 0

@dataclass
class Sp3OrbitRecord:
    t_gpst: float; sat_id: str
    x: float; y: float; z: float; clock: float
    x_sdev: float = 0.0; y_sdev: float = 0.0
    z_sdev: float = 0.0; clk_sdev: float = 0.0

@dataclass
class Sp3Orbit:
    source: str; version: str
    records: Dict[str, List[Sp3OrbitRecord]]
    interval: float = 900.0

@dataclass
class ClkRecord:
    t_gpst: float; sat_id: str
    clock_bias: float; clock_drift: float = 0.0; sigma: float = 0.0

@dataclass
class ClkData:
    source: str; records: Dict[str, List[ClkRecord]]; interval: float = 30.0

@dataclass
class AttitudeQuat:
    t_gpst: float; q0: float; q1: float; q2: float; q3: float
    sigma_deg: float = 0.0

@dataclass
class ErpRecord:
    mjd_utc: float; x_pole: float; y_pole: float
    ut1_utc: float; lod: float = 0.0
    dX: float = 0.0; dY: float = 0.0

@dataclass
class ErpData:
    source: str; records: List[ErpRecord]

@dataclass
class AntennaPco:
    sys: str; freq: str; north_mm: float; east_mm: float; up_mm: float

@dataclass
class AntennaPcv:
    sys: str; freq: str
    nadir_angles_deg: List[float]; pcv_mm: List[float]

@dataclass
class AntexData:
    version: str
    sat_pco: Dict[str, Dict[str, AntennaPco]]
    sat_pcv: Dict[str, Dict[str, AntennaPcv]]

@dataclass
class DcbRecord:
    t_start_mjd: float; t_end_mjd: float
    sat_id: str; type: str; bias: float; sigma: float = 0.0

@dataclass
class DcbData:
    source: str; records: List[DcbRecord]
```

---

#### 4.1.5 Layer 3：数据输入管理器（DataInputManager）

**核心职责**：调用Layer 1解析 → 确定主时间轴 → 插值对齐 → 输出AlignedDataBuffer

```python
# src/data/data_input_manager.py
class DataInputManager:
    def __init__(self, config: Dict):
        self.config = config
        self.aligned_buffer: Optional['AlignedDataBuffer'] = None

    def load_and_align(self, file_specs: List) -> 'AlignedDataBuffer':
        parsed = self._parse_all(file_specs)
        t_grid = self._build_time_grid(parsed['observations'])
        aligned = self._align_all(t_grid, parsed)
        self.aligned_buffer = AlignedDataBuffer(aligned, t_grid)
        return self.aligned_buffer

    def _parse_all(self, specs) -> Dict:
        result = {'observations':[], 'navigation':[], 'sp3_orbits':[],
                  'clk_clocks':[], 'attitude':None, 'erp':None,
                  'antex':None, 'dcb':None}
        for s in specs:
            p = ParserRegistry.get_parser(s.format_id)
            objs = list(p.parse(s.path))
            # 按format_id路由到对应字段（略）
        return result

    def _build_time_grid(self, obs_list) -> np.ndarray:
        return np.sort(np.array(list(set(e.t_gpst for e in obs_list))))

    def _align_all(self, t_grid, data) -> List:
        aligned = []
        for t in t_grid:
            so = {sid: self._interp_sp3(data['sp3_orbits'], sid, t)
                   for sp3 in data['sp3_orbits'] for sid in sp3.records}
            sc = {sid: self._interp_clk(data['clk_clocks'], sid, t)
                   for cl in data['clk_clocks'] for sid in cl.records}
            erp = self._interp_erp(data['erp'], t)
            att = self._interp_att(data['attitude'], t)
            dcb = self._lookup_dcb(data['dcb'], t)
            aligned.append(AlignedDataAtEpoch(
                t_gpst=t, observations={}, sat_orbits=so,
                sat_clocks=sc, attitude=att, erp=erp, dcbs=dcb))
        return aligned

    def _interp_sp3(self, sp3_list, sid, t):
        recs = sp3_list[0].records[sid]  # 简化：第一个SP3文件
        idx = np.searchsorted([r.t_gpst for r in recs], t)
        # Lagrange 7-point interpolation (略）
        return (0.0, 0.0, 0.0, 0.0)  # placeholder

    def _interp_clk(self, cl_list, sid, t):
        recs = cl_list[0].records[sid]
        idx = np.searchsorted([r.t_gpst for r in recs], t)
        r0, r1 = recs[idx-1], recs[idx]
        s = (t - r0.t_gpst) / (r1.t_gpst - r0.t_gpst)
        return r0.clock_bias * (1-s) + r1.clock_bias * s

    def _interp_erp(self, erp, t):
        t_mjd = erp.records[0].mjd_utc  # 简化
        return None  # placeholder

    def _interp_att(self, att_list, t):
        idx = np.searchsorted([a.t_gpst for a in att_list], t)
        return self._slerp(att_list[idx-1], att_list[idx], t)

    def _slerp(self, q0, q1, t):
        # Slerp implementation (略）
        return q0  # placeholder

    def _lookup_dcb(self, dcb, t):
        return {}  # placeholder
```

---

#### 4.1.6 插值方法详述

| 数据类型 | 插值方法 | 说明 |
|---------|---------|------|
| SP3轨道 | Lagrange 7点（SP3-d）或Hermite（SP3-c） | 高精度要求 |
| CLK钟差 | 线性插值 | 30s间隔，线性足够 |
| ERP参数 | 线性插值 | 天间隔，线性足够 |
| 姿态四元数 | **Slerp球面线性**（禁止线性！）| 四元数单位长度约束 |
| DCB偏差 | 分段常数/线性 | 天际常数，天内线性 |
| ANTEX天线 | 查表（不插值）| 按卫星PRN直接查 |

**SP3 Lagrange 7点插值**：
```
x(t) = Σᵢ xᵢ · Lᵢ(t),  Lᵢ(t) = Πₖ≠ᵢ (t-tₖ)/(tᵢ-tₖ)
```

**Slerp四元数插值**：
```
q(t) = [sin((1-s)Ω)·q₁ + sin(sΩ)·q₂] / sin(Ω)
s = (t-t₁)/(t₂-t₁),  Ω = acos(q₁·q₂)
```

---

#### 4.1.7 AlignedDataBuffer设计

```python
@dataclass
class AlignedDataAtEpoch:
    t_gpst: float
    observations: Dict[str, SatObservation]
    sat_orbits: Dict[str, 'SatOrbitAtEpoch']
    sat_clocks: Dict[str, float]
    attitude: Optional
    erp: Optional
    dcbs: Dict[str, float]

@dataclass
class SatOrbitAtEpoch:
    sat_id: str; t_gpst: float
    x_km: float; y_km: float; z_km: float
    clk_us: float
    x_dot_kms: float = 0.0; y_dot_kms: float = 0.0; z_dot_kms: float = 0.0

class AlignedDataBuffer:
    def __init__(self, epochs, t_grid):
        self.epochs = epochs; self.t_grid = t_grid
        self._epoch_dict = {e.t_gpst: e for e in epochs}

    def get_at(self, t_gpst: float) -> 'AlignedDataAtEpoch':
        return self._epoch_dict[t_gpst]

    def iterate_epochs(self):
        return iter(self.epochs)
```

---

#### 4.1.8 配置文件中的数据路径规格

```yaml
data_sources:
  observations:
    - format_id: "gracefo_obs"
      path: "/data/gracefo/L1B/GNV1B_2024001.txt"
      role: "primary"
  navigation:
    - format_id: "rinex_nav"
      path: "/data/gnss/brdc0010.24n"
      systems: ["G", "C"]
  precise_orbits:
    - format_id: "sp3"
      path: "/data/igs/igs22960.sp3"
      source: "IGS"
  precise_clocks:
    - format_id: "clk"
      path: "/data/igs/igs22960.clk"
      source: "IGS"
  attitude:
    - format_id: "sca1b"
      path: "/data/gracefo/SCA1B_2024001.txt"
  erp:
    - format_id: "iers_erp"
      path: "/data/iers/eopc04_2024.txt"
  antex:
    - format_id: "antex"
      path: "/data/igs/igs20.atx"
  dcb:
    - format_id: "code_dcb"
      path: "/data/code/P1P22401.DCB"
```

---

### 4.2 模块2：观测预处理器（ObservPrep）

**职责**：
1. 周跳检测与标记（GF组合 + MW组合）
2. 观测值质量控制（残差3σ剔除）
3. 形成观测方程（伪距、相位分别对应不同权重）
4. 天线PCO/PCV校正（利用姿态数据旋转PCO到SAT固定系）

**核心算法**：

**周跳检测**：
```
GF = L1 - L2 (cycles) → 平滑检验，跳变>0.05 cycles → 周跳
MW = (L1 - L2) - (C1 - C2)/λ_iono → 整数跳变 → 周跳
```

**天线PCO校正**（利用姿态四元数）：
```python
def apply_pco_correction(
    sat_pos_eci: np.ndarray,    # 卫星位置(ECI)
    sat_vel_eci: np.ndarray,    # 卫星速度(ECI)
    quaternion: np.ndarray,       # 姿态四元数 [q0,q1,q2,q3]
    pco_ecf: np.ndarray,         # PCO在卫星本体坐标系(mm)
    ant_excenter_ecf: np.ndarray # 天线偏心在卫星本体坐标系(mm)
) -> float:
    """
    计算天线相位中心偏移对测距的校正量
    返回: 校正距离 (m), 加到几何距离上
    """
    # 1. 四元数旋转: 本体系PCO → ECI系
    pco_eci = quat_rotate(quaternion, pco_ecf / 1000.0)  # mm→m
    
    # 2. 计算ECI系中的几何距离校正
    #    Δρ = (sat_pos_eci + pco_eci - rec_pos_eci)的方向点积pco_eci
    #    实际实现: Δρ = |r_sat_pco| - |r_sat| ≈ (r_sat·pco)/|r_sat|
    delta_rho = np.dot(sat_pos_eci + pco_eci, pco_eci) / np.linalg.norm(sat_pos_eci + pco_eci)
    return delta_rho
```

**PCV校正**：对每对收发卫星，根据高度角和方位角查表ANTEX，得到PCV校正值（mm），转为米。

### 4.3 模块3：力学模型（Dynamics）

**职责**：计算卫星在ECI系中受到的加速度。

**接口设计**（解耦关键）：
```python
class DynamicsModel(ABC):
    """力学模型抽象基类 — 所有动力学模型实现此接口"""
    
    @abstractmethod
    def compute_acceleration(
        self,
        t: float,                    # 当前时刻 (GPS秒)
        state: np.ndarray,            # [x,y,z, vx,vy,vz] (ECI)
        params: Dict[str, Any]       # 当前估计的参数值
    ) -> np.ndarray:                 # 返回 [ax,ay,az] (m/s²)
        ...
    
    @abstractmethod
    def get_parameter_names(self) -> List[str]:
        """返回此模型引入的待估参数名列表"""
        ...
    
    @abstractmethod
    def get_parameter_a_priori(self) -> Dict[str, Tuple[float, float]]:
        """返回(先验值, 先验标准差)字典"""
        ...
```

**具体实现类**：

```python
class TwoBodyDynamics(DynamicsModel):
    """二体问题 — 实时模式用"""
    def compute_acceleration(self, t, state, params):
        r = state[0:3]
        return -MU_EARTH * r / np.linalg.norm(r)**3

class J2Dynamics(DynamicsModel):
    """J2摄动 — 实时模式增强"""
    def compute_acceleration(self, t, state, params):
        r = state[0:3]
        norm_r = np.linalg.norm(r)
        factor = 1.5 * J2 * MU_EARTH * RE_EARTH**2 / norm_r**5
        ax = factor * r[0] * (5 * r[2]**2 / norm_r**2 - 1)
        ay = factor * r[1] * (5 * r[2]**2 / norm_r**2 - 1)
        az = factor * r[2] * (5 * r[2]**2 / norm_r**2 - 3)
        return np.array([ax, ay, az])

class GravityFieldDynamics(DynamicsModel):
    """地球重力场（球谐系数）— NRT/Delayed模式"""
    def __init__(self, max_degree: int, gravity_model: str):
        self.max_degree = max_degree
        self.C, self.S = load_gravity_model(gravity_model)  # 读取EIGEN-6S4等
        
    def compute_acceleration(self, t, state, params):
        # 球谐函数计算重力梯度 → 加速度
        # 使用完全规格化球谐系数
        return gravity_spherical_harmonic(state[0:3], self.C, self.S, self.max_degree)

class NBodyDynamics(DynamicsModel):
    """N体摄动（日月+大行星）— NRT/Delayed"""
    def __init__(self, de_ephemeris_path: str):
        self.ephem = load_de_ephemeris(de_ephemeris_path)  # DE440
        
    def compute_acceleration(self, t, state, params):
        # 计算日月位置（DE440）→ 第三者引力加速度
        sun_pos = self.ephem.sun_position(t)
        moon_pos = self.ephem.moon_position(t)
        a_sun = G * M_SUN * (sun_pos - r) / |sun_pos - r|^3 - G * M_SUN * sun_pos / |sun_pos|^3
        a_moon = ...  # 同理
        return a_sun + a_moon

class AtmosphericDragDynamics(DynamicsModel):
    """大气阻力 — NRT/Delayed"""
    def __init__(self, macro_model_path: str):
        self.macro = load_macro_model(macro_model_path)  # 卫星宏模型: 各面面积+光学系数
        
    def compute_acceleration(self, t, state, params):
        Cd = params["CD_drag_coeff"]  # 待估参数
        rho = nrlmsise00(t, state[0:3])  # 大气密度模型
        v_rel = state[3:6] - earth_rotation_velocity(state[0:3])
        drag_area = self.macro.get_drag_area(state[0:3], params)
        a_drag = -0.5 * Cd * rho * |v_rel| * v_rel * drag_area / satellite_mass
        return a_drag
    
    def get_parameter_names(self):
        return ["CD_drag_coeff"]

class SolarRadiationPressureDynamics(DynamicsModel):
    """太阳辐射压 — NRT/Delayed"""
    def compute_acceleration(self, t, state, params):
        Cr = params["CR_srp_coeff"]
        shadow = conical_shadow_model(t, state[0:3], sun_pos)  # 锥形地影
        srp_area = self.macro.get_srp_area(state[0:3], params)
        a_srp = Cr * shadow * solar_flux / c * srp_area / satellite_mass * sun_dir
        return a_srp

class EmpiricalAccelerationRTN(DynamicsModel):
    """RTN经验加速度 — Reduced-Dynamic核心"""
    def __init__(self, segment_length_sec: float):
        self.seg_len = segment_length_sec  # 6min 或 15min
        
    def compute_acceleration(self, t, state, params):
        # 确定当前时间段index
        seg_idx = int(t // self.seg_len)
        # RTN→ECI旋转矩阵
        R_rtn_to_eci = compute_rtn_rotation(state[0:3], state[3:6])
        # 从params中提取当前段的aR, aT, aN
        aR = params[f"emp_aR_seg{seg_idx}"]
        aT = params[f"emp_aT_seg{seg_idx}"]
        aN = params[f"emp_aN_seg{seg_idx}"]
        a_rtn = np.array([aR, aT, aN])
        return R_rtn_to_eci @ a_rtn
    
    def get_parameter_names(self):
        # 动态生成: 根据弧段长度和数据时长
        n_seg = int(arc_duration / self.seg_len)
        return [f"emp_a{R,T,N}_seg{i}" for i in range(n_seg) for R in ["R","T","N"]]
```

### 4.4 模块4：数值积分器（Integrator）

**职责**：对状态向量和变分方程进行数值积分。

**接口设计**：
```python
class Integrator(ABC):
    """数值积分器抽象基类"""
    
    @abstractmethod
    def integrate(
        self,
        t0: float,
        tf: float,
        y0: np.ndarray,           # 初始状态 [x,y,z, vx,vy,vz] + 变分(可选)
        dynamics: DynamicsModel,
        params: Dict[str, Any],
        step_output_times: List[float]  # 需要输出的时刻列表
    ) -> IntegrationResult:
        """
        从t0积分到tf，返回在step_output_times时刻的状态
        """
        ...
```

**具体实现**：

```python
class RK4Integrator(Integrator):
    """四阶Runge-Kutta — 实时模式用（简单快速）"""
    def integrate(self, t0, tf, y0, dynamics, params, step_output_times):
        # 标准RK4实现，固定步长
        dt = self.step_size  # 通常10s
        ...

class DOPRI8Integrator(Integrator):
    """Dormand-Prince 8阶 — 高精度，自适应步长"""
    def integrate(self, t0, tf, y0, dynamics, params, step_output_times):
        # DOPRI8(7)13方法，自适应步长，高精度
        # 同时积分变分方程（STM = State Transition Matrix）
        ...

class GaussJacksonIntegrator(Integrator):
    """
    Gauss-Jackson多步积分器 — GRACE-FO标准
    特点: 专为轨道力学优化，比RK方法快10-20倍
    """
    def integrate(self, t0, tf, y0, dynamics, params, step_output_times):
        # Gauss-Jackson 8阶预测-校正方法
        # 需要预先存储若干步的历史加速度
        # 输出为等间隔（通常60s或300s）
        ...
```

### 4.5 模块5：参数估计器（Estimator）

**职责**：根据观测数据和动力学模型，估计轨道状态向量和其他参数。

**接口设计**：
```python
class Estimator(ABC):
    """参数估计器抽象基类"""
    
    @abstractmethod
    def estimate(
        self,
        obs_data: List[ObsEpoch],
        dynamics: List[DynamicsModel],
        integrator: Integrator,
        initial_state: np.ndarray,
        initial_params: Dict[str, float],
        apriori_cov: np.ndarray
    ) -> EstimationResult:
        ...
```

**具体实现**：

```python
class LSQEstimator(Estimator):
    """
    最小二乘批处理估计器
    — 近实时/延时模式
    """
    def estimate(self, obs_data, dynamics, integrator, initial_state, initial_params, apriori_cov):
        max_iter = self.config["max_iterations"]  # 通常3-5
        convergence_threshold = 1e-8
        
        x = self._pack_state_params(initial_state, initial_params)
        P = apriori_cov
        
        for iteration in range(max_iter):
            # 1. 用当前x预测观测值（通过数值积分）
            pred_obs, H = self._compute_predictions_and_jacobian(
                x, obs_data, dynamics, integrator
            )
            # 2. 计算残差
            residuals = self._compute_residuals(obs_data, pred_obs)
            # 3. 构建法方程: (HᵀWH + W0⁻¹) Δx = HᵀW v
            N = H.T @ W @ H + W0_inv
            b = H.T @ W @ residuals + W0_inv @ (x_apriori - x)
            dx = solve(N, b)
            x = x + dx
            P = inv(N)
            # 4. 检查收敛
            if norm(dx) < convergence_threshold:
                break
        
        return EstimationResult(
            converged=True,
            state=self._unpack_state(x),
            params=self._unpack_params(x),
            covariance=P,
            residuals=residuals,
            iterations=iteration+1
        )
    
    def _compute_predictions_and_jacobian(self, x, obs_data, dynamics, integrator):
        """
        核心: 对每个观测历元，用当前轨道参数预测观测值
        并计算设计矩阵H（观测对参数的偏导数）
        
        使用数值积分+变分方程同时得到轨道和STM
        STM = ∂X(t)/∂X(t0) → 用于计算H矩阵中轨道参数的偏导数
        """
        # 用integrator同时积分轨道和变分方程
        integration_result = integrator.integrate_with_stm(
            t0, tf, x_state, dynamics, params, obs_times
        )
        # 对每个观测历元计算预测观测值和H行
        ...

class EKFEstimator(Estimator):
    """
    扩展卡尔曼滤波 — 实时模式
    """
    def estimate(self, obs_data, dynamics, integrator, initial_state, initial_params, apriori_cov):
        # 标准EKF: 预测步 + 更新步
        # 预测: X_k|k-1 = f(X_k-1|k-1), P_k|k-1 = F P Fᵀ + Q
        # 更新: K = P Hᵀ (H P Hᵀ + R)⁻¹, X_k|k = X_k|k-1 + K(y-h), P_k|k = (I-KH)P
        ...
```

### 4.6 模块6：模糊度固定器（AR_Module）

**职责**：将浮点模糊度固定为整数（PPP-AR），提升精度。

**接口设计**：
```python
class AmbiguityResolver(ABC):
    @abstractmethod
    def resolve(
        self,
        float_ambiguities: np.ndarray,    # 浮点模糊度 (cycles)
        covariance: np.ndarray,             # 模糊度协方差
        obs_data: List[ObsEpoch]
    ) -> AmbiguityResult:
        ...
```

**LAMBDA算法实现**：
```python
class LAMBDAAmbiguityResolver(AmbiguityResolver):
    """
    LAMBDA (Least-squares AMBiguity Decorrelation Adjustment)
    标准PPP-AR模糊度固定算法
    """
    def resolve(self, float_amb, cov, obs_data):
        # 1. Z变换（降相关）
        Z, D = integer_gaussian_transform(cov)
        y_tilde = Z.T @ float_amb
        
        # 2. 搜索整数向量（球状搜索）
        candidates = spherical_search(y_tilde, D, search_radius=3.0)
        
        # 3. 计算Ratio检验
        best = candidates[0]
        second = candidates[1] if len(candidates) > 1 else None
        ratio = second.norm_residual / best.norm_residual if second else float('inf')
        
        if ratio > self.config["ratio_threshold"]:  # 通常3.0
            return AmbiguityResult(
                fixed=True,
                integer_amb=best.z_int,
                ratio=ratio,
                fixed_ambiguity_vector=Z @ best.z_int
            )
        else:
            return AmbiguityResult(fixed=False, ratio=ratio)
```

### 4.7 模块7：质量评估器（QAModule）

**职责**：评估定轨结果质量，输出精度报告。

```python
class QualityAssessor:
    def evaluate_overlap(self, arc1_result, arc2_result, overlap_window_sec: float) -> OverlapStats:
        """重叠弧段检验：相邻弧段在重叠区的轨道差异RMS"""
        ...
    
    def evaluate_residuals(self, residuals: np.ndarray, obs_type: str) -> ResidualStats:
        """残差统计：均值、标准差、RMS、异常值比例"""
        ...
    
    def evaluate_sigma0(self, residuals, HWH_invHT: np.ndarray) -> float:
        """单位权重标准差 σ0 = sqrt(vᵀWv / (n - m))"""
        ...
```

---

## 4.X 卫星宏模型计算模块 (Satellite Macro Model Module)

> 宏模型（Macro Model）是自研动力学后端和 Orekit 后端的**共享输入模块**。
> 无论使用哪种后端，非保守力建模都需要同一份卫星几何-光学参数。
> 本节设计一个独立的 `MacroModel` 模块，负责加载配置、计算有效面积、生成各面板受力。

### 4.X.1 设计原则

```
                    ┌─────────────────────┐
                    │  config/macro_model │
                    │      .yaml          │
                    └────────┬────────────┘
                             │ 加载
                    ┌────────▼────────────┐
                    │   MacroModel 类     │
                    │  (纯Python, 独立)    │
                    │                     │
                    │ compute_eff_area()  │──► 用于自研 DynamicsModel
                    │ compute_drag_acc()  │──► 用于自研 DynamicsModel
                    │ compute_srp_acc()   │──► 用于自研 DynamicsModel
                    │ to_orekit_sc()      │──► 用于 OrekitAdapter
                    └─────────────────────┘
```

关键约束：
- **单一数据源**：宏模型参数只在 `config/macro_model.yaml` 一处定义，两端共享
- **后端无关**：`MacroModel` 类不 import 任何 Orekit 代码；Orekit 转换在适配层完成
- **可配置驱动**：所有面板数量、面积、光学系数均通过配置文件控制，不硬编码

### 4.X.2 宏模型配置文件 (`config/macro_model.yaml`)

将前述 JSON 宏模型升级为完整的 YAML 配置文件，包含卫星物理参数和所有面板定义：

```yaml
# config/macro_model.yaml
# ============================================================
# 卫星宏模型配置文件 — 用于大气阻力和太阳辐射压计算
# 本文件由 自研 DynamicsModel 和 OrekitAdapter 共同读取
# ============================================================

satellite:
  id: "SAR_SAT_01"
  name: "C-band SAR Satellite"
  description: >
    Sun-synchronous dawn-dusk orbit, ~500 km altitude,
    right-looking SAR with 25 m² antenna

  # ---- 质量参数 ----
  mass:
    dry_mass_kg: 1200.0                  # 干质量
    fuel_mass_kg: 80.0                   # 初始燃料质量
    total_mass_kg: 1280.0                # 发射时总质量 = dry + fuel
    mass_sigma_kg: 5.0                   # 质量不确定性 (1σ)
    # 燃料消耗模型（可选，用于长时间弧段的质量变化）
    fuel_consumption:
      enabled: false                     # 延时模式可启用
      rate_kg_per_day: 0.02              # 日均燃料消耗

  # ---- 质心参数 ----
  center_of_mass:
    offset_sbf_m: [0.05, 0.0, -0.15]     # 质心在星体系中的偏移 (m)
    uncertainty_m: [0.01, 0.01, 0.02]   # 偏移不确定性 (1σ)

  # ---- 面板定义 ----
  # 每个面板定义了面积、法向量（星体系）、光学系数
  # 法向量惯例：指向卫星外侧
  panels:
    - name: "Star Body +X (RAM/沿飞行前侧)"  # 大气阻力主导面
      area_m2: 6.25
      normal_vector_sbf: [1.0, 0.0, 0.0]      # 指向飞行方向
      optical:
        reflectivity_coeff: 0.30               # 镜面反射系数 ρ_s
        diffuse_coeff: 0.40                    # 漫反射系数 ρ_d
        absorption_coeff: 0.30                 # 吸收系数 = 1-ρ_s-ρ_d
      # 热控制材料：MLI 多层隔热，白色涂层
      material: "MLI with white coating"

    - name: "Star Body -X (Wake/沿飞行后侧)"
      area_m2: 6.25
      normal_vector_sbf: [-1.0, 0.0, 0.0]     # 指向尾流方向
      optical:
        reflectivity_coeff: 0.25
        diffuse_coeff: 0.45
        absorption_coeff: 0.30
      material: "MLI with white coating"
      note: "尾流面，阻力贡献近似为零"

    - name: "Star Body +Y (轨道面正法向/右侧视)"
      area_m2: 8.0
      normal_vector_sbf: [0.0, 1.0, 0.0]      # 指向轨道面正法向
      optical:
        reflectivity_coeff: 0.20
        diffuse_coeff: 0.50
        absorption_coeff: 0.30
      material: "MLI"

    - name: "Star Body -Y (轨道面负法向/左侧视)"
      area_m2: 8.0
      normal_vector_sbf: [0.0, -1.0, 0.0]
      optical:
        reflectivity_coeff: 0.20
        diffuse_coeff: 0.50
        absorption_coeff: 0.30
      material: "MLI"

    - name: "Star Body +Z (天顶方向)"
      area_m2: 8.0
      normal_vector_sbf: [0.0, 0.0, 1.0]      # 指向天顶
      optical:
        reflectivity_coeff: 0.15
        diffuse_coeff: 0.55
        absorption_coeff: 0.30
      material: "MLI with OSR radiator"
      note: "SRP主导面，晨昏轨道持续受照"

    - name: "Star Body -Z (天底/对地方向)"
      area_m2: 8.0
      normal_vector_sbf: [0.0, 0.0, -1.0]     # 指向地心
      optical:
        reflectivity_coeff: 0.20
        diffuse_coeff: 0.60
        absorption_coeff: 0.20
      material: "MLI"
      note: "SAR天线安装面"

    - name: "Solar Wing +Y (太阳翼受照面)"
      area_m2: 12.0
      normal_vector_sbf: [0.0, 1.0, 0.0]      # 受照面法向量；跟踪太阳时动态变化
      optical:
        reflectivity_coeff: 0.05               # 太阳能电池片镜面反射低
        diffuse_coeff: 0.15
        absorption_coeff: 0.80                 # 吸收率高达80%→光电转换
      material: "Triple-junction GaAs solar cells"
      note: "法向量在Y-Z面内跟踪太阳；需外部提供solar_wing_angle"

    - name: "Solar Wing -Y (太阳翼背面)"
      area_m2: 12.0
      normal_vector_sbf: [0.0, -1.0, 0.0]
      optical:
        reflectivity_coeff: 0.10
        diffuse_coeff: 0.30
        absorption_coeff: 0.60
      material: "Carbon fiber + OSR"
      note: "散热面，光学特性与受照面不同"

    - name: "SAR Antenna (合成孔径雷达天线面)"
      area_m2: 25.0
      normal_vector_sbf: [0.0, 0.0, -1.0]     # 天线波束指向天底方向
      optical:
        reflectivity_coeff: 0.08
        diffuse_coeff: 0.65
        absorption_coeff: 0.27
      material: "CFRP honeycomb + metallic coating"
      note: "最大面积面板；对阻力和SRP都敏感；法向量固定天底方向"

    - name: "SAR Antenna Back (SAR天线背面)"
      area_m2: 25.0
      normal_vector_sbf: [0.0, 0.0, 1.0]      # 天线背面指向天顶
      optical:
        reflectivity_coeff: 0.12
        diffuse_coeff: 0.58
        absorption_coeff: 0.30
      material: "CFRP honeycomb"
      note: "天线背面热控特性与反射面不同"

  # ---- 汇总参数（自动计算，仅供参考） ----
  summary:
    total_panel_count: 10
    total_panel_area_m2: 118.5
    mass_area_ratio_kg_per_m2: 10.80
    # 典型参考值，用于 IsotropicDrag 快速模式
    reference_drag_area_m2: 15.0
    reference_srp_area_m2: 35.0

  # ---- 太阳翼跟踪 ----
  solar_wing:
    tracking_mode: "sun-pointing"              # sun-pointing / fixed / body-fixed
    # sun-pointing: 太阳翼绕Y轴旋转跟踪太阳
    # fixed: 固定角度（如 launch_lock_angle）
    # body-fixed: 与星体固联
    max_rotation_rate_deg_per_sec: 0.5        # 最大旋转角速度
    # 如果无法获得实时太阳翼转角，可用近似：
    # β_angle_approximation: true  # 用β角近似太阳翼转角
```

### 4.X.3 MacroModel 核心类设计

```python
# pod_core/macro_model/macro_model.py

import numpy as np
from dataclasses import dataclass
from typing import List, Optional, Dict
import yaml

@dataclass
class PanelDef:
    """单个面板的物理定义"""
    name: str
    area_m2: float
    normal_vector_sbf: np.ndarray      # (3,) 星体系法向量
    reflectivity_coeff: float           # ρ_s 镜面反射
    diffuse_coeff: float               # ρ_d 漫反射
    absorption_coeff: float             # ρ_a = 1 - ρ_s - ρ_d
    material: str = ""
    note: str = ""

@dataclass
class SatellitePhysicalParams:
    """卫星物理参数"""
    satellite_id: str
    dry_mass_kg: float
    fuel_mass_kg: float
    total_mass_kg: float
    mass_sigma_kg: float
    center_of_mass_offset_sbf: np.ndarray  # (3,)
    com_uncertainty_m: np.ndarray          # (3,)
    reference_drag_area_m2: float
    reference_srp_area_m2: float

class MacroModel:
    """
    卫星宏模型 — 非保守力计算核心

    功能：
    1. 从 YAML 配置文件加载面板定义和卫星物理参数
    2. 计算给定姿态下对特定方向的有效截面积
    3. 计算大气阻力加速度矢量（自研后端用）
    4. 计算太阳辐射压加速度矢量（自研后端用）
    5. 导出 Orekit BoxAndSolarArraySpacecraft 对象（Orekit后端用）

    Usage:
        mm = MacroModel.from_yaml("config/macro_model.yaml")
        A_eff = mm.effective_area(attitude_quat, direction, "drag")
        a_drag = mm.compute_drag_acceleration(
            state_eci, attitude_quat, rho, v_sat_eci, Cd
        )
    """

    # ---- 物理常数 ----
    SPEED_OF_LIGHT = 299792458.0           # m/s
    SOLAR_CONSTANT = 1361.0                # W/m² at 1 AU
    SOLAR_PRESSURE_AT_1AU = SOLAR_CONSTANT / SPEED_OF_LIGHT  # ~4.54e-6 N/m²

    def __init__(self,
                 satellite_params: SatellitePhysicalParams,
                 panels: List[PanelDef],
                 solar_wing_config: Optional[Dict] = None):
        self._sat = satellite_params
        self._panels = panels
        self._solar_wing = solar_wing_config or {}

        # 预计算和校验
        self._validate_panels()
        self._n_panels = len(panels)

    # ================================================================
    # 工厂方法
    # ================================================================

    @classmethod
    def from_yaml(cls, yaml_path: str) -> "MacroModel":
        """从 YAML 配置文件加载宏模型"""
        with open(yaml_path, 'r') as f:
            data = yaml.safe_load(f)

        sat_cfg = data["satellite"]

        # 解析面板
        panels = []
        for p in sat_cfg["panels"]:
            panels.append(PanelDef(
                name=p["name"],
                area_m2=p["area_m2"],
                normal_vector_sbf=np.array(p["normal_vector_sbf"]),
                reflectivity_coeff=p["optical"]["reflectivity_coeff"],
                diffuse_coeff=p["optical"]["diffuse_coeff"],
                absorption_coeff=p["optical"]["absorption_coeff"],
                material=p.get("material", ""),
                note=p.get("note", ""),
            ))

        # 解析卫星参数
        sat_params = SatellitePhysicalParams(
            satellite_id=sat_cfg["id"],
            dry_mass_kg=sat_cfg["mass"]["dry_mass_kg"],
            fuel_mass_kg=sat_cfg["mass"]["fuel_mass_kg"],
            total_mass_kg=sat_cfg["mass"]["total_mass_kg"],
            mass_sigma_kg=sat_cfg["mass"]["mass_sigma_kg"],
            center_of_mass_offset_sbf=np.array(
                sat_cfg["center_of_mass"]["offset_sbf_m"]
            ),
            com_uncertainty_m=np.array(
                sat_cfg["center_of_mass"]["uncertainty_m"]
            ),
            reference_drag_area_m2=sat_cfg["summary"]["reference_drag_area_m2"],
            reference_srp_area_m2=sat_cfg["summary"]["reference_srp_area_m2"],
        )

        return cls(sat_params, panels, sat_cfg.get("solar_wing"))

    # ================================================================
    # 有效面积计算
    # ================================================================

    def effective_area(self,
                       attitude_quaternion: np.ndarray,  # [q0,q1,q2,q3] SBF→ECI
                       direction_eci: np.ndarray,         # 归一化方向 (ECI系)
                       mode: str = "drag",                # "drag" / "srp"
                       ) -> float:
        """
        计算给定姿态下对给定方向的有效截面积。

        算法：
          A_eff = Σ max(A_i * (n̂_i_ECI · direction), 0)
        即：各面板法向量转到ECI系后，取与目标方向的正投影面积之和。

        Args:
            attitude_quaternion: 星体→ECI 四元数 [q0, q1, q2, q3]
            direction_eci: 目标方向单位矢量 (ECI)
              - mode="drag": 使用 v̂_rel (相对大气来流方向)
              - mode="srp":  使用 -r̂_sun (太阳光入射反向)
            mode: "drag" 或 "srp"

        Returns:
            A_eff: 有效截面积 (m²)
        """
        # 四元数→旋转矩阵
        R = self._quaternion_to_matrix(attitude_quaternion)

        A_eff = 0.0
        for panel in self._panels:
            # 面板法向量从星体系转到ECI
            n_eci = R @ panel.normal_vector_sbf
            # 投影面积（仅计入正投影）
            cos_theta = max(np.dot(n_eci, direction_eci), 0.0)
            A_eff += panel.area_m2 * cos_theta

        return A_eff

    # ================================================================
    # 大气阻力加速度计算
    # ================================================================

    def compute_drag_acceleration(self,
                                  position_eci: np.ndarray,       # (3,) ECI位置 (m)
                                  velocity_eci: np.ndarray,       # (3,) ECI速度 (m/s)
                                  attitude_quaternion: np.ndarray, # [q0,q1,q2,q3]
                                  atmospheric_density: float,     # ρ (kg/m³)
                                  Cd: float,                      # 阻力系数
                                  ) -> np.ndarray:
        """
        计算大气阻力加速度矢量 (ECI系)。

        公式: a_drag = -0.5 * Cd * (A_eff/m) * ρ * |v_rel| * v̂_rel

        Args:
            position_eci: 卫星位置 (ECI, m)
            velocity_eci: 卫星速度 (ECI, m/s)
            attitude_quaternion: 姿态四元数
            atmospheric_density: 大气密度 ρ 来自 NRLMSISE-00 或 DTM-2013
            Cd: 阻力系数（待估参数）

        Returns:
            a_drag_eci: 阻力加速度矢量 (ECI, m/s²)
        """
        # 简化：假设大气与地球共旋
        # v_rel = v_sat - ω_earth × r_sat
        omega_earth = np.array([0.0, 0.0, 7.2921150e-5])  # rad/s
        v_atm = np.cross(omega_earth, position_eci)
        v_rel = velocity_eci - v_atm

        v_rel_norm = np.linalg.norm(v_rel)
        if v_rel_norm < 1e-12:
            return np.zeros(3)

        v_rel_unit = v_rel / v_rel_norm

        # 计算当前姿态下对来流方向的有效面积
        A_eff = self.effective_area(attitude_quaternion, v_rel_unit, mode="drag")

        # 阻力加速度
        a_magnitude = 0.5 * Cd * (A_eff / self._sat.total_mass_kg) * \
                      atmospheric_density * v_rel_norm**2

        return -a_magnitude * v_rel_unit

    # ================================================================
    # 太阳辐射压加速度计算
    # ================================================================

    def compute_srp_acceleration(self,
                                 position_eci: np.ndarray,        # (3,) ECI位置 (m)
                                 sun_direction_eci: np.ndarray,   # (3,) 太阳方向单位矢量 (ECI)
                                                                  # 指向太阳
                                 attitude_quaternion: np.ndarray,
                                 Cr: float,                       # 光压系数
                                 in_shadow: bool = False,
                                 ) -> np.ndarray:
        """
        计算太阳辐射压加速度矢量 (ECI系)，按面板分解。

        公式 (每面板):
          a_i = -(S₀/c) * (A_i * Cr / m) * cos(θ_i) *
                [(1-ρ_s)·r̂_sun + 2·(ρ_s/3 + ρ_d)·n̂_i]

        总加速度 = Σ a_i

        Args:
            position_eci: 卫星位置 (ECI)
            sun_direction_eci: 太阳方向单位矢量 (指向太阳)
            attitude_quaternion: 姿态四元数
            Cr: 光压系数（待估参数）
            in_shadow: 是否在地影中（True→返回零矢量）

        Returns:
            a_srp_eci: SRP加速度矢量 (ECI, m/s²)
        """
        if in_shadow:
            return np.zeros(3)

        R = self._quaternion_to_matrix(attitude_quaternion)

        # 光入射方向（从太阳指向卫星）：与太阳方向相反
        incidence_dir = -sun_direction_eci

        a_total = np.zeros(3)
        for panel in self._panels:
            n_eci = R @ panel.normal_vector_sbf
            cos_theta = max(np.dot(n_eci, incidence_dir), 0.0)

            if cos_theta < 1e-12:
                continue

            # 反射分量
            rho_s = panel.reflectivity_coeff
            rho_d = panel.diffuse_coeff

            # 光压公式 (Montenbruck & Gill §3.4)
            # a = -P * (A/m) * Cr * cosθ * [ (1 - ρ_s)*ŝ + 2*(ρ_s/3 + ρ_d)*n̂ ]
            P = self.SOLAR_PRESSURE_AT_1AU  # ~4.54e-6 N/m²
            area_mass_ratio = panel.area_m2 / self._sat.total_mass_kg

            # 各向同性分量（吸收 + 漫反射的均匀部分）
            a_isotropic = (1.0 - rho_s) * incidence_dir

            # 法向分量（镜面反射 + 漫反射的非均匀部分）
            a_normal = 2.0 * (rho_s / 3.0 + rho_d) * n_eci

            a_panel = -P * area_mass_ratio * Cr * cos_theta * (a_isotropic + a_normal)
            a_total += a_panel

        return a_total

    # ================================================================
    # 地影检测
    # ================================================================

    def check_shadow(self,
                     position_eci: np.ndarray,
                     sun_direction_eci: np.ndarray,
                     shadow_model: str = "conical",
                     ) -> float:
        """
        计算地影因子 ν ∈ [0, 1]

        0 = 全影 (umbra)
        (0,1) = 半影 (penumbra, 仅 conical 模型)
        1 = 全光照

        Args:
            position_eci: 卫星位置 (ECI, m)
            sun_direction_eci: 太阳方向单位矢量
            shadow_model: "cylindrical" (柱形, 无半影) / "conical" (锥形, 含半影)

        Returns:
            shadow_factor: 地影因子 [0, 1]
        """
        # 简化柱形地影模型
        R_EARTH = 6378137.0  # m

        # 卫星到太阳连线与地心距离
        r = np.linalg.norm(position_eci)
        sun_dir = sun_direction_eci

        # 卫星在地球阴影锥中的投影
        proj = np.dot(position_eci, sun_dir)

        if proj > 0:
            # 卫星在地球日照侧
            if shadow_model == "cylindrical":
                return 1.0
            else:
                # 锥形模型：计算到阴影轴的距离
                d_perp = np.sqrt(r**2 - proj**2)
                # 半影区判断
                R_SUN = 696340000.0  # 太阳半径 (m)
                AU = 149597870700.0  # 天文单位 (m)
                # 锥角
                alpha_umbra = np.arcsin((R_SUN - R_EARTH) / AU)
                alpha_penumbra = np.arcsin((R_SUN + R_EARTH) / AU)
                # ... 完整实现需要计算锥体交线
                # 此处返回简化结果
                if d_perp < R_EARTH:
                    return 0.0  # 全影
                else:
                    return 1.0  # 光照

        return 1.0  # 光照

    # ================================================================
    # Orekit 导出接口
    # ================================================================

    def to_orekit_config(self) -> dict:
        """
        导出为 Orekit适配层需要的配置字典。

        OrekitAdapter 使用此字典构造:
          - BoxAndSolarArraySpacecraft → DragForce / SolarRadiationPressure
          或
          - IsotropicDrag (简化模式)

        Returns:
            config dict with keys:
              - satellite_mass_kg
              - panels: list of {area, normal, reflectivity, diffuse, absorption}
              - reference_drag_area
              - reference_srp_area
        """
        return {
            "satellite_mass_kg": self._sat.total_mass_kg,
            "dry_mass_kg": self._sat.dry_mass_kg,
            "center_of_mass_offset_sbf": self._sat.center_of_mass_offset_sbf.tolist(),
            "panels": [
                {
                    "name": p.name,
                    "area_m2": p.area_m2,
                    "normal_vector_sbf": p.normal_vector_sbf.tolist(),
                    "reflectivity_coeff": p.reflectivity_coeff,
                    "diffuse_coeff": p.diffuse_coeff,
                    "absorption_coeff": p.absorption_coeff,
                }
                for p in self._panels
            ],
            "reference_drag_area_m2": self._sat.reference_drag_area_m2,
            "reference_srp_area_m2": self._sat.reference_srp_area_m2,
        }

    # ================================================================
    # 内部工具方法
    # ================================================================

    @staticmethod
    def _quaternion_to_matrix(q: np.ndarray) -> np.ndarray:
        """四元数 [q0,q1,q2,q3] → 3×3 旋转矩阵 (SBF→ECI)"""
        q0, q1, q2, q3 = q
        return np.array([
            [1 - 2*(q2**2 + q3**2), 2*(q1*q2 - q0*q3),     2*(q1*q3 + q0*q2)],
            [2*(q1*q2 + q0*q3),     1 - 2*(q1**2 + q3**2), 2*(q2*q3 - q0*q1)],
            [2*(q1*q3 - q0*q2),     2*(q2*q3 + q0*q1),     1 - 2*(q1**2 + q2**2)],
        ])

    def _validate_panels(self):
        """校验面板光学系数：ρ_s + ρ_d + ρ_a == 1.0（容差 1e-6）"""
        for i, p in enumerate(self._panels):
            total = p.reflectivity_coeff + p.diffuse_coeff + p.absorption_coeff
            if abs(total - 1.0) > 1e-6:
                raise ValueError(
                    f"Panel[{i}] '{p.name}': optical coeffs sum to "
                    f"{total:.6f} ≠ 1.0 (ρ_s={p.reflectivity_coeff}, "
                    f"ρ_d={p.diffuse_coeff}, ρ_a={p.absorption_coeff})"
                )

    # ---- 属性访问 ----
    @property
    def satellite_mass_kg(self) -> float:
        return self._sat.total_mass_kg

    @property
    def n_panels(self) -> int:
        return self._n_panels

    @property
    def panels(self) -> List[PanelDef]:
        return self._panels
```

### 4.X.4 MacroModel 在各后端的集成方式

#### 自研动力学后端

```python
# pod_core/dynamics/self_dynamics.py

class SelfDynamicsModel(DynamicsModel):
    def __init__(self, config: dict, macro_model: MacroModel):
        self._mm = macro_model
        self._gravity_model = load_gravity_field(config["gravity_field"])
        self._density_model = load_atmosphere_model(config["atmospheric_drag"])
        # ...

    def compute_acceleration(self, t, state_eci, params):
        pos = state_eci[0:3]
        vel = state_eci[3:6]

        # 重力
        a_grav = self._gravity_model.acceleration(pos)

        # 大气阻力 — 使用宏模型
        rho = self._density_model.get_density(pos, t)
        if params.get("estimate_Cd", False):
            Cd = params["Cd"]  # 从待估参数中取值
        else:
            Cd = self._config["atmospheric_drag"]["Cd_apriori"]

        a_drag = self._mm.compute_drag_acceleration(
            pos, vel,
            attitude_quaternion=self._get_attitude(t),
            atmospheric_density=rho,
            Cd=Cd,
        )

        # 太阳辐射压 — 使用宏模型
        sun_dir = self._get_sun_direction(t)
        is_shadow = self._mm.check_shadow(pos, sun_dir) < 0.01
        if params.get("estimate_Cr", False):
            Cr = params["Cr"]
        else:
            Cr = self._config["solar_radiation_pressure"]["Cr_apriori"]

        a_srp = self._mm.compute_srp_acceleration(
            pos, sun_dir,
            attitude_quaternion=self._get_attitude(t),
            Cr=Cr,
            in_shadow=is_shadow,
        )

        # ... 经验加速度等

        return a_grav + a_drag + a_srp + ...
```

#### Orekit 后端

```python
# OrekitAdapter 中宏模型的使用见 Section 5.4
```

---

## 5. Orekit适配层设计

### 5.1 设计目标

Orekit 是用 Java 编写的航天动力学库。系统用 Python 实现，需通过适配层使用 Orekit。
Orekit **同样需要卫星宏模型参数**来建模非保守力——它提供了 `BoxAndSolarArraySpacecraft` 和 `IsotropicDrag` 两个接口接收这些参数。

关键原则：
- **宏模型参数来自同一份 YAML 文件**（`config/macro_model.yaml`），自研和 Orekit 两端共享数据源
- Orekit 适配层负责将 `MacroModel` 对象转换为 Orekit Java 对象
- 力学模型配置（重力场阶数、密度模型选择、积分器类型等）由 `config/dynamics.yaml` 统一管理

### 5.2 Orekit 中的宏模型类映射

| Orekit Java 类 | 功能 | 对应我们的模块 | 需要的参数 |
|:---|------|:---|------|
| `BoxAndSolarArraySpacecraft` | 完整 Box-Wing 面板模型 | `MacroModel` (全部面板) | 每面板面积、法向量、ρ_s、ρ_d；太阳翼面积/法向量/转角 |
| `IsotropicDrag` | 简化各向同性阻力模型 | `MacroModel` (仅 reference_drag_area) | 固定截面积 + Cd |
| `DragForce` | 大气阻力力模型 | 使用上述 Spacecraft 对象 | Spacecraft + Atmosphere 模型 |
| `SolarRadiationPressure` | 太阳辐射压力模型 | 使用上述 Spacecraft 对象 | Spacecraft + 地影模型 |
| `DTM2000` / `NRLMSISE00` | 大气密度模型 | 对应我们的 `density_model` 配置 | F10.7, Ap, 时间, 位置 |
| `ConicalShadow` / `CylindricalShadow` | 地影模型 | 对应我们的 `shadow_model` 配置 | 太阳/地球半径、AU距离 |

### 5.3 Orekit Java 侧：NumericalPropagator 构造流程

以下 Java 伪代码展示 Orekit 如何消费宏模型参数（适配层通过 py4j 调用）。宏模型参数 `MacroModel.to_orekit_config()` 的字典被映射到 Java 对象：

```java
// ===== OrekitAdapter._create_orekit_propagator() 对应的Java逻辑 =====

// 1. 读取宏模型面板
BoxAndSolarArraySpacecraft spacecraft = new BoxAndSolarArraySpacecraft(
    1200.0,                          // satellite_mass_kg
    new Vector3D(6.25, ...),         // +X面板面积*法向量（含正负号）
    new Vector3D(6.25, ...),         // -X
    new Vector3D(8.0, ...),          // +Y
    new Vector3D(8.0, ...),          // -Y
    new Vector3D(8.0, ...),          // +Z
    new Vector3D(8.0, ...),          // -Z
    0.30, 0.40,                      // +X面: ρ_s, ρ_d
    0.25, 0.45,                      // -X面
    0.20, 0.50,                      // +Y面
    0.20, 0.50,                      // -Y面
    0.15, 0.55,                      // +Z面
    0.20, 0.60,                      // -Z面
    // 太阳翼
    new Vector3D(12.0, 0.0, 0.0),    // 太阳翼面积*法向量
    0.05, 0.15                       // 太阳翼 ρ_s, ρ_d
);

// 注：BoxAndSolarArraySpacecraft 仅支持6面板+太阳翼。
// 对于SAR天线等额外面板，需要扩展为自定义的
// ExtendedBoxAndSolarArraySpacecraft 或使用 Orekit 的
// GenericSpacecraft 接口（如果版本支持）。

// 2. 创建力模型并加入 propagator
NumericalPropagator propagator = new NumericalPropagator(integrator);

// 2a. 重力场
propagator.addForceModel(
    new HolmesFeatherstoneAttractionModel(
        FramesFactory.getITRF(IERSConventions.IERS_2010, true),
        GravityFieldFactory.getNormalizedProvider(150, 150)
    )
);

// 2b. 大气阻力 — 使用宏模型的 BoxAndSolarArraySpacecraft
if (dragEnabled) {
    Atmosphere atmosphere = new DTM2000(
        new MarshallSolarActivityFutureEstimation("Jan2000F10.txt"),
        CelestialBodyFactory.getSun()
    );
    DragForce drag = new DragForce(atmosphere, spacecraft);
    drag.setDragCoefficientGradientNames(
        new ParameterDriver("Cd", 2.2, 0.5, 1.5, 3.0).getName()
    );
    propagator.addForceModel(drag);
}

// 2c. 太阳辐射压 — 使用宏模型的 BoxAndSolarArraySpacecraft
if (srpEnabled) {
    SolarRadiationPressure srp = new SolarRadiationPressure(
        Constants.SUN_RADIUS,
        Constants.IAU_2012_NOMINAL_EARTH_EQUATORIAL_RADIUS,
        new ConicalShadow(Constants.SUN_RADIUS, Constants.IAU_2012_NOMINAL_EARTH_EQUATORIAL_RADIUS),
        spacecraft
    );
    srp.getCrParameterDriver().setValue(1.3);
    propagator.addForceModel(srp);
}

// 2d. 经验加速度 (RTN 分段常数)
if (empiricalEnabled) {
    ParametricAcceleration empirical = new ParametricAcceleration(
        FramesFactory.getEME2000(),
        new HarmonicAcceleration(
            Vector3D.ZERO,                     // 常值项 (通常设零)
            HarmonicAcceleration.Order.ORDER_0, // 仅常数项
            900.0                               // 分段长度 (秒)
        )
    );
    propagator.addForceModel(empirical);
}
```

### 5.4 Orekit 适配层 Python 接口（重新设计）

```python
# pod_core/orekit_adapter/orekit_dynamics.py

class OrekitDynamicsAdapter(DynamicsModel):
    """
    Orekit动力学适配层 — 完整版

    将外部配置文件 (dynamics.yaml + macro_model.yaml) 映射到
    Orekit NumericalPropagator 的完整力模型配置。

    构造流程：
      1. 加载 dynamics.yaml → 确定启用哪些力模型
      2. 加载 macro_model.yaml → 构造 BoxAndSolarArraySpacecraft
      3. 通过 py4j/Jpype 调用 Orekit Java API 构造 propagator
      4. 在 compute_acceleration() 中委托给 Orekit 计算
    """

    def __init__(self,
                 dynamics_config: dict,
                 macro_model: MacroModel,
                 data_dir: str = "./orekit-data"):
        """
        Args:
            dynamics_config: 来自 config/dynamics.yaml 的完整配置字典
            macro_model: 来自 config/macro_model.yaml 的 MacroModel 实例
            data_dir: Orekit数据包路径 (EOP, 重力场, 海潮等)
        """
        self._config = dynamics_config
        self._mm = macro_model
        self._data_dir = data_dir

        # 转换为 Orekit 格式
        self._orekit_spacecraft_config = macro_model.to_orekit_config()

        # 初始化 JVM
        if not OrekitJVM.is_started():
            OrekitJVM.start(data_dir)

        # 构造 Orekit NumericalPropagator
        self._propagator = self._build_propagator()

    def _build_propagator(self):
        """
        根据 dynamics.yaml 配置构造 Orekit NumericalPropagator

        配置到 Orekit 对象的映射：
        ┌──────────────────────┬───────────────────────────┐
        │ dynamics.yaml 段      │ Orekit 对象               │
        ├──────────────────────┼───────────────────────────┤
        │ gravity_field        │ HolmesFeatherstoneAttract. │
        │ n_body               │ ThirdBodyAttraction        │
        │ tides                │ SolidTides + OceanTides    │
        │ atmospheric_drag     │ DragForce + DTM2000/NRLMS  │
        │ solar_radiation_pres │ SolarRadiationPressure     │
        │ empirical_accel      │ ParametricAcceleration     │
        │ relativity           │ Relativity                 │
        │ integrator           │ DormandPrince853 / GaussJ  │
        └──────────────────────┴───────────────────────────┘
        """
        # 通过 py4j 调用 Java
        from py4j.java_gateway import JavaGateway
        gw = JavaGateway()

        # 积分器
        integrator = gw.jvm.org.orekit.propagation.numerical. \
            NumericalPropagator.tolerances(
                1e-13, 1e-13,
                gw.jvm.org.orekit.orbits.OrbitType.EQUINOCTIAL,
                gw.jvm.org.orekit.propagation.Propagator.DEFAULT_LAW
            )

        # 构造力模型列表
        force_models = gw.jvm.java.util.ArrayList()

        # 重力场
        if self._config.get("gravity_field", {}).get("enabled", True):
            gf_cfg = self._config["gravity_field"]
            max_degree = gf_cfg.get("max_degree", 150)
            gravity = gw.jvm.org.orekit.forces.gravity. \
                HolmesFeatherstoneAttractionModel(
                    gw.jvm.org.orekit.frames.FramesFactory.getITRF(
                        gw.jvm.org.orekit.utils.IERSConventions.IERS_2010, True
                    ),
                    gw.jvm.org.orekit.forces.gravity.potential. \
                        GravityFieldFactory.getNormalizedProvider(max_degree, 0)
                )
            force_models.add(gravity)

        # 大气阻力
        drag_cfg = self._config.get("atmospheric_drag", {})
        if drag_cfg.get("enabled", False):
            spacecraft = self._create_orekit_spacecraft(gw)
            atmosphere = self._create_atmosphere(gw, drag_cfg)
            drag = gw.jvm.org.orekit.forces.drag.DragForce(atmosphere, spacecraft)
            force_models.add(drag)

        # 太阳辐射压
        srp_cfg = self._config.get("solar_radiation_pressure", {})
        if srp_cfg.get("enabled", False):
            spacecraft = self._create_orekit_spacecraft(gw)
            shadow_model = srp_cfg.get("shadow_model", "conical")
            # ... 构造 SRP 力模型
            force_models.add(srp)

        # ... 经验加速度、N体、潮汐等

        propagator = gw.jvm.org.orekit.propagation.numerical. \
            NumericalPropagator(integrator)
        propagator.setForceModels(force_models)

        return propagator

    def _create_orekit_spacecraft(self, gw):
        """
        从 MacroModel 构造 Orekit BoxAndSolarArraySpacecraft

        映射：
          panels[i].area × normal_vector → Orekit Vector3D facet
          panels[i].optical.{ρ_s, ρ_d}   → Orekit double[] coefficients
        """
        cfg = self._orekit_spacecraft_config
        mass = cfg["satellite_mass_kg"]
        panels = cfg["panels"]

        # 构造6个本体面 + 太阳翼
        # BoxAndSolarArraySpacecraft 构造器参数：
        #   (mass, +X_area*n, -X_area*n, +Y_area*n, -Y_area*n,
        #    +Z_area*n, -Z_area*n,
        #    +X_ρs, +X_ρd, -X_ρs, -X_ρd, ..., 
        #    solar_array_area*n, solar_array_ρs, solar_array_ρd)
        # 注意：法向量带符号表示方向，面积取绝对值

        def panel_to_vec3d(panel, sign=1.0):
            a = panel["area_m2"]
            n = panel["normal_vector_sbf"]
            # Orekit: 面板贡献 = area * normal_direction
            return gw.jvm.org.hipparchus.geometry.euclidean.threed.Vector3D(
                sign * a * n[0], sign * a * n[1], sign * a * n[2]
            )

        # 找6个本体面
        body_panels = {
            "+X": None, "-X": None,
            "+Y": None, "-Y": None,
            "+Z": None, "-Z": None,
        }
        for p in panels:
            n = p["normal_vector_sbf"]
            if abs(n[0] - 1.0) < 0.01: body_panels["+X"] = p
            elif abs(n[0] + 1.0) < 0.01: body_panels["-X"] = p
            elif abs(n[1] - 1.0) < 0.01: body_panels["+Y"] = p
            elif abs(n[1] + 1.0) < 0.01: body_panels["-Y"] = p
            elif abs(n[2] - 1.0) < 0.01: body_panels["+Z"] = p
            elif abs(n[2] + 1.0) < 0.01: body_panels["-Z"] = p

        # 构造 BoxAndSolarArraySpacecraft
        sc = gw.jvm.org.orekit.forces.BoxAndSolarArraySpacecraft(
            mass,
            panel_to_vec3d(body_panels["+X"]),
            panel_to_vec3d(body_panels["-X"]),
            panel_to_vec3d(body_panels["+Y"]),
            panel_to_vec3d(body_panels["-Y"]),
            panel_to_vec3d(body_panels["+Z"]),
            panel_to_vec3d(body_panels["-Z"]),
            body_panels["+X"]["reflectivity_coeff"],
            body_panels["+X"]["diffuse_coeff"],
            body_panels["-X"]["reflectivity_coeff"],
            body_panels["-X"]["diffuse_coeff"],
            body_panels["+Y"]["reflectivity_coeff"],
            body_panels["+Y"]["diffuse_coeff"],
            body_panels["-Y"]["reflectivity_coeff"],
            body_panels["-Y"]["diffuse_coeff"],
            body_panels["+Z"]["reflectivity_coeff"],
            body_panels["+Z"]["diffuse_coeff"],
            body_panels["-Z"]["reflectivity_coeff"],
            body_panels["-Z"]["diffuse_coeff"],
            # 太阳翼
            12.0,  # solar array area
            gw.jvm.org.hipparchus.geometry.euclidean.threed.Vector3D(
                0.0, 1.0, 0.0  # solar array normal (跟踪太阳时动态)
            ),
            0.05, 0.15  # solar array ρ_s, ρ_d
        )

        return sc

    def _create_atmosphere(self, gw, drag_cfg):
        """根据配置创建大气密度模型 (Orekit Java)"""
        model_name = drag_cfg.get("density_model", "NRLMSISE-00")
        if model_name == "DTM-2013" or model_name == "DTM2000":
            return gw.jvm.org.orekit.forces.drag.atmosphere.DTM2000(
                gw.jvm.org.orekit.forces.drag.atmosphere.data. \
                    MarshallSolarActivityFutureEstimation(
                        "Jan2000F10.txt"
                    ),
                gw.jvm.org.orekit.bodies.CelestialBodyFactory.getSun()
            )
        else:
            # NRLMSISE-00
            return gw.jvm.org.orekit.forces.drag.atmosphere.NRLMSISE00(
                gw.jvm.org.orekit.forces.drag.atmosphere.data. \
                    MarshallSolarActivityFutureEstimation(
                        "Jan2000F10.txt"
                    ),
                gw.jvm.org.orekit.bodies.CelestialBodyFactory.getSun()
            )

    # ================================================================
    # DynamicsModel 接口实现
    # ================================================================

    def compute_acceleration(self, t: float, state: np.ndarray,
                             params: Dict) -> np.ndarray:
        """委托给 Orekit propagator 计算加速度"""
        # 将 Python state 转为 Orekit SpacecraftState
        orekit_state = self._python_to_orekit_state(t, state)

        # 更新 Cd/Cr 等参数
        if "Cd" in params:
            self._update_parameter("Cd", params["Cd"])
        if "Cr" in params:
            self._update_parameter("Cr", params["Cr"])

        # 调用 Orekit 计算
        derivatives = self._propagator.computeDerivatives(orekit_state)

        # 提取加速度 (ECI)
        acc = derivatives.getAdditionalDerivative("acceleration")
        return np.array([acc.getX(), acc.getY(), acc.getZ()])
```

### 5.5 Orekit 专用动力学配置文件 (`config/orekit_dynamics.yaml`)

虽然通用的 `config/dynamics.yaml`（Section 6）已涵盖基本动力学参数，但 Orekit 后端需要一些额外配置项（Java类路径、JVM参数、Orekit数据包路径等）。建议使用独立文件或主配置文件的 `orekit:` 段：

```yaml
# config/orekit_dynamics.yaml
# ============================================================
# Orekit 后端专用配置 — 覆盖/扩展 dynamics.yaml 中与 Orekit 相关的设置
# 此文件仅在 dynamics_backend: "orekit" 时被读取
# ============================================================

# ---- JVM 配置 ----
jvm:
  java_home: "auto"                    # "auto" 使用 JAVA_HOME 环境变量
  max_heap_mb: 4096                    # JVM 最大堆内存
  orekit_data_path: "./orekit-data"    # Orekit 数据包路径（重力场/EOP/海潮等）
  # orekit-data 下载: https://gitlab.orekit.org/orekit/orekit-data

# ---- 力模型开关 ----
force_models:
  # 重力场
  gravity:
    enabled: true
    model: "EIGEN-6S4"                     # 模型名（Orekit自动从数据包加载）
    max_degree: 150                         # 截断阶数
    max_order: 150                          # 截断次数
    use_tide_system: true                   # 使用永久潮汐系统 (zero_tide / tide_free)

  # N体摄动 — Orekit 使用 CelestialBodyFactory
  third_body:
    sun: true
    moon: true
    planets: ["Jupiter", "Venus"]          # 可选额外大行星

  # 固体潮
  solid_tides:
    enabled: true
    conventions: "IERS_2010"

  # 海潮
  ocean_tides:
    enabled: true
    model: "FES2014b"                      # Orekit 通过 OceanTides 类加载
    max_degree: 50

  # 大气阻力
  drag:
    enabled: true
    atmosphere_model: "DTM2000"            # DTM2000 / NRLMSISE00 / HarrisPriester
    spacecraft_model: "box_and_solar_array" # box_and_solar_array / isotropic
    estimate_Cd: true
    Cd_apriori: 2.2
    Cd_min: 1.5
    Cd_max: 3.0
    # 太阳活动数据
    solar_activity_file: "./orekit-data/SpaceWeather/Jan2000F10.txt"

  # 太阳辐射压
  srp:
    enabled: true
    shadow_model: "conical"                # conical / cylindrical / none
    spacecraft_model: "box_and_solar_array"
    estimate_Cr: true
    Cr_apriori: 1.3
    Cr_min: 0.8
    Cr_max: 2.0

  # 经验加速度 — Orekit ParametricAcceleration
  empirical:
    enabled: true
    direction: "RTN"                       # RTN / TNW / QSW
    segment_length_sec: 360                # 分段时长 (6min for delayed)
    harmonic_orders: [0]                   # [0] = 常数; [0,1] = 常数+一次谐波
    constrain_RTN_m_s2: [5e-10, 1e-9, 1e-9]

  # 相对论效应
  relativity:
    enabled: true
    # Orekit 使用 Relativity 力模型，自动计算 Schwarzschild 项

  # 极潮 / 海洋极潮
  pole_tide:
    enabled: true
  ocean_pole_tide:
    enabled: true

# ---- 轨道参数化 ----
orbit:
  type: "EQUINOCTIAL"                     # CARTESIAN / KEPLERIAN / EQUINOCTIAL / CIRCULAR
  # EQUINOCTIAL 避免 e=0/i=0 奇点，推荐

# ---- 数值积分器 ----
integrator:
  type: "DormandPrince853"                # DormandPrince853 / GaussJackson / GraggBulirschStoer
  min_step_sec: 0.001
  max_step_sec: 300.0
  position_tolerance_m: 1e-3              # 容差 (变步长时)
  # GaussJackson 固定步长时
  fixed_step_sec: 60.0

# ---- 坐标系 ----
frames:
  inertial: "EME2000"                      # EME2000 / ICRF / GCRF
  terrestrial: "ITRF"                      # ITRF-2014 / ITRF-2020
  conventions: "IERS_2010"

# ---- 参数估计 ----
estimation:
  # Orekit 的 ParameterDriver 列表
  estimated_parameters:
    - "initial_orbit_6d"                   # 初始轨道6参数
    - "Cd"                                 # 阻力系数
    - "Cr"                                 # 光压系数
    - "empirical_RTN"                      # 经验加速度参数组
```

### 5.6 配置文件层次关系

```
用户配置层:
  ┌─────────────────────┐
  │  pod_config.yaml    │  ← 顶层配置：选择后端、模式、GNSS系统
  │  dynamics_backend:  │
  │    "self" / "orekit"│
  └──────┬──────────────┘
         │
    ┌────▼───────────────┐
    │ dynamics.yaml      │  ← 力学模型参数（与后端无关的通用参数）
    │ (Section 6.1)      │     重力场阶数、Cd/Cr先验值、积分步长等
    └──────┬──────────────┘
         │
    ┌────▼───────────────┐
    │ macro_model.yaml   │  ← 卫星几何-光学参数（两端共享）
    │ (Section 4.X.2)    │     面板面积、法向量、ρ_s、ρ_d、质量
    └──────┬──────────────┘
         │
    ┌────▼───────────────┐
    │ orekit_dynamics.   │  ← Orekit 专用扩展配置（仅在Orekit后端时加载）
    │ yaml (Section 5.5) │     JVM参数、Orekit类名、Orekit数据包路径
    └────────────────────┘

两种后端的加载流程：

  [自研后端]
    load("macro_model.yaml") → MacroModel 实例
    load("dynamics.yaml")    → 传给 SelfDynamicsModel.__init__()
    → 自研代码直接调用 mm.compute_drag_acceleration() 等

  [Orekit后端]
    load("macro_model.yaml") → MacroModel 实例
    load("dynamics.yaml")    → 传给 OrekitDynamicsAdapter.__init__()
    load("orekit_dynamics.yaml") → Orekit专用配置覆盖
    → adapter 调用 mm.to_orekit_config() 转为 Java 对象
    → 通过 py4j 构造 Orekit NumericalPropagator
```

### 5.7 配置切换

```yaml
# pod_config.yaml — 顶层只需要一个开关
global:
  dynamics_backend: "self"    # "self" = 纯Python自研 (默认)
                              # "orekit" = Orekit Java后端

# 当 dynamics_backend == "self" 时：
#   - 加载 dynamics.yaml → SelfDynamicsModel
#   - 加载 macro_model.yaml → MacroModel
#   - 忽略 orekit_dynamics.yaml

# 当 dynamics_backend == "orekit" 时：
#   - 加载 dynamics.yaml + orekit_dynamics.yaml → OrekitDynamicsAdapter
#   - 加载 macro_model.yaml → MacroModel
#   - macro_model.to_orekit_config() → Orekit BoxAndSolarArraySpacecraft
```

### 5.8 Orekit 的 BoxAndSolarArraySpacecraft 局限性

Orekit 内置的 `BoxAndSolarArraySpacecraft` 仅支持 **6面星本体 + 1个太阳翼**，对于SAR卫星有以下不足：

| 局限 | 影响 | 解决方案 |
|------|------|---------|
| 仅6个本体面 | 无法建模SAR天线正面+背面 | ① 将SAR天线面积合并到-Z面（光学系数等效）② 自定义扩展 `GenericBoxAndSolarArraySpacecraft` |
| 仅1个太阳翼 | 无法区分+Solar Wing和-Solar Wing | 将两面光学系数加权平均 |
| 太阳翼固定法向量 | 无法动态跟踪太阳 | 通过外部传入 `solar_wing_angle` 并实时更新 |

对于精度要求最高的延时模式，建议使用 **方案②**：扩展 Orekit 的 Spacecraft 接口以支持任意数量面板，或使用 `GenericSpacecraft` + 自定义 ForceModel。简化方案（方案①）在近实时模式下足以满足精度需求。

---

## 6. 配置系统设计

### 6.1 配置文件结构（YAML）

```yaml
# SAR_POD_Config.yaml
# ============================================================
# 遥感SAR卫星精密定轨系统配置文件
# ============================================================

# ---- 全局设置 ----
global:
  satellite_id: "SAR_SAT_01"
  satellite_mass_kg: 1200.0
  processing_mode: "nrt"          # "realtime" / "nrt" / "delayed"
  dynamics_backend: "self"        # "self" / "orekit"
  gnss_systems: ["GPS", "BDS"]  # ["GPS"] / ["BDS"] / ["GPS", "BDS"]
  use_joint_pod: true            # GPS+BDS联合定轨

# ---- 数据路径 ----
data_paths:
  obs_file: "data/rinex/obs/GRaceFO_2024015_0000.nav"
  nav_file: "data/rinex/nav/brdc0010.24n"
  sp3_file: "data/precise/igs22641.sp3"   # 近实时/延时模式
  clk_file: "data/precise/igs22641.clk"
  antex_file: "data/antenna/igs20.atx"
  dcb_file: "data/bias/P1P2110100.DCB"
  erp_file: "data/eop/iers14_C04.erp"
  attitude_file: "data/attitude/gracefo_att_2024015.txt"
  macro_model_file: "data/satellite/sar_macro_model.json"
  de_ephemeris: "data/ephemeris/de440.bsp"

# ---- 观测预处理 ----
obs_preprocessing:
  sampling_interval_sec: 30       # 降采样间隔（近实时=30s，延时=10s）
  elevation_cutoff_deg: 10.0      # 高度角截止（近实时=10°，延时=5°）
  cycle_slip_gf_threshold_cycles: 0.05
  cycle_slip_mw_threshold_cycles: 0.1
  outlier_rejection_sigma: 3.0     # 残差3σ剔除
  use_phase: true                  # 是否使用载波相位
  use_pseudorange: true            # 是否使用伪距
  phase_weight_factor: 1.0        # 相位权重相对伪距的倍数

# ---- 力学模型 ----
dynamics:
  # 地球重力场
  gravity_field:
    model: "EIGEN-6S4"
    max_degree: 80                 # 实时=2, 近实时=50-80, 延时=150
    use_tidal_variations: true     # 时变重力C20-C40
    tidal_model: "FES2014b"
  
  # N体摄动
  n_body:
    use_sun: true
    use_moon: true
    use_planets: false             # 大行星（近实时=false，延时=true）
    ephemeris: "DE440"
  
  # 潮汐
  tides:
    solid_earth_tides: true
    ocean_tides: true
    ocean_tide_model: "FES2014b"  # 近实时=FES2004, 延时=FES2014b
    pole_tide: true
    ocean_pole_tide: true
  
  # 大气阻力
  atmospheric_drag:
    enabled: true
    density_model: "NRLMSISE-00"  # 或 "DTM-2013"
    estimate_Cd: true              # 是否估计阻力系数Cd
    Cd_apriori: 2.2
    Cd_apriori_sigma: 0.5         # 先验标准差
    Cd_process_noise: 1e-8        # 过程噪声 (m²/s³)
    macro_model_file: "data/satellite/sar_macro_model.json"
  
  # 太阳辐射压
  solar_radiation_pressure:
    enabled: true
    shadow_model: "conical"        # "cylindrical" 或 "conical"
    estimate_Cr: true
    Cr_apriori: 1.3
    Cr_apriori_sigma: 0.3
    Cr_process_noise: 1e-8
  
  # 经验加速度 (RTN)
  empirical_acceleration:
    enabled: true                  # 实时=false，近实时/延时=true
    segment_length_sec: 900        # 15min（近实时）或360（6min，延时）
    constrain_RTN: [1e-9, 1e-9, 1e-9]  # RTN方向先验过程噪声 (m/s²)
    adaptive_segmentation: true     # 地影区自动缩短分段

### 6.X 动力学参数可配置表 (Dynamics Parameter Reference)

> 下表列出所有可配置的动力学参数，含默认值、有效范围、单位和三种处理模式的推荐值。
> 参数通过 `config/dynamics.yaml` 或 `config/pod_config.yaml` 的 `dynamics` 段配置。
> 标记 ★ 的参数为**待估参数**（在最小二乘/滤波中解算），其余为固定配置。

#### 6.X.1 重力场参数

| 参数名 | 键路径 | 类型 | 默认值 | 有效范围 | 单位 | 实时 | 近实时 | 延时 | 说明 |
|--------|--------|------|--------|---------|------|:----:|:-----:|:----:|------|
| 重力场模型 | `gravity_field.model` | string | `EIGEN-6S4` | EGM2008, EIGEN-6S4, GOCO06s, XGM2019e | — | J2 only | EIGEN-6S4 | EIGEN-6S4 | 静态重力场球谐模型 |
| 最大阶数 | `gravity_field.max_degree` | int | 80 | 2–2190 | — | 2 | 50–80 | 150 | 球谐展开截断阶数；实时只保留J2 |
| 时变重力 | `gravity_field.use_tidal_variations` | bool | true | true/false | — | false | true | true | C20/C30/C40时变项（来自IERS规范） |
| 海潮模型 | `gravity_field.tidal_model` | string | `FES2014b` | FES2004, FES2014b, GOT4.10c | — | — | FES2004 | FES2014b | 海潮负荷形变模型 |

##### 重力场模型选择策略

| 模型 | 完整阶数 | 50×50精度 | 数据基础 | 推荐场景 |
|------|:-------:|----------|---------|---------|
| **EIGEN-6S4** | 2190 | ~3 cm (GEOID) | GRACE + GOCE + LAGEOS + 地面 | 延时、科研级POD |
| **GOCO06s** | 300 | ~2 cm (GEOID) | GRACE + GOCE + SLR | 延时（强调时变精度） |
| **EGM2008** | 2190 | ~5 cm (GEOID) | GRACE + 地面 + 测高 | 兼容性最好，广泛使用 |
| **XGM2019e** | 760 | ~3 cm (GEOID) | GOCO06s + 地面 + 测高 | 近实时推荐 |

#### 6.X.2 N体摄动 & 潮汐参数

| 参数名 | 键路径 | 类型 | 默认值 | 有效范围 | 单位 | 实时 | 近实时 | 延时 | 说明 |
|--------|--------|------|--------|---------|------|:----:|:-----:|:----:|------|
| 太阳引力 | `n_body.use_sun` | bool | true | true/false | — | true | true | true | 太阳第三体引力摄动 |
| 月球引力 | `n_body.use_moon` | bool | true | true/false | — | true | true | true | 月球第三体引力摄动 |
| 大行星 | `n_body.use_planets` | bool | false | true/false | — | false | false | true | Jupiter/Venus/Mars/Saturn等大行星 |
| 行星历表 | `n_body.ephemeris` | string | `DE440` | DE430, DE440, INPOP | — | — | DE440 | DE440 | JPL行星历表 |
| 固体潮 | `tides.solid_earth_tides` | bool | true | true/false | — | false | true | true | IERS 2010固体潮模型 |
| 海潮 | `tides.ocean_tides` | bool | true | true/false | — | false | true | true | 海潮负荷形变 |
| 海潮模型 | `tides.ocean_tide_model` | string | `FES2014b` | FES2004, FES2014b | — | — | FES2004 | FES2014b | 海潮球谐或网格模型 |
| 极潮 | `tides.pole_tide` | bool | true | true/false | — | false | false | true | 极移固体潮 |
| 海洋极潮 | `tides.ocean_pole_tide` | bool | true | true/false | — | false | false | true | 海洋极潮负荷 |

#### 6.X.3 大气阻力参数

| 参数名 | 键路径 | 类型 | 默认值 | 有效范围 | 单位 | 实时 | 近实时 | 延时 | 说明 |
|--------|--------|------|--------|---------|------|:----:|:-----:|:----:|------|
| 启用阻力 | `atmospheric_drag.enabled` | bool | true | true/false | — | false | true | true | 是否启用大气阻力建模 |
| 大气密度模型 | `atmospheric_drag.density_model` | string | `NRLMSISE-00` | NRLMSISE-00, DTM-2013, JB2008 | — | — | NRLMSISE-00 | DTM-2013 | 热层大气密度模型；DTM-2013精度更高 |
| ★ 估计Cd | `atmospheric_drag.estimate_Cd` | bool | true | true/false | — | — | true | true | 是否将Cd作为待估参数 |
| ★ Cd先验值 | `atmospheric_drag.Cd_apriori` | float | 2.2 | 1.5–3.0 | — | — | 2.2 | 2.2 | 阻力系数先验值；SAR卫星典型范围2.0-2.5 |
| ★ Cd先验σ | `atmospheric_drag.Cd_apriori_sigma` | float | 0.5 | 0.1–1.0 | — | — | 0.5 | 0.3 | Cd先验标准差 |
| ★ Cd过程噪声 | `atmospheric_drag.Cd_process_noise` | float | 1e-8 | 1e-10–1e-6 | m²/s³ | — | 1e-8 | 1e-8 | Cd随机游走过程噪声（EKF模式） |
| 宏模型文件 | `atmospheric_drag.macro_model_file` | path | `data/satellite/sar_macro_model.json` | 任意路径 | — | — | ✓ | ✓ | 卫星Box-Wing宏模型JSON路径 |
| F10.7太阳通量 | `atmospheric_drag.f107_source` | string | `auto` | auto, file, constant | sfu | — | auto | file | 太阳10.7cm射电通量（驱动大气密度模型） |
| 地磁指数 | `atmospheric_drag.ap_source` | string | `auto` | auto, file, constant | — | — | auto | file | Ap地磁活动指数 |

##### 大气密度模型对比

| 模型 | 精度 (1σ) | 适用高度 | 输入参数 | 推荐场景 |
|------|:--------:|---------|---------|---------|
| **NRLMSISE-00** | ~15-20% | 0–1000 km | F10.7, Ap, DOY, UT, Lat, Lon, Alt | 通用，NRT首选 |
| **DTM-2013** | ~10-15% | 120–1500 km | F10.7, Kp, DOY, UT, Lat, Lon, Alt | 高精度需求，延时推荐 |
| **JB2008** | ~15-18% | 120–2500 km | F10.7, S10, M10, Y10, Dst | 太阳活动高年；需要更多太阳指数输入 |

#### 6.X.4 太阳辐射压 (SRP) 参数

| 参数名 | 键路径 | 类型 | 默认值 | 有效范围 | 单位 | 实时 | 近实时 | 延时 | 说明 |
|--------|--------|------|--------|---------|------|:----:|:-----:|:----:|------|
| 启用SRP | `solar_radiation_pressure.enabled` | bool | true | true/false | — | false | true | true | 是否启用SRP建模 |
| 地影模型 | `solar_radiation_pressure.shadow_model` | string | `conical` | cylindrical, conical | — | — | conical | conical | 地影模型；conical=锥形（含半影），cylindrical=柱形 |
| ★ 估计Cr | `solar_radiation_pressure.estimate_Cr` | bool | true | true/false | — | — | true | true | 是否将Cr作为待估参数 |
| ★ Cr先验值 | `solar_radiation_pressure.Cr_apriori` | float | 1.3 | 0.8–2.0 | — | — | 1.3 | 1.3 | SRP系数先验值 |
| ★ Cr先验σ | `solar_radiation_pressure.Cr_apriori_sigma` | float | 0.3 | 0.1–0.5 | — | — | 0.3 | 0.2 | Cr先验标准差 |
| ★ Cr过程噪声 | `solar_radiation_pressure.Cr_process_noise` | float | 1e-8 | 1e-10–1e-6 | m²/s³ | — | 1e-8 | 1e-8 | Cr随机游走过程噪声 |
| 太阳常数 | `solar_radiation_pressure.solar_constant` | float | 1361.0 | 1360–1362 | W/m² | — | 1361.0 | 1361.0 | 1 AU处太阳辐照度；一般不需修改 |
| 反照率建模 | `solar_radiation_pressure.albedo_enabled` | bool | false | true/false | — | false | false | true | 地球反照率辐射压（延时可选，量级 ~10⁻⁹ m/s²） |
| IR辐射建模 | `solar_radiation_pressure.ir_enabled` | bool | false | true/false | — | false | false | false | 地球红外辐射压（量级极低，一般关闭） |

#### 6.X.5 经验加速度 (Empirical Acceleration — Reduced-Dynamic核心)

| 参数名 | 键路径 | 类型 | 默认值 | 有效范围 | 单位 | 实时 | 近实时 | 延时 | 说明 |
|--------|--------|------|--------|---------|------|:----:|:-----:|:----:|------|
| 启用经验加速度 | `empirical_acceleration.enabled` | bool | true | true/false | — | false | true | true | Reduced-Dynamic模式核心开关 |
| ★ 分段长度 | `empirical_acceleration.segment_length_sec` | int | 900 | 180–1800 | s | — | 900 (15min) | 360 (6min) | 经验加速度分段常数时长；延时更短→更多参数→更高精度 |
| ★ R方向约束 | `empirical_acceleration.constrain_RTN[0]` | float | 1e-9 | 1e-10–1e-7 | m/s² | — | 1e-9 | 5e-10 | 径向经验加速度先验约束 |
| ★ T方向约束 | `empirical_acceleration.constrain_RTN[1]` | float | 1e-9 | 1e-10–1e-7 | m/s² | — | 1e-9 | 1e-9 | 切向经验加速度先验约束（T方向最弱约束） |
| ★ N方向约束 | `empirical_acceleration.constrain_RTN[2]` | float | 1e-9 | 1e-10–1e-7 | m/s² | — | 1e-9 | 1e-9 | 法向经验加速度先验约束 |
| 自适应分段 | `empirical_acceleration.adaptive_segmentation` | bool | true | true/false | — | — | true | true | 地影区自动缩短分段长度（通常半分段） |

##### 分段长度与参数数量关系

| 弧段长度 | 分段时长 | 分段数 | RTN参数总数 (×3) | 时间分辨率 | 存储/计算开销 | 推荐场景 |
|---------|:-------:|:-----:|:---------------:|:---------:|:-----------:|---------|
| 2 h (NRT) | 15 min | 8 | 24 | 低 | 低 | **近实时**（2min处理时限） |
| 2 h (NRT) | 10 min | 12 | 36 | 中 | 中 | 近实时（精度优先） |
| 24 h (Delayed) | 6 min | 240 | 720 | 高 | 高 | **延时**（5min处理时限） |
| 24 h (Delayed) | 15 min | 96 | 288 | 中 | 中 | 延时（快速模式） |
| 24 h (Delayed) | 3 min | 480 | 1440 | 最高 | 最高 | 延时（科研级，超限） |

##### 经验加速度约束策略说明

```
约束公式（用于最小二乘正则化或EKF过程噪声）：

  E[ΔaR²] = σ_R² · Δt    (径向约束最强—轨道力学径向约束天然强)
  E[ΔaT²] = σ_T² · Δt    (切向约束最弱—沿轨方向摄动最复杂)
  E[ΔaN²] = σ_N² · Δt    (法向约束中等)

其中 Δt = segment_length_sec

典型约束关系: σ_R ≈ 0.5·σ_N ≈ 0.25·σ_T
即 constrain_RTN = [5e-10, 2e-9, 1e-9] (延时)
                  [1e-9,  4e-9, 2e-9] (近实时)
```

#### 6.X.6 数值积分参数

| 参数名 | 键路径 | 类型 | 默认值 | 有效范围 | 单位 | 实时 | 近实时 | 延时 | 说明 |
|--------|--------|------|--------|---------|------|:----:|:-----:|:----:|------|
| 积分器类型 | `integrator.type` | string | `GaussJackson` | RK4, DOPRI8, GaussJackson | — | RK4 | GaussJackson | GaussJackson | Gauss-Jackson多步法精度最高 |
| 积分步长 | `integrator.step_size_sec` | float | 60.0 | 10–300 | s | 30 | 60 | 60 | 固定步长积分器步长 |
| 相对容差 | `integrator.integ_rel_tol` | float | 1e-13 | 1e-10–1e-16 | — | — | — | 1e-13 | 变步长积分器相对容差 |
| 绝对容差 | `integrator.integ_abs_tol` | float | 1e-13 | 1e-10–1e-16 | — | — | — | 1e-13 | 变步长积分器绝对容差 |
| 变步长 | `integrator.var_step` | bool | false | true/false | — | true | false | false | true=DOPRI8变步长, false=GaussJackson固定步长 |
| 轨道类型 | `integrator.orbit_type` | string | `cartesian` | cartesian, equinoctial, keplerian | — | cartesian | equinoctial | equinoctial | 轨道参数化方式；equinoctial避免奇点 |

#### 6.X.7 相对论效应参数

| 参数名 | 键路径 | 类型 | 默认值 | 有效范围 | 单位 | 实时 | 近实时 | 延时 | 说明 |
|--------|--------|------|--------|---------|------|:----:|:-----:|:----:|------|
| 广义相对论 | `relativity.enabled` | bool | true | true/false | — | false | true | true | Schwarzschild项（1/c²量级） |
| Lense-Thirring | `relativity.lense_thirring` | bool | false | true/false | — | false | false | false | 惯性系拖曳效应（1/c³量级；对LEO ~mm级，通常忽略） |

#### 6.X.8 总待估参数维度汇总

以典型24h弧段为例，各模式的待估参数数量：

| 参数类别 | 每项数量 | 实时 (EKF) | 近实时 (2h NRT) | 延时 (24h Final) |
|---------|:-------:|:----------:|:---------------:|:----------------:|
| ★ 初始轨道状态 | 6 | 6 | 6 | 6 |
| ★ Cd (阻力系数) | 1 | — | 1 | 1 |
| ★ Cr (光压系数) | 1 | — | 1 | 1 |
| ★ 接收机钟差 | 1/历元 | ~720 | ~720 | ~8640 |
| ★ 载波相位模糊度 | 1/sat/弧段 | ~8 | ~8 | ~8 |
| ★ 经验加速度 (RTN×N) | 3N | — | 24 (15min×8) | 720 (6min×240) |
| **总参数** | — | **~734** | **~760** | **~9376** |

> 延时模式参数数量最多（~1万），但得益于更长的数据弧段（24h vs 2h），观测方程条件数更好。

# ---- 数值积分 ----
integrator:
  type: "GaussJackson"            # "RK4" / "DOPRI8" / "GaussJackson"
  step_size_sec: 60.0             # 积分步长
  integ_rel_tol: 1e-13
  integ_abs_tol: 1e-13
  var_step: false                 # 固定步长(GaussJackson) 或 变步长(DOPRI8)

# ---- 参数估计 ----
estimator:
  type: "LSQ"                     # "LSQ" (近实时/延时) 或 "EKF" (实时)
  max_iterations: 4
  convergence_threshold: 1e-8
  lambda_reg: 1e-6               # Levenberg-Marquardt阻尼（可选）
  
  # 待估参数列表（根据模式和配置动态生成）
  estimated_parameters:
    - "initial_state"              # 初始位置速度 (6个参数)
    - "receiver_clock_bias"         # 接收机钟差 (每历元1个参数 或 随机游走)
    - "carrier_phase_bias_L1"      # 相位模糊度 (每卫星1个参数)
    - "carrier_phase_bias_L2"
    - "CD_drag_coeff"              # 大气阻力系数
    - "CR_srp_coeff"              # 光压系数
    - "emp_aR_seg*"                # 经验加速度 (动态数量)
  
  # 权重
  obs_weights:
    pseudorange_L1_sigma_m: 3.0   # 伪距L1权重（标准差米）
    pseudorange_L2_sigma_m: 3.0
    phase_L1_sigma_cycles: 0.02   # 相位权重（标准差周）
    phase_L2_sigma_cycles: 0.02

# ---- 模糊度固定 (PPP-AR) ----
ambiguity_resolution:
  enabled: true                   # 近实时=false（浮点解），延时=true
  method: "LAMBDA"              # 唯一支持的方法
  ratio_threshold: 3.0           # Ratio检验阈值
  max_amb_search: 99             # 最大搜索模糊度数量
  fix_mode: "partial"             # "partial"(部分固定) 或 "full"(全固定)

# ---- 输出 ----
output:
  sp3_output_file: "output/orbit/sp3/gr1a0010.26l"
  residual_file: "output/residuals/gr1a0010_residuals.csv"
  covariance_file: "output/covariance/gr1a0010_cov.bin"
  log_file: "output/logs/gr1a0010.log"
  report_file: "output/reports/gr1a0010_report.json"
  write_interval_sec: 900         # SP3输出间隔（秒）

# ---- 性能/时效 ----
performance:
  parallel_arc_segments: 4        # 并行弧段数（近实时=4，延时=12）
  num_threads: 4                  # 并行线程数
  use_gpu_acceleration: false     # 未来扩展：GPU加速力学模型
  max_processing_time_sec: 120    # 近实时=120s，延时=300s
```

---

## 7. 文件存储与数据管理

> **设计原则：能不用数据库就不用。** POD系统是以文件为中心的科学计算软件，核心数据是GNSS观测文件、星历文件、轨道结果文件等。这些数据的天然存储形式是文件系统，而非数据库。保持简单，减少组件，提高可靠性。

### 7.1 为什么不用数据库？

POD 系统本质上是一个**以文件为中心的科学计算软件**。输入端是 GNSS 观测文件（RINEX）、精密星历（SP3）、钟差文件（CLK）等标准格式文件；输出端是轨道结果文件（orb）、残差文件等。这些数据的天然存储形式是文件系统，引入数据库反而增加了不必要的复杂度。

| 考量维度 | 数据库方案 | 文件系统方案 | 结论 |
|---------|-----------|------------|------|
| **科学计算的天然数据流** | 文件→解析→入库→查询→出库→计算 | 文件→直接读取→计算→直接写出 | 文件系统更直接 |
| **依赖复杂度** | 需增加 SQLite/PostgreSQL 依赖 | 零外部依赖（Python 标准库） | 文件系统胜出 |
| **可复现性** | 需要从数据库导出数据才能复现 | `runs/` 目录自带完整快照 | 文件系统更可靠 |
| **调试便利性** | 需 SQL 查询或 GUI 工具查看数据 | `ls`/`cat`/文本编辑器直接查看 | 文件系统更直观 |
| **部署复杂度** | 需初始化数据库、管理迁移 | 只需目录结构，git clone 即用 | 文件系统更简单 |
| **版本管理** | 数据库 dump 文件，diff 困难 | 标准 Git 版本控制，diff/merge 自然 | 文件系统更灵活 |
| **团队协作** | 数据库需额外同步机制 | Git 直接管理配置和索引，数据文件通过共享存储 | 文件系统更轻量 |

**核心原则：**

> POD 软件的价值在于科学计算（轨道确定算法），而非数据管理系统。每增加一个外部依赖，就增加一个潜在的故障点。对于 MVP 阶段和中等规模数据处理（日处理数十弧段），文件系统（目录结构 + JSON 索引）完全足够。

### 7.2 目录结构规范

软件遵循"约定优于配置"原则，固定目录结构如下：

```
sar_pod/
├── config/                         # 配置文件（Git 管理）
│   ├── pod_config.yaml             # 顶层配置：后端选择、GNSS系统、处理模式
│   ├── dynamics.yaml               # 力学模型参数
│   └── satellites/
│       ├── sentinel1.yaml          # 卫星配置：质量/面积/姿态
│       └── gf3.yaml
│
├── data/                           # 本地数据仓库（不纳入 Git）
│   ├── rinex/                      # GNSS 观测文件
│   │   └── {YYYY}/
│   │       └── {DOY}/
│   │           └── {satellite}/
│   ├── sp3/                        # 精密星历
│   ├── clk/                        # 精密钟差
│   ├── erp/                        # 地球定向参数
│   ├── ant/                        # 天线相位中心文件 (ANTEX)
│   └── dcbs/                       # 码偏差文件 (DCB)
│
├── runs/                           # 处理运行记录（本节核心）
│   └── {YYYYMMDD_HHMMSS}/          # 每次运行的不可变记录
│       ├── inputs/                 # 本次使用的输入文件副本（符号链接）
│       ├── configs/                # 本次使用的配置文件副本
│       ├── outputs/                # 本次产出的轨道文件、残差文件
│       ├── logs/                   # 本次运行日志
│       └── run_info.json           # 运行元数据：时间、卫星、弧段、状态
│
├── results/                        # 发布级结果汇总
│   └── {satellite}/
│       └── {YYYY}/
│           └── {DOY}/
│               └── orbits/         # 最终轨道产品
│
├── data_index.json                 # 全局数据索引（替代数据库查询）
│
├── src/                            # 源代码
│   ├── pod_engine/
│   ├── io/
│   ├── data_manager/
│   └── utils/
│
└── tests/
```

**关键约定：**
- `config/` 纳入 Git，保证配置可追溯
- `data/` 和 `runs/` 不纳入 Git（通过 `.gitignore`），数据量太大
- `results/` 中的最终产品纳入 Git LFS 或共享存储
- 所有路径使用 `pathlib.Path`，跨平台兼容

### 7.3 文件命名约定

统一采用以下命名规范，使文件名能自描述数据属性：

| 文件类型 | 命名模板 | 示例 | 说明 |
|---------|---------|------|------|
| RINEX 观测 | `{SAT}_{YYYY}{DOY}{HH}{MM}_{TYPE}.rnx` | `S1A_20251531200_OBS.rnx` | `TYPE` = OBS/NAV |
| SP3 精密星历 | `{AC}_{WWWW}{DOW}_{HHU}_{ORB}.sp3` | `GFZ0MGXRAP_20251530000_01D_05M_ORB.SP3` | AC=分析中心 |
| CLK 钟差 | `{AC}_{WWWW}{DOW}_{HHU}_{CLK}.clk` | `GFZ0MGXRAP_20251530000_01D_30S_CLK.CLK` | |
| ERP 参数 | `{AC}_{WWWW}{DOW}_{HHU}_{ERP}.erp` | `COD0MGXRAP_20251820000_01D_01D_ERP.ERP` | |
| ANTEX 天线 | `{AC}_{WWWW}_{HHU}_{ANT}.atx` | `IGS20_2025000.atx` | |
| 轨道输出 | `{SAT}_{YYYY}{DOY}_{MODE}.orb` | `S1A_2025153_NRT.orb` | MODE=NRT/STD/FINAL |
| 残差输出 | `{SAT}_{YYYY}{DOY}_{MODE}.res` | `S1A_2025153_NRT.res` | |

**解析工具：** 提供 `FilenameParser` 类，支持从文件名自动提取卫星、时间、类型等元数据：

```python
from sar_pod.io import FilenameParser

parser = FilenameParser("S1A_20251531200_OBS.rnx")
assert parser.satellite == "S1A"
assert parser.year == 2025
assert parser.doy == 153
assert parser.data_type == "RINEX_OBS"
```

### 7.4 `runs/` 目录设计：可复现性的核心

每次 POD 处理运行创建一个**不可变时间戳目录**，包含该次运行的全部上下文：

```
runs/
├── 20260608_003528/                # 2026-06-08 00:35:28 启动的运行
│   ├── inputs/
│   │   ├── S1A_20261591200_OBS.rnx      # → 符号链接到 data/rinex/.../
│   │   ├── GFZ0MGXRAP_...SP3            # → 符号链接到 data/sp3/.../
│   │   └── ...
│   ├── configs/
│   │   ├── pod_config.yaml              # 本次使用的配置副本
│   │   ├── dynamics.yaml
│   │   └── satellites/
│   │       └── sentinel1.yaml
│   ├── outputs/
│   │   ├── S1A_2026159_NRT.orb          # 轨道解算结果
│   │   ├── S1A_2026159_NRT.res          # 残差文件
│   │   └── residuals_plot.png           # 残差可视化
│   ├── logs/
│   │   ├── processing.log               # 完整处理日志
│   │   └── timing.csv                   # 各步骤耗时统计
│   └── run_info.json                    # 运行元数据
```

**`run_info.json` 格式：**

```json
{
  "run_id": "20260608_003528",
  "satellite": "S1A",
  "mode": "NRT",
  "arc": {"start": "2026-06-08T00:00:00", "end": "2026-06-08T02:00:00"},
  "gnss_systems": ["G", "E"],
  "pod_backend": "PRIDE",
  "config_hash": "a1b2c3d4e5f6...",
  "input_files": {
    "rinex": {"file": "S1A_20261591200_OBS.rnx", "sha256": "abc123..."},
    "sp3": {"file": "GFZ0MGXRAP_20261590000_01D_05M_ORB.SP3", "sha256": "def456..."},
    "clk": {"file": "GFZ0MGXRAP_20261590000_01D_30S_CLK.CLK", "sha256": "789abc..."}
  },
  "status": "SUCCESS",
  "start_time": "2026-06-08T00:35:28",
  "end_time": "2026-06-08T00:37:15",
  "duration_sec": 107,
  "quality": {
    "rms_phase_residual_mm": 5.3,
    "overlap_3d_rms_cm": null
  }
}
```

**优势：**
1. **完全可复现**：要复现某次运行，只需 `cd runs/20260608_003528`，所有输入和配置都在
2. **自然去重**：同一卫星连续弧段每次运行独立目录，不会相互覆盖
3. **自动存档**：`runs/` 按时间自然排序，`ls runs/ | head -20` 即可查看最近运行
4. **零查询开销**：不需要数据库就知道某次运行用了什么数据 — 直接看目录内容
5. **故障排查**：出错时直接 `cat runs/20260608_003528/logs/processing.log`

### 7.5 索引文件：替代数据库查询

当需要在大量运行记录中按时间范围查找时，遍历 `runs/` 下所有 `run_info.json` 效率较低。引入轻量级 JSON 索引文件：

**`data_index.json` 结构：**

```json
{
  "version": "1.0",
  "last_updated": "2026-06-08T00:37:15",
  "runs": [
    {
      "run_id": "20260608_003528",
      "run_dir": "runs/20260608_003528",
      "satellite": "S1A",
      "mode": "NRT",
      "arc_start": "2026-06-08T00:00:00",
      "arc_end": "2026-06-08T02:00:00",
      "status": "SUCCESS",
      "rms_residual_mm": 5.3
    },
    {
      "run_id": "20260608_020015",
      "run_dir": "runs/20260608_020015",
      "satellite": "S1A",
      "mode": "NRT",
      "arc_start": "2026-06-08T02:00:00",
      "arc_end": "2026-06-08T04:00:00",
      "status": "SUCCESS",
      "rms_residual_mm": 4.8
    }
  ]
}
```

**查询工具：** 提供 `RunIndex` 类，支持基于内存排序和过滤的快速查询：

```python
from sar_pod.data_manager import RunIndex

index = RunIndex.load("data_index.json")

# 按时间范围查询
runs = index.query(
    satellite="S1A",
    mode="NRT",
    time_range=("2026-06-08T00:00:00", "2026-06-08T04:00:00")
)

# 按状态筛选失败运行
failed = index.query(status="FAILED")
```

**性能分析：** 日处理 48 弧段（NRT 模式），一年约 17,500 条记录。JSON 文件大小约 5 MB，解析和过滤在 Python 中 < 0.1 秒 — 完全不需要数据库的 B-Tree 索引。

### 7.6 配置版本管理

利用 Git 管理 `config/` 目录，每次修改自动记录。运行 `run_info.json` 中的 `config_hash` 字段记录该次运行使用的配置版本：

```python
import hashlib
import subprocess

def get_config_hash(config_dir: str) -> str:
    """Calculate content hash of config/ directory for run reproducibility."""
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=config_dir,
        capture_output=True, text=True
    )
    if result.returncode == 0:
        return result.stdout.strip()  # Git commit hash
    # Fallback: compute file tree hash
    return compute_file_tree_hash(config_dir)
```

**工作流：**
1. 修改 `config/satellites/gf3.yaml` → `git commit -m "update GF-3 mass to 2400kg"`
2. 运行 POD → `run_info.json` 记录 `config_hash: "abc123def"`
3. 需要复现时 → `git checkout abc123def` 恢复当时配置

### 7.8 与 Section 4.1 数据对齐架构的衔接

Section 4.1 中的数据对齐架构 `DataInputManager` 与本节文件管理方案无缝配合：

```
                    ┌──────────────────────────┐
                    │    DataInputManager      │
                    │  (Section 4.1)           │
                    │                          │
                    │  load_and_align()        │
                    └──────────┬───────────────┘
                               │
              ┌────────────────┼────────────────┐
              │                │                │
              ▼                ▼                ▼
    ┌─────────────┐  ┌─────────────┐  ┌─────────────┐
    │ RINEX 文件   │  │ SP3/CLK 文件 │  │ ERP/ANT 文件 │
    │ data/rinex/  │  │ data/sp3/   │  │ data/erp/   │
    └─────────────┘  └─────────────┘  └─────────────┘
              │                │                │
              └────────────────┼────────────────┘
                               │
                               ▼
                    ┌──────────────────────────┐
                    │   AlignedDataBuffer      │
                    │  (内存中的对齐数据)        │
                    └──────────┬───────────────┘
                               │
                               ▼
                    ┌──────────────────────────┐
                    │     POD Engine           │
                    └──────────┬───────────────┘
                               │
                               ▼
                    ┌──────────────────────────┐
                    │   runs/{run_id}/          │
                    │   ├── inputs/  (符号链接)  │
                    │   ├── configs/ (副本)      │
                    │   ├── outputs/ (结果)      │
                    │   └── logs/    (日志)      │
                    └──────────┬───────────────┘
                               │
                               ▼
                    ┌──────────────────────────┐
                    │   data_index.json        │
                    │   (更新运行记录索引)        │
                    └──────────────────────────┘
```

**数据对齐流程：**

1. `DataInputManager.load_and_align()` 通过 `FilenameParser` 定位 `data/` 下的文件
2. 对齐完成后，`AlignedDataBuffer` 持有内存中的对齐数据（不经过数据库中转）
3. POD 引擎计算完成后，将本次运行的所有上下文保存到 `runs/{run_id}/`
4. 更新 `data_index.json`，追加本次运行记录

---

## 8. 软件目录结构设计

```
sar_pod/
├── config/                         # 配置文件（Git 管理）
│   ├── pod_config.yaml             # 顶层配置（选择后端、模式、GNSS系统）
│   ├── dynamics.yaml               # 力学模型参数（重力场/阻力/SRP/经验加速度）
│   ├── satellites/
│   │   ├── sentinel1.yaml        # 卫星配置：质量/面积/姿态
│   │   └── gf3.yaml
│   ├── config_realtime.yaml       # 实时模式覆盖配置
│   ├── config_nrt.yaml           # 近实时模式覆盖配置
│   └── config_delayed.yaml       # 延时模式覆盖配置
│
├── data/                           # 本地数据仓库（不纳入 Git）
│   ├── rinex/                      # GNSS 观测文件
│   │   └── {YYYY}/
│   │       └── {DOY}/
│   │           └── {satellite}/
│   ├── sp3/                        # 精密星历
│   ├── clk/                        # 精密钟差
│   ├── erp/                        # 地球定向参数
│   ├── ant/                        # 天线相位中心文件 (ANTEX)
│   └── dcbs/                       # 码偏差文件 (DCB)
│
├── runs/                           # 处理运行记录（不可变）
│   └── {YYYYMMDD_HHMMSS}/          # 每次运行的不可变记录
│       ├── inputs/                 # 本次使用的输入文件副本（符号链接）
│       ├── configs/                # 本次使用的配置文件副本
│       ├── outputs/                # 本次产出的轨道文件、残差文件
│       ├── logs/                   # 本次运行日志
│       └── run_info.json           # 运行元数据
│
├── results/                        # 发布级结果汇总
│   └── {satellite}/
│       └── {YYYY}/
│           └── {DOY}/
│               └── orbits/         # 最终轨道产品
│
├── data_index.json                 # 全局数据索引（替代数据库查询）
│
├── src/                            # 源代码
│   └── sar_pod/
│       ├── __init__.py
│       ├── main.py                  # 主程序入口
│       │
│       ├── data_loader/             # 数据加载层
│       │   ├── __init__.py
│       │   ├── rinex_reader.py     # RINEX观测/导航文件读取
│       │   ├── sp3_reader.py       # SP3精密轨道读取
│       │   ├── clk_reader.py        # CLK精密钟差读取
│       │   ├── antex_reader.py      # ANTEX天线模型读取
│       │   ├── dcb_reader.py        # DCB偏差读取
│       │   ├── erp_reader.py        # ERP地球定向参数读取
│       │   └── attitude_reader.py   # 姿态四元数读取
│       │
│       ├── macro_model/              # 卫星宏模型模块（独立共享模块）
│       │   ├── __init__.py
│       │   ├── macro_model.py        # MacroModel 核心类（加载/计算/导出）
│       │   ├── panel_def.py          # PanelDef + SatellitePhysicalParams 数据结构
│       │   └── shadow_check.py       # 地影检测（cylindrical / conical）
│       │
│       ├── data_manager/              # 数据管理模块（替代 database/）
│       │   ├── __init__.py
│       │   ├── run_index.py          # RunIndex 类：JSON 索引加载/查询/更新
│       │   ├── run_recorder.py       # RunRecorder 类：创建 runs/{run_id}/ 目录结构
│       │   ├── filename_parser.py    # FilenameParser 类：文件名→元数据解析
│       │   └── config_version.py     # 配置版本管理：Git hash + 文件哈希
│       │
│       ├── obs_preprocessing/       # 观测预处理层
│       │   ├── __init__.py
│       │   ├── cycle_slip_detector.py  # 周跳检测(GF+MW)
│       │   ├── outlier_rejector.py     # 异常值剔除(3σ)
│       │   ├── elevation_filter.py     # 高度角滤波
│       │   ├── pco_pcv_corrector.py   # 天线PCO/PCV校正(含姿态)
│       │   ├── ionospheric_corrector.py # 电离层校正(无电离层组合)
│       │   └── obs_combination.py      # 观测值组合(GF/MW/无电离层)
│       │
│       ├── dynamics/               # 力学模型层
│       │   ├── __init__.py
│       │   ├── base_dynamics.py    # DynamicsModel抽象基类
│       │   ├── two_body.py         # 二体问题
│       │   ├── j2_dynamics.py      # J2摄动
│       │   ├── gravity_field.py     # 地球重力场(球谐)
│       │   ├── n_body.py           # N体摄动
│       │   ├── solid_tide.py       # 固体潮
│       │   ├── ocean_tide.py       # 海潮
│       │   ├── atmospheric_drag.py # 大气阻力
│       │   ├── srp.py             # 太阳辐射压
│       │   ├── empirical_accel.py  # RTN经验加速度
│       │   └── relativistic.py     # 相对论效应
│       │
│       ├── integrator/             # 数值积分层
│       │   ├── __init__.py
│       │   ├── base_integrator.py  # Integrator抽象基类
│       │   ├── rk4.py             # RK4积分器
│       │   ├── dopri8.py          # DOPRI8自适应积分器
│       │   ├── gauss_jackson.py   # Gauss-Jackson多步积分器
│       │   └── stm_propagator.py  # 状态转移矩阵传播
│       │
│       ├── estimator/              # 参数估计层
│       │   ├── __init__.py
│       │   ├── base_estimator.py   # Estimator抽象基类
│       │   ├── lsq_estimator.py    # LSQ批处理估计器
│       │   ├── ekf_estimator.py    # EKF扩展卡尔曼滤波
│       │   ├── design_matrix.py    # 设计矩阵H计算
│       │   └── weight_matrix.py    # 权重矩阵W计算
│       │
│       ├── ambiguity/              # 模糊度固定层
│       │   ├── __init__.py
│       │   ├── base_ar.py          # AmbiguityResolver抽象基类
│       │   ├── lambda_search.py    # LAMBDA搜索算法
│       │   ├── decorrelation.py    # Z变换降相关
│       │   └── ratio_test.py       # Ratio检验
│       │
│       ├── quality_assessment/     # 质量评估层
│       │   ├── __init__.py
│       │   ├── overlap_test.py     # 重叠弧段检验
│       │   ├── residual_analysis.py # 残差分析
│       │   ├── sigma0_calc.py      # σ0计算
│       │   └── report_generator.py # 精度报告生成
│       │
│       ├── orekit_adapter/         # Orekit适配层(可选)
│       │   ├── __init__.py
│       │   ├── orekit_loader.py    # Orekit JVM启动
│       │   ├── orekit_dynamics.py  # Orekit动力学适配
│       │   ├── orekit_integrator.py # Orekit积分器适配
│       │   └── orekit_coords.py    # Orekit坐标转换适配
│       │
│       ├── output/                 # 输出层
│       │   ├── __init__.py
│       │   ├── sp3_writer.py       # SP3格式轨道输出
│       │   ├── csv_writer.py       # CSV残差输出
│       │   ├── covariance_writer.py # 协方差矩阵输出
│       │   └── log_writer.py       # 日志输出
│       │
│       └── utils/                  # 工具函数
│           ├── __init__.py
│           ├── time_utils.py        # 时间系统转换(GPS/UTC/UT1)
│           ├── coord_transforms.py  # 坐标转换(ECI/ECEF/RTN)
│           ├── quaternion_math.py   # 四元数运算
│           ├── constants.py         # 物理常数(MU_EARTH, C, etc.)
│           └── config_parser.py     # YAML配置解析
│
├── tests/                      # 单元测试
│   ├── test_data_loader/
│   ├── test_obs_preprocessing/
│   ├── test_dynamics/
│   ├── test_integrator/
│   ├── test_estimator/
│   └── test_ambiguity/
│
├── docs/                       # 文档
│   ├── algorithm_design.md
│   ├── api_reference.md
│   └── user_manual.md
│
├── requirements.txt
├── setup.py
└── README.md
```

---

## 9. 主程序流程设计

### 9.1 主流程伪代码

```python
# main.py - 主程序入口

def main(config_path: str):
    # 1. 加载配置
    config = load_config(config_path)
    
    # 2. 根据处理模式选择配置子集
    mode = config["global"]["processing_mode"]
    if mode == "realtime":
        cfg = config["realtime_overrides"]  # 简化配置
    elif mode == "nrt":
        cfg = config["nrt_overrides"]
    else:  # delayed
        cfg = config["delayed_overrides"]
    
    # 3. 加载数据
    data_loader = DataLoader(config["data_paths"])
    obs_data = data_loader.load_observations()
    attitude_data = data_loader.load_attitude()
    if mode != "realtime":
        precise_orbits = data_loader.load_sp3()
        precise_clocks = data_loader.load_clk()
    nav_data = data_loader.load_navigation()  # 广播星历（所有模式都需要）
    
    # 4. 观测预处理
    prep = ObservPrep(config["obs_preprocessing"])
    prep.set_attitude(attitude_data)
    preprocessed_obs = prep.process(obs_data)
    
    # 5. 构建动力学模型列表
    dynamics_list = build_dynamics_models(config["dynamics"], mode)
    
    # 6. 选择积分器
    integrator = build_integrator(config["integrator"], mode)
    
    # 7. 并行处理弧段
    arc_segments = split_arcs(preprocessed_obs, config["performance"]["parallel_arc_segments"])
    
    with ThreadPoolExecutor(max_workers=config["performance"]["num_threads"]) as executor:
        futures = []
        for arc in arc_segments:
            future = executor.submit(
                process_single_arc,
                arc, dynamics_list, integrator, config, mode
            )
            futures.append(future)
        arc_results = [f.result() for f in futures]
    
    # 8. 拼接弧段 + 重叠弧段检验
    final_orbit = merge_arcs(arc_results)
    qa = QualityAssessor(config["output"])
    qa_report = qa.evaluate(final_orbit, arc_results)
    
    # 9. 输出
    writer = OutputWriter(config["output"])
    writer.write_sp3(final_orbit)
    writer.write_report(qa_report)
    
    return final_orbit, qa_report

def process_single_arc(arc_data, dynamics_list, integrator, config, mode):
    """处理单个弧段"""
    # 初值：用广播星历或前一个弧段结果
    initial_state = get_initial_state(arc_data, config)
    
    # 构建估计器
    if mode == "realtime":
        estimator = EKFEstimator(config["estimator"])
    else:
        estimator = LSQEstimator(config["estimator"])
    
    # 执行估计
    result = estimator.estimate(
        obs_data=arc_data,
        dynamics=dynamics_list,
        integrator=integrator,
        initial_state=initial_state,
        initial_params=get_initial_params(config),
        apriori_cov=get_apriori_cov(config)
    )
    
    # 若启用AR，进行模糊度固定
    if config["ambiguity_resolution"]["enabled"] and result.converged:
        ar = LAMBDAmbiguityResolver(config["ambiguity_resolution"])
        ar_result = ar.resolve(result.float_ambiguities, result.covariance, arc_data)
        if ar_result.fixed:
            # 用固定模糊度重新估计
            result = estimator.reestimate_with_fixed_ambiguities(
                result, ar_result.fixed_ambiguity_vector
            )
    
    return result
```

### 9.2 近实时模式处理时序（2小时数据，≤2分钟）

```
T=0s     数据到达
T=0-5s   数据加载 + 格式解析
T=5-15s  观测预处理（周跳检测+质量控制）
T=15-20s 弧段分割（2h → 4×30min弧段）
T=20-120s 4个弧段并行处理（每个~25s）
          ├─ 初值设置 (广播星历)
          ├─ LSQ迭代1: 积分+设计矩阵 (~10s)
          ├─ LSQ迭代2: (~8s)
          ├─ LSQ迭代3: (~7s)
          └─ 收敛，输出弧段结果
T=120-130s 弧段拼接 + 重叠检验
T=130-135s 输出SP3 + 报告
T=135s     完成 ✓
```

### 9.3 延时模式处理时序（24小时数据，≤5分钟）

```
T=0s      数据到达
T=0-10s   数据加载
T=10-30s  观测预处理（不降采样，全量处理）
T=30-40s  弧段分割（24h → 12×2h弧段）
T=40-280s 12弧段并行（3线程，4轮，每轮~60s）
          每弧段: LSQ 3-4次迭代 + LAMBDA AR
T=280-290s 弧段拼接
T=290-300s 重叠弧段检验 + 报告
T=300s     完成 ✓
```

---

## 10. GPS与北斗联合定轨设计

### 10.1 数据层

```python
# 观测数据结构支持多GNSS
@dataclass
class ObsEpoch:
    time: GPSTime
    satellites: Dict[str, SatObs]  # key: "G01", "C01", "E01"等
    # SatObs包含该卫星的所有观测值

# 联合定轨时，GPS和北斗观测值同时参与估计
# 区别在于：
#   - 相位波长不同 (GPS L1=0.19029m, BDS B1=0.19204m)
#   - 模糊度参数分开（每卫星系统独立）
#   - DCB校正不同（GPS用CODE DCB, BDS用CAS BSX）
#   - 天线PCO/PCV模型不同（ANTEX中G系和C系分开）
```

### 10.2 模糊度参数设计

```
GPS L1模糊度: N_G01_L1, N_G02_L1, ...
GPS L2模糊度: N_G01_L2, N_G02_L2, ...
BDS B1模糊度: N_C01_B1, N_C01_B1, ...
BDS B3模糊度: N_C01_B3, N_C01_B3, ...
```

### 10.3 系统间偏差（ISB）参数

联合定轨需估计GPS与BDS间的系统间偏差：
```
ISB_GPS_BDS: 1个参数（所有BDS观测共享）
```

---

## 11. 关键算法详细说明

### 11.1 状态向量定义

**实时模式（EKF）**：12-20个参数
```
X = [x, y, z, vx, vy, vz,      # 轨道状态 (6)
     c*bias_rec,                  # 接收机钟差 (1)
     Cd, Cr]                      # 阻力/光压系数 (2)
     + 经验加速度(可选)
```

**近实时模式（LSQ）**：~50-200个参数（2h弧段）
```
X = [x0, y0, z0, vx0, vy0, vz0,  # 初始状态 (6)
     c*bias_rec_0, c*bias_rec_1,..., # 接收机钟差(每历元或随机游走)
     N_G01_L1, N_G01_L2, ...,      # 相位模糊度 (每卫星×每频点)
     Cd, Cr,                        # 阻力/光压系数
     aR_0, aT_0, aN_0, ...,       # 经验加速度 (2h/15min×3 = 24个)]
```

**延时模式（LSQ+AR）**：~500-2000个参数（24h弧段）
```
X = [同上，但弧段更长，经验加速度更多(24h/6min×3=720个)]
```

### 11.2 设计矩阵H的计算

观测方程：$y = h(X) + v$

设计矩阵H = ∂h/∂X，分块结构：

```
H = [H_orbit | H_clock | H_ambiguity | H_dynamics_params]
```

- `H_orbit`: 几何距离对轨道状态的偏导数 = ∂ρ/∂X = (r_sat - r_rec)/|r_sat - r_rec| · ∂r_sat/∂X
  - ∂r_sat/∂X 通过STM（状态转移矩阵）计算
- `H_clock`: 对钟差参数的偏导数 = 1（伪距）或0（相位，模糊度吸收）
- `H_ambiguity`: 对模糊度参数的偏导数 = λ（相位）或0（伪距）
- `H_dynamics_params`: 对Cd/Cr/经验加速度的偏导数，通过变分方程计算

### 11.3 变分方程（Variational Equations）

状态转移矩阵STM的定义：
$$\Phi(t, t_0) = \frac{\partial X(t)}{\partial X(t_0)}$$

变分方程：
$$\dot{\Phi} = A(t) \Phi, \quad \Phi(t_0) = I_6$$

其中 $A(t) = \frac{\partial f(t, X(t))}{\partial X}$（加速度对状态的偏导数矩阵）

**实现**：在数值积分时，同时积分状态向量X和STM（42维 = 6状态 + 36 STM）。

---

## 12. 实施路线图

### Phase 1: 核心框架 + 实时模式（2周）

- [ ] 项目骨架搭建（目录结构、配置系统）
- [ ] 数据加载器：RINEX观测+导航文件读取
- [ ] 观测预处理：周跳检测+质量控制
- [ ] 力学模型：二体+J2
- [ ] 积分器：RK4
- [ ] 估计器：EKF
- [ ] 输出：SP3格式
- [ ] 单元测试

### Phase 2: 近实时模式（3周）

- [ ] 数据加载器：SP3/CLK/ANTEX读取
- [ ] 观测预处理：PCO/PCV校正（含姿态）
- [ ] 力学模型：完整Reduced-Dynamic（重力场50×50+N体+潮汐+阻力+光压+经验加速度）
- [ ] 积分器：Gauss-Jackson
- [ ] 估计器：LSQ批处理
- [ ] 并行弧段处理
- [ ] 性能优化（满足2分钟约束）
- [ ] 质量评估：重叠弧段检验

### Phase 3: 延时模式 + PPP-AR（3周）

- [ ] 力学模型：完整150×150重力场+时变
- [ ] 模糊度固定：LAMBDA算法
- [ ] 估计器：LSQ+AR重新估计
- [ ] 性能优化（满足5分钟约束）
- [ ] GPS+北斗联合定轨
- [ ] DCB校正

### Phase 4: Orekit适配层（2周）

- [ ] Orekit JVM集成（py4j）
- [ ] OrekitDynamicsAdapter
- [ ] OrekitIntegratorAdapter
- [ ] 配置切换测试

### Phase 5: 验证与文档（2周）

- [ ] 用GRACE-FO L1B数据验证
- [ ] 精度评估（与JPL官方轨道比较）
- [ ] 用户文档
- [ ] API文档

---

## 13. 验证方案（多源星载数据 + 多方法验证）

> 本章系统梳理可用于遥感SAR卫星POD算法验证的**所有公开星载GNSS数据集**，以及**多种独立验证方法论**。
> 目标是建立一个分级验证体系：从简单比对到独立几何验证，层层深入。

---

### 13.1 可用于POD验证的公开星载GNSS数据集

以下数据集均可免费获取，按**与SAR卫星的相似度**排序（最相似的排最前）。

---

#### 12.1.1 Sentinel-1A / Sentinel-1B（最推荐 — C波段SAR卫星）

| 项目 | 详情 |
|------|------|
| **任务性质** | ESA哥白尼计划，C波段合成孔径雷达（SAR）卫星 |
| **轨道高度** | ~693 km，太阳同步轨道（SSO），LTAN ~18:00 |
| **星载GNSS** | GPS L1/L2（NovAtel OEMV-1G接收机），1Hz采样 |
| **数据格式** | RINEX 2.11 / 3.0x OBS + NAV |
| **精密轨道参考** | ESA官方科学轨道（~1-2cm精度，SLR验证） |
| **数据下载** | **ESA Copernicus Open Access Hub** — https://scihub.copernicus.eu/ |
| **GNSS原始数据** | **ESA GNSS RINEX Archive** — https://navipedia.esa.int/GNSS_Data/ |
| **备选入口** | **EUMETSAT** — https://www.eumetsat.int/sentinel-1-gnss-data |
| **精密轨道下载** | **ESA Sentinel-1 POD Products** — https://sentinel.esa.int/web/sentinel/missions/sentinel-1/products |
| **文件命名示例** | `S1A_OPER_AUX_POEORB_OPOD_20240101T000000_V20240102T000000_001.rmc.nc` |
| **更新频率** | 最终精密轨道约T+18天发布 |
| **推荐用途** | **主要验证数据集** — 与你的目标SAR卫星最相似（都是C波段SAR，高度~700km） |
| **注意事项** | Sentinel-1A/B的GNSS数据需要向ESA申请；精密轨道可直接从ESA网站下载NC/XML格式 |

**精密轨道精度（ESA产品）**：
```
Sentinel-1A Final Precise Orbit (POEORB):
  3D RMS vs SLR: ~1.5-2.0 cm
  径向精度: ~0.5-1.0 cm
  切向精度: ~1.0-1.5 cm
  法向精度: ~1.0-1.5 cm
```

---

#### 12.1.2 TerraSAR-X / TanDEM-X（高精度验证 — X波段SAR）

| 项目 | 详情 |
|------|------|
| **任务性质** | DLR/ASL X波段SAR卫星（双星编队） |
| **轨道高度** | ~514 km，太阳同步轨道 |
| **星载GNSS** | GPS L1/L2（Astrium/Thales接收机），1Hz，含GLONASS |
| **数据格式** | RINEX 2.11 OBS/NAV；DAS产品（GFZ格式） |
| **精密轨道参考** | DLR/AIUB联合解（~1cm精度） |
| **数据下载** | **DLR Oberpfaffenhofen** — https://www.dlr.de/eoc/terrasar-x/ |
| **公开数据入口** | **GFZ ISDC** — https://isdc.gfz-potsdam.de/index.php?id=32&L=1 |
| **精密轨道下载** | **AIUB CODTS** — ftp://ftp.aiub.unibe.ch/CODE/YYYY/ （CODTS产品含TerraSAR-X轨道） |
| **备选** | **GFZ POD Service** — https://doi.org/10.5880/GFZ.1.1.TDX.2023 |
| **推荐用途** | **高精度验证** — TerraSAR-X的精密轨道精度是公开数据中最高的之一（~1cm 3D） |
| **注意事项** | 部分数据需要向DLR申请许可；GFZ ISDC提供部分公开数据 |

**精密轨道精度（AIUB/DLR产品）**：
```
TerraSAR-X Final Orbit (AIUB solution):
  3D RMS vs SLR: ~1.0-1.5 cm
  径向精度: ~0.3-0.5 cm (最高精度)
```

---

#### 12.1.3 GRACE-FO（已覆盖 — 补充细节）

| 项目 | 详情 |
|------|------|
| **轨道高度** | ~490 km（最低，对大气阻力最敏感） |
| **星载GNSS** | GPS L1/L2（BlackJack接收机），10s降采样L1B |
| **精密轨道参考** | JPL/GFZ/AIUB三线一致解（~1-2cm） |
| **数据下载** | NASA PODAAC（已在前文详述） |
| **推荐用途** | **基准验证** — 最完整、最长期的数据集，适合算法调优 |
| **独特优势** | KBR（K波段测距）提供卫星间精密距离，可用于验证POD在编队飞行场景 |

---

#### 12.1.4 Sentinel-3A / Sentinel-3B（多仪器卫星，GPS+SLR）

| 项目 | 详情 |
|------|------|
| **任务性质** | ESA哥白尼计划，海洋/陆地观测（多光谱+SAR） |
| **轨道高度** | ~814 km，SSO |
| **星载GNSS** | GPS L1/L2（1Hz，RINEX格式） |
| **其他定轨手段** | SLR（激光反射镜阵列）、DORIS |
| **精密轨道参考** | EUMETSAT/ESA联合解（~2cm） |
| **数据下载** | **EUMETSAT** — https://www.eumetsat.int/sentinel-3-gnss-data |
| **精密轨道** | **EUMETSAT NTC Products** — https://navigator.eumetsat.int/product/ |
| **推荐用途** | **SLR独立验证** — Sentinel-3有SLR数据，可做GNSS vs SLR交叉验证 |

---

#### 12.1.5 Jason-3 / Jason-4（高度计卫星，GPS+DORIS+SLR，最高精度验证）

| 项目 | 详情 |
|------|------|
| **任务性质** | CNES/NASA高度计卫星（海洋测高） |
| **轨道高度** | ~1336 km（高于大部分LEO，对大气阻力不敏感） |
| **星载GNSS** | GPS + DORIS + LRA（激光反射镜） |
| **精密轨道参考** | CNES/NASA/ESA联合解（**~0.5-1.0cm 3D**，全球最高精度） |
| **数据下载** | **AVISO+** — https://www.aviso.altimetry.fr/en/data/products/jason-3.html |
| **精密轨道** | **CNES GDR Products** — https://data.nodc.noaa.gov/coris/data/FAO_GDR/ |
| **推荐用途** | **黄金标准验证** — Jason-3的精密轨道是全球定轨精度的黄金标准（~0.5cm精度） |
| **注意事项** | Jason-3高度过高（1336km），动力学环境与SAR卫星（~700km）差异较大；但作为精度验证的"真值"非常有说服力 |

---

#### 12.1.6 CryoSat-2（极轨卫星，GPS+DORIS+SLR）

| 项目 | 详情 |
|------|------|
| **任务性质** | ESA极地冰盖监测卫星 |
| **轨道高度** | ~717 km，高偏心轨道（轨道倾角92°，覆盖两极） |
| **星载GNSS** | GPS L1/L2（1Hz） |
| **其他定轨手段** | DORIS、SLR |
| **精密轨道参考** | ESA/AIUB解（~1-2cm） |
| **数据下载** | **ESA CryoSat-2 Data** — https://earth.esa.int/eogateway/missions/cryosat |
| **精密轨道** | **ESA Precise Orbit** — https://sentinel.esa.int/web/sentinel/missions/cryosat-2/products |
| **推荐用途** | **极区轨道验证** — CryoSat-2的极轨特性对经验加速度模型（N方向）是很好的考验 |

---

#### 12.1.7 MetOp-A/B/C（EUMETSAT气象卫星，GPS+GRAS）

| 项目 | 详情 |
|------|------|
| **任务性质** | EUMETSAT极轨气象卫星 |
| **轨道高度** | ~817 km，SSO |
| **星载GNSS** | GPS L1/L2（GRAS掩星接收机），1Hz |
| **精密轨道参考** | EUMETSAT（~2-3cm） |
| **数据下载** | **EUMETSAT** — https://www.eumetsat.int/metop-gnss-data |
| **推荐用途** | **NRT模式验证** — MetOp数据更新频繁，适合验证近实时处理 |

---

#### 12.1.8 FengYun-3 / FY-3（中国气象卫星，GPS+BD）

| 项目 | 详情 |
|------|------|
| **任务性质** | 中国气象局（CMA）极轨气象卫星 |
| **轨道高度** | ~836 km，SSO |
| **星载GNSS** | GPS L1/L2 **+ 北斗 B1/B3**（首批搭载北斗的LEO卫星） |
| **精密轨道参考** | 中国IGS分析中心（iGMAS）产品 |
| **数据下载** | **国家卫星气象中心** — http://www.nsmc.org.cn/nsmc/cn/home/index.html |
| **备选** | **iGMAS** — http://www.igmas.org/product/ |
| **推荐用途** | **北斗联合定轨验证** — FY-3是公开数据中少有的同时有GPS+北斗星载观测的数据集 |

---

#### 12.1.9 数据集对比总表

| 卫星 | 高度(km) | GPS | 北斗 | 精密轨道精度 | SLR | DORIS | 获取难度 | 推荐指数 |
|------|----------|:----:|:----:|:------------:|:---:|:------:|:-------:|:-------:|
| **Sentinel-1A/B** | 693 | ✓ | ✗ | ~1.5cm | ✗ | ✗ | 中 | ⭐⭐⭐⭐⭐ |
| **TerraSAR-X** | 514 | ✓ | ✗ | ~1.0cm | ✓ | ✗ | 高 | ⭐⭐⭐⭐⭐ |
| **GRACE-FO** | 490 | ✓ | ✗ | ~1.5cm | ✓ | ✗ | 低 | ⭐⭐⭐⭐ |
| **Sentinel-3A/B** | 814 | ✓ | ✗ | ~2.0cm | ✓ | ✓ | 中 | ⭐⭐⭐⭐ |
| **Jason-3** | 1336 | ✓ | ✗ | ~0.5cm | ✓ | ✓ | 中 | ⭐⭐⭐（高度过高） |
| **CryoSat-2** | 717 | ✓ | ✗ | ~1.5cm | ✓ | ✓ | 中 | ⭐⭐⭐⭐ |
| **MetOp-A/B/C** | 817 | ✓ | ✗ | ~2.5cm | ✗ | ✓ | 中 | ⭐⭐⭐ |
| **FengYun-3** | 836 | ✓ | ✓ | ~3.0cm | ✗ | ✗ | 高 | ⭐⭐⭐（北斗验证） |

---

### 13.2 精度验证方法论（多方法交叉验证）

> 单一验证方法可能存在系统性偏差。建议采用**三级验证体系**：
> 1. **一级验证**：与官方精密轨道直接比对（最直接）
> 2. **二级验证**：独立几何观测验证（SLR/DORIS，不依赖GNSS）
> 3. **三级验证**：内部一致性检验（重叠弧段、重收敛检验）

---

#### 12.2.1 一级验证：与官方精密轨道直接比对

**原理**：将你的POD结果与官方分析中心发布的精密轨道（参考真值）逐历元比对，计算3D RMS。

**步骤**：
```
1. 下载官方精密轨道（SP3格式）
2. 将你的轨道结果也输出为SP3格式
3. 时间对齐（内插官方轨道到你的历元）
4. 坐标系对齐（确保都在ITRF框架下）
5. 计算差异向量：Δr(t) = r_est(t) - r_ref(t)
6. 统计：3D RMS = sqrt( mean(||Δr||²) )
```

**关键注意**：
- 官方轨道SP3文件的参考框架可能是IGS20/ITRF2014/ITRF2020，需要与你的结果框架一致
- 使用**非重叠区间**的历元进行统计（重叠部分可能人为降低RMS）
- 分别报告**R/T/N方向**的RMS（径向精度通常最高，法向最低）

**参考轨道质量分级**（用于选择比对对象）：
```
金牌: Jason-3 CNES (0.5cm), TerraSAR-X AIUB (1.0cm)
银牌: GRACE-FO JPL (1.5cm), Sentinel-1A ESA (1.5cm)
铜牌: Sentinel-3 EUMETSAT (2.0cm), MetOp EUMETSAT (2.5cm)
```

---

#### 12.2.2 二级验证：SLR（卫星激光测距）独立几何验证

**原理**：SLR地面站向卫星发射激光脉冲，测量往返时间，得到mm级精度的几何距离。SLR完全独立于GNSS，是验证GNSS POD结果的**黄金标准**。

**哪些卫星有SLR？**
- Jason-3 ✅（LRA激光反射镜阵列）
- Sentinel-3A/B ✅
- CryoSat-2 ✅
- TerraSAR-X ✅
- GRACE-FO ✅
- Sentinel-1A/B ❌（无LRA）

**验证步骤**：
```
1. 获取SLR观测数据
   下载: https://ilrs.gsfc.nasa.gov/data_and_products/data.html
   格式: CRD (Consolidated Ranging Data) 或传统格式

2. 用你的POD轨道计算理论SLR距离
   ρ_SLR(t) = || r_stat_ECI(t) - r_sat_ECI(t) || + 相对论校正

3. 与SLR观测距离比对
   Residual = ρ_obs - ρ_calc
   统计: O-C RMS (Observation minus Calculation)

4. 合格标准:
   SLR O-C RMS < 2-3 cm  → 你的POD轨道精度可信
   SLR O-C RMS > 5 cm   → 可能存在系统误差
```

**SLR数据下载**：
| 数据产品 | 下载地址 | 格式 |
|---------|---------|------|
| ILRS快速SLR | https://ilrs.gsfc.nasa.gov/data_and_products/rapid_products.html | CRD |
| ILRS最终SLR | https://ilrs.gsfc.nasa.gov/data_and_products/data.html | CRD/传统 |
| CSR SLRF2014 | https://csr.utexas.edu/slr/toric/products.html | 测站坐标 |

---

#### 12.2.3 二级验证：DORIS独立验证

**原理**：DORIS（Doppler Orbitography by Radiopositioning Integrated on Satellite）是法国CNES开发的星载多普勒测速系统，独立于GNSS。

**哪些卫星有DORIS？**
- Jason-3 ✅
- Sentinel-3A/B ✅
- CryoSat-2 ✅
- HY-2（中国）✅
- SPOT系列（部分）✅

**验证方法**：与SLR类似，用你的POD轨道计算理论DORIS距离，与观测值比对。

**DORIS数据下载**：
- **AVISO+ (CNES)** — https://www.aviso.altimetry.fr/en/data/products/jason-3/doris.html
- **IDS (International DORIS Service)** — https://ids-doris.org/data-products.html

---

#### 12.2.4 三级验证：内部一致性检验

**方法A：重叠弧段检验（Overlap Arc Test）**

```
将24h数据分成：
  弧段1: [T, T+24h]
  弧段2: [T+18h, T+42h]   ← 与弧段1有6h重叠

分别用弧段1和弧段2独立定轨，
比较重叠区域（T+18h ~ T+24h）的轨道差异。

合格标准: 重叠区域3D RMS < 2 cm
```

**方法B：重收敛检验（Re-convergence Test）**

```
用不同初值（广播星历 vs 运动学轨道 vs 随机初值）分别定轨，
比较收敛后的轨道结果。

合格标准: 不同初值得到的轨道3D差异 < 1 cm
```

**方法C：模糊度固定率检验**

```
统计LAMBDA固定解的成功率：
  AR成功率 = N_fixed / N_total
  
合格标准:
  延时模式: AR成功率 > 85%
  近实时模式: AR成功率 > 70%（若不启用AR则不适用）
```

---

#### 12.2.5 四级验证：预测轨道精度检验（针对实时模式）

**原理**：实时模式输出的轨道需要外推（预测）一段时间。检验预测轨道与事后精密轨道的差异。

**步骤**：
```
1. 实时模式在T时刻输出轨道（含预测部分）
2. 用事后精密轨道作为参考
3. 比较预测部分（T ~ T+Δt）的轨道差异
4. 绘制预测误差随时间增长的曲线

合格标准（经验值）:
  Δt=15min: 预测误差 < 5 cm
  Δt=1h:   预测误差 < 20 cm
  Δt=3h:   预测误差 < 50 cm
```

---

### 13.3 各数据集精密轨道参考获取详细指南

#### Sentinel-1A/B 精密轨道

| 产品类型 | 延迟 | 精度 | 下载地址 | 文件格式 |
|---------|------|------|---------|---------|
| **POEORB**（精确轨道） | T+18-27d | ~1.5cm | https://sentinel.esa.int/web/sentinel/missions/sentinel-1/products | NC/XML |
| **RESORB**（快速轨道） | T+3d | ~5cm | 同上 | NC/XML |

**文件解析**（Python）：
```python
import xarray as xr
ds = xr.open_dataset("S1A_OPER_AUX_POEORB_OPOD_*.nc")
time = ds["time"][:]      # 时间戳
x = ds["X_ECF"][:]       # ECEF X坐标 (m)
y = ds["Y_ECF"][:]
z = ds["Z_ECF"][:]
```

---

#### TerraSAR-X 精密轨道

| 产品类型 | 来源 | 精度 | 下载 |
|---------|------|------|------|
| **AIUB CODTS** | AIUB FTP | ~1.0cm | ftp://ftp.aiub.unibe.ch/CODE/YYYY/CODTS... |
| **DLR GFZ** | GFZ ISDC | ~1.2cm | https://isdc.gfz-potsdam.de/ |

---

#### GRACE-FO 精密轨道（补充）

| 分析中心 | 产品 | 精度 | 下载 |
|---------|------|------|------|
| JPL | RL04 | ~1.5cm | https://podaac.jpl.nasa.gov/dataset/GRACEFO_L2_CSpriteVariations_RL04 |
| GFZ | RL04 | ~1.5cm | https://isdc.gfz-potsdam.de/grace-fo/ |
| AIUB | RL04 | ~1.2cm | ftp://ftp.aiub.unibe.ch/CODE/YYYY/ |

---

#### Jason-3 精密轨道（黄金标准）

| 产品 | 分析中心 | 精度 | 下载 |
|------|---------|------|------|
| GDR-F (Final) | CNES/NASA | ~0.5cm | https://www.aviso.altimetry.fr/en/data/products/jason-3.html |
| GDR-INT (Intermediate) | CNES | ~1.0cm | 同上 |

---

### 13.4 验证指标与合格标准汇总

| 验证层级 | 指标 | 近实时目标 | 延时目标 | 实时目标 |
|---------|------|----------|---------|---------|
| **一级**（轨道比对） | 3D RMS | < 5cm | < 2cm | < 20cm |
|  | 径向RMS | < 2cm | < 0.5cm | < 10cm |
|  | 切向RMS | < 3cm | < 1.0cm | < 15cm |
|  | 法向RMS | < 4cm | < 1.5cm | < 18cm |
| **二级**（SLR O-C） | SLR O-C RMS | < 4cm | < 2cm | N/A |
| **三级**（重叠弧段） | 重叠6h 3D RMS | < 3cm | < 1.5cm | N/A |
| **四级**（预测误差） | 15min预测 | N/A | N/A | < 5cm |
|  | 1h预测 | N/A | N/A | < 20cm |
| **质量**（模糊度） | AR固定率 | N/A | > 85% | N/A |
|  | Sigma0相位 | < 0.015周 | < 0.008周 | < 0.02周 |

---

### 13.5 推荐验证流程图

```
                    ┌─────────────────────────────────┐
                    │   选择验证数据集                  │
                    │ (Sentinel-1A / TerraSAR-X /    │
                    │  GRACE-FO / Jason-3)           │
                    └──────────┬──────────────────────┘
                               │
                    ┌──────────▼──────────────────────┐
                    │   下载数据                       │
                    │  - 星载GNSS观测 (RINEX OBS)     │
                    │  - 精密轨道参考 (SP3)            │
                    │  - (可选) SLR数据 (CRD)          │
                    └──────────┬──────────────────────┘
                               │
              ┌────────────────▼────────────────┐
              │   运行你的POD算法                 │
              │   (实时/近实时/延时三模式)        │
              └────────┬─────────────────────┬──┘
                       │                     │
            ┌──────────▼──────────┐   ┌─────▼──────────────┐
            │ 一级验证：轨道比对    │   │ 二级验证：SLR O-C   │
            │ 3D RMS, R/T/N RMS   │   │ (若有SLR数据)       │
            └──────────┬──────────┘   └─────┬──────────────┘
                       │                     │
                       └──────────┬──────────┘
                                  │
                    ┌─────────────▼──────────────┐
                    │   三级验证：内部一致性        │
                    │   重叠弧段 + AR固定率        │
                    └─────────────┬──────────────┘
                                  │
                    ┌─────────────▼──────────────┐
                    │   输出验证报告                │
                    │   (精度表格 + 残差图)        │
                    └─────────────────────────────┘
```

---

### 13.6 自动化验证框架建议

为支持持续集成（CI），建议构建自动化验证框架：

```python
# validate_pod.py - 自动化验证脚本框架

class PODValidator:
    def __init__(self, config):
        self.config = config
        self.results = {}
    
    def run_validation(self, dataset_name, mode):
        """对指定数据集运行完整验证"""
        # 1. 下载/加载数据
        data = self.load_dataset(dataset_name, mode)
        
        # 2. 运行POD
        pod_result = self.run_pod(data, mode)
        
        # 3. 一级验证：轨道比对
        ref_orbits = self.load_reference_orbits(dataset_name)
        rms_3d, rms_r, rms_t, rms_n = self.compare_orbits(
            pod_result.orbit_sp3, ref_orbits.sp3
        )
        
        # 4. 二级验证：SLR（如果可用）
        if dataset_name in self.SLR_CAPABLE:
            slr_residuals = self.validate_slr(
                pod_result.orbit_sp3, data.slr_crd
            )
            slr_rms = np.std(slr_residuals)
        else:
            slr_rms = None
        
        # 5. 三级验证：重叠弧段
        overlap_rms = self.overlap_test(pod_result.arc_results)
        
        # 6. 生成报告
        report = self.generate_report({
            "dataset": dataset_name,
            "mode": mode,
            "rms_3d": rms_3d,
            "rms_r": rms_r, "rms_t": rms_t, "rms_n": rms_n,
            "slr_rms": slr_rms,
            "overlap_rms": overlap_rms,
            "ar_fix_rate": pod_result.ar_fix_rate,
            "sigma0_phase": pod_result.sigma0_phase,
        })
        
        return report
    
    def load_dataset(self, name, mode):
        """加载指定数据集"""
        if name == "sentinel1a":
            return self.load_sentinel1a(mode)
        elif name == "terrasarx":
            return self.load_terrasarx(mode)
        # ...
```

---

---

## 14. 附录：重要常数与公式

### 14.1 物理常数

```python
MU_EARTH = 3.986004418e14    # m³/s² (WGS84)
RE_EARTH = 6378137.0          # m (WGS84椭球长半轴)
J2 = 1.0826359e-3            # 无量纲
C = 299792458.0               # m/s (光速)
F = -2.0 * C / 0.19029e-5   # 电离层因子 (频间)
```

### 14.2 关键公式

**几何距离**：
$$\rho = |\vec{r}_{sat}^{ECI}(t_{sat}) - \vec{r}_{rec}^{ECI}(t_{rec})|$$

**相对论钟差校正**（接收机端）：
$$\Delta t_{rel} = -2\frac{\vec{r}_{sat} \cdot \vec{v}_{sat}}{c^2}$$

**相位观测方程**：
$$\lambda \cdot L = \rho + c(dt_{rec} - dt_{sat}) + T + I + \lambda \cdot N + \epsilon$$

**伪距观测方程**：
$$P = \rho + c(dt_{rec} - dt_{sat}) + T + I + \epsilon$$

---

*文档结束*
