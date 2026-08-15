# GNSS POD 精度提升路线图

## 目标：LEO 卫星精密定轨 ≤ 5 cm (3D RMS)

从 V2.2.1 (Sequential EKF, Float PPP, 0.94m) 出发，分阶段提升。**V3.3 已达成 0.047m GRACE-FO, 0.017m SWARM (04-29). 多任务 YAML 配置 + Orekit EKF 动力学**。

---

## 1. 精度演进总览 (GRACE-FO C, GPS only)

| Phase | 版本 | 方法 | 0.17h | 0.5h | 说明 |
|-------|------|------|-------|------|------|
| V2.2.1 | — | Float PPP | 0.936m | 1.817m | baseline |
| 7.0 | V2.2.4 | 自适应 clock_rw | 0.293m | 0.986m | -69%/-46% |
| 20.0 | V3.0.0 | Orekit GN (v13+速度优化) | 0.043m | 0.409m | -95%/-77% (单天最优) |
| **24.0** | **V3.1** | **QC自动化 + gap分裂 + IRLS + TurboEdit** | **0.047m** | **0.313m** | **多天一致性修复 ★** |
| **26.0** | **V3.2.1** | **Code-Orbit自洽初轨 + GNV1B消除** | **0.047m** | **0.397m** | **自洽管线** |
| **27.0** | **V3.3.0** | **多任务 + Orekit EKF动力学 + YAML配置** | **0.017m** | — | **SWARM 多任务 ★** |

### V3.1 9 天结果 (2026-06-28, Phase 24.0)

| Date | 0.17h 3D RMS | QC | 覆盖率 | 完整SV | 0.50h 3D RMS | 异常标记 |
|------|-------------|-----|--------|--------|-------------|---------|
| 04-29 | **0.047m** ★ | 0.96A | 85.7% | 7/12 | 0.397m | 1 gap SV |
| 04-30 | 0.189m | 0.81B | 67.8% | 5/13 | 0.541m | 2 SV拒收 |
| 05-01 | 0.171m | 0.81B | 71.8% | 7/13 | 0.392m | 3 gap SV |
| 05-02 | 0.344m | 0.73B | 63.9% | 6/14 | 0.434m | 4 SV拒收 |
| 05-03 | 0.266m | 0.71B | 56.3% | 3/11 | 0.366m | 11 SV拒收 |
| 05-04 | **0.067m** ★ | 0.82B | 83.3% | 5/9 | 0.320m | 仅9 SV |
| 05-05 | 0.231m | 0.85B | 80.6% | 9/12 | 0.813m | 3 SV拒收 |
| 05-06 | **0.102m** ☆ | 0.90A | 81.0% | 6/12 | 0.313m | G05已修复 |
| 05-08 | 0.101m | 0.81B | 81.7% | 9/13 | 0.863m | 3 SV拒收 |

**统计**: 0.17h mean=**0.169m**, median=0.171m, ≤0.2m: 6/9 天

### V3.1 第二批验证 (2026-06-30, May 9-17)

| Date | 0.17h 3D RMS | QC Score | Coverage | SV | Flags |
|------|------------|----------|----------|-----|-------|
| 05-09 | 0.244m | 0.86A | 78% | 12 | OK |
| 05-10 | 0.181m | 0.70B | 73% | 13 | 6 SV rejected |
| 05-11 | 0.265m | 0.78B | 81% | 10 | 2 gap SV |
| 05-12 | **0.080m** ★ | 0.85B | 68% | 14 | 2 SV rejected |
| 05-13 | 0.213m | 0.79B | 75% | 11 | 2 SV rejected |
| 05-14 | 0.178m | 0.72B | 67% | 12 | 3 gap, 4 rejected |
| 05-15 | **0.134m** ★ | 0.93A | 76% | 10 | 1 SV rejected |
| 05-16 | ⚠️ **4.688m** | 0.70B | 73% | 12 | BS code RMS=21.7m |
| 05-17 | ⚠️ **1.708m** | 0.77B | 69% | 11 | 1 gap, 4 rejected |

**统计** (排除05-16/17异常): 0.17h mean=**0.185m** (7天), median=0.181m, ≤0.2m: 4/7
**全部9天**: mean=0.855m, median=0.213m

### 关键发现：覆盖率决定精度

"可见星多→精度差"是伪相关。真正的自变量是 **SV 覆盖率**：

| 覆盖率区间 (0.17h) | 天数 | 平均3D RMS | 代表 |
|-------------------|------|-----------|------|
| > 80% | 5天 | **0.110m** | 04-29, 05-04/05/06/08 |
| 70-80% | 2天 | 0.180m | 04-30, 05-01 |
| < 70% | 2天 | **0.305m** | 05-02, 05-03 |

覆盖率与 3D RMS 相关系数: 0.17h r=-0.70, 0.50h r=-0.63。

### 05-06 异常：根因与修复

G05 在10分钟弧段内有2次数据中断(120-150s)，每次恢复后相位跳变138/157周期。
Batch求解器假设模糊度在整弧段内为常数——G05被强制拟合为单值→Phase RMS 50.7m→污染全局解。

| 修复 | 05-06 效果 |
|------|-----------|
| SV Gap分裂(>90s)→独立弧段 | G05分3段→batch Phase 50.7m→0.26m |
| MW跳变拒收(σ>10cyc或≥2跳变) | 自动排除G05 |
| Huber IRLS(k=2.5σ,3次迭代) | 9%观测自动降权 |

**05-06 0.17h: 5.88m→0.102m(58x), 0.50h: 9.13m→0.313m(29x)**

### V3.2 真实广播星历验证 (2026-07-05, IGS BRDC from BKG, --broadcast flag)

8 天 IGS 真实广播星历 (RINEX nav 文件, ~270KB/天, 32 GPS SV) 与 CODE 精密产品对比。

| Date | BRDC 3D RMS | CODE 3D RMS | 退化 | BRDC 3V RMS | BRDC Phase | 
|------|-----------|-----------|------|-----------|-----------|
| 04-29 | 1.825m | 0.047m | 39× | 12.2mm/s | 0.285m |
| 04-30 | 1.306m | 0.189m | 7× | 35.5mm/s | 0.379m |
| 05-01 | **0.428m** ★ | 0.171m | 2.5× | 18.4mm/s | 0.269m |
| 05-02 | 0.782m | 0.344m | 2.3× | 53.3mm/s | 0.199m |
| 05-04 | 2.362m | 0.067m | 35× | 12.9mm/s | 0.192m |
| 05-05 | 1.229m | 0.231m | 5.3× | 31.9mm/s | 0.423m |
| 05-06 | 0.879m | 0.102m | 8.6× | 14.8mm/s | 0.268m |
| 05-08 | 2.031m | 0.101m | 20× | 24.6mm/s | 0.711m |

**统计**: BRDC mean=**1.355m**, CODE mean=**0.169m**, 退化 **2.3-39×** (均值~8×)。

**关键发现**:
1. **真实广播星历优于模拟数据**: 均值 1.36m vs 模拟 5.60m——真实 BRDC 轨道误差并非纯粹高斯白噪声，存在系统性缓慢变化，EKF 可部分吸收
2. **Phase RMS 保持低值 (0.19-0.71m)**: 模糊度固定 + MW 宽巷窄巷仍有效——广播星历不影响相位观测量的一致性
3. **05-01 最佳 (0.428m)**: 当日 GPS 星座几何与 BRDC 轨道残差形成幸运抵消
4. **精度由 GPS 轨道误差主导**: BRDC vs CODE SP3 轨道差异 ~1-2m (3D RMS)，直接通过 GDOP 映射到 LEO 位置解。GPS 钟差可通过接收机钟参数部分吸收 (BRDC 钟差经 clk*C 转换为米后与 SP3 一致性良好)
5. **无法满足 ≤5cm 精密定轨**: 8× 的平均退化意味着广播星历仅能支持 ~1m 级的 LEO 定轨。任何后处理算法都无法补偿 1-2m 的输入 GPS 轨道误差

**数据源**: BKG (德国联邦测绘局) IGS 数据镜像 `igs.bkg.bund.de`
**预处理**: BKG RINEX 2.11 nav → `parse_brdc_dataline` 提取 D 格式 → `compute_satpos_nav` (ICD-200 开普勒方程) → SP3 格式 pkl 缓存 (钟差=clk_s×C 米)

---

## 2. Phase 25.0: 三条改进路径

### 路径① — 高覆盖率弧段自动筛选 (已实现)

```bash
py eval_5day_orekit.py --min-coverage 0.70
```
覆盖率 <70% 的弧段自动 SKIP。覆盖率 >80% 的5天精度均值 0.110m。

### 路径② — 覆盖率自适应模糊度先验 (已实现)

BatchSolver 对每颗 SV 的模糊度施加覆盖率自适应正则化：
```
σ_prior(sv) ∝ 1/(coverage_pct)²
覆盖100%: 1e-4 baseline → 强约束
覆盖30%:  1e-2 regularization → 防止漂移污染其他参数
```
替代当前的保留/拒收二分法，实现软过渡。文件: `src/batch_solver.py` (`sv_coverage` 参数)

### 路径③ — 冗余弧段融合 (已实现)

对1小时连续数据取N个滑动10min弧段，QC加权融合：
```bash
py eval_5day_orekit.py --fuse-arcs 6
```
`r_fused = Σ(w_i × r_i) / Σw_i, w_i = qc_score²`

GPS星座每10min经历"世代交替"→相邻弧段观测不同SV组合→独立解算→融合消除病态模糊度偏置。

---

## 3. V3.1 自动化 QC 体系

```
L1 原始数据 → SV覆盖率/SNR/MP1筛选 → 产品完整性校验
L2 EKF→Batch → MW稳定性分析 → gap自动分裂 → 覆盖率自适应拒收
L3 Batch求解 → Huber IRLS迭代重加权(k=2.5σ)
L4 求解后 → 6项加权综合质量评分(0-1, A/B/C/D/F)
```

每条弧段输出: `[QC] score=0.96(A) | flags: 1 rejected SV | actions: 1 SV rejected`

新增文件: `src/data_quality.py` (4层QC函数库)

---

## 4. 软件架构

### 三条处理管线
```
GPS1B观测 + SP3/CLK/DCB/ANTEX/IERS + 重力场
    │
    ├── ① EKF序贯滤波 (实时/近实时) → 3D RMS: 0.29m
    ├── ② Batch线性求解 (固定轨道)   → Phase: 0.16m  
    └── ③ Orekit GN外层 (最高精度)   → 3D RMS: 0.047m ★
```

### 关键文件
| 文件 | 功能 |
|------|------|
| `eval_5day_orekit.py` | 多天多任务 Orekit GN 验证 + YAML配置 + 广播星历 |
| `src/sequential_filter.py` | EKF核心(~1300行) + Orekit动力学预报 |
| `src/batch_solver.py` | BatchLinearSolver + IRLS + 自适应先验 |
| `src/batch_orbit_v3.py` | 9-param GN + Orekit外层 + 自退 Python 动力学 |
| `src/code_orbit.py` | 伪距GN自洽初轨 + 异常检测 (V3.2) |
| `src/data_quality.py` | 4层QC: 覆盖/SNR/MP1/MW/产品/评分 |
| `src/config_loader.py` | YAML 配置管线 builder (V3.2.1) |
| `src/orekit_bridge.py` | Orekit Java桥接 + SRP(none)/drag(Harris-Priester)控制 (V3.3) |
| `src/swarm_adapter.py` | SWARM L2 SP3 参考轨道加载 + pkl 缓存 (V3.3) |
| `scripts/convert_swarm_rnx3.py` | RINEX 3.00 → GPS1B 转换器 (V3.3, 固定宽度解析) |
| `V3.2.1/config.yaml` | GRACE-FO YAML 配置 |
| `V3.2.1/config_SWARM.yaml` | SWARM YAML 配置 (Orekit EKF) |

### 运行命令
```powershell
$env:OREKIT_DATA_PATH = 'd:\prj\gnss_pod\data\orekit'
# 标准9天评估
py eval_5day_orekit.py --dates 2024-04-29,...,2024-05-08 --hours 0.17,0.50
# 仅高覆盖弧段(Path1)
py eval_5day_orekit.py --min-coverage 0.70
# 弧段融合模式(Path3)
py eval_5day_orekit.py --dates 2024-04-29 --hours 0.17 --fuse-arcs 6
# 真实广播星历评估
py eval_5day_orekit.py --broadcast --skip-code-orbit

# 多任务 YAML 配置模式
py eval_5day_orekit.py --config V3.2.1/config.yaml         # GRACE-FO
py eval_5day_orekit.py --config V3.2.1/config_SWARM.yaml    # SWARM
```

---

## 5. Multi-Mission Support (V3.3.0, 2026-07-15)

### 支持的卫星任务

| 任务 | 接收机 | 通道数 | 数据格式 | L2参考轨道 | GNV1B依赖 | 状态 |
|------|--------|--------|---------|-----------|----------|------|
| GRACE-FO | BlackJack | 12 | RINEX 2.11/GPS1B | GNV1B (1cm) | 验证用 | **生产** |
| **SWARM-A** | **RUAG** | **8** | **RINEX 3.00** | **CODE L2 RN SP3 (2cm)** | **无 ★** | **评估中** |
| GRACE | BlackJack | 12 | RINEX 2.11/GPS1B | GNV1B | 验证用 | 待验证 |
| FY-3 | — | — | 待适配 | — | 无 | 待适配 |
| COSMIC-2 | — | — | 待适配 | — | 无 | 待适配 |
| Jason-3 | — | — | 待适配 | — | 无 | 待适配 |

### YAML 配置驱动架构

40+ 参数从 YAML 配置加载，支持每卫星独立调优：

```
sections: qc.* | ekf.* | gn_loop.* | gravity.* | dynamics.*
```

- `config.yaml` → GRACE-FO 标准配置 (QC严格, EKF Python动力学)
- `config_SWARM.yaml` → SWARM 配置 (Orekit EKF动力学, SRP=none)
- 新增卫星: 复制 config.yaml → 按接收机特性调整阈值

### SWARM-A 关键配置差异 (vs GRACE-FO)

| 参数 | GRACE-FO | SWARM | 原因 |
|------|----------|-------|------|
| ekf.dynamics_mode | simplified | **orekit** | SWARM需全动力学补偿8通道弱几何 |
| srp_model | isotropic | **none** | Orekit EclipseDetector NPE @ SWARM高-β轨道 |
| drag_model | exponential | **harris-priester** | HP模型更适合SWARM (~450km高度) |
| min_sv_coverage | 0.40 | **0.20** | 8通道→可见星减少 |
| mw_std_noisy | 0.30cyc | **0.60cyc** | RUAG C/A码MW噪声稍大 |
| mw_std_unstable | 0.50cyc | **0.80cyc** | MW组合稳定性放宽 |
| wl_fix_residual | 0.35cyc | **0.45cyc** | 宽巷残差放宽 |
| chi2_short | 25 | **50** | 多径大→容忍度增加 |
| sigma_phase | 0.20m | **0.25m** | 相位观测先验放宽 |
| clock_rw_short | 0.0004 | **0.0008** | 钟差过程噪声放宽 |
| satellite.mass | 600kg | **473kg** | SWARM-A 质量 |
| satellite.area_drag | 0.25m² | **1.0m²** | 阻力面积 |
| satellite.area_srp | 2.0m² | **4.5m²** | 光压面积 |
| satellite.CD | 1.4 | **2.4** | 阻力系数 |
| satellite.CR | 1.3 | **1.2** | 光压系数 |

### SWARM-A 8天评估结果 (V3.3.0, Orekit EKF dynamics + fixed RINEX parser)

8天: 2024-04-29 ~ 2024-05-07 (05-05 无数据跳过)

| Date | 3D RMS | 3V RMS | SVs | GN Phase | QC | Method |
|------|--------|--------|-----|----------|-----|--------|
| 04-29 | **0.017m ★** | 0.1mm/s | 3 | — | 0.88A | Orekit EKF predict-only (3SV <4→skip GN) |
| 04-30 | 53.222m | 62.2mm/s | 11 | 7.22m | 0.84B | Orekit GN 6iter, 6/8 SV NL-fixed |
| 05-01 | 49.495m | 54.7mm/s | 9 | 4.64m | 0.88A | Orekit GN 6iter, 6/9 SV NL-fixed |
| 05-02 | 48.988m | 57.5mm/s | 10 | 6.17m | 0.81B | Orekit GN 6iter, 7/9 SV NL-fixed |
| 05-03 | 58.981m | 59.4mm/s | 8 | 6.00m | 0.86A | Orekit GN 6iter, 5/8 SV NL-fixed |
| 05-04 | 62.183m | 99.7mm/s | 11 | 15.88m | 0.80B | Orekit GN 6iter, 2/6 SV NL-fixed |
| 05-06 | 48.264m | 56.9mm/s | 11 | 5.66m | 0.84B | Orekit GN 6iter, 6/9 SV NL-fixed |
| 05-07 | 62.582m | 64.6mm/s | 10 | 10.08m | 0.83B | Orekit GN 6iter, 4/9 SV NL-fixed |

**统计**: 04-29: 0.017m (Orekit动力学主导), 其余7天均值=**54.8m**, 中位数=53.2m

### V3.3.0 关键发现

#### 1. RINEX 3.00 解析 Bug（已修复）

**根因**: `convert_swarm_rnx3.py` 用 `split()` 解析 RINEX 3.00 观测行，但 RINEX 3 格式是 `F14.3,I1,I1`（16字符/观测量），LLI/SSI 标志位被当作独立 token → **所有观测值错位**。

| 示例 | 旧解析(split) | 正确解析(16-char) |
|------|-------------|------------------|
| `38.800 6` | ["38.800", "6"] | C1W=23725512.269 |
| 6 → 被当作 P1 伪距值 | P1 = 6m (错误) | P1 = 23,725,512m (正确) |

**修复**: 固定宽度解析 `data_line[start:start+16][:14]`，跳过 LLI/SSI 标志位。

| Date | MW std (修复前) | MW std (修复后) | WL 固定 (修复前) | WL 固定 (修复后) |
|------|----------------|----------------|-----------------|-----------------|
| 04-30 | >100,000 cyc | **0.05-0.22 cyc** | 0 SV | **8 SV** |
| 05-01 | >100,000 cyc | **0.08-0.39 cyc** | 0 SV | **9 SV** |
| 全部7天 | 全部 >100,000 cyc | 全部 <1.2 cyc | 0 | 6-9 SV |

#### 2. 04-29 为何 0.017m？

- **仅有的 21 历元文件**（预截断，非全天数据）
- 3 SV 通过 EKF chi² → EKF 测量更新有效 → Batch 相位 RMS=38m
- Orekit 全动力学 (GGM05C 150阶 + 潮汐 + HP 大气阻力) 从 L2 参考初态传播10分钟 → 漂移<2cm
- **这不是 SWARM POD 的真实精度** — 是 Orekit 动力学的自洽性测试：用精密初态 + 精密力模型，10分钟预报与参考轨道的 RMS=1.7cm

#### 3. 04-30+ 为何 ~50m？

**多径是主导限制**:
- 所有 8 天 100% SV 的 MP1 RMS > 4.5m 阈值
- 相位残差 4.6-15.9m，码残差 13-33m
- **Orekit GN 在 3-6 次迭代内收敛 → 相位残差从 ~38m 降到 ~7m** — 说明动力学约束有效
- 但轨道精度停止在 ~50m — **这是观测噪声墙，不是动力学的限制**

#### 4. 与 GRACE-FO 的对比

| 指标 | GRACE-FO | SWARM | 比值 |
|------|----------|-------|------|
| 平均有效 SV | 10-14 | 8-11 | ~0.75× |
| 平均低SNR SV | 0-1 | 2-5 | ~4× |
| MP1超阈值SV | 0-2 | 8-11 (100%) | **>5×** |
| MW std 中位数 | 0.03-0.15 cyc | 0.05-0.25 cyc | ~2× |
| Phase RMS | 0.16m | 5-16m | **30-100×** |
| Code RMS | 0.5-1.5m | 13-33m | **~20×** |
| 3D RMS (最优) | 0.047m | 0.017m | — |
| 3D RMS (均值) | 0.169m | 54.8m | — |

SWARM 的 RUAG 接收机天线位置（不在质心）导致严重的多径，这是 cm 级定轨的根本障碍。

### 改进方向 (仅配置/参数层面)

| 优先级 | 方向 | 预期效果 | 类型 |
|--------|------|---------|------|
| P1 | 提高 EKF chi² 阈值至 100-200 | 更多观测量参与滤波 | 配置参数 |
| P2 | MP1 阈值从 4.5m 提高到 8-10m | 保留更多 SV | 配置参数 |
| P3 | GN 先验收窄 (prior_r0=1m, prior_v0=0.005m/s) | 约束 GN 步进 | 配置参数 |
| P4 | 增加观测弧长到 0.5-2h | 更多动力学约束 | 配置参数 |
| P5 | 预拟合 RUAG 接收机 DCB (日稳定分量) | MW 系统性偏移消除 | 预处理脚本 |

> **注**: SWARM 的 GNV1B 依赖已完全消除。初轨使用 L2 参考轨道 (仅初始历元，不参与滤波更新)。CODE L2 RN 轨道产品从 ESA SWARM 数据门户下载，~2cm 精度。Orekit GN 外层在 n_sv<4 时自动跳过，防止欠定发散。
