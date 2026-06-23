# GRACE-FO PPP 精密定轨系统

**GRACE-FO Precise Orbit Determination using PPP (Precise Point Positioning)**

验证结果（2024年5月全月，31天）：
- **E（东西）RMS**: 1.96 ± 0.99 cm
- **N（南北）RMS**: 1.66 ± 0.78 cm
- **U（垂直）RMS**: 11.46 ± 1.93 cm
- **3D RMS**: 11.80 ± 1.96 cm（范围 8.67 ~ 17.65 cm）
- 成功率：89%

## 核心算法

| 改正项 | 说明 |
|--------|------|
| **光行时改正** | 迭代 2 次，考虑地球自转 |
| **Sagnac 效应** | `(Ω_E/C) · (x_sat·y_rcv - y_sat·x_rcv)` |
| **相对论时钟** | Einstein 项 + 二阶 Doppler - 13m 重力红移 |
| **天线 PCV** | GPS 卫星天线相位中心变化（高度角依赖，最多 ~5mm） |
| **对流层延迟** | Saastamoinen 模型 + 随机游走湿分量 |
| **电离层-free** | L1/L2 双频组合消除一阶电离层延迟 |

## 文件结构

```
gnss_pod/
├── src/
│   ├── ppp.py           # 核心 PPP 求解器（已验证）
│   ├── sp3_loader.py    # SP3 精密星历加载
│   └── ...
├── run_ppp.py           # 主程序（单日）
├── ppp_may2024.py       # 全月批量处理
├── gen_may_report.py    # HTML 报告生成
├── output/
│   ├── ppp_may_2024-05-*_4h.csv   # 每日结果
│   └── report_may2024.html         # 全月报告
└── data/                # GRACE-FO L1B 数据（需下载）
```

## 使用方法

### 1. 下载 GRACE-FO 数据

```bash
# 编辑 ISDC_URL 为实际下载地址
python3 ppp_may2024.py  # 自动下载并处理 5 月全月
```

### 2. 单日处理

```bash
python3 run_ppp.py --year 2024 --month 5 --day 1 --hours 4 --interval 30
```

### 3. 生成报告

```bash
python3 gen_may_report.py
```

## 推送到 GitLab

```bash
git remote add origin https://gitlab.com/dhs1314/gnss_pod.git
git push -u origin master
```

（首次推送需要 GitLab Personal Access Token，密码方式已废弃）

## 验证方法

- **参考轨道**: ISDC/GFZ GNV1B（~2 cm 精度）
- **过滤条件**: 3D 误差 < 1m（去除极轨发散历元）
- **采样间隔**: 30 秒
- **每段时长**: 4 小时
- **数据来源**: GRACE-FO L1B RL04（ISDC）

## 参考

- ISDC: https://isdc-data.gfz.de/grace-fo/
- GRACE-FO POD: GNV1B 产品（ECEF 坐标，1 Hz）
- GPS 星历: 广播星历（简化圆轨道，演示用）
