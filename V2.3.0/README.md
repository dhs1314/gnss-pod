# V2.3.0 — GRACE-FO PPP-AR POD

## Quick Start

```powershell
powershell -File run_v230.ps1 -Mode 017     # 0.17h EKF (0.293m)
powershell -File run_v230.ps1 -Mode batch   # 0.17h + Batch solver
powershell -File run_v230.ps1 -Mode all     # 0.17h + 0.5h
powershell -File run_v230.ps1 -Mode fullday # 全天线验证
```

## Accuracy

| Arc | EKF Orbit | Phase RMS (EKF) | Phase RMS (Batch) |
|-----|-----------|-----------------|-------------------|
| 0.17h | 0.293m | 0.276m | **0.160m** (-42%) |
| 0.5h | 0.986m | 1.207m | **0.211m** (-83%) |

## New in V2.3.0 (vs V2.2.4)

- **Framework v3 Batch Solver**: Fixed EKF orbit + global clock/amb LSQ
- **Arc-based AR**: Continuous tracking arc ambiguity resolution
- **Satellite Config DB**: 11 LEO satellite physical parameters
- **CLK1B parser**: GRACE-FO receiver clock offset (USO, 0.03ns precision)
- **CODE OSB parser**: SINEX BIAS format for undifferenced NL fixing
- **Multi-date validation**: Batch processing for multiple dates

## Requirements

- Python 3.12.6, numpy 2.2.6, scipy 1.17.1, matplotlib 3.10.9
- Data directory at `../data/` (one level up from V2.3.0/)
- CODE precision products (SP3, CLK, DCB) for processing date

## Directory

```
V2.3.0/
├── run_sequential_pod.py   ← Main entry point
├── run_v230.ps1            ← One-click run script
├── validate_v224.py        ← Full-day assessment
├── fullday_batch_v3.py     ← Batch solver full-day
├── eval_day2.py            ← 6-hour assessment
├── VERSION.md              ← Full version documentation
├── README.md               ← This file
├── src/                    ← 35 Python source files
└── references/             ← Improvement roadmap
```

## Run Commands

### EKF (best 0.293m at 0.17h)
```powershell
py -3.12 run_sequential_pod.py \
  --date 2024-04-29 --hours 0.17 --interval 30 --grace-id C \
  --dynamics-mode simplified \
  --sp3-file ../data/CODE/2024/COD0OPSFIN_20241200000_01D_05M_ORB.SP3 \
  --clk-file ../data/CODE/2024/COD0OPSFIN_20241200000_01D_30S_CLK.CLK \
  --dcb-file ../data/CODE/2024/P1P22404.DCB \
  --antex-file ../data/igs14.atx \
  --iers-c04 ../data/IERS/eopc04_IAU2000.txt \
  --enable-phase-windup --enable-relativity \
  --ar-min-epochs 6 --gravity-nmax 90
```

### Batch Solver Phase RMS Assessment
```powershell
# Add --batch-lsq-v2 to any of the above commands
```

### Arc-based AR (0.246m)
```powershell
# Add --arc-ar
```

### CODE OSB (experimental)
```powershell
# Add --osb-file ../data/CODE/2024/COD0OPSFIN_20241200000_01D_01D_OSB.BIA
```

### Orekit Dynamics (optional)
```powershell
# Set JAVA_HOME + OREKIT_DATA_PATH, add --dynamics-mode orekit
```
