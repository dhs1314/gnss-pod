# V3.1.0 — GRACE-FO PPP-AR POD with Automated QC

**Precision Orbit Determination for GRACE-FO using GNSS carrier-phase ambiguity resolution.**

## Quick Start

```powershell
$env:OREKIT_DATA_PATH = "d:/prj/gnss_pod/data/orekit"
$env:JAVA_HOME = "C:/Program Files/JetBrains/PyCharm Community Edition 2024.3.5/jbr"

# 9-day Orekit GN validation with automated QC
py eval_5day_orekit.py --dates 2024-04-29,...,2024-05-08 --hours 0.17,0.50
```

## Accuracy (2024-04-29 ~ 05-08, 9 days)

| Arc | Best | Mean | Median |
|-----|------|------|--------|
| 0.17h | **0.047m** | **0.169m** | 0.171m |
| 0.50h | **0.313m** | **0.493m** | 0.397m |

**Key**: 05-06 anomaly auto-fixed from 5.88m → 0.102m (58×). Each arc scored (0-1) with A/B/C/D/F grade.

## What's New vs V3.0.0

- **Automated QC** (`src/data_quality.py`): 4-layer data quality with per-arc scoring
- **SV Gap Detection + Auto-Split**: Corrupted SVs automatically segmented into independent ambiguity arcs
- **Robust Estimation**: Huber IRLS (k=2.5σ) + coverage-adaptive ambiguity priors
- **Arc MW + OSB NL Fixing**: Full-arc Melbourne-Wübbena averaging with CODE satellite biases

## Architecture

Three pipelines: EKF Sequential, Batch Linear Solver, Orekit GN Outer Loop.
See VERSION.md for full documentation and changelog.

## Requirements

- Python 3.12+ with orekit-jpype 13.1.5, numpy, matplotlib, jpype1
- Java Runtime 8+
- Orekit data files in data/orekit/

## Data

Required data products (download separately):
- GPS1B, GNV1B from JPL PO.DAAC
- CODE SP3/CLK/DCB from AIUB FTP
- IGS ANTEX, IERS C04, GGM05C.gfc from respective archives

## License

Research code. See VERSION.md for full changelog.
