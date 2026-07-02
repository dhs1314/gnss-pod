# V3.2.0 -- GRACE-FO PPP-AR POD: Self-Consistent + Automated QC

**Precision Orbit Determination without external orbit products.**

## Quick Start

```powershell
$env:OREKIT_DATA_PATH = "d:/prj/gnss_pod/data/orekit"
$env:JAVA_HOME = "C:/Program Files/JetBrains/PyCharm Community Edition 2024.3.5/jbr"

# 9-day validation (self-consistent, no GNV1B required)
py eval_5day_orekit.py --dates 2024-04-29,...,2024-05-08 --hours 0.17,0.50
```

## Accuracy (9 days, 2024-04-29 ~ 05-08)

| Arc | Best | Mean | Median |
|-----|------|------|--------|
| 0.17h | **0.047m** | **0.169m** | 0.171m |
| 0.50h | **0.313m** | **0.493m** | 0.397m |

Each arc scored 0-1 with A/B/C/D/F grade. Best arc: 0.96A.

## What's New vs V3.1.0

- **Self-Consistent Initial Orbit** (`src/code_orbit.py`): Pseudorange GN eliminates GNV1B dependency
- **Code-Level Outlier Detection**: 5m coarse detection before precise pipeline
- **4-Step Pipeline**: Code WLS -> Code GN -> EKF -> Batch AR + Orekit GN
- **Broadcast Ephemeris Simulation**: `--broadcast` evaluates accuracy limits

## Requirements

- Python 3.12+ with orekit-jpype 13.1.5, numpy, matplotlib, jpype1
- Java Runtime 8+
- Orekit data files in data/orekit/

## License

Research code. See VERSION.md for full changelog.
