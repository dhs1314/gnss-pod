# V3.0.0 — GRACE-FO PPP-AR POD

**Precision Orbit Determination for GRACE-FO using GNSS carrier-phase ambiguity resolution.**

## Quick Start

```powershell
# Set environment
$env:OREKIT_DATA_PATH = "d:/prj/gnss_pod/data/orekit"
$env:JAVA_HOME = "C:/Program Files/JetBrains/PyCharm Community Edition 2024.3.5/jbr"

# Run EKF + Batch solver (Phase RMS report)
py run_sequential_pod.py --date 2024-04-29 --hours 0.17 --interval 30 --grace-id C `
  --gravity-nmax 150 --gravity-model GGM05C `
  --sp3-file data/CODE/2024/COD0OPSFIN_20241200000_01D_05M_ORB.SP3 `
  --clk-file data/CODE/2024/COD0OPSFIN_20241200000_01D_30S_CLK.CLK `
  --dcb-file data/CODE/2024/P1P22404.DCB --dcb-p1c1-file data/CODE/2024/P1C12404.DCB `
  --antex-file data/igs14.atx --iers-c04 data/IERS/eopc04_IAU2000.txt `
  --enable-phase-windup --enable-relativity --ar-min-epochs 6 `
  --batch-lsq-v2

# 5-day Orekit GN validation
py eval_5day_orekit.py
```

## Accuracy (2024-04-29, GRACE-FO C)

| Method | 0.17h | 0.50h |
|--------|-------|-------|
| EKF Sequential Filter | 0.293m | 0.986m |
| Orekit GN Outer Loop | **0.043m** | **0.409m** |

## Architecture

Three pipelines: EKF Sequential Filter, Batch Linear Solver, Orekit GN Outer Loop.
See VERSION.md for full documentation.

## Requirements

- Python 3.12+ with orekit-jpype 13.1.5, numpy, matplotlib, jpype1
- Java Runtime 8+
- Orekit data files in data/orekit/

## Data

Required data products not included (download separately):
- GPS1B, GNV1B from JPL PO.DAAC
- CODE SP3/CLK/DCB from AIUB FTP
- IGS ANTEX, IERS C04 from respective archives
- GGM05C.gfc from ICGEM

## License

Research code. See VERSION.md for full changelog.
