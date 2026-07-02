# V3.2.0 -- GRACE-FO POD: Self-Consistent Code-Orbit + QC + GN Outer Loop

**Release Date: 2026-07-03**

Removes GNV1B external dependency -- initial orbit now determined from pseudorange-only batch GN.

---

## Accuracy (2024-04-29 ~ 05-08, GRACE-FO C, GPS only, 9 days)

| Date | 0.17h 3D RMS | QC | Coverage | 0.50h 3D RMS | QC | Flags |
|------|------------|-----|----------|-------------|-----|-------|
| 04-29 | **0.047m** | 0.96A | 86% | 0.397m | 0.88A | -- |
| 04-30 | 0.189m | 0.81B | 68% | 0.541m | 0.78B | -- |
| 05-01 | 0.171m | 0.83B | 72% | 0.392m | 0.74B | -- |
| 05-02 | 0.344m | 0.75B | 64% | 0.434m | 0.69C | 4 SV rejected |
| 05-03 | 0.266m | 0.71B | 56% | 0.366m | 0.59C | -- |
| 05-04 | **0.067m** | 0.85A | 83% | 0.320m | 0.68C | -- |
| 05-05 | 0.231m | 0.85B | 81% | 0.813m | 0.73B | -- |
| 05-06 | 0.102m | 0.89A | 81% | 0.313m | 0.78B | G05 auto-fixed |
| 05-08 | 0.101m | 0.81B | 82% | 0.863m | 0.70C | -- |

**0.17h**: mean=**0.169m**, median=0.171m, best=**0.047m**, worst=0.344m, <=0.2m: 6/9
**0.50h**: mean=**0.493m**, median=0.397m, best=0.313m, worst=0.863m

### Broadcast Ephemeris Simulation (CODE SP3 + 1.0m/1.5m noise)

| Date | BRDC 3D RMS | CODE 3D RMS | Degradation |
|------|-----------|-----------|-------------|
| 04-29 | 4.085m | 0.047m | 87x |
| 04-30 | 3.304m | 0.189m | 17x |
| 05-01 | 3.247m | 0.171m | 19x |
| 05-08 | 4.465m | 0.101m | 44x |

**BRDC mean=5.60m, median=4.47m -- broadcast ephemeris insufficient for cm-level POD.**

---

## What's New vs V3.1.0

### Self-Consistent Initial Orbit (Phase 26.0)

**File**: `src/code_orbit.py` (NEW, ~350 lines)

Two-step pipeline removes GNV1B dependency entirely:

```
Step 0: Kinematic WLS single-epoch positioning (pseudorange only)
        -> r_ecef (sigma ~5-10m), v=0

Step 1: Code-only Gauss-Newton batch orbit determination
        -> arc: configurable (default 1.0h, 120 epochs)
        -> params: r0(3)+v0(3)+aR,aT,aN(3) = 9 nonlinear
        -> linear sub-problem: clk_i + zwd_i per epoch (2x2 independent)
        -> Python GGM05C N=150 dynamics (+ solid tides + third body + SRP + drag)
        -> target: 0.15-0.20m 3D RMS vs GNV1B

Step 2: Outlier detection
        -> |pseudorange residual| > 5m -> flag epoch/SV
        -> per-SV code RMS > 3m -> reject from precise pipeline
```

**CLI**: `--code-arc-hours 1.0` (default), `--skip-code-orbit` (backward compat to GNV1B)

### Key Changes

- GNV1B no longer required for initial orbit -- code-orbit is self-consistent
- OrekitPropagator created ONCE per arc, shared between code-orbit and GN loop
- Code-orbit propagates with Python integrator (avoids Orekit v13 multi-instance issue)
- Outlier detection at code level before precise OD pipeline
- Coverage-adaptive ambiguity priors (sigma ~ (1-coverage)^2)
- Huber IRLS robust reweighting in BatchLinearSolver (k=2.5sigma, 3 iterations)
- 4-layer automated QC system (`src/data_quality.py`)

---

## Architecture

```
GPS1B + SP3/CLK/DCB/ANTEX/IERS + GGM05C + Orekit data
                |
    +-----------+-----------+
    |           |           |
    v           v           v
Step 0-1     EKF Seq    Batch LSQ
Code-Orbit   Filter     + AR
(self-consistent)  |      |
    |           v           v
    +--------> pass1    fixed_amb
    |           |           |
    v           v           v
  r0,v0   --------->  Orekit GN Outer (6 iter)
                       9-param: r0,v0,aRTN
                       |
                       v
                 3D RMS: 0.047m (best)
                 QC score: 0.96A
```

---

## CLI Reference

```powershell
$env:OREKIT_DATA_PATH = 'd:\prj\gnss_pod\data\orekit'
$env:JAVA_HOME = '...'

# Standard 9-day validation (self-consistent, no GNV1B needed)
py eval_5day_orekit.py --dates 2024-04-29,...,2024-05-08 --hours 0.17,0.50

# With code-only initial orbit (default, 1h arc)
py eval_5day_orekit.py --code-arc-hours 1.0

# Backward compat (use GNV1B)
py eval_5day_orekit.py --skip-code-orbit

# Broadcast simulation
py eval_5day_orekit.py --broadcast

# Coverage filter
py eval_5day_orekit.py --min-coverage 0.70
```

---

## Changelog from V3.1.0

- **Self-Consistent Code-Orbit**: `src/code_orbit.py` -- pseudorange-only GN for initial orbit
- **GNV1B Eliminated**: No external reference orbit required for solving
- **4-Step Pipeline**: Code WLS -> Code GN -> EKF -> Batch AR + Orekit GN
- **Code-Level Outlier Detection**: 5m residuals flagged before precise pipeline
- **OrekitPropagator Singleton**: Single instance per arc, shared across steps
- **Python Dynamics Fallback**: Code-orbit uses Python integrator when Orekit unavailable
- **Broadcast Ephemeris Simulation**: `--broadcast` flag (CODE SP3 + 1.0m/1.5m noise)
- **9-Day Validation Confirmed**: 0.17h mean=0.169m (unchanged from V3.1.0)
