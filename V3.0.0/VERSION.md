# V3.0.0 — GRACE-FO PPP-AR POD with Orekit GN Outer Loop

## Release Date: 2026-06-23

---

## Accuracy (2024-04-29, GRACE-FO C, GPS only)

| Arc | EKF 3D RMS | EKF Phase | Batch Phase | **Orekit GN 3D RMS** | GN Phase |
|-----|-----------|-----------|-------------|----------------------|----------|
| 0.17h | 0.293m | 0.276m | 0.160m | **0.043m** | 0.151m |
| 0.50h | 0.986m | 1.207m | 0.211m | **0.409m** | 0.233m |

### 5-Day Validation (2024-04-29 ~ 05-03)

| Date | 0.17h | 0.50h | SVs(0.17/0.5) |
|------|-------|-------|----|
| 04-29 | **0.042m** | 0.409m | 12/16 |
| 04-30 | 0.201m | 0.539m | 13/18 |
| 05-01 | 0.169m | 0.357m | 13/19 |
| 05-02 | 0.372m | 0.450m | 14/19 |
| 05-03 | 0.288m | 0.735m | 11/19 |

0.17h: mean=0.214m, median=0.201m, best=**0.042m**
0.50h: mean=0.498m, median=0.450m, best=**0.357m**

---

## Architecture

### Three Processing Pipelines

```
                         GPS1B Obs + SP3 + CLK + DCB + ANTEX + IERS
                                        |
              ┌─────────────────────────┼─────────────────────────┐
              ▼                         ▼                         ▼
    ┌──────────────────┐    ┌──────────────────────┐    ┌──────────────────────┐
    │ ① EKF Sequential │    │ ② BatchLinearSolver  │    │ ③ Orekit GN Outer    │
    │    Filter        │    │    (on EKF orbit)     │    │    Loop              │
    └────────┬─────────┘    └──────────┬───────────┘    └──────────┬───────────┘
             │                        │                           │
             ▼                        ▼                           ▼
      3D RMS: 0.293m           Phase: 0.160m             3D RMS: 0.043m ★
      Phase:  0.276m           (orbit unchanged)         Phase:  0.151m
```

### ① EKF Sequential Filter (V2.2.4 baseline)

**File**: `src/sequential_filter.py`

State vector (ECI): `[r(3), v(3), aR, aT, aN, zwd, clk, amb_1..amb_Nsv]`

Per-epoch processing:
- `predict()`: Orbit propagation (Python GGM05C Nmax=150 integrator, 30s steps with 10s sub-steps)
- `process_epoch()`: GPS code/phase measurement update
- MW combination accumulation → WL integer fixing (>=6 epochs, sigma < 0.35 cyc)
- SD NL bootstrapping

Key parameters (adaptive):
- `elev_exp_phase=1.0`, `elev_exp_code=` 1.0 (0.17h) / 0.70 (>=0.3h)
- `clock_rw=` 0.0004 (0.17h) / 0.001 (>=0.3h)
- `chi2_threshold=` 25 (0.17h) / 100 (>=0.3h)
- `ar_min_epochs=6`

Dynamics: Python GGM05C Nmax=150 spherical harmonics + J2 + third-body (Sun, Moon) + drag + SRP

### ② BatchLinearSolver (Framework v3)

**File**: `src/batch_solver.py`

Given a fixed EKF orbit trajectory, solves ALL epochs simultaneously:
```
  P_if = geo_full + clk_i + zwd_i * mf + noise
  L_if = geo_full + clk_i + zwd_i * mf + amb_k + noise
```

Parameters: clk(N_epoch) + zwd(N_epoch) + amb(N_sv), solved by single normal equation system.

Eliminates EKF sequential cold-start problem → phase RMS improvement -42% (0.17h) to -83% (0.5h).

Optional features:
- OSB NL fixing (`--osb-file`): CODE satellite biases for undifferenced N1 integers
- CLK1B prior: epoch-to-epoch clock difference constraints (not applicable for JPL data)

### ③ Orekit GN Outer Loop (Phase 20.0, ★ Produces 3D RMS results)

**Files**: `src/batch_orbit_v3.py`, `src/orekit_bridge.py`

Two-step Gauss-Newton optimization of orbit initial state:

**Step 1**: BatchLinearSolver fixes IF ambiguities on initial orbit (Orekit pure propagation, ~0.14m vs GNV1B)

**Step 2**: GN iterations (6 rounds):
- Orekit continuous arc propagation (Nmax=150 gravity + solid tides + ocean tides + third-body + drag + SRP + relativity)
- Compute clock-differenced residuals
- Build analytical Jacobian: `J = ∂res/∂(r0, v0, aR, aT, aN)`
- Solve normal equations with priors: `H·dx = -g`
- Line-search backtracking (α = 1.0 → 0.125) to guarantee cost reduction

Orekit force models:
```java
HolmesFeatherstoneAttractionModel (GGM05C Nmax=150)
+ SolidTides (IERS 2010)
+ OceanTides (FES2004, degree 50)
+ ThirdBodyAttraction (Sun, Moon)
+ SolarRadiationPressure (IsotropicRadiationSingleCoefficient)
+ DragForce (IsotropicDrag + SimpleExponentialAtmosphere)
+ Relativity (Schwarzschild)
```

Speed optimizations (13.7x vs initial implementation):
- `propagate_continuous_arc()`: Single Orekit propagator reused across epochs
- Phi pre-computation cache: STM changes <0.5% for r0 adjustments <100m
- `_precompute_sat_eci()`: Satellite ECI positions computed once for all GN iterations

GN time: ~45s (0.17h), ~145s (0.5h)

---

## Data Products

| Product | Source | Path | Purpose |
|---------|--------|------|---------|
| GPS1B .rnx/.pkl | JPL/GRACE-FO L1B | data/gracefo/{year}/{date}/ | GPS L1/L2/P1/P2 raw obs |
| GNV1B .txt | JPL/GRACE-FO L1B | data/gracefo/{year}/{date}/ | Reference orbit (~2cm) |
| CODE Final SP3 | CODE/AIUB | data/CODE/{year}/ | GPS precise orbit (2.5cm) |
| CODE 30s CLK | CODE | data/CODE/{year}/ | GPS precise satellite clocks |
| CODE P1P2/P1C1 DCB | CODE | data/CODE/{year}/ | Code biases |
| IGS14.atx | IGS | data/ | GPS satellite antenna PCO |
| IERS C04 ERP | IERS | data/IERS/ | Earth orientation parameters |
| GGM05C.gfc | ICGEM | data/gravity/ | Static gravity field Nmax=150 |
| EIGEN-6C4_N200.gfc | ICGEM | data/gravity/ | Alternative gravity model Nmax=200 |
| Orekit data | orekit.org | data/orekit/ | EOP, CSSI weather, DE-440 ephemerides |

---

## CLI Reference (run_sequential_pod.py)

### EKF Only (V2.2.4 baseline)
```powershell
py run_sequential_pod.py --date 2024-04-29 --hours 0.17 --interval 30 --grace-id C \
  --gravity-nmax 150 \
  --sp3-file data/CODE/2024/COD0OPSFIN_20241200000_01D_05M_ORB.SP3 \
  --clk-file data/CODE/2024/COD0OPSFIN_20241200000_01D_30S_CLK.CLK \
  --dcb-file data/CODE/2024/P1P22404.DCB \
  --dcb-p1c1-file data/CODE/2024/P1C12404.DCB \
  --antex-file data/igs14.atx \
  --iers-c04 data/IERS/eopc04_IAU2000.txt \
  --enable-phase-windup --enable-relativity --ar-min-epochs 6
```

### EKF + Batch Solver (Phase RMS report)
```powershell
# Add --batch-lsq-v2 to the above command
```

### Arc-based Ambiguity Resolution
```powershell
# Add --arc-ar to the EKF command
```

---

## Key New Files vs V2.3.0

| File | Description |
|------|-------------|
| `src/orekit_bridge.py` | Orekit v13.1.5 interface: propagate_continuous_arc, DragForce, SRP |
| `src/batch_orbit_v3.py` | Gauss-Newton outer loop with Orekit + piecewise RTN + line-search |
| `src/batch_solver.py` | BatchLinearSolver with OSB NL fixing + kinematic mode + CLK1B |
| `eval_5day_orekit.py` | 5-day Orekit GN validation script |
| `test_orekit_gn.py` | Single-arc Orekit GN test harness |
| `download_orekit_data.py` | Orekit data package downloader (~20 MB) |

---

## Requirements

- Python 3.12+
- orekit-jpype 13.1.5
- numpy, matplotlib, jpype1
- Java Runtime (JRE 8+)
- OREKIT_DATA_PATH pointing to data/orekit/

## Installation

```powershell
pip install orekit-jpype numpy matplotlib jpype1
$env:OREKIT_DATA_PATH = "d:/prj/gnss_pod/data/orekit"
$env:JAVA_HOME = "C:/Program Files/..."  # path to JRE
```

---

## Changelog from V2.3.0

- **Orekit GN Outer Loop**: 9-parameter (r0, v0, aR, aT, aN) Gauss-Newton with line-search
- **Orekit v13 API Fix**: DragForce wrapper, IsotropicRadiationSingleCoefficient SRP
- **Continuous Arc Propagation**: Single Orekit propagator reused across epochs (2x speedup)
- **Phi Pre-computation Cache**: Python STM cached across GN iterations (10x speedup)
- **sat_eci Pre-computation**: Satellite positions transformed once per GN solve (2x speedup)
- **Piecewise RTN**: Per-segment empirical accelerations (3-min segments)
- **Cd/CR Estimation**: Analytic Jacobian via force model linearity
- **EIGEN-6C4 Gravity Model**: Downloaded, truncated to Nmax=200, integrated with CLI switch
- **Orekit Data Package**: CSSI weather, SOLFSMY, DE-440 ephemerides downloaded
- **HarrisPriester + NRLMSISE00 Drag Models**: Integrated with Orekit v13 API
- **Kinematic Batch Solver**: Per-epoch ECEF position estimation (experimental)
- **OSB NL Fixing**: CODE satellite bias undifferenced integer resolution
- **CLK1B Clock Prior**: USO clock difference constraints
- **5-Day Validation Pipeline**: `eval_5day_orekit.py`
- **GitHub Push**: Repository at https://github.com/dhs1314/gnss-pod

Total speedup vs V2.3.0 GN loop: **13.7x** (624s → 45.6s for 0.17h arc)
