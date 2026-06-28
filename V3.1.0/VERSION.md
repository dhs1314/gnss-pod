# V3.1.0 — GRACE-FO PPP-AR POD with Automated QC & Multi-Day Consistency

## Release Date: 2026-06-29

---

## Accuracy (2024-04-29 ~ 05-08, GRACE-FO C, GPS only, 9 days)

| Date | 0.17h 3D RMS | QC Score | Coverage | 0.50h 3D RMS | Flags |
|------|------------|----------|----------|-------------|-------|
| 04-29 | **0.047m** | 0.96A | 86% | 0.397m | 1 SV rejected |
| 04-30 | 0.189m | 0.81B | 68% | 0.541m | 2 SV rejected |
| 05-01 | 0.171m | 0.83B | 72% | 0.392m | clean |
| 05-02 | 0.344m | 0.75B | 64% | 0.434m | 4 SV rejected |
| 05-03 | 0.266m | 0.71B | 56% | 0.366m | 1 gap SV |
| 05-04 | **0.067m** | 0.85A | 83% | 0.320m | 1 SV rejected |
| 05-05 | 0.231m | 0.85B | 81% | 0.813m | 3 SV rejected |
| 05-06 | **0.102m** ★ | 0.89A | 81% | 0.313m | G05 auto-fixed |
| 05-08 | 0.101m | 0.81B | 82% | 0.863m | 3 SV rejected |

**0.17h**: mean=**0.169m**, median=0.171m, best=**0.047m**, worst=0.344m, ≤0.2m: 6/9
**0.50h**: mean=**0.493m**, median=0.397m, best=**0.313m**, worst=0.863m

### Key breakthrough: 05-06 anomaly auto-fixed

G05 had 2 data gaps with 138/157-cycle phase jumps. V3.1 automatically:
1. Detected the gaps (epoch diffs > 90s)
2. Split G05 into 3 independent ambiguity arcs
3. MW jump screening rejected the corrupted segments
4. Huber IRLS down-weighted remaining outliers

**Result: 5.88m → 0.102m (58× improvement)** on 0.17h

---

## What's New vs V3.0.0

### Automated Data Quality System (Phase 24.0)

**File**: `src/data_quality.py` (new, ~565 lines)

4-layer QC architecture:
```
L1 Raw Data   → SV coverage / SNR / MP1 multipath / SP3+CLK integrity
L2 EKF→Batch  → MW stability analysis / auto gap-split / coverage-adaptive reject
L3 Batch LSQ  → Huber IRLS iterative reweighting (k=2.5σ, 3 iterations)
L4 Post-solve → 6-factor weighted quality score (0-1, A/B/C/D/F grade)
```

Each arc outputs: `[QC] score=0.96(A) | flags: 1 rejected SV | actions: 1 SV rejected`

### Robust Estimation

- **SV Gap Detection**: Epoch jumps > 2 (90s) trigger automatic segmentation
  - G05 → G05_S0, G05_S1, G05_S2 (3 independent ambiguity parameters)
- **Huber IRLS** in BatchLinearSolver: `robust_reweight=True`
  - k=2.5σ threshold, 3 iterations, automatic outlier down-weighting
- **MW Jump Screening**: σ_MW > 10 cyc or ≥2 jumps → SV rejected
- **Coverage-Adaptive Ambiguity Priors**: σ ∝ (1-coverage)² for partial-coverage SVs

### Pipeline Improvements

- TurboEdit cycle-slip detection enabled by default (`use_cycle_slip=True`)
- Arc-level MW + OSB NL fixing in BatchOrbitLSQv3.solve()
- Coverage statistics in result dict (avg_cov, n_full_cov, n_partial, n_low_cov)
- CLI: `--min-coverage 0.70` (auto-skip low-coverage arcs)
- CLI: `--fuse-arcs N` (sliding-window arc fusion, experimental)

### Evaluation & Plotting

- 3-panel side-by-side plot (3D RMS curves + Phase RMS + Coverage vs Accuracy scatter)
- QC scores annotated on data points
- Bubble size proportional to SV count in scatter plot
- Date labels on coverage-accuracy points

---

## Architecture (unchanged from V3.0.0)

### Three Processing Pipelines

```
                       GPS1B Obs + SP3 + CLK + DCB + ANTEX + IERS
                                      │
            ┌─────────────────────────┼─────────────────────────┐
            ▼                         ▼                         ▼
  ┌──────────────────┐    ┌──────────────────────┐    ┌──────────────────────┐
  │ ① EKF Sequential │    │ ② BatchLinearSolver  │    │ ③ Orekit GN Outer    │
  │    Filter        │    │    (on EKF orbit)     │    │    Loop              │
  └────────┬─────────┘    └──────────┬───────────┘    └──────────┬───────────┘
           │                        │                           │
           ▼                        ▼                           ▼
    3D RMS: 0.293m           Phase: 0.160m             3D RMS: 0.047m ★
    Phase:  0.276m           (orbit unchanged)         Phase:  0.151m
```

Each pipeline now benefits from the automated QC layer and robust estimation.

---

## CLI Reference

### Standard Multi-Day Validation
```powershell
$env:OREKIT_DATA_PATH = 'd:\prj\gnss_pod\data\orekit'
$env:JAVA_HOME = 'C:\Program Files\JetBrains\PyCharm Community Edition 2024.3.5\jbr'

py eval_5day_orekit.py --dates 2024-04-29,...,2024-05-08 --hours 0.17,0.50

# Coverage-filtered (Path ①): only arcs with >70% SV coverage
py eval_5day_orekit.py --min-coverage 0.70

# Arc fusion (Path ③, experimental): 6 sliding 10-min arcs per hour
py eval_5day_orekit.py --dates 2024-04-29 --hours 0.17 --fuse-arcs 6
```

### EKF + Batch (unchanged from V3.0.0)
```powershell
py run_sequential_pod.py --date 2024-04-29 --hours 0.17 --interval 30 --grace-id C \
  --gravity-nmax 150 --gravity-model GGM05C \
  --sp3-file data/CODE/2024/COD0OPSFIN_20241200000_01D_05M_ORB.SP3 \
  --clk-file data/CODE/2024/COD0OPSFIN_20241200000_01D_30S_CLK.CLK \
  --dcb-file data/CODE/2024/P1P22404.DCB --dcb-p1c1-file data/CODE/2024/P1C12404.DCB \
  --antex-file data/igs14.atx --iers-c04 data/IERS/eopc04_IAU2000.txt \
  --enable-phase-windup --enable-relativity --ar-min-epochs 6 --batch-lsq-v2
```

---

## Key Files

| File | Description | New in V3.1 |
|------|-------------|-------------|
| `src/data_quality.py` | 4-layer QC: coverage/SNR/MP1/MW/product/scoring | **NEW** |
| `src/batch_solver.py` | BatchLinearSolver + IRLS + coverage-adaptive priors | **UPDATED** |
| `src/batch_orbit_v3.py` | 9-param GN + Orekit + arcMW + OSB NL fixing | **UPDATED** |
| `eval_5day_orekit.py` | Multi-day Orekit GN validation + QC + plotting | **UPDATED** |
| `src/sequential_filter.py` | EKF core (unchanged) | — |
| `src/orekit_bridge.py` | Orekit v13 interface (unchanged) | — |

---

## Changelog from V3.0.0

- **Automated QC System**: `src/data_quality.py` — 4-layer data quality framework
- **SV Gap Detection**: Auto-detect epoch gaps >90s → split into independent ambiguity arcs
- **Huber IRLS**: Robust reweighting in BatchLinearSolver (k=2.5σ, 3 iterations)
- **Coverage-Adaptive Priors**: σ_amb ∝ (1-coverage)² for soft prior on partial-coverage SVs
- **TurboEdit Enabled**: `use_cycle_slip=True` by default in validation
- **Arc MW + OSB NL**: Full arc-level MW averaging + OSB non-differenced NL fixing in GN outer loop
- **Quality Score**: 6-factor weighted 0-1 score per arc with A/B/C/D/F grading
- **Multi-Day Plot**: 3-panel side-by-side (3D RMS + Phase RMS + Coverage-vs-Accuracy)
- **CLI Options**: `--min-coverage`, `--fuse-arcs` for arc quality filtering and fusion
- **05-06 Anomaly Fixed**: Automatic gap detection + segmentation → 58× accuracy improvement
- **9-Day Validation**: Confirmed 0.169m mean accuracy across 9 days with full QC reporting

---

## Requirements

- Python 3.12+
- orekit-jpype 13.1.5
- numpy, matplotlib, jpype1
- Java Runtime (JRE 8+)
- OREKIT_DATA_PATH pointing to data/orekit/
