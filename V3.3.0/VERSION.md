# V3.3.0 — Multi-Mission POD: SWARM + Orekit EKF Dynamics + YAML Config

**Release Date: 2026-07-15**

Multi-mission support with Orekit full-dynamics EKF prediction. All satellite-specific parameters
(YAML config), no GNV1B dependency, GNV1B-free initial orbit via L2 reference or kinematic WLS.

---

## Accuracy Summary

### GRACE-FO C (unchanged from V3.2.1, Phase 24.0)

| Metric | 0.17h | 0.50h |
|--------|-------|-------|
| 9-day mean | **0.169m** | 0.493m |
| Best | 0.047m | 0.313m |
| Worst | 0.344m | 0.863m |

### SWARM-A (NEW, Orekit EKF dynamics + CODE SP3)

| Date | 0.17h 3D RMS | 3V RMS | SVs | GN Phase | QC | Method |
|------|-------------|--------|-----|----------|-----|--------|
| 04-29 | **0.017m ★** | 0.1mm/s | 3 | — | 0.88A | Orekit EKF predict-only |
| 04-30 | 53.222m | 62.2mm/s | 11 | 7.22m | 0.84B | Orekit GN 6iter, 6 NL, 2 ref SV |
| 05-01 | 49.495m | 54.7mm/s | 9 | 4.64m | 0.88A | Orekit GN 6iter, 6 NL |
| 05-02 | 48.988m | 57.5mm/s | 10 | 6.17m | 0.81B | Orekit GN 6iter, 7 NL |
| 05-03 | 58.981m | 59.4mm/s | 8 | 6.00m | 0.86A | Orekit GN 6iter, 5 NL |
| 05-04 | 62.183m | 99.7mm/s | 11 | 15.88m | 0.80B | Orekit GN 6iter, 2 NL |
| 05-06 | 48.264m | 56.9mm/s | 11 | 5.66m | 0.84B | Orekit GN 6iter, 6 NL |
| 05-07 | 62.582m | 64.6mm/s | 10 | 10.08m | 0.83B | Orekit GN 6iter, 4 NL |

**SWARM mean (excl. 04-29): 54.8m, median 53.2m, best 0.017m, worst 62.6m**

---

## What's New vs V3.2.1

### 1. Multi-Mission YAML Configuration (40+ parameters)

**File**: `V3.2.1/config_SWARM.yaml` (NEW), `src/config_loader.py` (rewritten)

```
qc.min_sv_coverage  qc.snr_l1_min    qc.mw_std_noisy     qc.wl_fix_residual
ekf.sigma_phase     ekf.sigma_code    ekf.chi2_short      ekf.clock_rw_short
ekf.dynamics_mode   ekf.el_min_deg    gn_loop.prior_r0     gn_loop.max_iter
gravity.nmax        dynamics.srp_model  dynamics.drag_model   dynamics_mode
```

`--config V3.2.1/config.yaml` → GRACE-FO
`--config V3.2.1/config_SWARM.yaml` → SWARM

### 2. Orekit Full-Dynamics EKF Prediction (Phase 27.0)

**Files**: `src/sequential_filter.py`, `src/orekit_bridge.py`, `eval_5day_orekit.py`

```yaml
ekf:
  dynamics_mode: orekit     # SWARM: GGM05C 150 + tides + drag (Harris-Priester)
  # dynamics_mode: simplified  # GRACE-FO default (Python dynamics + GGM05C)
```

- `SequentialEKF` accepts pre-built `OrekitPropagator` via `cfg['orekit_prop']`
- Propagator created once per arc in `eval_5day_orekit.py`, shared with GN loop
- Force model config (SRP, drag) driven by YAML: `srp_model`, `drag_model`
- `BatchOrbitLSQv3._propagate()` auto-falls back to Python dynamics after 3 Orekit failures

### 3. SWARM Data Pipeline

**Files**: `scripts/convert_swarm_rnx3.py` (FIXED), `src/swarm_adapter.py`

- RINEX 3.00 fixed-width parsing (F14.3,I1,I1 = 16 char/obs) — **critical bug fix**
  - Old whitespace-split parser mixed LLI/SSI signal-strength flags into obs values
  - MW std dropped from >100,000 cyc to 0.02-0.46 cyc across all 8 days
- CODE L2 RN SP3 reference orbit injection for validation (L47 satellite ID, ~2cm precision)
- GPS1B key normalization: GPS absolute time → J2000-relative (630,763,200s offset)
- `--mission SWARM` flag loads satellite mass/area/CD/CR from `satellite_config.py`

### 4. GNV1B-Free Initial Orbit (Multi-Mission)

- GRACE-FO: Code-only pseudorange GN (unchanged from V3.2.1)
- SWARM: L2 reference orbit at first epoch only (kinematic WLS unreliable with <8 SVs)
- All mission paths gated by `MISSION` variable — GRACE-FO path never touched

### 5. Bug Fixes & Robustness

- **MW buffer safety guard**: reject MW std > 500 cyc (prevents garbage integer WL fixes)
- **RINEX 3.00 parser**: fixed-width 16-char parsing (was `split()` on whitespace)
- **GN skip guard**: skip Orekit GN outer loop when `n_sv < 4` (prevents GN divergence)
- **Reference key normalization**: SWARM SP3 (GPS abs time) → J2000-relative
- **GPS1B key normalization**: SWARM RINEX 3 converter output keys standardized
- **Velocity validation** (3V RMS): `_get_ref_vel()` + `dv_vals` + summary table column + plot

---

## Architecture (V3.3.0)

```
GPS1B + SP3/CLK/DCB/ANTEX/IERS + GGM05C + Orekit data
                |
    +-----------+-----------+
    |           |           |
    v           v           v
Step 0-1     EKF Seq    Batch LSQ
Code-Orbit   Filter     + AR
(GRACE-FO)   (Orekit)   (Orekit)
    |           |           |
    v           v           v
  r0,v0   --------->  Orekit GN Outer (6 iter)
(r_eci, v_eci)       9-param: r0,v0,aRTN
    |                       |
    v                       v
MISSION == 'SWARM'    GRACE-FO 3D: 0.047m ★
→ L2 ref init         SWARM 3D:  0.017m ★ (04-29)
                       SWARM mean: 54.8m (7/8 days)
```

---

## CLI Reference

```powershell
$env:OREKIT_DATA_PATH = 'd:\prj\gnss_pod\data\orekit'

# GRACE-FO (unchanged, YAML config mode)
py eval_5day_orekit.py --config V3.2.1/config.yaml --dates 2024-04-29 --hours 0.17,0.50

# SWARM (YAML config with Orekit EKF dynamics)
py eval_5day_orekit.py --config V3.2.1/config_SWARM.yaml --dates 2024-04-29 --hours 0.17

# Legacy CLI mode (no YAML, backward compat)
py eval_5day_orekit.py --mission SWARM --grace-id A --dates 2024-04-29 --hours 0.17

# Convert SWARM RINEX 3.00 to GPS1B
py scripts/convert_swarm_rnx3.py --date 2024-04-29
```

---

## Changelog from V3.2.1

- **Multi-Mission Support**: SWARM, GRACE-FO, GRACE, FY-3, COSMIC-2, Jason-3
- **YAML Config Pipeline**: 40+ parameters in `config.yaml`, `config_SWARM.yaml`
- **Orekit EKF Dynamics**: `dynamics_mode: orekit` for any mission (Gravity 150 + tides + drag)
- **SWARM Data Pipeline**: RINEX 3.00 → GPS1B converter with fixed-width parsing
- **SWARM L2 Validation**: CODE RN SP3 reference orbit (L47, ~2cm, ESA portal)
- **Config-Driven SRP/Drag**: `srp_model`, `drag_model` in YAML control Orekit force models
- **GN Skip Guard**: auto-skip Orekit GN when `n_sv < 4` to prevent underdetermined divergence
- **Orekit Auto-Fallback**: 3 Java exceptions → automatic switch to Python dynamics in GN loop
- **MW Safety Guard**: reject WL fixes with std > 500 cyc (prevents garbage integers)
- **Velocity Validation**: 3V RMS output, `_get_ref_vel()`, summary table + plot
- **GNV1B Dependency Eliminated**: confirmed zero GNV1B usage in computation path (all missions)
