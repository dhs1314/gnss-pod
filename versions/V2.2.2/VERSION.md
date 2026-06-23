# V2.2.2 — Sequential EKF + WL-AR PPP

**Date**: 2026-06-10

## Key Features

### Phase 3.0: PPP Ambiguity Resolution (WL-only)
- Self-calibrated receiver WL bias estimation (median fractional part across tracked SVs)
- Melbourne-Wübbena combination with per-SV accumulation and fixing
- WL-fixed ambiguity state correction and covariance tightening (P_amb → 0.10 m²)
- SV pruning for multi-hour arcs (remove SVs unseen for prune_timeout seconds)

### Phase 2.3: Measurement Corrections
- Satellite antenna PCO correction (ANTEX IGS14.atx)
- Phase wind-up correction (Wu et al. 1993)
- Relativistic Shapiro delay correction
- DCB correction for IF code observations (α × P1-C1 only)
- IERS C04 Earth Orientation Parameters (via astropy)
- Cycle-slip detection (TurboEdit MW+GF) — disabled by default

### Phase 2.2: Sequential EKF
- STM-based state and covariance propagation
- Gravity: GGM05C Nmax=90
- Empirical RTN accelerations (Gauss-Markov)
- Troposphere: Saastamoinen ZHD + GMF mapping, ZWD random walk
- Clock: random walk (TCXO, σ_rw=0.032 m/√s)
- Scalar Kalman updates with chi² innovation testing

## State Vector (ECI)
```
[r(3), v(3), aR, aT, aN, zwd, clk, amb_1..amb_Nsv]
```

## Verified Accuracy (2024-04-29, GRACE-FO C)

| Arc   | Float PPP (Phase 2.3) | WL-AR (Phase 3.0) | Improvement |
|-------|-----------------------|-------------------|-------------|
| 0.17h | 0.936m                | **0.803m**        | 14%         |
| 0.5h  | 1.816m                | **1.005m**        | 45%         |
| 2h    | 6.63m                 | **6.06m**         | 9%          |

WL fix rate: 92-94% (0.35-cyc threshold, min 5-10 epochs)
Receiver WL bias: +0.22 ± 0.01 cyc

## Limitations
- **Dynamics**: GGM05C-only gravity limits arc length to ~1h. 2h+ arcs diverge due to unmodeled solid tides, ocean tides, and higher-degree gravity.
- **NL ambiguity**: WL-only fixing leaves N1 float (σ ≈ 0.3m/SV), coupling into position.
- **Cycle slips**: TurboEdit GF threshold too tight for GRACE data — keep disabled for best results.
- **sigma_acc**: 1e-3 is optimal for GGM05C-only. Lower values cause filter divergence.

## Usage
```bash
# 0.5h arc with all corrections + WL-AR
py -3.12 run_sequential_pod.py \
  --date 2024-04-29 --hours 0.5 --interval 30 --grace-id C \
  --sp3-file data/CODE/2024/COD0OPSFIN_20241200000_01D_05M_ORB.SP3 \
  --clk-file data/CODE/2024/COD0OPSFIN_20241200000_01D_30S_CLK.CLK \
  --antex-file data/igs14.atx \
  --dcb-file data/CODE/2024/P1P22404.DCB \
  --iers-c04 data/IERS/eopc04_IAU2000.txt \
  --ar-min-epochs 6 --enable-phase-windup --enable-relativity
```

## Files
```
V2.2.2/
  run_sequential_pod.py        # Main entry: sequential EKF POD
  run_kf_ppp.py                # Kinematic KF (no dynamics)
  run_batch_pod.py             # Batch LSQ POD
  run_gps1b_rnx_ppp.py         # Single-epoch PPP
  VERSION.md                   # This file
  src/
    sequential_filter.py       # SequentialEKF + MWBuffer (WL-AR)
    precision_products.py      # DCB, IERS C04, ANTEX readers
    ambiguity.py               # MW/NL/IF ambiguity functions
    measurement_corrections.py # Wind-up, relativity, PCO
    cycle_slip.py              # TurboEdit MW+GF detector
    orbit_integrator.py        # RK4/DP8 integrator with STM
    orbit_dynamics.py          # GGM05C gravity + 3rd-body + SRP + drag
    gravity_model.py           # ICGEM GFC reader + spherical harmonic synthesis
    solid_tides.py             # Solid Earth tides (IERS 2010)
    relativity_orbit.py        # Relativistic orbit corrections
    coordinates.py             # ECEF↔ECI transforms
    troposphere.py             # Saastamoinen + GMF
    batch_estimator.py         # Batch LSQ estimator
    batch_lsq.py               # LSQ system builder
    gps1b_loader.py            # GPS1B ASCII loader
    gps1b_rnx_loader.py        # GPS1B RINEX loader
    gps1a_loader.py            # GPS1A ASCII loader
    sp3_loader.py              # SP3 precise ephemeris reader
    sp3_ftp.py                 # SP3 FTP downloader
    sp3_probe.py               # SP3 format probe
    ppp.py                     # Single-epoch PPP processor
    fetch_data.py              # GRACE-FO data downloader
    crin2rin.py                # CRINEX → RINEX converter
    third_body.py              # Sun/Moon/planet positions (JPL DE)
    srp.py                     # Solar radiation pressure
    empirical.py               # Empirical acceleration models
    orekit_bridge.py           # Orekit dynamics integration (experimental)
    plotting.py                # Result visualization
```
