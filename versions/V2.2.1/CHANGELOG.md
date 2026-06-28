# CHANGELOG — v2.2.1 (2026-06-07)

## New Feature: Sequential EKF Reduced-Dynamic POD

新增 `src/sequential_filter.py` 和 `run_sequential_pod.py`，实现逐历元序贯扩展卡尔曼滤波精密定轨。

### Architecture
- **State vector (ECI)**: [r(3), v(3), aR, aT, aN, zwd, clk, amb_1..amb_Nsv]
- **Dynamics**: GGM05C gravity (Nmax=90) + empirical RTN accelerations (Gauss-Markov, τ=600s)
- **Measurements**: GPS1B ionosphere-free code + phase, scalar sequential update (Joseph form)
- **Troposphere**: Saastamoinen ZHD + GMF mapping (`src/troposphere.py`)
- **SV biases**: Zero-meaned per-SV code bias with on-the-fly clock-removed estimation for new SVs

### Key Fixes (during tuning)
1. **Clock process noise** (14x improvement at 2h): Q_clk reduced from `1.0*dt` to `0.001*dt` (σ_rw=0.032 m/√s, stable TCXO)
2. **Phase sigma** (1.8x improvement at 0.17h): Increased from 0.01 to 0.20 to prevent ambiguity lock-in after first corrupted phase update
3. **sv_bias zero-meaning**: 60-epoch window with constellation-mean subtraction; on-the-fly uses state.clk for receiver clock removal

### Accuracy vs Baselines

| Arc   | Seq-EKF V2.2.1 | KF (no-clk_r) | Batch LSQ |
|-------|---------------|---------------|-----------|
| 0.17h | 0.78m         | ~0.64m        | 3.16m     |
| 0.5h  | 1.16m         | 0.64m         | 10-17m    |
| 2h    | 6.63m         | ~9m           | unstable  |

### Known Limitations
- 2h arcs degrade after ~1h when SV count exceeds ~27 (on-the-fly bias quality degrades with accumulated position error)
- Measurement-driven approach (sigma_acc=1e-3) is correct for GGM05C-only dynamics; tighter constraints harm per-epoch empirical estimation

### Files Added
- `src/sequential_filter.py` — SequentialEKF class
- `src/troposphere.py` — Saastamoinen ZHD + GMF troposphere model
- `run_sequential_pod.py` — Main sequential POD script

### Existing Files (carried from v2.1)
- `src/orbit_integrator.py`, `src/orbit_dynamics.py`, `src/coordinates.py`
- `src/gravity_model.py`, `src/empirical.py`, `src/srp.py`, `src/third_body.py`
- `src/sp3_loader.py`, `src/gps1b_rnx_loader.py`, `src/gps1b_loader.py`
- `src/batch_estimator.py`
- `run_batch_pod.py`, `run_kf_ppp.py`
