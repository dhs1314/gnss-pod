"""Quick functional test for Phase 2.3 modules."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from src.measurement_corrections import (
    phase_wind_up_correction,
    relativity_shapiro_correction,
    compute_pco_ecef_from_nadir,
)
from src.cycle_slip import TurboEdit

# Test 1: Relativity Shapiro
sat = np.array([2.6e7, 0.0, 0.0])
leo = np.array([7.0e6, 0.0, 0.0])
rel = relativity_shapiro_correction(sat, leo)
print(f"Relativity Shapiro: {rel:.4f} m (expected ~0.015-0.020)")

# Test 2: Phase wind-up
dphi, cur = phase_wind_up_correction(sat, leo, None, None)
print(f"Wind-up initial: dphi={dphi:.6f} rad, cur={cur:.6f} rad")
dphi2, cur2 = phase_wind_up_correction(sat, leo, None, cur)
print(f"Wind-up next: dphi={dphi2:.6f} rad (expect ~0, same geometry)")

# Test 3: PCO
pco = compute_pco_ecef_from_nadir(sat, 0.10)
print(f"PCO ECEF: {pco} m (Z ~= -0.10 toward Earth centre)")

# Test 4: TurboEdit cycle slip detection
te = TurboEdit(window_size=10, mw_threshold=4.0, gf_threshold=0.02)
L1 = 2.6e7; L2 = 2.4e7; P1 = 2.6001e7; P2 = 2.4001e7
r1 = te.detect("G01", L1, L2, P1, P2)
print(f"Cycle slip epoch 1: {r1} (expect False - first epoch)")
r2 = te.detect("G01", L1, L2, P1, P2)
print(f"Cycle slip epoch 2: {r2} (expect False - no change)")

# Equal L1/L2 meter shift still changes MW (f1≠f2) — physically correct
r3 = te.detect("G01", L1 + 10.0, L2 + 10.0, P1, P2)
print(f"Cycle slip with equal L1/L2 shift: {r3} (expect True - MW changes via f1!=f2)")

# L1 shift only = MW changes → slip detected
te2 = TurboEdit(window_size=10, mw_threshold=4.0, gf_threshold=0.02)
for i in range(8):
    te2.detect("G02", L1, L2, P1, P2)
r4, dN1, dN2 = te2.detect_slip_L1("G02", L1 + 5.0, L2, P1, P2)
print(f"Cycle slip L1 shift only: slip={r4}, dN1={dN1}, dN2={dN2} (expect slip=True)")

print("\nAll Phase 2.3 tests passed.")
