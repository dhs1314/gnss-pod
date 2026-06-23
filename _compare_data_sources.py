#!/usr/bin/env python3
"""Compare RINEX vs ASCII GPS1B data for the same epoch"""
import sys; sys.path.insert(0, 'src')
import pickle
from gps1b_rnx_loader import load_gps1b_rnx

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
F1_SQ, F2_SQ = F1*F1, F2*F2
ALPHA = F1_SQ / (F1_SQ - F2_SQ)
BETA  = -F2_SQ / (F1_SQ - F2_SQ)
LAM1 = C / F1
LAM2 = C / F2

print(f"ALPHA={ALPHA:.10f}, BETA={BETA:.10f}")
print(f"LAM1={LAM1:.6f}m, LAM2={LAM2:.6f}m")

# Load both
rnx_data = load_gps1b_rnx('data/GPS1B_2024-04-29_C_04.rnx')
ascii_data = pickle.load(open('data/gracefo/2024/2024-04-29/GPS1B_2024-04-29_C_04.pkl', 'rb'))

# Find common epoch
rnx_sods = sorted(rnx_data.keys())
ascii_sods = sorted(ascii_data.keys())

sod = rnx_sods[0]  # 767620800
print(f"\nEpoch sod={sod}")
print(f"RINEX has? {sod in rnx_data}")
print(f"ASCII has? {sod in ascii_data}")

# Compare common SVs
rnx_ep = rnx_data[sod]
ascii_ep = ascii_data[sod]

common_svs = sorted(set(rnx_ep.keys()) & set(ascii_ep.keys()))
print(f"Common SVs: {common_svs}")

for sv in common_svs[:5]:
    r = rnx_ep[sv]
    a = ascii_ep[sv]

    # Recompute P_if from RINEX raw values
    rnx_P_if_computed = ALPHA * r['P1'] + BETA * r['P2']
    rnx_L_if_computed = ALPHA * (r['L1_cyc'] * LAM1) + BETA * (r['L2_cyc'] * LAM2)

    print(f"\n{sv}:")
    print(f"  RINEX: P1={r['P1']:.4f} P2={r['P2']:.4f} C1={r['C1']:.4f}")
    print(f"  RINEX: L1_cyc={r['L1_cyc']:.0f} L2_cyc={r['L2_cyc']:.0f}")
    print(f"  RINEX: L1_m={r['L1_cyc']*LAM1:.4f} L2_m={r['L2_cyc']*LAM2:.4f}")
    print(f"  RINEX: P_if stored={r['P_if']:.4f}  computed={rnx_P_if_computed:.4f}  diff={r['P_if'] - rnx_P_if_computed:.4f}")
    print(f"  RINEX: L_if stored={r['L_if']:.4f}  computed={rnx_L_if_computed:.4f}  diff={r['L_if'] - rnx_L_if_computed:.4f}")
    print(f"  ASCII: P1={a['P1']:.4f} P2={a['P2']:.4f}")
    print(f"  ASCII: L1_phase={a['L1_phase']:.4f} L2_phase={a['L2_phase']:.4f}")
    print(f"  ASCII: P_if={a['P_if']:.4f} L_if={a['L_if']:.4f}")

    # ASCII IF computed from ASCII raw values
    ascii_P_if = ALPHA * a['P1'] + BETA * a['P2']
    ascii_L_if = ALPHA * a['L1'] + BETA * a['L2']
    print(f"  ASCII: P_if computed={ascii_P_if:.4f} (vs stored={a['P_if']:.4f}) diff={ascii_P_if - a['P_if']:.4f}")
    print(f"  ASCII: L_if computed={ascii_L_if:.4f} (vs stored={a['L_if']:.4f}) diff={ascii_L_if - a['L_if']:.4f}")

    # Compare L_if between sources
    print(f"  COMPARE:")
    print(f"    L_if: RINEX={r['L_if']:.4f} ASCII={a['L_if']:.4f} diff={r['L_if']-a['L_if']:.4f}m")
    print(f"    P_if: RINEX={r['P_if']:.4f} ASCII={a['P_if']:.4f} diff={r['P_if']-a['P_if']:.4f}m")
    print(f"    P1:   RINEX={r['P1']:.4f} ASCII={a['P1']:.4f} diff={r['P1']-a['P1']:.4f}m")
    print(f"    P2:   RINEX={r['P2']:.4f} ASCII={a['P2']:.4f} diff={r['P2']-a['P2']:.4f}m")
    print(f"    B_if: RINEX={r['L_if']-r['P_if']:.4f} ASCII={a['L_if']-a['P_if']:.4f}")
