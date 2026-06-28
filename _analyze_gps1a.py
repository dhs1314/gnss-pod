#!/usr/bin/env python3
"""Analyze GPS1A vs GPS1B RINEX phase data — check integer ambiguity preservation"""
import tarfile

C = 299792458.0
F1, F2 = 1575.42e6, 1227.60e6
LAM1 = C / F1
LAM2 = C / F2

def load_rnx_epochs(tgz, member, max_epochs=3):
    t = tarfile.open(tgz)
    data = t.extractfile(member).read().decode("ascii", errors="replace")
    t.close()

    idx = data.find("END OF HEADER")
    lines = data[idx + 50:].split("\n")

    epochs = []
    cur_epoch = None
    cur_svs = []

    for line in lines:
        line = line.strip()
        if not line:
            continue
        # Epoch header: starts with year
        if len(line) > 30 and line[0:2] == "24":
            if cur_epoch and cur_svs:
                epochs.append((cur_epoch, cur_svs))
                if len(epochs) >= max_epochs:
                    break
            cur_epoch = line
            cur_svs = []
        elif cur_epoch and len(line) > 40:
            # Data line: parse L1, L2, C1, P1, P2
            parts = line.split()
            if len(parts) >= 5:
                cur_svs.append([float(x) for x in parts[:5]])
    if cur_epoch and cur_svs and len(epochs) < max_epochs:
        epochs.append((cur_epoch, cur_svs))
    return epochs

for level in ["1A", "1B"]:
    tgz = f"data/gracefo_{level}_2024-04-29_RL04.ascii.noLRI.tgz"
    member = f"GPS{level}_2024-04-29_C_04.rnx"
    epochs = load_rnx_epochs(tgz, member, max_epochs=2)

    print(f"=== GPS{level} RINEX ===")
    for ep_line, svs in epochs:
        # Parse epoch header
        parts = ep_line.split()
        yr, mo, dy, hr, mn = int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])
        sec = float(parts[5])
        n_sv = int(parts[7]) if len(parts) > 7 else len(svs)
        print(f"  Epoch: {yr}-{mo:02d}-{dy:02d} {hr:02d}:{mn:02d}:{sec:06.3f}, {n_sv} SVs")

        for i, vals in enumerate(svs):
            L1_cyc, L2_cyc, C1_m, P1_m, P2_m = vals
            L1_m = L1_cyc * LAM1
            L2_m = L2_cyc * LAM2

            # The key test: L1 - P1/lam in cycles
            diff_cyc = L1_cyc - P1_m / LAM1
            # Also check IF combination
            alpha = F1**2 / (F1**2 - F2**2)
            beta = -F2**2 / (F1**2 - F2**2)
            L_if_m = alpha * L1_m + beta * L2_m
            P_if_m = alpha * P1_m + beta * P2_m
            B_if = L_if_m - P_if_m

            print(f"    SV{i:2d}: L1={L1_cyc:15.5f}cyc  P1={P1_m:12.3f}m  "
                  f"L1-P1/lam={diff_cyc:10.1f}cyc  B_if={B_if:8.3f}m")
            if i >= 3:
                break
        print()
    print()
