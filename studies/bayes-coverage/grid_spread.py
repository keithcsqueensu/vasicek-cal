# SPDX-License-Identifier: Apache-2.0
"""S-9: where is the parity grid too coarse for a grid posterior? From pinned results only.

A grid posterior puts its mass on grid points. Where the estimator's spread across replicates is
comparable to the grid spacing, the posterior's quantiles are set by the grid rather than the data.
The spread is the pinned standard deviation of logit(estimate) over the recovery replicates
(tests/golden/recovery/summary.csv, pd_sd_u_hex and rho_sd_u_hex), and the spacing that of the parity
grid's logit axes (61 points over PD [1e-4, 0.2], 41 over rho [1e-3, 0.5]; D-115).

Prints, per parameter, the groups B-D scenarios whose spread is under half a spacing (coarse) and the
number with at least 1.5 spacings (adequate), with the adequate ones' pinned profile verdicts.

    python studies/bayes-coverage/grid_spread.py
"""

import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def logit(v):
    return math.log(v / (1 - v))


def main():
    with open(ROOT / "tests" / "golden" / "recovery" / "summary.csv", newline="") as f:
        rec = {int(r["scenario"]): r for r in csv.DictReader(line for line in f if not line.startswith("#"))}
    spacing = {"pd": (logit(0.2) - logit(1e-4)) / 60, "rho": (logit(0.5) - logit(1e-3)) / 40}
    bcd = [s for s in sorted(rec) if rec[s]["pd_verdict"] != "DEFERRED"]
    for name, sp in spacing.items():
        ratio = {s: float.fromhex(rec[s][f"{name}_sd_u_hex"]) / sp for s in bcd}
        coarse = [s for s in bcd if ratio[s] < 0.5]
        adequate = [s for s in bcd if ratio[s] >= 1.5]
        verdicts = {}
        for s in adequate:
            v = rec[s][f"{name}_profile_verdict"]
            verdicts[v] = verdicts.get(v, 0) + 1
        print(f"{name}: spacing {sp:.4f} in logit; spread/spacing {min(ratio.values()):.2f}-{max(ratio.values()):.2f} "
              f"over B-D")
        print(f"  coarse (under 0.5): {', '.join(map(str, coarse)) or 'none'}")
        print(f"  adequate (1.5 or more): {len(adequate)} scenarios: {', '.join(map(str, adequate))}")
        print(f"  their pinned profile verdicts: {verdicts}")


if __name__ == "__main__":
    main()
