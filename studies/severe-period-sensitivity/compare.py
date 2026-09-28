# SPDX-License-Identifier: Apache-2.0
"""S-34: score predictions K1-K10 (PREDICTION.md) against the full-matrix summary.

Inputs: SUMMARY, the summary.csv of study_severe_period over the full matrix; and, from the
repository, the pinned recovery summary (tests/golden/recovery/summary.csv) for the groups of D-136
and the pinned q profile coverage. Prints one Markdown row per prediction, then a consistency check:
the run's original-panel q coverage against the pinned one. Standard library only. Committed before
the full-matrix results existed; where a prediction's wording admits two readings the stricter is
used and stated.

    python studies/severe-period-sensitivity/compare.py SUMMARY
"""

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
VARIANTS = ["1x1-in-100", "2x1-in-100", "1x1-in-1000", "2x1-in-1000"]
# The large-n reference (PREDICTION.md's table, from large_n_reference.py): median shifts in SE units.
REF_Q = {20: [0.81, 1.46, 1.40, 2.45], 40: [0.59, 1.11, 1.02, 1.91], 100: [0.38, 0.73, 0.67, 1.29]}
REF_RHO = {20: [0.68, 1.20, 1.28, 2.24], 40: [0.48, 0.91, 0.92, 1.73], 100: [0.31, 0.60, 0.60, 1.16]}


def read_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def main():
    h = float.fromhex
    rows = read_rows(sys.argv[1])
    sp = {(int(r["scenario"]), r["variant"]): r for r in rows}
    rec = {int(r["scenario"]): r for r in read_rows(ROOT / "tests" / "golden" / "recovery" / "summary.csv")}
    group = {s: "A" if r["pd_verdict"] == "DEFERRED" else {"20": "B", "40": "C", "100": "D"}[r["periods"]]
             for s, r in rec.items()}
    ids = {g: [s for s in sorted(rec) if group[s] == g] for g in "ABCD"}
    bcd = ids["B"] + ids["C"] + ids["D"]
    periods = {s: int(rec[s]["periods"]) for s in rec}
    npd = {s: round(float(rec[s]["pd"]) * int(rec[s]["obligors"]), 6) for s in rec}
    R = int(rows[0]["replicates"])
    out = []

    def row(n, predicted, result, held):
        out.append(f"| {n} | {predicted} | {result} | {'held' if held else 'not held'} |")

    def shift(s, v, a):
        return h(sp[(s, v)][f"shift_{a}_se_hex"])

    def rel(s, v, a):
        return h(sp[(s, v)][f"rel_{a}_hex"])

    # K1: stricter reading: the SE-unit median and the relative median both positive.
    fails = {a: [(s, v) for s in bcd for v in VARIANTS if not (shift(s, v, a) > 0 and rel(s, v, a) > 0)]
             for a in ("pd", "rho", "q")}
    row("K1", "B-D: median shifts of PD, rho, q positive in every scenario-variant (200 each; SE-unit and relative "
        "medians both)", "; ".join(f"{a}: {200 - len(f)} of 200" for a, f in fails.items()),
        all(not f for f in fails.values()))

    # K2: cells (PD, rho, n) with all three T outside group A.
    bad2, cells = [], 0
    for s20 in (s for s in rec if periods[s] == 20):
        cell = [s20, s20 + 3, s20 + 6]
        if any(group[s] == "A" for s in cell):
            continue
        cells += 1
        for v in VARIANTS:
            q = [shift(s, v, "q") for s in cell]
            if not (q[0] > q[1] > q[2]):
                bad2.append(f"{s20}/{v}")
    row("K2", "B-D: q shift (SE) falls with T in every cell outside A, each variant",
        f"{cells} cells; not falling: {', '.join(bad2) or 'none'}", not bad2)

    # K3: two periods against one, same severity.
    ratios = [shift(s, VARIANTS[i + 1], "q") / shift(s, VARIANTS[i], "q") for s in bcd for i in (0, 2)]
    in3 = sum(1.55 <= r <= 2.05 for r in ratios)
    row("K3", "B-D: two/one q shift ratio 1.55-2.05 in at least 90 of 100 scenario-severity pairs",
        f"{in3} of {len(ratios)} (range {min(ratios):.2f}-{max(ratios):.2f})", in3 >= 90)

    # K4 and K6: n*PD >= 100 (n = 10^4, PD >= 1%); K5: n*PD <= 10.
    big = [s for s in bcd if npd[s] >= 100]
    small = [s for s in bcd if npd[s] <= 10]

    def near_reference(name, a, ref):
        ratio = [shift(s, v, a) / ref[periods[s]][i] for s in big for i, v in enumerate(VARIANTS)]
        within = sum(abs(x - 1) <= 0.25 for x in ratio)
        row(name, f"n*PD >= 100 ({len(big)} scenarios): {a} shift within ±25% of the reference in at least 58 of "
            f"{len(ratio)}", f"{within} of {len(ratio)} (ratio to reference {min(ratio):.2f}-{max(ratio):.2f})",
            within >= 58)

    near_reference("K4", "q", REF_Q)
    below = [shift(s, v, "q") < REF_Q[periods[s]][i] for s in small for i, v in enumerate(VARIANTS)]
    row("K5", f"n*PD <= 10 ({len(small)} scenarios): q shift below the reference in at least 69 of {len(below)}",
        f"{sum(below)} of {len(below)}", sum(below) >= 69)
    near_reference("K6", "rho", REF_RHO)

    # K7: exceedance of the original q interval.
    def above(s, v):
        return int(sp[(s, v)]["q_above_old_upper"]) / R

    a1 = [s for s in bcd if above(s, VARIANTS[0]) > 0.05]
    a2 = sum(above(s, VARIANTS[3]) >= 0.40 for s in ids["B"])
    a3 = [f"{s}/{v}" for s in ids["D"] for v in VARIANTS if above(s, v) > 0.10]
    b_range = [above(s, VARIANTS[3]) for s in ids["B"]]
    row("K7", "q-hat above the old upper end: 1x1-in-100 at most 5% in every B-D; 2x1-in-1000 at least 40% in at "
        "least 8 of 12 B; every D variant at most 10%",
        f"1x1-in-100 over 5%: {', '.join(map(str, a1)) or 'none'}; 2x1-in-1000 at 40% or more in B: {a2} of 12 "
        f"(range {min(b_range):.3f}-{max(b_range):.3f}); D over 10%: {', '.join(a3) or 'none'}",
        not a1 and a2 >= 8 and not a3)

    # K8: coverage of the unchanged truth against the pinned q profile coverage.
    def cov(s, v):
        return int(sp[(s, v)]["q_covered"])

    def pinned(s):
        return int(rec[s]["q_profile_covered"])

    c1 = sum(cov(s, VARIANTS[0]) >= pinned(s) for s in bcd)
    c2 = sum(cov(s, VARIANTS[3]) <= pinned(s) - 0.05 * R for s in ids["B"])
    c3 = sum(cov(s, VARIANTS[3]) < pinned(s) for s in ids["D"])
    row("K8", "q coverage: 1x1-in-100 at least pinned in at least 40 of 50; 2x1-in-1000 at least 0.05 below pinned "
        "in at least 8 of 12 B, and below pinned in at least 15 of 21 D",
        f"{c1} of 50; {c2} of 12; {c3} of 21", c1 >= 40 and c2 >= 8 and c3 >= 15)

    # K9 and K10: relative size.
    k9 = [rel(s, VARIANTS[2], "q") for s in ids["B"]]
    row("K9", "T = 20: 1x1-in-1000 raises q-hat by a median of at least 15% in at least 10 of 12 B",
        f"{sum(v >= 0.15 for v in k9)} of 12 (range {min(k9):+.3f} to {max(k9):+.3f})", sum(v >= 0.15 for v in k9) >= 10)
    k10 = {s: [rel(s, v, "q") for v in VARIANTS] for s in (0, 3, 6)}
    row("K10", "Scenarios 0, 3, 6: q-hat's median relative shift at most 0 under all four variants",
        "; ".join(f"{s}: " + ", ".join(f"{x:+.3f}" for x in v) for s, v in k10.items()),
        all(x <= 0 for v in k10.values() for x in v))

    print("| # | Prediction | Result | Held |")
    print("|---|---|---|---|")
    print("\n".join(out))
    print()
    mism = [s for s in rec if (s, "none") in sp and cov(s, "none") != pinned(s)]
    print(f"Consistency: original-panel q coverage equals the pinned coverage in {len(rec) - len(mism)} of {len(rec)} "
          f"scenarios{'; differing: ' + ', '.join(map(str, mism)) if mism else ''}.")


if __name__ == "__main__":
    main()
