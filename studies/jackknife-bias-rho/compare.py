# SPDX-License-Identifier: Apache-2.0
"""The shared jackknife run: score predictions J1-J23 (PREDICTION.md and its two addenda).

Inputs:
  FULL      summary.csv of study_jackknife_run over the full matrix (J1-J18);
  POLISHED  summary.csv of study_jackknife_run --polished over the study subset (J19-J23), and the
            registered arms on the same subset for comparison;
  and, from the repository: the pinned recovery summary (tests/golden/recovery/summary.csv) for the
  groups of D-136, the parity profile, percentile and S-23 verdicts, and the reviewed bootstrap list
  (tests/recovery/recovery.hpp) for J9's boundary-breakdown verdicts.
Prints one Markdown row per prediction. Standard library only. Committed before the full-matrix
results were seen; where a prediction's wording admits two readings the stricter is used and stated.

    python studies/jackknife-bias-rho/compare.py FULL_SUMMARY POLISHED_SUMMARY
"""

import csv
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
R = 1000
BAND = (0.95 - 3.29 * math.sqrt(0.95 * 0.05 / R), 0.95 + 3.29 * math.sqrt(0.95 * 0.05 / R))
SUBSET_INFORMATIVE = [29, 37, 68, 49, 43, 51]


def read(path):
    with open(path, newline="") as f:
        return {int(r["scenario"]): r for r in csv.DictReader(line for line in f if not line.startswith("#"))}


def cls(c):
    return "below" if c < BAND[0] else "above" if c > BAND[1] else "pass"


def main():
    jk, pol = read(sys.argv[1]), read(sys.argv[2])
    rec = read(ROOT / "tests" / "golden" / "recovery" / "summary.csv")
    group = {s: "A" if r["pd_verdict"] == "DEFERRED" else {"20": "B", "40": "C", "100": "D"}[r["periods"]]
             for s, r in rec.items()}
    ids = {g: [s for s in sorted(rec) if group[s] == g] for g in "ABCD"}
    bcd = ids["B"] + ids["C"] + ids["D"]
    h = float.fromhex
    rows = []

    def row(n, predicted, result, held):
        rows.append(f"| {n} | {predicted} | {result} | {'held' if held else 'not held'} |")

    def cov(r, col):
        return int(r[col]) / R

    # --- S-3 ---
    n1 = sum(abs(h(jk[s]["rho_tilde_bias_hex"])) < abs(h(jk[s]["rho_bias_hex"])) for s in bcd)
    row("J1", "|bias(ρ̃)| < |bias(ρ̂)| in at least 40 of 50 (B-D)", f"{n1} of 50", n1 >= 40)
    rb = sum(h(jk[s]["rho_tilde_rmse_hex"]) > h(jk[s]["rho_rmse_hex"]) for s in ids["B"])
    rd = sum(abs(h(jk[s]["rho_tilde_rmse_hex"]) / h(jk[s]["rho_rmse_hex"]) - 1) <= 0.03 for s in ids["D"])
    row("J2", "RMSE(ρ̃) > RMSE(ρ̂) in at least 9 of 12 (B); within ±3% in at least 15 of 21 (D)",
        f"B: {rb} of 12; D: {rd} of 21", rb >= 9 and rd >= 15)
    diffs = {s: cov(jk[s], "shifted_rho_covered") - cov(jk[s], "parity_rho_covered") for s in bcd}
    worst = max(diffs.values(), key=abs)
    higher = sum(d > 0 for d in diffs.values())
    row("J3", "B-D: shifted coverage within 0.015 of parity in every scenario, higher in at least 30 of 50",
        f"largest |difference| {abs(worst):.3f} ({sum(abs(d) > 0.015 for d in diffs.values())} beyond 0.015); "
        f"higher in {higher}", all(abs(d) <= 0.015 for d in diffs.values()) and higher >= 30)
    moved = [s for s in (29, 55, 74, 72) if jk[s]["shifted_class"] == "PASS"]
    row("J4", "At most 2 of the ρ small-T findings (29, 55, 74, 72) move into the band",
        f"moved: {', '.join(map(str, moved)) or 'none'}", len(moved) <= 2)
    broken = [s for s in bcd if rec[s]["rho_profile_verdict"] == "PASS" and jk[s]["shifted_class"] != "PASS"]
    row("J5", "At most 2 B-D ρ profile PASS verdicts leave the band with the shift",
        f"{len(broken)}: {', '.join(map(str, broken)) or 'none'}", len(broken) <= 2)
    a5 = sum(int(jk[s]["clamped"]) >= 0.05 * R for s in ids["A"])
    cd1 = [s for s in ids["C"] + ids["D"] if int(jk[s]["clamped"]) >= 0.01 * R]
    row("J6", "ρ̃ clamped in at least 5% of replicates in at least 10 of 31 (A); under 1% in every C and D",
        f"A: {a5} of 31; C/D at 1% or more: {', '.join(map(str, cd1)) or 'none'}", a5 >= 10 and not cd1)

    # --- S-5 ---
    up = sum(cov(jk[s], "rho_bca_covered") > cov(jk[s], "rho_pct_covered") for s in bcd)
    bb = sum(jk[s]["rho_bca_class"] == "BELOW" for s in ids["B"])
    dp = sum(jk[s]["rho_bca_class"] == "PASS" for s in ids["D"])
    row("J7", "ρ B-D: BCa > percentile in at least 45 of 50; BCa below in at least 8 of 12 (B); at least 14 of 21 PASS (D)",
        f"{up} of 50; B below {bb}; D PASS {dp}", up >= 45 and bb >= 8 and dp >= 14)
    p8 = sum(jk[s]["pd_bca_class"] == "PASS" for s in bcd)
    row("J8", "PD B-D: BCa in the band in at least 35 of 50", f"{p8} of 50", p8 >= 35)
    src = (ROOT / "tests" / "recovery" / "recovery.hpp").read_text(encoding="utf-8")
    boundary = [(int(s), int(p)) for s, p in re.findall(r"\{(\d+), ([01]), BootstrapDiagnosis::BoundaryBreakdown\}", src)]
    stay = sum(jk[s][("pd" if p == 0 else "rho") + "_bca_class"] == "BELOW" for s, p in boundary)
    nc = sum(max(int(jk[s]["pd_bca_not_computed"]), int(jk[s]["rho_bca_not_computed"])) >= 0.02 * R for s in rec)
    row("J9", f"Of the {len(boundary)} boundary-breakdown verdicts at least 45 stay below; BCa not computed in at least "
        "2% of replicates in at least 10 scenarios", f"{stay} stay below; not computed at 2% or more in {nc} scenarios",
        stay >= 45 and nc >= 10)
    below = sum(jk[s][p + "_bca_class"] == "BELOW" for s in rec for p in ("pd", "rho"))
    above = sum(jk[s][p + "_bca_class"] == "ABOVE" for s in rec for p in ("pd", "rho"))
    row("J10", "Totals over 162: BCa below 60-100; above at most 5", f"below {below}; above {above}",
        60 <= below <= 100 and above <= 5)

    # --- S-21 ---
    ranges = {20: (0.45, 0.85), 40: (0.4, 0.75), 100: (0.3, 0.65)}
    out11 = [s for s in bcd if not ranges[int(rec[s]["periods"])][0] <= h(jk[s]["infl_rho_median_hex"])
             <= ranges[int(rec[s]["periods"])][1]]
    cells_bad = []
    for s20 in (s for s in rec if rec[s]["periods"] == "20"):
        cell = [s20, s20 + 3, s20 + 6]
        if any(group[s] == "A" for s in cell):
            continue
        v = [h(jk[s]["infl_rho_median_hex"]) for s in cell]
        if not (v[0] > v[1] > v[2]):
            cells_bad.append(s20)
    row("J11", "B-D median largest |Δρ| in range by T (0.45-0.85, 0.4-0.75, 0.3-0.65) in every scenario, falling with "
        "T in every cell outside A", f"outside range: {', '.join(map(str, out11)) or 'none'}; cells not falling "
        f"(by T = 20 scenario): {', '.join(map(str, cells_bad)) or 'none'}", not out11 and not cells_bad)
    share = {s: int(jk[s]["extreme_z_top"]) / max(int(jk[s]["extreme_z_assessed"]), 1) for s in bcd}
    low = [s for s, v in share.items() if v < 0.60]
    row("J12", "B-D: the most influential period has the most extreme Z in at least 60% of replicates, every scenario",
        f"range {min(share.values()):.2f}-{max(share.values()):.2f}; below 60%: {len(low)}", not low)
    ratios = {s: h(jk[s]["pair_ratio_median_hex"]) for s in bcd}
    row("J13", "B-D: median leave-two/leave-one ratio 1.3-1.9 in every scenario",
        f"range {min(ratios.values()):.2f}-{max(ratios.values()):.2f}", all(1.3 <= v <= 1.9 for v in ratios.values()))
    gaps = {s: h(jk[s]["refit_gap_max_hex"]) for s in bcd}
    clean = {s: h(jk[s]["refit_gap_clean_max_hex"]) for s in bcd}
    row("J14", "Refined vs exact delete-one estimates under 0.05 SE in every checked replicate outside A",
        f"largest {max(gaps.values()):.3f} SE (scenario {max(gaps, key=lambda s: gaps[s])}); over unflagged delete-one fits "
        f"{max(clean.values()):.3f}; scenarios at 0.05 or more: {sum(v >= 0.05 for v in gaps.values())}",
        all(v < 0.05 for v in gaps.values()))

    # --- J15-J18: the q Wald sub-arms against S-23 ---
    s23_below = [s for s in rec if rec[s]["q_wald_verdict"] == "KNOWN FINDING"
                 and int(rec[s]["q_wald_covered"]) / int(rec[s]["unflagged"]) < 0.95]

    def m(s):
        return int(rec[s]["unflagged"])

    better = sum(int(jk[s]["qa_covered"]) / m(s) > int(rec[s]["q_wald_covered"]) / m(s) for s in s23_below)
    a_in = sum(jk[s]["qa_class"] == "PASS" for s in s23_below)
    row("J15", f"(a) above S-23's Wald in at least 14 of {len(s23_below)}; (a) in the band in 5-12 of them",
        f"above in {better}; in the band in {a_in}", better >= 14 and 5 <= a_in <= 12)
    b_ge = sum(int(jk[s]["qb_covered"]) >= int(jk[s]["qa_covered"]) for s in s23_below)
    b_in = sum(jk[s]["qb_class"] == "PASS" for s in s23_below)
    row("J16", "(b) covers at least as well as (a) in at least 14 of the 17; in the band at least as often",
        f"(b) ≥ (a) in {b_ge}; in the band: (b) {b_in}, (a) {a_in}", b_ge >= 14 and b_in >= a_in)
    lower = sum(h(jk[s]["corr_err_se_jack_hex"]) < h(jk[s]["corr_err_se_delta_hex"]) for s in s23_below)
    row("J17", "error-SE correlation lower for SE_J than SE_Δ in at least 12 of the 17", f"{lower} of 17", lower >= 12)
    s23_pass = [s for s in rec if rec[s]["q_wald_verdict"] == "PASS"]
    leave_a = sum(jk[s]["qa_class"] != "PASS" for s in s23_pass)
    leave_b = sum(jk[s]["qb_class"] != "PASS" for s in s23_pass)
    assessed = [s for s in rec if rec[s]["q_wald_verdict"] != "DEFERRED"]
    b_above = sum(jk[s]["qb_class"] == "ABOVE" for s in assessed)
    row("J18", f"of S-23's {len(s23_pass)} Wald PASS at most 3 leave under (a) and at most 3 under (b); (b) above the "
        f"band in at most 5 of {len(assessed)}", f"(a) {leave_a}; (b) {leave_b}; (b) above {b_above}",
        leave_a <= 3 and leave_b <= 3 and b_above <= 5)

    # --- J19-J23: the polished arm on the subset ---
    def c(s, col):
        return cov(pol[s], col)

    j19 = c(29, "shifted_p_covered") >= 0.91 and h(pol[29]["rho_tilde_p_rmse_hex"]) <= 0.8 * h(pol[29]["rho_tilde_rmse_hex"])
    row("J19", "29: polished shifted coverage at least 0.91; ρ̃_p's RMSE at least 20% below the registered ρ̃'s",
        f"{c(29, 'shifted_p_covered'):.3f}; RMSE {h(pol[29]['rho_tilde_p_rmse_hex']):.4f} vs "
        f"{h(pol[29]['rho_tilde_rmse_hex']):.4f}", j19)
    d20 = {s: c(s, "shifted_p_covered") - c(s, "shifted_rho_covered") for s in (37, 68, 49, 43, 51)}
    row("J20", "37, 68, 49, 43, 51: polished and registered shifted coverages within 0.01",
        "; ".join(f"{s}: {d:+.3f}" for s, d in d20.items()), all(abs(d) <= 0.01 for d in d20.values()))
    d21 = {s: c(s, "shifted_p_covered") - c(s, "parity_rho_covered") for s in (4, 7)}
    row("J21", "4, 7: polished shifted coverage still at least 0.03 below parity",
        "; ".join(f"{s}: {d:+.3f}" for s, d in d21.items()), all(d <= -0.03 for d in d21.values()))
    d22 = {s: c(s, "shifted_p_covered") - c(s, "parity_rho_covered") for s in SUBSET_INFORMATIVE}
    n22 = sum(abs(d) <= 0.015 for d in d22.values())
    row("J22", "At least 5 of the 6 informative: polished shifted coverage within 0.015 of parity",
        f"{n22} of 6 (" + "; ".join(f"{s}: {d:+.3f}" for s, d in d22.items()) + ")", n22 >= 5)

    def qa(s, col):
        return int(pol[s][col]) / int(pol[s]["unflagged"])

    ok23 = all(qa(s, "qa_p_covered") >= qa(s, "qa_covered") - 0.005 for s in SUBSET_INFORMATIVE) and \
        qa(29, "qa_p_covered") > qa(29, "qa_covered")
    row("J23", "q (a): polished at least registered - 0.005 in each of the 6, higher at 29",
        "; ".join(f"{s}: {qa(s, 'qa_p_covered'):.3f} vs {qa(s, 'qa_covered'):.3f}" for s in SUBSET_INFORMATIVE), ok23)

    print("| # | Prediction | Result | Held |\n|---|---|---|---|")
    print("\n".join(rows))


if __name__ == "__main__":
    main()
