# SPDX-License-Identifier: Apache-2.0
"""S-23: score the pre-registered predictions (PREDICTION.md, P1-P11) against the pinned run.

Reads tests/golden/recovery/summary.csv (written by recovery_harness --write) and prints one
Markdown row per prediction: what was predicted, what happened, held or not. Standard library
only. The rules below were committed while the full-matrix run was in progress, before any of
its results were seen. Where PREDICTION.md's wording admits two readings, the stricter one is
used and stated in the row.

Groups as in D-136: A = PD's Wald verdict DEFERRED; otherwise B, C, D by T = 20, 40, 100.

    python studies/derived-quantity-intervals/compare.py
"""

import csv
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
R = 1000
BAND = (0.95 - 3.29 * math.sqrt(0.95 * 0.05 / R), 0.95 + 3.29 * math.sqrt(0.95 * 0.05 / R))


def load():
    with open(ROOT / "tests" / "golden" / "recovery" / "summary.csv", newline="") as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    out = {}
    for r in rows:
        s = {k: r[k] for k in r}
        s["id"] = int(r["scenario"])
        s["T"] = int(r["periods"])
        s["group"] = "A" if r["pd_verdict"] == "DEFERRED" else {20: "B", 40: "C", 100: "D"}[s["T"]]
        s["cov"] = int(r["q_profile_covered"]) / R
        s["below"] = s["cov"] < BAND[0]
        s["above"] = s["cov"] > BAND[1]
        s["pass"] = not (s["below"] or s["above"])
        s["pd_cov"] = int(r["pd_profile_covered"]) / R
        s["rho_cov"] = int(r["rho_profile_covered"]) / R
        s["mre"] = float.fromhex(r["q_median_rel_err_hex"])
        s["logw"] = float.fromhex(r["q_profile_log_width_median_hex"])
        out[s["id"]] = s
    return out


def ids(xs):
    return ", ".join(str(x) for x in sorted(xs)) or "none"


def main():
    sc = load()
    group = {g: [s for s in sc.values() if s["group"] == g] for g in "ABCD"}
    rows = []

    def row(name, predicted, result, held):
        rows.append(f"| {name} | {predicted} | {result} | {'held' if held else 'not held'} |")

    d = group["D"]
    out_d = [s["id"] for s in d if not s["pass"]]
    row("P1", "Group D: all 21 PASS; at most 1 out of the band", f"{len(d) - len(out_d)} of {len(d)} PASS; "
        f"out of the band: {ids(out_d)}", len(d) == 21 and len(out_d) <= 1)

    c = group["C"]
    below_c = {s["id"] for s in c if s["below"]}
    allowed_c = {31, 32, 41, 49, 58, 66, 68, 75}
    pass_c = sum(s["pass"] for s in c)
    row("P2", "Group C: at least 14 of 17 PASS; at most 3 below, all among 31, 32, 41, 49, 58, 66, 68, 75",
        f"{pass_c} of {len(c)} PASS; below: {ids(below_c)}",
        len(c) == 17 and pass_c >= 14 and len(below_c) <= 3 and below_c <= allowed_c)

    b = group["B"]
    below_b = {s["id"] for s in b if s["below"]}
    covs_b = [s["cov"] for s in b]
    row("P3", "Group B: coverage 0.910-0.950 throughout; 2-6 of 12 below; 29, 55, 74 among them",
        f"coverage {min(covs_b):.3f}-{max(covs_b):.3f}; below: {ids(below_b)}",
        len(b) == 12 and all(0.910 <= x <= 0.950 for x in covs_b) and 2 <= len(below_b) <= 6
        and {29, 55, 74} <= below_b)

    a = group["A"]
    below_a = {s["id"] for s in a if s["below"]}
    above_a = [s for s in a if s["above"]]
    allowed_a = {0, 1, 3, 4, 6, 9, 12, 15, 18, 21, 27, 30, 33}
    box_ok = all(int(s["q_profile_box_limited"]) >= 0.25 * R
                 and int(s["q_profile_box_limited_lo"]) > 0.5 * int(s["q_profile_box_limited"]) for s in above_a)
    s72 = sc[72]
    ok72 = s72["below"] or s72["cov"] <= BAND[0] + 0.005
    row("P4", "Group A: at most 3 below, 72 below or within 0.005 of the lower edge; 8-16 above, all among "
        "0, 1, 3, 4, 6, 9, 12, 15, 18, 21, 27, 30, 33, each box-limited in at least 25% of replicates, mostly "
        "at the lower end",
        f"below: {ids(below_a)}; 72 at {s72['cov']:.3f}; above: {ids(s['id'] for s in above_a)}; box-limited "
        f"condition {'met' if box_ok else 'not met'}",
        len(below_a) <= 3 and ok72 and 8 <= len(above_a) <= 16
        and {s["id"] for s in above_a} <= allowed_a and box_ok)

    n_pass = sum(s["pass"] for s in sc.values())
    n_above = sum(s["above"] for s in sc.values())
    n_below = sum(s["below"] for s in sc.values())
    row("P5", "Totals: PASS 55-72; above 8-18; below 3-10", f"PASS {n_pass}; above {n_above}; below {n_below}",
        55 <= n_pass <= 72 and 8 <= n_above <= 18 and 3 <= n_below <= 10)

    below_all = [s for s in sc.values() if s["below"]]
    frac = {s["id"]: int(s["q_profile_miss_below"]) / (R - int(s["q_profile_covered"])) for s in below_all}
    row("P6", "Every below-band scenario: more than 2/3 of non-covering intervals entirely below the true q",
        "; ".join(f"{k}: {v:.2f}" for k, v in sorted(frac.items())) or "no below-band scenario",
        all(v > 2.0 / 3.0 for v in frac.values()))

    inf = b + c + d
    track = [s for s in inf if min(s["pd_cov"], s["rho_cov"]) - 0.015 <= s["cov"] <= max(s["pd_cov"], s["rho_cov"])
             + 0.015]
    row("P7", "Groups B-D: q's coverage within [min(PD, rho) - 0.015, max(PD, rho) + 0.015] in at least 40 of 50",
        f"{len(track)} of {len(inf)}", len(inf) == 50 and len(track) >= 40)

    in_range = [s for s in b if -0.10 <= s["mre"] <= -0.01]
    cells, shrink = 0, 0
    ratios = []
    for s20 in sc.values():
        if s20["T"] != 20:
            continue
        s100 = sc[s20["id"] + 6]  # same (PD, rho, n), T = 100
        if s20["group"] == "A" or s100["group"] == "A":
            continue
        cells += 1
        shrink += abs(s100["mre"]) < abs(s20["mre"])
        ratios.append(s100["logw"] / s20["logw"])
    row("P8", "Group B: median(q-hat/q - 1) between -10% and -1% in at least 9 of 12; its magnitude smaller at "
        "T = 100 than at T = 20 in at least 80% of cells with neither in group A",
        f"{len(in_range)} of 12 in [-10%, -1%] (range {min(s['mre'] for s in b):+.3f} to "
        f"{max(s['mre'] for s in b):+.3f}); smaller at T = 100 in {shrink} of {cells} cells",
        len(in_range) >= 9 and cells > 0 and shrink >= 0.8 * cells)

    widths_b = [math.exp(s["logw"]) for s in b]
    row("P9", "Group B: median width ratio (hi/lo) 1.5-5 in every scenario; median log-ratio at T = 100 0.35-0.60 "
        "of that at T = 20 in every cell where both are informative",
        f"group B ratios {min(widths_b):.2f}-{max(widths_b):.2f}; T = 100 / T = 20 log-ratios "
        f"{min(ratios):.3f}-{max(ratios):.3f} (median {statistics.median(ratios):.3f}) over {len(ratios)} cells",
        all(1.5 <= w <= 5.0 for w in widths_b) and all(0.35 <= x <= 0.60 for x in ratios))

    assessed = [s for s in sc.values() if s["q_wald_verdict"] != "DEFERRED"]
    w_below = []
    for s in assessed:
        m = R - int(s["flagged"])
        lo = 0.95 - 3.29 * math.sqrt(0.95 * 0.05 / m)
        if int(s["q_wald_covered"]) / m < lo:
            w_below.append(s)
    p_below_assessed = sum(s["below"] for s in assessed)
    in_b = sum(s["group"] == "B" for s in w_below)
    miss_lo = sum(int(s["q_wald_miss_below"]) for s in w_below)
    miss_hi = sum(int(s["q_wald_miss_above"]) for s in w_below)
    row("P10", "Wald: among the assessed scenarios, at least as many below the band as the profile, concentrated "
        "in B (more than half), misses mostly below the truth",
        f"{len(assessed)} assessed; Wald below {len(w_below)} ({in_b} in B) vs profile below {p_below_assessed}; "
        f"misses below/above {miss_lo}/{miss_hi}",
        len(w_below) >= p_below_assessed and in_b > len(w_below) / 2 and miss_lo > miss_hi)

    boot = [int(s["q_boot_covered"]) / R for s in sc.values()]
    b_below = sum(x < BAND[0] for x in boot)
    b_above = sum(x > BAND[1] for x in boot)
    b_pass = len(boot) - b_below - b_above
    row("P11", "Bootstrap: below the band in at least 55 of 81; none above; at most 20 PASS",
        f"below {b_below}; above {b_above}; PASS {b_pass}", b_below >= 55 and b_above == 0 and b_pass <= 20)

    surprise = [s["id"] for s in inf if s["cov"] <= min(s["pd_cov"], s["rho_cov"]) - 0.02]
    print("| # | Prediction | Result | Held |\n|---|---|---|---|")
    print("\n".join(rows))
    print(f"\nSurprise check (profile coverage for q 0.02 or more below the lower of PD's and rho's, "
          f"groups B-D): {ids(surprise)}.")
    print(f"Groups: A {len(a)}, B {len(b)}, C {len(c)}, D {len(d)}.")


if __name__ == "__main__":
    main()
