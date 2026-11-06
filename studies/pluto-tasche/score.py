# SPDX-License-Identifier: Apache-2.0
"""S-6: score the run against PREDICTION.md (P1-P6), after the checks it requires.

Written and merged before the run, tested only on the synthetic fixture of --self-test. The readings below are
fixed here, before any result exists.

Inputs: the run's rows (study_pluto_tasche --out, S-6's seed 0x533036504C55544F) and reference.csv (reference.py,
the registered exact numbers). Nothing else: the flag bits it names are tests/engine_flags.json's, which the
self-test checks.

Checks before scoring; any failure stops the scoring (exit 1):
  1. the rows are the registered design: the seed in the header, each scenario's PD, rho, T and n equal to
     reference.csv's, and every (scenario, replicate) with the five gammas;
  2. the tool's bounds are the reference's: per panel, each Pluto-Tasche bound covers exactly when D >= D_min
     (reference.csv's d_min for that scenario, method and gamma). The two compute the bound independently, so this
     checks both;
  3. each bound is a function of D: within a scenario and gamma, panels with the same total have the same bound.

Readings:
  - a cell is (scenario, method, gamma); the methods are independent, oracle and fixed (reference.csv's rho_0.12);
    tier 1 is the 18 scenarios with n PD = 1, the other tier the 9 with n PD = 0.1;
  - coverage is covered panels over R; pooled coverage is summed over the cells named (R equal, so the mean);
  - "shortfall" is gamma minus coverage, "excess" the fixed arm's coverage minus the oracle's, both pooled per
    (rho, gamma) over the tier's 6 scenarios. The registered numbers are the exact pools rounded to 3 decimals
    (the self-test checks they are reference.csv's); within +/-0.02 is |simulated - registered| <= 0.02;
  - a cell's simulated median ratio is the median of bound / PD over its R panels (the mean of the middle two);
  - P5 (iii) and (iv): the stated numbers are exact (reference.csv's, as the self-test checks). Their simulated
    counterpart, the median over the tier's scenarios of the cells' simulated medians, is held if it lies in
    [median of the cells' ratio_lo, median of their ratio_hi]: the median is monotone in each argument, so the cells'
    own bands carry over;
  - P6's laws, ordered: independent <= (b, 0.02) <= (b, 0.12) <= (b, 0.24). A panel carries the independent bound,
    the oracle (its own rho) and the fixed arm (0.12), so per rho the panel's order is: independent <= both;
    at rho 0.02 oracle <= fixed; at 0.24 fixed <= oracle; at 0.12 they are equal. "Follows the order" is >= for the
    higher law's coverage or median ratio;
  - a range is inclusive, with 1e-12 of slack for floating point; "below" is strict;
  - a prediction is held only on the full matrix (27 scenarios x 1,000 replicates x 5 gammas).
Reported without a prediction (PREDICTION.md): coverage at n PD = 0.1 beside the exact values; the profile and
Jeffreys upper limits (coverage, median ratio, and paired against the oracle and the independent bound, with D-176's
z = (b - c) / sqrt(b + c) over the discordant panels); not-computed and truncated profile ends; the 5th and 95th
percentiles of bound / PD per cell (in --summary-out); tier-1 coverage by T.

    python studies/pluto-tasche/score.py ROWS.csv [--summary-out FILE]
    python studies/pluto-tasche/score.py --self-test
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import io
import json
import math
import re
import statistics
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REFERENCE = HERE / "reference.csv"

SEED = "0x533036504C55544F"
GAMMAS = (0.5, 0.75, 0.9, 0.95, 0.99)
METHODS = ("independent", "oracle", "fixed")
REF_METHOD = {"independent": "independent", "oracle": "oracle", "fixed": "rho_0.12"}
RHOS = (0.02, 0.12, 0.24)
REPLICATES = 1000
Z_BAND = 3.29  # D-131
SLACK = 1e-12
SIZE_TOL = 0.02

# Registered numbers (PREDICTION.md), per gamma 0.5 / 0.75 / 0.9 / 0.95 / 0.99.
P2_SHORTFALL = {0.12: (0.016, 0.098, 0.103, 0.100, 0.054), 0.24: (0.053, 0.184, 0.217, 0.217, 0.154)}
P2_BELOW_MIN = 55
P3_EXCESS_002 = (0.059, 0.111, 0.063, 0.036, 0.008)
P3_SHORTFALL_024 = (0.025, 0.089, 0.102, 0.081, 0.040)
P3_EXACT_BELOW_024 = 27
P3_BELOW_MIN_024 = 26
P4_IN_BAND_MIN = 268
P5_IN_BAND_MIN = 403
P5_TIER1_095 = {"oracle": (1.34, 1.52, 1.79), "independent": (1.30, 1.27, 1.23)}
P5_TIER01_095 = {"oracle": (2.32, 2.55, 3.02), "independent": (2.29, 2.29, 2.29)}
P5_TIER01_MAX_099 = (18, 6.81)  # the oracle's largest exact median ratio at n PD = 0.1, gamma 0.99

# Flag bits (tests/engine_flags.json; the self-test checks them).
UPPER_TRUNCATED = 2  # kIntervalUpperTruncated
NOT_COMPUTED = 4  # kIntervalNotComputed
POSTERIOR_REFUSED = 2  # kPosteriorRefused
POSTERIOR_NUMERIC = 8  # kPosteriorNumeric

COLUMNS = ["scenario", "replicate", "pd_true", "rho_true", "periods", "obligors", "gamma", "d_total", "independent",
           "oracle", "fixed", "profile", "jeffreys", "cov_independent", "cov_oracle", "cov_fixed", "cov_profile",
           "cov_jeffreys", "fit_flags", "profile_flags", "jeffreys_flags"]


def num(s: str) -> float:
    return math.nan if "nan" in s.lower() else float(s)


def read_reference(path: Path = REFERENCE) -> dict:
    """{(scenario, method, gamma): row}, methods as the tool names them."""
    back = {v: k for k, v in REF_METHOD.items()}
    out = {}
    with open(path, newline="", encoding="utf-8") as f:
        for raw in csv.DictReader(line for line in f if not line.startswith("#")):
            r = {k: num(v) for k, v in raw.items() if k != "method"}
            out[(int(r["scenario"]), back[raw["method"]], r["gamma"])] = r
    return out


def scenarios(ref: dict) -> dict:
    """{scenario: (pd, rho, T, n)}."""
    return {k[0]: (r["pd"], r["rho"], int(r["periods"]), int(r["obligors"])) for k, r in ref.items()}


def tier1(sc: tuple) -> bool:
    return abs(sc[3] * sc[0] - 1.0) < 1e-9


def gi(g: float) -> int:
    return GAMMAS.index(g)


def median(x) -> float:
    return statistics.median(list(x))


def quantile7(x: list[float], p: float) -> float:
    s = sorted(x)
    h = (len(s) - 1) * p
    lo = math.floor(h)
    return s[lo] + (h - lo) * (s[min(lo + 1, len(s) - 1)] - s[lo])


def in_range(v: float, lo: float, hi: float) -> bool:
    return lo - SLACK <= v <= hi + SLACK


class Run:
    """The rows, indexed: per (scenario, gamma) a list over replicates of the parsed row."""

    def __init__(self, path: Path):
        self.header_lines = []
        self.rows = defaultdict(dict)  # (scenario, gamma) -> {replicate: row}
        with open(path, newline="", encoding="utf-8") as f:
            lines = []
            for line in f:
                if line.startswith("#"):
                    self.header_lines.append(line.rstrip("\n"))
                else:
                    lines.append(line)
        reader = csv.DictReader(lines)
        self.columns = reader.fieldnames
        for r in reader:
            row = {k: num(r[k]) for k in ("pd_true", "rho_true", "gamma", *METHODS, "profile", "jeffreys")}
            for k in ("scenario", "replicate", "periods", "obligors", "d_total", "fit_flags", "profile_flags",
                      "jeffreys_flags", *(f"cov_{m}" for m in (*METHODS, "profile", "jeffreys"))):
                row[k] = int(r[k])
            self.rows[(row["scenario"], row["gamma"])][row["replicate"]] = row

    def cell(self, sid: int, g: float) -> list[dict]:
        d = self.rows.get((sid, g), {})
        return [d[r] for r in sorted(d)]


def checks(run: Run, ref: dict) -> tuple[list[str], bool]:
    """The three checks before scoring: (failures, full matrix)."""
    fail = []
    if run.columns != COLUMNS:
        fail.append(f"check 1: columns {run.columns} are not the tool's")
        return fail, False
    if not any(SEED.lower() in h.lower() for h in run.header_lines):
        fail.append(f"check 1: the header does not name the seed {SEED}")
    sc = scenarios(ref)
    seen = {sid for sid, _ in run.rows}
    for sid in sorted(seen - set(sc)):
        fail.append(f"check 1: scenario {sid} is not one of the registered 27")
    full = seen == set(sc)
    reps = None
    for sid in sorted(seen & set(sc)):
        pd, rho, t, n = sc[sid]
        per_gamma = [set(run.rows.get((sid, g), {})) for g in GAMMAS]
        if any(p != per_gamma[0] for p in per_gamma):
            fail.append(f"check 1: scenario {sid} does not have every gamma for every replicate")
        if per_gamma[0] != set(range(len(per_gamma[0]))):
            fail.append(f"check 1: scenario {sid}'s replicates are not 0..R-1")
        reps = len(per_gamma[0]) if reps is None else reps
        if len(per_gamma[0]) != reps:
            fail.append(f"check 1: scenario {sid} has {len(per_gamma[0])} replicates, others {reps}")
        for g in GAMMAS:
            for row in run.cell(sid, g):
                if (row["pd_true"], row["rho_true"], row["periods"], row["obligors"]) != (pd, rho, t, n):
                    fail.append(f"check 1: scenario {sid} replicate {row['replicate']} is not the registered cell")
                    break
            totals = {row["replicate"]: row["d_total"] for row in run.cell(sid, g)}
            if g != GAMMAS[0] and totals != {r["replicate"]: r["d_total"] for r in run.cell(sid, GAMMAS[0])}:
                fail.append(f"check 1: scenario {sid}'s totals differ between gammas")
    full = full and reps == REPLICATES
    if fail:
        return fail, full
    # Check 2: coverage exactly when D >= D_min. Check 3: the bound is a function of D.
    mism = defaultdict(int)
    nonfunc = 0
    for (sid, g), reps_ in run.rows.items():
        by_d = {}
        for row in reps_.values():
            for m in METHODS:
                want = int(row["d_total"] >= ref[(sid, m, g)]["d_min"])
                if row[f"cov_{m}"] != want or int(row[m] >= row["pd_true"]) != want:
                    mism[(sid, m, g)] += 1
            key = row["d_total"]
            vals = tuple(row[m] for m in METHODS)
            if by_d.setdefault(key, vals) != vals:
                nonfunc += 1
    for (sid, m, g), k in sorted(mism.items()):
        fail.append(f"check 2: scenario {sid} {m} gamma {g}: {k} panels whose coverage is not D >= D_min "
                    f"({int(ref[(sid, m, g)]['d_min'])})")
    if nonfunc:
        fail.append(f"check 3: {nonfunc} panels whose bound differs from another panel's with the same total")
    return fail, full


def coverage(cell: list[dict], key: str) -> float:
    return sum(r[f"cov_{key}"] for r in cell) / len(cell)


def median_ratio(cell: list[dict], key: str) -> float:
    vals = [r[key] / r["pd_true"] for r in cell if not math.isnan(r[key])]
    return median(vals) if vals else math.nan


def score(run: Run, ref: dict, full: bool, summary_out: Path | None) -> list[tuple[str, str, str, bool]]:
    sc = scenarios(ref)
    t1 = [s for s in sorted(sc) if tier1(sc[s])]
    t0 = [s for s in sorted(sc) if not tier1(sc[s])]
    by_rho = {tier: {rho: [s for s in ids if sc[s][1] == rho] for rho in RHOS} for tier, ids in (("1", t1), ("0.1", t0))}
    cov = {}
    med = {}
    for sid in sc:
        for g in GAMMAS:
            cell = run.cell(sid, g)
            for m in (*METHODS, "profile", "jeffreys"):
                cov[(sid, m, g)] = coverage(cell, m)
                med[(sid, m, g)] = median_ratio(cell, m)
    results = []

    # P1
    exact_ok = sum(ref[(s, "oracle", g)]["coverage"] >= g for s in sc for g in GAMMAS)
    low = [(s, g, cov[(s, "oracle", g)]) for s in t1 for g in GAMMAS
           if not cov[(s, "oracle", g)] >= g - Z_BAND * math.sqrt(g * (1 - g) / REPLICATES) - SLACK]
    held = exact_ok == 135 and not low
    detail = f"exact >= gamma in {exact_ok} of 135; simulated above gamma - 3.29 SE in {90 - len(low)} of 90"
    if low:
        detail += "; below: " + ", ".join(f"{s} @ {g} {c:.3f}" for s, g, c in low)
    results.append(("P1", "oracle covers at >= gamma", detail, held))

    # P2
    cells = [(s, g) for rho in (0.12, 0.24) for s in by_rho["1"][rho] for g in GAMMAS]
    exact_below = sum(ref[(s, "independent", g)]["coverage"] < g for s, g in cells)
    sim_below = sum(cov[(s, "independent", g)] < g for s, g in cells)
    parts, size_ok = [], True
    for rho in (0.12, 0.24):
        vals = []
        for k, g in enumerate(GAMMAS):
            short = g - statistics.mean(cov[(s, "independent", g)] for s in by_rho["1"][rho])
            ok = abs(short - P2_SHORTFALL[rho][k]) <= SIZE_TOL + SLACK
            size_ok &= ok
            vals.append(f"{short:.3f}{'' if ok else ' OUT'}")
        parts.append(f"rho {rho}: " + " / ".join(vals))
    held = exact_below == 60 and size_ok and sim_below >= P2_BELOW_MIN
    results.append(("P2", "independent bound undercovers at rho >= 0.12",
                    (f"exact below gamma in {exact_below} of 60; shortfall {'; '.join(parts)}; simulated below gamma "
                     f"in {sim_below} of 60 (>= {P2_BELOW_MIN})"), held))

    # P3
    p3 = []
    ok_all = True
    # (i) rho 0.02
    ids = by_rho["1"][0.02]
    dmin_ok = all(ref[(s, "fixed", g)]["d_min"] <= ref[(s, "oracle", g)]["d_min"] for s in ids for g in GAMMAS)
    paired = sum(1 for s in ids for g in GAMMAS for r in run.cell(s, g) if r["cov_oracle"] and not r["cov_fixed"])
    ge = sum(cov[(s, "fixed", g)] >= cov[(s, "oracle", g)] for s in ids for g in GAMMAS)
    vals, size_ok = [], True
    for k, g in enumerate(GAMMAS):
        ex = statistics.mean(cov[(s, "fixed", g)] - cov[(s, "oracle", g)] for s in ids)
        ok = abs(ex - P3_EXCESS_002[k]) <= SIZE_TOL + SLACK
        size_ok &= ok
        vals.append(f"{ex:.3f}{'' if ok else ' OUT'}")
    ok_i = dmin_ok and paired == 0 and ge == 30 and size_ok
    ok_all &= ok_i
    p3.append(f"(i) D_min <= oracle's {'in all 30' if dmin_ok else 'NOT in all 30'}; panels the oracle covers and "
              f"fixed does not: {paired}; fixed >= oracle in {ge} of 30; excess {' / '.join(vals)}")
    # (ii) rho 0.24
    ids = by_rho["1"][0.24]
    dmin_ok = all(ref[(s, "fixed", g)]["d_min"] >= ref[(s, "oracle", g)]["d_min"] for s in ids for g in GAMMAS)
    le = sum(cov[(s, "fixed", g)] <= cov[(s, "oracle", g)] for s in ids for g in GAMMAS)
    exact_below = sum(ref[(s, "fixed", g)]["coverage"] < g for s in ids for g in GAMMAS)
    sim_below = sum(cov[(s, "fixed", g)] < g for s in ids for g in GAMMAS)
    vals, size_ok = [], True
    for k, g in enumerate(GAMMAS):
        short = g - statistics.mean(cov[(s, "fixed", g)] for s in ids)
        ok = abs(short - P3_SHORTFALL_024[k]) <= SIZE_TOL + SLACK
        size_ok &= ok
        vals.append(f"{short:.3f}{'' if ok else ' OUT'}")
    ok_ii = (dmin_ok and le == 30 and exact_below == P3_EXACT_BELOW_024 and size_ok
             and sim_below >= P3_BELOW_MIN_024)
    ok_all &= ok_ii
    p3.append(f"(ii) D_min >= oracle's {'in all 30' if dmin_ok else 'NOT in all 30'}; fixed <= oracle in {le} of 30; "
              f"exact below gamma in {exact_below} of 30; shortfall {' / '.join(vals)}; simulated below gamma in "
              f"{sim_below} of 30 (>= {P3_BELOW_MIN_024})")
    # (iii) rho 0.12, both tiers
    ids = by_rho["1"][0.12] + by_rho["0.1"][0.12]
    differ = sum(1 for s in ids for g in GAMMAS for r in run.cell(s, g) if r["fixed"] != r["oracle"])
    ok_all &= differ == 0
    p3.append(f"(iii) panels where fixed and oracle differ at rho 0.12: {differ}")
    results.append(("P3", "the rho_a = 0.12 arm's errors", "; ".join(p3), ok_all))

    # P4
    out = [(s, m, g, cov[(s, m, g)]) for s in t1 for m in METHODS for g in GAMMAS
           if not in_range(cov[(s, m, g)], ref[(s, m, g)]["band_lo"], ref[(s, m, g)]["band_hi"])]
    detail = f"in band in {270 - len(out)} of 270 (>= {P4_IN_BAND_MIN})"
    if out:
        detail += "; out: " + ", ".join(f"{s} {m} @ {g} {c:.3f} [{ref[(s, m, g)]['band_lo']:.3f}, "
                                        f"{ref[(s, m, g)]['band_hi']:.3f}]" for s, m, g, c in out)
    results.append(("P4", "simulated coverage matches the exact (tier 1)", detail, 270 - len(out) >= P4_IN_BAND_MIN))

    # P5
    out = [(s, m, g, med[(s, m, g)]) for s in sc for m in METHODS for g in GAMMAS
           if not in_range(med[(s, m, g)], ref[(s, m, g)]["ratio_lo"], ref[(s, m, g)]["ratio_hi"])]
    ok_i = 405 - len(out) >= P5_IN_BAND_MIN
    rising = all(ref[(s, m, GAMMAS[k])]["median_ratio"] < ref[(s, m, GAMMAS[k + 1])]["median_ratio"]
                 for s in sc for m in METHODS for k in range(4))
    tiers_ok = True
    tier_txt = []
    for tier, registered in (("1", P5_TIER1_095), ("0.1", P5_TIER01_095)):
        for m, regs in registered.items():
            vals = []
            for rho, reg in zip(RHOS, regs, strict=True):
                ids = by_rho[tier][rho]
                exact = median(ref[(s, m, 0.95)]["median_ratio"] for s in ids)
                sim = median(med[(s, m, 0.95)] for s in ids)
                lo = median(ref[(s, m, 0.95)]["ratio_lo"] for s in ids)
                hi = median(ref[(s, m, 0.95)]["ratio_hi"] for s in ids)
                ok = round(exact, 2) == reg and in_range(sim, lo, hi)
                tiers_ok &= ok
                vals.append(f"{sim:.2f}{'' if ok else ' OUT'}")
            tier_txt.append(f"n PD {tier} {m} {' / '.join(vals)}")
    mx = max(t0, key=lambda s: ref[(s, "oracle", 0.99)]["median_ratio"])
    max_ok = (mx, round(ref[(mx, "oracle", 0.99)]["median_ratio"], 2)) == P5_TIER01_MAX_099
    detail = (f"(i) in band in {405 - len(out)} of 405 (>= {P5_IN_BAND_MIN}); (ii) exact medians rise with gamma: "
              f"{'yes' if rising else 'NO'}; (iii)/(iv) at gamma 0.95, simulated median over scenarios (rho 0.02 / "
              f"0.12 / 0.24): {'; '.join(tier_txt)}; largest at n PD 0.1, gamma 0.99: scenario {mx}")
    if out:
        detail += "; out: " + ", ".join(f"{s} {m} @ {g} {v:.3f} [{ref[(s, m, g)]['ratio_lo']:.3f}, "
                                        f"{ref[(s, m, g)]['ratio_hi']:.3f}]" for s, m, g, v in out[:20])
    results.append(("P5", "bound/true ratio by gamma", detail, ok_i and rising and tiers_ok and max_ok))

    # P6
    def pairs(rho):  # (lower, upper) per panel
        base = [("independent", "oracle"), ("independent", "fixed")]
        return base + ([("oracle", "fixed")] if rho <= 0.12 else [("fixed", "oracle")])

    broken = defaultdict(int)
    for s in sc:
        for g in GAMMAS:
            for r in run.cell(s, g):
                for lo, hi in pairs(sc[s][1]):
                    if not r[lo] <= r[hi]:
                        broken[(s, lo, hi)] += 1
    cov_bad = [(s, lo, hi, g) for s in t1 for g in GAMMAS for lo, hi in pairs(sc[s][1])
               if not cov[(s, hi, g)] >= cov[(s, lo, g)]]
    med_bad = [(s, lo, hi, g) for s in sc for g in GAMMAS for lo, hi in pairs(sc[s][1])
               if not med[(s, hi, g)] >= med[(s, lo, g)]]
    detail = (f"panels out of order: {sum(broken.values())}; (i) tier-1 coverage out of order in {len(cov_bad)} of "
              f"270 cell pairs; (ii) median ratio out of order in {len(med_bad)} of 405")
    if broken:
        detail += "; " + ", ".join(f"{s} {lo} > {hi}: {k}" for (s, lo, hi), k in sorted(broken.items()))
    results.append(("P6", "ordering of the bounds", detail, not broken and not cov_bad and not med_bad))

    print_results(results, full)
    reported(run, ref, sc, t0, t1, by_rho, cov, med)
    if summary_out:
        write_summary(summary_out, run, ref, sc, cov, med)
    return results


def print_results(results, full: bool) -> None:
    print("## Predictions\n")
    if not full:
        print("NOT THE FULL MATRIX: verdicts are computed but none is held.\n")
    print("| # | Prediction | Result | Verdict |")
    print("|---|---|---|---|")
    for k, what, detail, held in results:
        print(f"| {k} | {what} | {detail} | {'held' if held and full else 'NOT HELD'} |")
    n = sum(h and full for *_, h in results)
    print(f"\n{n} of {len(results)} held.\n")


def z_paired(a: list[int], b: list[int]) -> tuple[int, int, float]:
    """D-176: discordant panels a-only (b10) and b-only (b01), and z = (b10 - b01) / sqrt(b10 + b01)."""
    b10 = sum(1 for x, y in zip(a, b, strict=True) if x and not y)
    b01 = sum(1 for x, y in zip(a, b, strict=True) if y and not x)
    return b10, b01, ((b10 - b01) / math.sqrt(b10 + b01) if b10 + b01 else 0.0)


def reported(run, ref, sc, t0, t1, by_rho, cov, med) -> None:
    print("## Reported, not predicted\n")
    print("### Coverage at n PD = 0.1 (simulated, exact), mean over the 3 scenarios per rho\n")
    print("| rho | gamma | independent | oracle | fixed | cells out of band |")
    print("|---|---|---|---|---|---|")
    for rho in RHOS:
        ids = by_rho["0.1"][rho]
        for g in GAMMAS:
            cells = []
            nout = 0
            for m in METHODS:
                sim = statistics.mean(cov[(s, m, g)] for s in ids)
                ex = statistics.mean(ref[(s, m, g)]["coverage"] for s in ids)
                cells.append(f"{sim:.3f} ({ex:.3f})")
                nout += sum(not in_range(cov[(s, m, g)], ref[(s, m, g)]["band_lo"], ref[(s, m, g)]["band_hi"])
                            for s in ids)
            print(f"| {rho} | {g} | {' | '.join(cells)} | {nout} of 9 |")
    print("\n### Profile and Jeffreys upper limits: pooled coverage and the median over scenarios of the median ratio\n")
    print("| n PD | rho | gamma | profile | Jeffreys | oracle | profile vs oracle (b10, b01, z) | "
          "Jeffreys vs oracle (b10, b01, z) | Jeffreys vs independent (b10, b01, z) |")
    print("|---|---|---|---|---|---|---|---|---|")
    for tier in ("1", "0.1"):
        for rho in RHOS:
            ids = by_rho[tier][rho]
            for g in GAMMAS:
                def pooled(m, g=g, ids=ids):
                    return statistics.mean(cov[(s, m, g)] for s in ids)

                def mr(m, g=g, ids=ids):
                    v = [med[(s, m, g)] for s in ids if not math.isnan(med[(s, m, g)])]
                    return f"{median(v):.2f}" if v else "nan"

                col = {m: [r[f"cov_{m}"] for s in ids for r in run.cell(s, g)]
                       for m in ("profile", "jeffreys", "oracle", "independent")}
                zs = [z_paired(col[a], col[b]) for a, b in (("profile", "oracle"), ("jeffreys", "oracle"),
                                                            ("jeffreys", "independent"))]
                print(f"| {tier} | {rho} | {g} | {pooled('profile'):.3f}, {mr('profile')} | "
                      f"{pooled('jeffreys'):.3f}, {mr('jeffreys')} | {pooled('oracle'):.3f}, {mr('oracle')} | "
                      + " | ".join(f"{b10}, {b01}, {z:+.1f}" for b10, b01, z in zs) + " |")
    print("\n### Ends: profile not computed or truncated above, and Jeffreys refused or numeric (panels, per scenario)\n")
    print("| scenario | n PD | profile not computed | profile upper truncated (gamma 0.99) | Jeffreys none |")
    print("|---|---|---|---|---|")
    for s in sorted(sc):
        c = run.cell(s, 0.99)
        nc = sum(1 for r in c if r["profile_flags"] & NOT_COMPUTED)
        tr = sum(1 for r in c if r["profile_flags"] & UPPER_TRUNCATED)
        jn = sum(1 for r in c if r["jeffreys_flags"] & (POSTERIOR_REFUSED | POSTERIOR_NUMERIC))
        if nc or tr or jn:
            print(f"| {s} | {'1' if s in t1 else '0.1'} | {nc} | {tr} | {jn} |")
    print("\n### Tier-1 coverage by T (pooled over rho and PD, simulated)\n")
    print("| T | gamma | independent | oracle | fixed |")
    print("|---|---|---|---|---|")
    for t in sorted({sc[s][2] for s in t1}):
        ids = [s for s in t1 if sc[s][2] == t]
        for g in GAMMAS:
            print(f"| {t} | {g} | " + " | ".join(f"{statistics.mean(cov[(s, m, g)] for s in ids):.3f}"
                                                for m in METHODS) + " |")
    print()


def write_summary(path: Path, run, ref, sc, cov, med) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["scenario", "pd", "rho", "periods", "obligors", "method", "gamma", "coverage", "exact_coverage",
                     "band_lo", "band_hi", "median_ratio", "exact_median_ratio", "ratio_lo", "ratio_hi", "ratio_p05",
                     "ratio_p95", "missing"])
        for s in sorted(sc):
            pd, rho, t, n = sc[s]
            for m in (*METHODS, "profile", "jeffreys"):
                for g in GAMMAS:
                    c = run.cell(s, g)
                    vals = [r[m] / r["pd_true"] for r in c if not math.isnan(r[m])]
                    rr = ref.get((s, m, g))
                    w.writerow([s, pd, rho, t, n, m, g, f"{cov[(s, m, g)]:.6g}",
                                *(f"{rr[k]:.6g}" if rr else "" for k in ("coverage", "band_lo", "band_hi")),
                                f"{med[(s, m, g)]:.6g}",
                                *(f"{rr[k]:.6g}" if rr else "" for k in ("median_ratio", "ratio_lo", "ratio_hi")),
                                f"{quantile7(vals, 0.05):.6g}" if vals else "nan",
                                f"{quantile7(vals, 0.95):.6g}" if vals else "nan", len(c) - len(vals)])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="S-6: score the run against PREDICTION.md (P1-P6)")
    ap.add_argument("rows", nargs="?")
    ap.add_argument("--summary-out")
    ap.add_argument("--reference", default=str(REFERENCE))
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args(argv)
    if a.self_test:
        return self_test()
    if not a.rows:
        ap.error("ROWS.csv is required")
    ref = read_reference(Path(a.reference))
    run = Run(Path(a.rows))
    fail, full = checks(run, ref)
    print("## Checks before scoring\n")
    for k in (1, 2, 3):
        mine = [x for x in fail if x.startswith(f"check {k}")]
        print(f"- check {k}: {len(mine)} failures" + ("" if not mine else ":\n  - " + "\n  - ".join(mine[:50])))
    print()
    if fail:
        print("Scoring stopped: a check failed.")
        return 1
    score(run, ref, full, Path(a.summary_out) if a.summary_out else None)
    return 0


# --- self-test: a synthetic fixture built from reference.csv ------------------------------------------------------


def fixture(path: Path, ref: dict, tweak=None) -> None:
    """Rows that match the reference: each scenario's totals are the true total's quantiles at (r + 0.5) / R, read
    off the survival function at the knots the reference gives (D_min with its coverage, and D_med at just over 1/2);
    each bound is logit-linear in D through the reference's own points (D_min - 1, D_min, D_lo, D_med, D_hi), then
    put in P6's order panel by panel. tweak(row) may edit a row (a dict of strings) before it is written."""
    sc = scenarios(ref)
    with open(path, "w", newline="", encoding="utf-8") as f:
        f.write(f"# S-6 (D-214, D-215): synthetic fixture; seed {SEED}; dgp_kernel scalar\n")
        w = csv.writer(f, lineterminator="\n")
        w.writerow(COLUMNS)
        for s in sorted(sc):
            pd, rho, t, n = sc[s]
            knots = {0: 1.0}
            for m in METHODS:
                for g in GAMMAS:
                    r = ref[(s, m, g)]
                    knots[int(r["d_min"])] = max(knots.get(int(r["d_min"]), 0.0), r["coverage"])
            dmed = int(ref[(s, "oracle", 0.5)]["d_med"])
            knots[dmed] = max(knots.get(dmed, 0.0), 0.5 + 0.6 / REPLICATES)
            ds = sorted(knots)
            for k in range(len(ds) - 2, -1, -1):  # survival non-increasing in D
                knots[ds[k]] = max(knots[ds[k]], knots[ds[k + 1]])
            totals = [max(d for d in ds if knots[d] >= (r + 0.5) / REPLICATES - 1e-15) for r in range(REPLICATES)]

            def bound(m, g, d, pd=pd, s=s):
                r = ref[(s, m, g)]
                pts = {int(r["d_min"]): r["margin_at_logit"]}
                if r["d_min"] > 0:
                    pts[int(r["d_min"]) - 1] = r["margin_below_logit"]
                lg = math.log(pd / (1 - pd))
                for dk, rk in (("d_lo", "ratio_lo"), ("d_med", "median_ratio"), ("d_hi", "ratio_hi")):
                    pts[int(r[dk])] = math.log(pd * r[rk] / (1 - pd * r[rk])) - lg
                xs = sorted(pts)
                if len(xs) == 1:
                    u = pts[xs[0]] + 0.1 * (d - xs[0])
                else:
                    j = max(0, min(len(xs) - 2, max(i for i in range(len(xs)) if xs[i] <= d) if d >= xs[0] else 0))
                    x0, x1 = xs[j], xs[j + 1]
                    u = pts[x0] + (pts[x1] - pts[x0]) * (d - x0) / (x1 - x0)
                return 1.0 / (1.0 + math.exp(-(lg + u)))

            for r, d in enumerate(totals):
                for g in GAMMAS:
                    b = {m: bound(m, g, d) for m in METHODS}
                    b["oracle"] = max(b["oracle"], b["independent"])
                    b["fixed"] = max(b["fixed"], b["independent"])
                    if rho < 0.12:
                        b["fixed"] = max(b["fixed"], b["oracle"])
                    elif rho > 0.12:
                        b["oracle"] = max(b["oracle"], b["fixed"])
                    else:
                        b["fixed"] = b["oracle"]
                    prof = b["oracle"] * 1.05
                    jef = b["independent"] * 1.1
                    row = {"scenario": s, "replicate": r, "pd_true": repr(pd), "rho_true": repr(rho), "periods": t,
                           "obligors": n, "gamma": repr(g), "d_total": d, **{m: repr(b[m]) for m in METHODS},
                           "profile": repr(prof), "jeffreys": repr(jef),
                           **{f"cov_{m}": int(b[m] >= pd) for m in METHODS}, "cov_profile": int(prof >= pd),
                           "cov_jeffreys": int(jef >= pd), "fit_flags": 0, "profile_flags": 0, "jeffreys_flags": 0}
                    if tweak:
                        tweak(row)
                    w.writerow([row[c] for c in COLUMNS])


def self_test() -> int:
    failures = []

    def expect(name, cond):
        if not cond:
            failures.append(name)

    ref = read_reference()
    sc = scenarios(ref)
    t1 = [s for s in sc if tier1(sc[s])]
    expect("27 scenarios, 18 in tier 1", len(sc) == 27 and len(t1) == 18)
    expect("405 reference rows", len(ref) == 405)
    # The registered numbers are reference.csv's, rounded (PREDICTION.md).
    for rho, vals in P2_SHORTFALL.items():
        for g, v in zip(GAMMAS, vals, strict=True):
            ex = g - statistics.mean(ref[(s, "independent", g)]["coverage"] for s in t1 if sc[s][1] == rho)
            expect(f"P2 shortfall rho {rho} gamma {g} rounds to {v}", round(ex, 3) == v)
    for g, v in zip(GAMMAS, P3_EXCESS_002, strict=True):
        ex = statistics.mean(ref[(s, "fixed", g)]["coverage"] - ref[(s, "oracle", g)]["coverage"]
                             for s in t1 if sc[s][1] == 0.02)
        expect(f"P3 excess gamma {g} rounds to {v}", round(ex, 3) == v)
    for g, v in zip(GAMMAS, P3_SHORTFALL_024, strict=True):
        ex = g - statistics.mean(ref[(s, "fixed", g)]["coverage"] for s in t1 if sc[s][1] == 0.24)
        expect(f"P3 shortfall gamma {g} rounds to {v}", round(ex, 3) == v)
    # The flag bits, as the engine defines them.
    flags = json.loads((ROOT / "tests" / "engine_flags.json").read_text(encoding="utf-8"))["interval"]
    expect("kIntervalNotComputed's value", flags["kIntervalNotComputed"]["value"] == NOT_COMPUTED)
    expect("kIntervalUpperTruncated's value", flags["kIntervalUpperTruncated"]["value"] == UPPER_TRUNCATED)
    posterior = (ROOT / "engine" / "posterior.hpp").read_text(encoding="utf-8")
    expect("kPosteriorRefused is bit 1", re.search(r"kPosteriorRefused = 1u << 1,", posterior) is not None)
    expect("kPosteriorNumeric is bit 3", re.search(r"kPosteriorNumeric = 1u << 3,", posterior) is not None)

    def run(argv):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = main(argv)
        return rc, out.getvalue() + err.getvalue()

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        good = tmp / "rows.csv"
        fixture(good, ref)
        summ = tmp / "summary.csv"
        rc, text = run([str(good), "--summary-out", str(summ)])
        expect(f"clean fixture scores (rc {rc})", rc == 0)
        expect("checks pass", all(f"check {k}: 0 failures" in text for k in (1, 2, 3)))
        expect("full matrix", "NOT THE FULL MATRIX" not in text)
        for k in range(1, 7):
            expect(f"P{k} held on the clean fixture", re.search(rf"\| P{k} \|.*\| held \|", text) is not None)
        expect("summary written", summ.exists() and len(summ.read_text(encoding="utf-8").splitlines()) == 1 + 27 * 5 * 5)
        expect("reported sections", "### Coverage at n PD = 0.1" in text and "### Profile and Jeffreys" in text)

        # Check 2 stops: a coverage flag that is not D >= D_min.
        bad = tmp / "bad.csv"
        fixture(bad, ref, tweak=lambda r: r.update(cov_oracle=1 - r["cov_oracle"])
                if r["scenario"] == 4 and r["replicate"] == 7 and r["gamma"] == repr(0.95) else None)
        rc, text = run([str(bad)])
        expect("check 2 stops", rc == 1 and "check 2: scenario 4 oracle gamma 0.95: 1 panels" in text)

        # Check 1 stops: a scenario's rho that is not the registered one.
        fixture(bad, ref, tweak=lambda r: r.update(rho_true="0.13") if r["scenario"] == 1 and r["replicate"] == 0
                else None)
        rc, text = run([str(bad)])
        expect("check 1 stops", rc == 1 and "check 1: scenario 1 replicate 0 is not the registered cell" in text)

        # Not the full matrix: one scenario's last replicate dropped everywhere is caught by check 1; a short run
        # (every scenario with fewer replicates) is scored but not held.
        short = tmp / "short.csv"
        with open(good, encoding="utf-8") as f, open(short, "w", encoding="utf-8", newline="") as o:
            for line in f:
                if line[0].isdigit() and int(line.split(",")[1]) >= 500:
                    continue
                o.write(line)
        rc, text = run([str(short)])
        expect("a short run is not held", rc == 0 and "NOT THE FULL MATRIX" in text and "0 of 6 held" in text)

        # P6: a total on which the independent bound is above both factor bounds is NOT HELD. The total is at or
        # above every D_min, so the coverage flags (check 2) do not change, and every panel with it moves (check 3).
        d_hi = max(int(ref[(25, m, 0.9)]["d_min"]) for m in METHODS)

        def swap(r):
            if r["scenario"] == 25 and r["gamma"] == repr(0.9) and r["d_total"] == d_hi:
                r["independent"] = repr(max(float(r["oracle"]), float(r["fixed"])) * 1.001)

        p6 = tmp / "p6.csv"
        fixture(p6, ref, tweak=swap)
        rc, text = run([str(p6)])
        expect("P6 broken by a swapped total", rc == 0 and re.search(r"\| P6 \|.*\| NOT HELD \|", text) is not None
               and "25 independent > oracle" in text)

        # P4: tier-1 cells whose band excludes the simulated coverage, through a modified reference. Two out is
        # held (>= 268 of 270), three is not.
        lines = REFERENCE.read_text(encoding="utf-8").splitlines(keepends=True)
        head = [ln for ln in lines if ln.startswith("#")]
        body = list(csv.DictReader(ln for ln in lines if not ln.startswith("#")))
        targets = [(s, "oracle", 0.95) for s in sorted(t1)]
        for nbad, verdict in ((2, "held"), (3, "NOT HELD")):
            moved = set(targets[:nbad])
            refp = tmp / f"ref{nbad}.csv"
            with open(refp, "w", newline="", encoding="utf-8") as o:
                o.writelines(head)
                w = csv.DictWriter(o, fieldnames=list(body[0]), lineterminator="\n")
                w.writeheader()
                for row in body:
                    key = (int(row["scenario"]), {"rho_0.12": "fixed"}.get(row["method"], row["method"]),
                           float(row["gamma"]))
                    if key in moved:
                        row = {**row, "band_lo": "2", "band_hi": "3"}
                    w.writerow(row)
            rc, text = run([str(good), "--reference", str(refp)])
            m = re.search(r"\| P4 \|.*\| (held|NOT HELD) \|", text)
            expect(f"P4 with {nbad} cells out is {verdict}", rc == 0 and m is not None and m.group(1) == verdict)

    for f in failures:
        print("FAIL:", f)
    print(f"self-test: {len(failures)} failures")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
