# SPDX-License-Identifier: Apache-2.0
"""S-35's M6 arm: score the run against PREDICTION.md's M6 addendum (L16-L28), after the checks it requires.

Written during the arm's run without reading any of its output (the owner's rules, D-203): tested only on the
synthetic fixture of --self-test, and merged before the run's last segment started. The readings below are fixed
here, before any result exists (PREDICTION.md, "Addendum, 2026-10-11").

Inputs: the arm's rows (study_lrpd_serial_m6 --out, S-35's seed) and committed files only:
  information_m6.csv   the credible flags (credible_phi at the cell's true phi)
  summary.csv          S-35's full run, per phi and scenario (compare.py --summary-out): the binomial MLE's PD profile
                       coverage (L18, L19) and its median relative error, the check below; its SHA-256 is checked
                       against MANIFEST.json first
  tests/engine_flags.json   the flag bits

Checks:
  - every panel is fitted with a finite log-likelihood. A panel not fitted fails the checks, as in the validation's
    W1, but every prediction is still reported and the failed cells are named (owner, D-203); its intervals miss,
    and the phi-hat and PD-hat statistics are over the fitted panels. The exit status is then 1;
  - stopping the scoring: a fitted panel's log-likelihood at the phi = 0 scan point is at least the binomial MLE's
    minus 1e-8, and its maximum at least the phi = 0 point's minus 1e-10 (PREDICTION.md, "Checks before scoring");
  - stopping the scoring: the binomial fits are S-35's: per phi and scenario, the median relative error of binom_pd
    (every panel) equals summary.csv's pd_median_rel at its six significant digits (so "binomial" means the same
    fits on the same panels).
Reported, not checks: the placement and doubled-rule flag counts, and the comparator's panel counts.

Readings:
  - pooled coverage or share is the mean of the per-scenario values (R equal), as compare.py pools;
  - a profile interval covers when computed and lo <= truth <= hi; not computed misses; every replicate counts;
  - the comparator's percentile interval covers when it has panels (boot_used > 0) and lo <= PD <= hi;
  - phi-hat is at a bound when |phi-hat| = 0.95 exactly (the search clamps phi to [-0.95, 0.95]);
  - a range is inclusive, compared with 1e-12 of slack for floating point; "under", "below" and "above" are strict;
  - a prediction with a range is held only inside it (L23: 0.55-0.80, not just above 0.5);
  - L21: median relative error is median(PD-hat) / PD - 1 (compare.py); "closer to 0" is |M6| < |binomial|;
    "differs by less than 0.01" is |M6 - binomial| < 0.01, the two relative errors;
  - L25's qualitative parts hold at each of phi = 0, 0.3 and 0.6; its strict part at phi = 0 and 0.6;
  - L24 and L25 get one verdict for the strict part (group D, credible) and one for the qualitative parts.
Deviation from the registration (D-203): the q interval, registered as reported without a prediction, is not
computed by study_lrpd_serial_m6, so it cannot be reported.

    python studies/lrpd-serial-correlation/score_m6.py ROWS.csv [--summary-out FILE]
    python studies/lrpd-serial-correlation/score_m6.py --self-test

Standard library only.
"""

from __future__ import annotations

import argparse
import ast
import contextlib
import csv
import hashlib
import io
import json
import math
import re
import statistics as st
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
INFO = HERE / "information_m6.csv"
SUMMARY = HERE / "summary.csv"
SEED = "0x4D31385245434F56"  # S-35's panels: the recovery seed (tests/recovery/recovery.hpp kSeed)
PHIS = (0.0, 0.3, 0.6)
SCENARIOS = 81
REPLICATES = 1000
GROUP_A = {0, 1, 2, 3, 4, 6, 7, 9, 10, 12, 15, 18, 19, 20, 21, 22, 23, 24, 27, 28, 30, 33, 36, 39, 45, 46, 48, 54,
           57, 60, 72}  # as compare.py (the self-test checks they agree)
SUBSET = (29, 37, 72, 4, 68, 49, 7, 43, 51)
PHI_CREDIBLE = {0.0: (43,), 0.3: (43,), 0.6: (68, 49, 43, 51)}  # L26 as registered
PHI_SHORT = (29, 37, 72)  # L27: the T = 20 subset scenarios
BOOT_REPLICATES, BOOT_PHI, BOOT_PANELS = 100, 0.6, 199
PHI_BOUND = 0.95
EPS = 1e-12
_FLAGS = json.loads((ROOT / "tests" / "engine_flags.json").read_text(encoding="utf-8"))
FLAG_QUAD = _FLAGS["estimate"]["kFlagQuadratureUnconverged"]["value"]  # the doubled-rule check
FLAG_PLACEMENT = _FLAGS["estimate"]["kFlagPlacementCheck"]["value"]
NOT_COMPUTED = _FLAGS["interval"]["kIntervalNotComputed"]["value"]


def group(sid: int, periods: int) -> str:
    return "A" if sid in GROUP_A else {20: "B", 40: "C", 100: "D"}[periods]


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(line for line in fh if not line.startswith("#")))


def inside(x: float, lo: float, hi: float) -> bool:
    return lo - EPS <= x <= hi + EPS


def mean(xs):
    xs = list(xs)
    return sum(xs) / len(xs) if xs else float("nan")


def median(xs):
    return st.median(xs) if xs else float("nan")  # a cell with no fitted panel


def summary_sha256_ok() -> tuple[bool, str]:
    """summary.csv against its MANIFEST.json entry, hashed as committed (LF; a Windows checkout may hold CRLF)."""
    entry = json.loads((HERE / "MANIFEST.json").read_text(encoding="utf-8"))["files"]["summary.csv"]
    m = re.search(r"sha256\s+([0-9a-f]{64})", entry)
    if m is None:
        return False, "no sha256 in MANIFEST.json's entry"
    want = m.group(1)
    got = hashlib.sha256(SUMMARY.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    return got == want, got


def covers(row: dict, name: str, truth: float) -> bool:
    if int(row[f"prof_flags_{name}"]) & NOT_COMPUTED:
        return False
    return float(row[f"prof_lo_{name}"]) <= truth <= float(row[f"prof_hi_{name}"])


class Cell:
    """What one (phi, scenario) contributes, kept small: the rows are streamed, not held."""

    def __init__(self, pd: float, rho: float, periods: int, obligors: int):
        self.pd, self.rho, self.periods, self.obligors = pd, rho, periods, obligors
        self.seen: list[int] = []  # every replicate, fitted or not
        self.reps: list[int] = []
        self.pd_hat: list[float] = []
        self.binom_pd: list[float] = []
        self.phi_hat: list[float] = []
        self.cov = {"pd": [], "rho": [], "phi": []}
        self.phi_computed = 0
        self.boot = []  # (replicate, percentile covers, profile covers) on the comparator's replicates
        self.boot_panels = []  # (used, excluded) where the comparator ran
        self.maxima: list[int] = []
        self.gap: list[float] = []
        self.flag_placement = 0
        self.flag_quad = 0
        self.not_fitted = 0


def load(path: str):
    """Stream the rows: (cells, row count, check failures that stop the scoring, panels not fitted)."""
    cells: dict[tuple[float, int], Cell] = {}
    bad: list[str] = []
    unfit: list[tuple[float, int, int]] = []
    n = 0
    with open(path, newline="") as fh:
        for r in csv.DictReader(line for line in fh if not line.startswith("#")):
            n += 1
            key = (float(r["phi"]), int(r["scenario"]))
            c = cells.get(key)
            if c is None:
                c = cells[key] = Cell(float(r["pd_true"]), float(r["rho_true"]), int(r["periods"]), int(r["obligors"]))
            rep = int(r["replicate"])
            c.seen.append(rep)
            c.binom_pd.append(float(r["binom_pd"]))
            boot_slot = key[0] == BOOT_PHI and key[1] in SUBSET and rep < BOOT_REPLICATES
            where = f"phi {r['phi']} scenario {r['scenario']} replicate {r['replicate']}"
            if r["fitted"] != "1" or not math.isfinite(float(r["loglik"])):
                unfit.append((key[0], key[1], rep))
                c.not_fitted += 1
                for name in ("pd", "rho", "phi"):
                    c.cov[name].append(False)  # every replicate counts: an unfitted panel's intervals miss
                if boot_slot:
                    c.boot.append((rep, False, False))
                continue
            if float(r["phi0_loglik"]) < float(r["binom_loglik"]) - 1e-8:
                bad.append(f"{where}: the phi = 0 point below the binomial MLE")
            elif float(r["loglik"]) < float(r["phi0_loglik"]) - 1e-10:
                bad.append(f"{where}: the maximum below the phi = 0 point")
            c.reps.append(rep)
            c.pd_hat.append(float(r["pd"]))
            c.phi_hat.append(float(r["phi_hat"]))
            c.cov["pd"].append(covers(r, "pd", c.pd))
            c.cov["rho"].append(covers(r, "rho", c.rho))
            phi_done = not int(r["prof_flags_phi"]) & NOT_COMPUTED
            c.phi_computed += phi_done
            c.cov["phi"].append(phi_done and covers(r, "phi", key[0]))
            used, excluded = int(r["boot_used"]), int(r["boot_excluded"])
            if used + excluded > 0:
                c.boot_panels.append((used, excluded))
            if boot_slot:
                pct = used > 0 and float(r["boot_lo_pd"]) <= c.pd <= float(r["boot_hi_pd"])
                c.boot.append((rep, pct, c.cov["pd"][-1]))
            c.maxima.append(int(r["scan_local_maxima"]))
            c.gap.append(float(r["second_maximum_gap"]))
            flags = int(r["flags"])
            c.flag_placement += bool(flags & FLAG_PLACEMENT)
            c.flag_quad += bool(flags & FLAG_QUAD)
    return cells, n, bad, unfit


def verdict(ok: bool, complete: bool) -> str:
    return ("held" if ok else "NOT HELD") if complete else ("dry run: " + ("held" if ok else "not held"))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("rows", nargs="?")
    ap.add_argument("--summary-out")
    ap.add_argument("--dry-run", action="store_true",
                    help="any seed and shape, to check the plumbing; nothing printed is a result")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args(argv)
    if a.self_test:
        return self_test()
    if not a.rows:
        ap.error("ROWS.csv is required")
    with open(a.rows) as fh:
        head = fh.readline()
    if a.dry_run:
        print("DRY RUN: plumbing check only; nothing below is a result\n")
    elif SEED.lower() not in head.lower():
        print(f"{a.rows}: not the arm's rows (no seed {SEED} in its header)", file=sys.stderr)
        return 2

    ok_hash, got = summary_sha256_ok()
    if not ok_hash:
        print(f"summary.csv's SHA-256 {got} is not MANIFEST.json's", file=sys.stderr)
        return 2
    info = {(int(r["scenario"]), float(r["phi"])): r["credible_phi"] == "1" for r in read_csv(INFO)}
    s35 = {(float(r["phi"]), int(r["scenario"])): r for r in read_csv(SUMMARY)}

    cells, n_rows, bad, unfit = load(a.rows)
    scen = sorted({s for (_, s) in cells})
    sizes = {len(c.seen) for c in cells.values()}
    full = (len(scen) == SCENARIOS and set(cells) == {(p, s) for p in PHIS for s in scen}
            and all(sorted(c.seen) == list(range(REPLICATES)) for c in cells.values()))
    print(f"rows: {n_rows}; {len(scen)} scenarios x {len({p for (p, _) in cells})} phi; replicates per cell "
          f"{sorted(sizes)}; {'full matrix' if full else 'NOT THE FULL MATRIX'}")
    if not full and not a.dry_run:
        print("the rows are not 81 scenarios x 3 phi x 1,000 replicates: not scored", file=sys.stderr)
        return 1

    # ---- checks -------------------------------------------------------------------------------------
    mismatch = []
    for (phi, s), c in sorted(cells.items()):
        want = s35.get((phi, s))
        got_rel = f"{st.median(c.binom_pd) / c.pd - 1:.6g}" if c.binom_pd else "none"
        if want is None or got_rel != want["pd_median_rel"]:
            mismatch.append(f"phi {phi} scenario {s}: binomial median relative error {got_rel} against S-35's "
                       f"{want['pd_median_rel'] if want else 'missing'}")
    if a.dry_run:  # smoke rows are other panels: reported only
        print(f"dry run: the binomial fits differ from S-35's in {len(mismatch)} cells")
    else:
        bad += mismatch
    print(f"checks: {len(bad)} failures that stop the scoring (the log-likelihood checks on {n_rows - len(unfit)} "
          f"fitted panels; the binomial fits against S-35's summary in {len(cells)} cells); panels not fitted "
          f"{len(unfit)}")
    for b in bad[:20]:
        print("  " + b)
    if bad:
        print("CHECKS FAILED: not scored", file=sys.stderr)
        return 1
    if unfit:
        failed = sorted({(p, s) for p, s, _ in unfit})
        print(f"CHECKS FAILED: {len(unfit)} panels not fitted, in {len(failed)} cells; every prediction is still "
              "reported (owner, D-203): their intervals miss, phi-hat and PD-hat statistics are over fitted panels")
        for p, s in failed:
            print(f"  failed cell: phi {p} scenario {s}: {cells[(p, s)].not_fitted} of {len(cells[(p, s)].seen)} "
                  f"not fitted (replicates {', '.join(str(r) for q, t, r in unfit if (q, t) == (p, s))})")
    n_place = sum(c.flag_placement for c in cells.values())
    n_quad = sum(c.flag_quad for c in cells.values())
    print(f"reported: placement check flagged {n_place} of {n_rows} fits, doubled rule {n_quad}; either is a result, "
          "not a filter")
    boot_jobs = [bp for c in cells.values() for bp in c.boot_panels]
    print(f"reported: comparator intervals {len(boot_jobs)} (registered: {len(SUBSET) * BOOT_REPLICATES}); panels "
          f"refitted with finite estimates {sum(u for u, _ in boot_jobs)}, not {sum(e for _, e in boot_jobs)} "
          f"(of {BOOT_PANELS} per interval)")

    # ---- per-cell summaries ---------------------------------------------------------------------------
    summ = {}
    for (phi, s), c in cells.items():
        x = {
            "group": group(s, c.periods), "periods": c.periods, "obligors": c.obligors,
            "credible_phi": int(info[(s, phi)]),
            "prof_pd": mean(c.cov["pd"]), "prof_rho": mean(c.cov["rho"]),
            "prof_phi": mean(c.cov["phi"]) if c.phi_computed else float("nan"),
            "phi_neg": mean(v < 0 for v in c.phi_hat), "phi_bound": mean(abs(v) == PHI_BOUND for v in c.phi_hat),
            "phi_mean": mean(c.phi_hat), "phi_median": median(c.phi_hat),
            "pd_median_rel": median(c.pd_hat) / c.pd - 1, "binom_median_rel": st.median(c.binom_pd) / c.pd - 1,
            "binom_prof_pd": float(s35[(phi, s)]["prof_pd"]),
            "boot_pd": mean(p for _, p, _ in c.boot) if c.boot else float("nan"),
            "boot_prof_pd": mean(q for _, _, q in c.boot) if c.boot else float("nan"),
            "local_maxima_mean": mean(c.maxima), "multiple_maxima": mean(m > 1 for m in c.maxima),
            "flag_placement": c.flag_placement, "flag_quad": c.flag_quad, "not_fitted": c.not_fitted,
        }
        summ[(phi, s)] = x

    def pooled(phi, key, sel):
        return mean(summ[(phi, s)][key] for s in scen if sel(s, summ[(phi, s)]))

    def in_group(*gs):
        return lambda s, x: x["group"] in gs

    def cred_d(s, x):
        return x["group"] == "D" and x["credible_phi"] == 1

    rows_out = []
    strict, qual = [], []

    def report(label, prediction, result, ok, kind="strict"):
        v = verdict(ok, full)
        (strict if kind == "strict" else qual).append(ok)
        rows_out.append((label, prediction, result, v))
        print(f"| {label} | {prediction} | {result} | {v}{' (qualitative)' if kind != 'strict' else ''} |")

    print()
    print("| # | Prediction | Result | Verdict |")
    print("|---|---|---|---|")
    # L16
    d = {phi: pooled(phi, "prof_pd", in_group("D")) for phi in (0.3, 0.6)}
    report("L16", "PD profile coverage, group D, in 0.92-0.96 at phi 0.3 and 0.6",
           f"phi 0.3 {d[0.3]:.4f}; phi 0.6 {d[0.6]:.4f}", all(inside(v, 0.92, 0.96) for v in d.values()))
    # L17
    want17 = {(0.6, "B"): (0.82, 0.92), (0.6, "C"): (0.87, 0.95), (0.3, "B"): (0.88, 0.95), (0.3, "C"): (0.90, 0.96)}
    got17 = {k: pooled(k[0], "prof_pd", in_group(k[1])) for k in want17}
    report("L17", "PD profile coverage: phi 0.6 B 0.82-0.92, C 0.87-0.95; phi 0.3 B 0.88-0.95, C 0.90-0.96",
           "; ".join(f"phi {p} {g} {got17[(p, g)]:.4f}" for p, g in want17),
           all(inside(got17[k], *want17[k]) for k in want17))
    # L18, L19
    gain = {(phi, g): pooled(phi, "prof_pd", in_group(g)) - pooled(phi, "binom_prof_pd", in_group(g))
            for phi in (0.0, 0.6) for g in "BCD"}
    report("L18", "phi 0.6: PD profile coverage minus the binomial MLE's >= 0.10 in each of B, C, D",
           "; ".join(f"{g} {gain[(0.6, g)]:+.4f}" for g in "BCD"), all(gain[(0.6, g)] >= 0.10 - EPS for g in "BCD"))
    report("L19", "phi 0: PD profile coverage minus the binomial MLE's in -0.03..+0.01 in each of B, C, D",
           "; ".join(f"{g} {gain[(0.0, g)]:+.4f}" for g in "BCD"),
           all(inside(gain[(0.0, g)], -0.03, 0.01) for g in "BCD"))
    # L20
    r20 = {g: pooled(0.6, "prof_rho", in_group(g)) for g in "BCD"}
    report("L20", "phi 0.6: rho profile coverage in 0.88-0.96 in each of B, C, D",
           "; ".join(f"{g} {r20[g]:.4f}" for g in "BCD"), all(inside(v, 0.88, 0.96) for v in r20.values()))
    # L21
    inf = [s for s in scen if summ[(0.6, s)]["group"] != "A"]
    closer = sum(abs(summ[(0.6, s)]["pd_median_rel"]) < abs(summ[(0.6, s)]["binom_median_rel"]) for s in inf)
    near = sum(abs(summ[(0.6, s)]["pd_median_rel"] - summ[(0.6, s)]["binom_median_rel"]) < 0.01 for s in inf)
    report("L21", "phi 0.6: median relative error of PD-hat closer to 0 than the binomial's in 19-31 of 50; "
           "differs by < 0.01 in >= 40 of 50", f"closer in {closer} of {len(inf)}; within 0.01 in {near}",
           19 <= closer <= 31 and near >= 40)
    # L22
    n22 = sum(cred_d(s, summ[(0.0, s)]) for s in scen)
    s22 = pooled(0.0, "phi_neg", cred_d)
    report("L22", "phi 0: share of phi-hat < 0 over group D's credible scenarios in 0.50-0.60",
           f"{s22:.4f} ({n22} scenarios)", inside(s22, 0.50, 0.60))
    # L23
    a20 = lambda s, x: x["group"] == "A" and x["periods"] == 20
    b23, a23 = pooled(0.0, "phi_neg", in_group("B")), pooled(0.0, "phi_neg", a20)
    report("L23", "phi 0: share of phi-hat < 0 over group B in 0.55-0.80; group A's T = 20 higher than B's",
           f"B {b23:.4f}; A T = 20 {a23:.4f}", inside(b23, 0.55, 0.80) and a23 > b23, "qualitative")
    # L24
    d24, b24 = pooled(0.6, "phi_neg", cred_d), pooled(0.6, "phi_neg", in_group("B"))
    n24 = sum(cred_d(s, summ[(0.6, s)]) for s in scen)
    report("L24", "phi 0.6: share of phi-hat < 0 under 0.5% in group D's credible scenarios",
           f"{d24:.4%} ({n24} scenarios)", d24 < 0.005)
    report("L24 (B)", "phi 0.6: share of phi-hat < 0 in group B in 1%-15%", f"{b24:.4%}", inside(b24, 0.01, 0.15),
           "qualitative")
    # L25
    d25 = {phi: pooled(phi, "phi_bound", cred_d) for phi in (0.0, 0.6)}
    b25 = {phi: pooled(phi, "phi_bound", in_group("B")) for phi in PHIS}
    a25 = {phi: pooled(phi, "phi_bound", a20) for phi in PHIS}
    report("L25", "phi-hat at +/-0.95 under 0.5% in group D's credible scenarios at phi 0 and 0.6",
           "; ".join(f"phi {p} {v:.4%} ({sum(cred_d(s, summ[(p, s)]) for s in scen)} scenarios)" for p, v in d25.items()),
           all(v < 0.005 for v in d25.values()))
    report("L25 (B, A)", "phi-hat at +/-0.95: group B in 1%-20% at each phi; group A's T = 20 at least 10% at each phi",
           "; ".join(f"phi {p}: B {b25[p]:.2%}, A T = 20 {a25[p]:.2%}" for p in PHIS),
           all(inside(b25[p], 0.01, 0.20) and a25[p] >= 0.10 - EPS for p in PHIS), "qualitative")
    # L26
    got26 = {(p, s): summ[(p, s)]["prof_phi"] for p, ss in PHI_CREDIBLE.items() for s in ss}
    listed = {p: sorted(s for s in SUBSET if summ[(p, s)]["credible_phi"]) for p in PHIS}
    agree = all(listed[p] == sorted(PHI_CREDIBLE[p]) for p in PHIS)
    report("L26", "phi interval coverage >= 0.92 where phi is credible (subset)",
           "; ".join(f"{s} at {p} {v:.4f}" for (p, s), v in got26.items())
           + ("" if agree else f" [information_m6.csv's credible subset scenarios: {listed}]"),
           all(v >= 0.92 - EPS for v in got26.values()))
    # L27
    got27 = {s: summ[(0.6, s)]["prof_phi"] for s in PHI_SHORT}
    report("L27", "phi 0.6: phi interval coverage below 0.93 in >= 2 of 29, 37, 72",
           "; ".join(f"{s} {v:.4f}" for s, v in got27.items()), sum(v < 0.93 for v in got27.values()) >= 2,
           "qualitative")
    # L28
    diff28 = {s: summ[(0.6, s)]["boot_pd"] - summ[(0.6, s)]["boot_prof_pd"] for s in SUBSET}
    within = sum(abs(v) <= 0.05 + EPS for v in diff28.values())
    report("L28", "comparator: percentile PD coverage within +/-0.05 of the profile's (replicates 0-99, phi 0.6) "
           "in >= 7 of 9", f"{within} of 9; " + "; ".join(
               f"{s} {summ[(0.6, s)]['boot_pd']:.2f} vs {summ[(0.6, s)]['boot_prof_pd']:.2f}" for s in SUBSET),
           within >= 7)

    print(f"\nstrict: {sum(strict)} of {len(strict)} held; qualitative: {sum(qual)} of {len(qual)} held"
          + (f"; CHECKS FAILED: {len(unfit)} panels not fitted (named above)" if unfit else ""))
    print("\nReported without a prediction, per phi (informative groups B-D unless stated):")
    for phi in PHIS:
        sel = in_group("B", "C", "D")
        print(f"  phi {phi}: phi-hat mean {pooled(phi, 'phi_mean', sel):+.4f}, median of medians "
              f"{st.median(summ[(phi, s)]['phi_median'] for s in scen if sel(s, summ[(phi, s)])):+.4f}; "
              f"share of fits with more than one local maximum along the scan "
              f"{pooled(phi, 'multiple_maxima', lambda s, x: True):.4f} (all scenarios); flagged: placement "
              f"{sum(summ[(phi, s)]['flag_placement'] for s in scen)}, doubled rule "
              f"{sum(summ[(phi, s)]['flag_quad'] for s in scen)}")
    gaps = [g for c in cells.values() for g in c.gap if math.isfinite(g)]
    if gaps:
        print(f"  gap to the second maximum where there is one: {len(gaps)} fits, median {st.median(gaps):.4g}")
    print("  the q interval: NOT REPORTED, a deviation from the registration (D-203): study_lrpd_serial_m6 does not "
          "compute it")
    print("  per-scenario phi-hat mean and median: --summary-out")

    if a.summary_out:
        keys = list(next(iter(summ.values())))
        with open(a.summary_out, "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(["phi", "scenario", *keys])
            for (phi, s), x in sorted(summ.items()):
                w.writerow([phi, s, *[f"{x[k]:.6g}" if isinstance(x[k], float) else x[k] for k in keys]])
    return 1 if unfit else 0


# ---- self-test: a synthetic fixture only, never the arm's rows ---------------------------------------------

HEADER = ("phi,scenario,replicate,pd_true,rho_true,periods,obligors,binom_pd,binom_rho,binom_loglik,binom_flags,"
          "fitted,pd,rho,phi_hat,se_pd,se_rho,se_phi,loglik,phi0_pd,phi0_rho,phi0_loglik,doubled_check,"
          "placement_check,scan_local_maxima,scan_unconverged,second_maximum_gap,evaluations,scores,flags,"
          "prof_lo_pd,prof_hi_pd,prof_flags_pd,prof_lo_rho,prof_hi_rho,prof_flags_rho,prof_lo_phi,prof_hi_phi,"
          "prof_flags_phi,prof_residual,boot_lo_pd,boot_hi_pd,boot_lo_rho,boot_hi_rho,boot_used,boot_excluded")


def fixture(path: Path, replicates: int, seed: str = SEED, tweak=None) -> None:
    """Rows with known answers. Per replicate r: phi-hat -0.1 for even r and +0.1 for odd, -0.95 where r % 100 == 0;
    PD and rho intervals cover unless r % 20 == 0 (0.95), phi's unless r % 25 == 0 (0.96), the comparator's unless
    r % 10 == 0 (0.90); PD-hat's median relative error S-35's binomial one plus 0.002."""
    meta = {int(r["scenario"]): r for r in read_csv(INFO) if float(r["phi"]) == 0.0}
    s35 = {(float(r["phi"]), int(r["scenario"])): r for r in read_csv(SUMMARY)}
    nan = "nan"
    with open(path, "w", newline="") as fh:
        fh.write(f"# S-35 M6 arm (synthetic fixture); seed {seed}; dgp_kernel none\n# profile flags: fixture\n")
        fh.write(HEADER + "\n")
        for phi in PHIS:
            for s in range(SCENARIOS):
                m = meta[s]
                pd, rho = float(m["pd"]), float(m["rho"])
                rel = float(s35[(phi, s)]["pd_median_rel"])
                bpd = pd * (1 + rel)
                for r in range(replicates):
                    ph = -PHI_BOUND if r % 100 == 0 else (-0.1 if r % 2 == 0 else 0.1)
                    cpd, cphi = r % 20 != 0, r % 25 != 0
                    lo_pd, hi_pd = (pd * 0.9, pd * 1.1) if cpd else (pd * 1.1, pd * 1.2)
                    lo_rho, hi_rho = (rho * 0.9, rho * 1.1) if cpd else (rho * 1.1, rho * 1.2)
                    phi_on = s in SUBSET
                    lo_phi, hi_phi = ((phi - 0.1, phi + 0.1) if cphi else (phi + 0.1, phi + 0.2)) if phi_on else (nan, nan)
                    boot = phi == BOOT_PHI and s in SUBSET and r < BOOT_REPLICATES
                    blo, bhi = ((pd * 0.9, pd * 1.1) if r % 10 != 0 else (pd * 1.1, pd * 1.2)) if boot else (nan, nan)
                    v = {
                        "phi": phi, "scenario": s, "replicate": r, "pd_true": pd, "rho_true": rho,
                        "periods": m["periods"], "obligors": m["obligors"], "binom_pd": bpd, "binom_rho": rho,
                        "binom_loglik": -100.0, "binom_flags": 0, "fitted": 1, "pd": pd * (1 + rel + 0.002),
                        "rho": rho, "phi_hat": ph, "se_pd": 0.001, "se_rho": 0.01, "se_phi": 0.1, "loglik": -99.0,
                        "phi0_pd": bpd, "phi0_rho": rho, "phi0_loglik": -100.0, "doubled_check": 0,
                        "placement_check": 0, "scan_local_maxima": 1 + (r % 4 == 0), "scan_unconverged": 0,
                        "second_maximum_gap": 0.5 if r % 4 == 0 else nan, "evaluations": 100, "scores": 10,
                        "flags": FLAG_PLACEMENT if r == 1 else 0,
                        "prof_lo_pd": lo_pd, "prof_hi_pd": hi_pd, "prof_flags_pd": 0,
                        "prof_lo_rho": lo_rho, "prof_hi_rho": hi_rho, "prof_flags_rho": 0,
                        "prof_lo_phi": lo_phi, "prof_hi_phi": hi_phi, "prof_flags_phi": 0 if phi_on else NOT_COMPUTED,
                        "prof_residual": 0, "boot_lo_pd": blo, "boot_hi_pd": bhi, "boot_lo_rho": nan,
                        "boot_hi_rho": nan, "boot_used": BOOT_PANELS if boot else 0, "boot_excluded": 0,
                    }
                    if tweak:
                        tweak(v)
                    fh.write(",".join(str(v[k]) for k in HEADER.split(",")) + "\n")


def self_test() -> int:
    global REPLICATES
    failures = []

    def expect(name, cond):
        if not cond:
            failures.append(name)

    tree = ast.parse((HERE / "compare.py").read_text(encoding="utf-8"))
    theirs = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign)
                  and any(isinstance(t, ast.Name) and t.id == "GROUP_A" for t in n.targets))
    expect("GROUP_A as compare.py's", theirs == GROUP_A)
    expect("summary.csv's hash", summary_sha256_ok()[0])

    def run(argv):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = main(argv)
        return rc, out.getvalue() + err.getvalue()

    saved = REPLICATES
    REPLICATES = 100  # the fixture's matrix; the comparator's replicates 0-99 are all of it
    try:
        with tempfile.TemporaryDirectory() as tmp:
            good = Path(tmp) / "rows.csv"
            fixture(good, REPLICATES)
            summ_path = Path(tmp) / "summary.csv"
            rc, text = run([str(good), "--summary-out", str(summ_path)])
            expect(f"clean fixture scores (rc {rc})", rc == 0)
            expect("checks pass", "checks: 0 failures" in text)
            expect("placement flag counted", f"placement check flagged {len(PHIS) * SCENARIOS} of" in text)
            expect("L16 0.95 held", "| L16 |" in text and "phi 0.3 0.9500; phi 0.6 0.9500 | held |" in text)
            expect("L17 0.95 at B phi 0.6 not held", "phi 0.6 B 0.9500" in text and
                   re.search(r"\| L17 \|.*\| NOT HELD \|", text) is not None)
            expect("L21 all within 0.01", "within 0.01 in 50" in text)
            expect("L22 share 0.5", re.search(r"\| L22 \|.*\| 0\.5000 \(", text) is not None)
            expect("L25 bound share 1%", "phi 0.0 1.0000% (" in text)
            expect("L26 0.96", "43 at 0.6 0.9600" in text)
            expect("L28 0.90 vs 0.95 within 0.05", "9 of 9;" in text and "29 0.90 vs 0.95" in text)
            expect("summary written", summ_path.exists() and len(read_csv(summ_path)) == len(PHIS) * SCENARIOS)

            bad = Path(tmp) / "bad.csv"
            fixture(bad, REPLICATES, tweak=lambda v: v.update(loglik=-101.0) if v["replicate"] == 7 else None)
            rc, text = run([str(bad)])
            expect("a maximum below the phi = 0 point stops", rc == 1 and "the maximum below the phi = 0 point" in text)
            unfit = Path(tmp) / "unfit.csv"
            fixture(unfit, REPLICATES, tweak=lambda v: v.update(fitted=0) if v["replicate"] == 3 else None)
            rc, text = run([str(unfit)])
            expect("a panel not fitted fails the checks", rc == 1 and "CHECKS FAILED: 243 panels not fitted" in text)
            expect("its cells are named", "failed cell: phi 0.0 scenario 0: 1 of 100 not fitted (replicates 3)" in text)
            expect("every prediction still reported", all(f"| L{k} |" in text for k in range(16, 29)))
            expect("an unfitted panel's interval misses", "phi 0.3 0.9400; phi 0.6 0.9400" in text)
            expect("the q interval a deviation", "NOT REPORTED, a deviation from the registration" in text)
            other = Path(tmp) / "other.csv"
            fixture(other, REPLICATES, tweak=lambda v: v.update(binom_pd=v["binom_pd"] * 1.01))
            rc, text = run([str(other)])
            expect("binomial fits not S-35's stop", rc == 1 and "binomial median relative error" in text)
            seed = Path(tmp) / "seed.csv"
            fixture(seed, 2, seed="0x7A11DA7E")
            rc, _ = run([str(seed)])
            expect("another seed refused", rc == 2)
            rc, text = run([str(seed), "--dry-run"])
            expect("dry run on a small fixture", "DRY RUN" in text and "NOT THE FULL MATRIX" in text)
    finally:
        REPLICATES = saved
    for f in failures:
        print("self-test FAILED:", f)
    print(f"self-test: {'ok' if not failures else f'{len(failures)} failures'}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
