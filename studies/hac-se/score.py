# SPDX-License-Identifier: Apache-2.0
"""S-48: score the run against PREDICTION.md (H1-H6), after the checks it requires.

Written and merged before the run, tested only on the synthetic fixture of --self-test. The readings below are
fixed here, before any result exists.

Inputs: the run's rows (study_hac_se --out, S-35's seed) and committed files only:
  reference.csv                                     the registered reference numbers (reference.py)
  ../lrpd-serial-correlation/summary.csv            S-35's full run, per phi and scenario: the check below; its
                                                    SHA-256 is checked against that study's MANIFEST.json first
  ../lrpd-serial-correlation/summary_m6.csv         the M6 arm's profile PD coverage (H5's comparator); hash checked
  tests/engine_flags.json                           the flag bits

Checks before scoring (PREDICTION.md):
  1. the fits are S-35's: per phi and scenario, the flagged count, and the iid Wald and profile coverage of PD and
     rho, equal summary.csv's to its printed digits (six significant, compare.py's format). If not, the scoring stops;
  2. the HAC SE is missing only where it should be: its count on unflagged fits is reported, per lag, and every such
     fit is named. Reported, not a stop.

Readings:
  - an unflagged fit is one with an iid Wald interval (indicator >= 0). The HAC interval's coverage is over the
    unflagged fits; a fit without a HAC SE misses. Pooled coverage is covered fits over unflagged fits, summed over
    a group's scenarios;
  - "paired" (H3, H4): HAC minus iid Wald, each pooled over the same unflagged fits;
  - the reference (H1, H2, H6) is reference.csv's, pooled per group as registered: weighted by its pinned
    unflagged counts, at full precision (the table's three decimals are its rounding; the self-test checks it);
  - H6's ratio is HAC SE / iid SE in logit PD per fit, over the unflagged fits that have both SEs (finite and
    positive); its root mean square is pooled over the group's fits; within +/-10% is |RMS / reference - 1| <= 0.10;
  - H5's band edge is the registered 0.9273; the AR(1) MLE's pooled profile PD coverage is the mean over the group's
    scenarios of summary_m6.csv's prof_pd at phi = 0.6 (R equal), which the self-test checks against D-208's
    0.910, 0.929 and 0.942;
  - a range is inclusive, with 1e-12 of slack for floating point; "below", "above" and orderings are strict;
  - a prediction is held only on the full matrix (81 scenarios x 1,000 replicates x phi 0, 0.3, 0.6).
Reported without a prediction (PREDICTION.md): the lag sensitivity (coverage and SE ratio at L = 1 and 2 x rule),
group A, the per-scenario verdicts against the recovery band (3.29 standard errors of 0.95 on the scenario's
unflagged count), rho's SE ratio, the side of the misses, the unflagged fits without a HAC SE, and the fixed-b
interval's coverage.

    python studies/hac-se/score.py ROWS.csv [--summary-out FILE]
    python studies/hac-se/score.py --self-test

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
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
S35 = ROOT / "studies" / "lrpd-serial-correlation"
REFERENCE = HERE / "reference.csv"
SEED = "0x4D31385245434F56"  # S-35's panels: the recovery seed (tests/recovery/recovery.hpp kSeed)
PHIS = (0.0, 0.3, 0.6)
SCENARIOS = 81
REPLICATES = 1000
GROUP_A = {0, 1, 2, 3, 4, 6, 7, 9, 10, 12, 15, 18, 19, 20, 21, 22, 23, 24, 27, 28, 30, 33, 36, 39, 45, 46, 48, 54,
           57, 60, 72}  # as compare.py (the self-test checks they agree)
PERIODS = (20, 40, 100)
GROUPS = ("B", "C", "D")
BAND_LO = 0.9273  # H5, as registered
BAND_Z = 3.29  # TOL_RECOVERY_COVERAGE_BAND_Z
M6_REGISTERED = {"B": 0.910, "C": 0.929, "D": 0.942}  # D-208, as PREDICTION.md cites them
LAGS = ("", "_l1", "_l2x")
EPS = 1e-12
NOT_COMPUTED = json.loads((ROOT / "tests" / "engine_flags.json").read_text(encoding="utf-8"))[
    "interval"]["kIntervalNotComputed"]["value"]


def group(sid: int, periods: int) -> str:
    return "A" if sid in GROUP_A else {20: "B", 40: "C", 100: "D"}[periods]


def scenario_periods(sid: int) -> int:
    return PERIODS[(sid // 3) % 3]  # tests/recovery/recovery.hpp scenario()


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(line for line in fh if not line.startswith("#")))


def manifest_sha256_ok(name: str) -> tuple[bool, str]:
    """A committed S-35 file against its MANIFEST.json entry, hashed as committed (LF; a Windows checkout may hold
    CRLF). summary.csv is the main run's file; summary_m6.csv the M6 arm's."""
    man = json.loads((S35 / "MANIFEST.json").read_text(encoding="utf-8"))
    files = man["files"] if name == "summary.csv" else man["m6_arm"]["files"]
    m = re.search(r"sha256\s+([0-9a-f]{64})", str(files.get(name, "")))
    if m is None:
        return False, f"no sha256 in MANIFEST.json's entry for {name}"
    got = hashlib.sha256((S35 / name).read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    return got == m.group(1), got


def fmt6(x: float) -> str:
    return f"{x:.6g}"


def cov(ind):
    """Coverage among indicators: 1 covers; 0 misses; -1 excluded (flagged)."""
    kept = [x for x in ind if x >= 0]
    return sum(kept) / len(kept) if kept else float("nan")


def cov_all(ind):
    return sum(x == 1 for x in ind) / len(ind) if ind else float("nan")


def logit(v: float) -> float:
    return math.log(v) - math.log1p(-v)


class Report:
    def __init__(self):
        self.verdicts = {}

    def row(self, label, prediction, result, ok):
        self.verdicts[label] = ok
        print(f"| {label} | {prediction} | {result} | {ok} |")


def held(ok: bool, complete: bool) -> str:
    return ("held" if ok else "NOT HELD") if complete else ("would hold" if ok else "would not hold")


def reference_pools():
    """reference.csv pooled per group and phi, as registered: weights the pinned unflagged counts."""
    ref = read_csv(REFERENCE)
    out = {}
    for g in "ABCD":
        for phi in PHIS:
            sel = [r for r in ref if r["group"] == g and float(r["phi"]) == phi and int(r["unflagged"]) > 0]
            w = sum(int(r["unflagged"]) for r in sel)
            for name in ("pd", "rho"):
                for tag in LAGS:
                    out[(g, phi, f"{name}_hac{tag}")] = sum(float(r[f"{name}_hac_pred{tag}"]) * int(r["unflagged"])
                                                            for r in sel) / w
                    out[(g, phi, f"{name}_rms{tag}")] = math.sqrt(
                        sum(float(r[f"{name}_se_ratio{tag}"]) ** 2 * int(r["unflagged"]) for r in sel) / w)
    return out


def m6_pools():
    rows = read_csv(S35 / "summary_m6.csv")
    out = {}
    for g in GROUPS:
        vals = [float(r["prof_pd"]) for r in rows if float(r["phi"]) == 0.6 and r["group"] == g]
        out[g] = sum(vals) / len(vals)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="S-48: score the run against PREDICTION.md (H1-H6)")
    ap.add_argument("rows", nargs="?")
    ap.add_argument("--summary-out")
    ap.add_argument("--s35-summary", help=argparse.SUPPRESS)  # the self-test's fixture summary
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args(argv)
    if a.self_test:
        return self_test()
    if not a.rows:
        ap.error("ROWS.csv is required")

    with open(a.rows, newline="") as fh:
        head = fh.readline()
    if f"seed {SEED}" not in head:
        print(f"{a.rows}: not S-48's rows (no seed {SEED} in its header)", file=sys.stderr)
        return 2
    print(head.strip())
    if a.s35_summary is None:
        for name in ("summary.csv", "summary_m6.csv"):
            ok, got = manifest_sha256_ok(name)
            if not ok:
                print(f"{name}: SHA-256 {got} does not match S-35's MANIFEST.json", file=sys.stderr)
                return 2
    s35 = {(float(r["phi"]), int(r["scenario"])): r for r in read_csv(a.s35_summary or S35 / "summary.csv")}

    rows = read_csv(a.rows)
    by = {}
    for r in rows:
        by.setdefault((float(r["phi"]), int(r["scenario"])), []).append(r)
    scen = sorted({s for _, s in by})
    phis = sorted({p for p, _ in by})
    R = min(len(v) for v in by.values())
    complete = (tuple(phis) == PHIS and scen == list(range(SCENARIOS)) and R == REPLICATES
                and all(len(v) == REPLICATES for v in by.values()))
    print(f"{len(rows)} rows: {len(scen)} scenarios x {len(phis)} phi x {R} replicates"
          f"{'' if complete else '; NOT THE FULL MATRIX: verdicts are marked would (not) hold'}")

    # ---- check 1: the fits are S-35's ------------------------------------------------------------------------------
    bad = []
    for (phi, s), rs in sorted(by.items()):
        p = s35.get((phi, s))
        if p is None:
            bad.append(f"phi {phi} scenario {s}: not in S-35's summary")
            continue
        flagged = sum(int(r["wald_pd"]) == -1 for r in rs)
        if flagged != int(p["flagged"]):
            bad.append(f"phi {phi} scenario {s}: flagged {flagged} against S-35's {p['flagged']}")
        for name in ("pd", "rho"):
            got = {f"wald_{name}": fmt6(cov([int(r[f"wald_{name}"]) for r in rs])),
                   f"prof_{name}": fmt6(cov_all([int(r[f"prof_{name}"]) for r in rs]))}
            for k, v in got.items():
                if v != p[k]:
                    bad.append(f"phi {phi} scenario {s}: {k} {v} against S-35's {p[k]}")
    print(f"check 1 (the fits are S-35's): {len(bad)} failures")
    for b in bad[:20]:
        print("  " + b)
    if bad:
        print("CHECKS FAILED: not scored", file=sys.stderr)
        return 1

    # ---- check 2: unflagged fits without a HAC SE ------------------------------------------------------------------
    for tag in LAGS:
        missing = [(phi, s, int(r["replicate"])) for (phi, s), rs in sorted(by.items()) for r in rs
                   if int(r["wald_pd"]) >= 0 and int(r[f"hac_flags{tag}"]) & NOT_COMPUTED]
        lag = {"": "the rule", "_l1": "L = 1", "_l2x": "2 x rule"}[tag]
        print(f"check 2 (unflagged fits without a HAC SE, {lag}): {len(missing)}")
        for phi, s, r in missing:
            print(f"  phi {phi} scenario {s} replicate {r}")
        wrong = [(phi, s) for (phi, s), rs in by.items() for r in rs
                 if int(r["wald_pd"]) == -1 and not int(r[f"hac_flags{tag}"]) & NOT_COMPUTED]
        if wrong:
            print(f"  NOTE: {len(wrong)} flagged fits carry a HAC SE (the engine should skip them)")

    grp = {s: group(s, int(by[(phis[0], s)][0]["periods"])) for s in scen}
    ref = reference_pools()
    m6 = m6_pools()

    def unflagged(rs):
        return [r for r in rs if int(r["wald_pd"]) >= 0]

    def pool(phi, g, key):
        """Covered over unflagged fits, summed over the group's scenarios."""
        fits = [r for s in scen if grp[s] == g for r in unflagged(by.get((phi, s), []))]
        return sum(int(r[key]) == 1 for r in fits) / len(fits) if fits else float("nan")

    def rms(phi, g, name, tag=""):
        ratios = [float(r[f"hac_se_{name}_u{tag}"]) / float(r[f"se_{name}_u"])
                  for s in scen if grp[s] == g for r in unflagged(by.get((phi, s), []))
                  if not int(r[f"hac_flags{tag}"]) & NOT_COMPUTED
                  and all(math.isfinite(v) and v > 0 for v in (float(r[f"hac_se_{name}_u{tag}"]),
                                                              float(r[f"se_{name}_u"])))]
        return (math.sqrt(sum(x * x for x in ratios) / len(ratios)) if ratios else float("nan")), len(ratios)

    print()
    print("| # | Prediction | Result | Held |")
    print("|---|---|---|---|")
    rep = Report()
    # H1, H2
    for label, name in (("H1", "pd"), ("H2", "rho")):
        cells = []
        for g in GROUPS:
            for phi in PHIS:
                got, want = pool(phi, g, f"hac_{name}"), ref[(g, phi, f"{name}_hac")]
                cells.append((g, phi, got, want, abs(got - want) <= 0.04 + EPS))
        n_ok = sum(c[4] for c in cells)
        res = "; ".join(f"{g} {phi} {got:.3f} ({want:.3f}{'' if ok else ', OUT'})" for g, phi, got, want, ok in cells)
        rep.row(label, f"{name.upper()} HAC coverage within 0.04 of the reference in >= 8 of 9",
                f"{n_ok} of 9: {res}", held(n_ok >= 8, complete))
    # H3
    gain = {g: pool(0.6, g, "hac_pd") - pool(0.6, g, "wald_pd") for g in GROUPS}
    ok = gain["D"] > gain["C"] > gain["B"] and gain["D"] >= 0.08 - EPS
    rep.row("H3", "phi 0.6 paired PD gain D > C > B, D >= +0.08",
            "; ".join(f"{g} {gain[g]:+.3f}" for g in GROUPS), held(ok, complete))
    # H4
    short = {g: pool(0.0, g, "wald_pd") - pool(0.0, g, "hac_pd") for g in GROUPS}
    ok = all(short[g] > 0 for g in GROUPS) and short["B"] > short["C"] > short["D"] and \
        0.03 - EPS <= short["B"] <= 0.10 + EPS
    rep.row("H4", "phi 0 PD shortfall (Wald - HAC) > 0 in each, B > C > D, B in 0.03-0.10",
            "; ".join(f"{g} {short[g]:+.3f}" for g in GROUPS), held(ok, complete))
    # H5
    h5 = {g: pool(0.6, g, "hac_pd") for g in GROUPS}
    ok = all(h5[g] < BAND_LO and h5[g] < m6[g] for g in GROUPS)
    rep.row("H5", f"phi 0.6 PD HAC below {BAND_LO} and below the AR(1) MLE's profile in each",
            "; ".join(f"{g} {h5[g]:.3f} (M6 {m6[g]:.3f})" for g in GROUPS), held(ok, complete))
    # H6
    cells = []
    for phi in (0.0, 0.6):
        for g in GROUPS:
            got, n = rms(phi, g, "pd")
            want = ref[(g, phi, "pd_rms")]
            cells.append((phi, g, got, want, n, abs(got / want - 1) <= 0.10 + EPS))
    ok = all(c[5] for c in cells)
    rep.row("H6", "RMS HAC/iid SE ratio (logit PD) within +/-10% of the reference, phi 0 and 0.6",
            "; ".join(f"{phi} {g} {got:.3f} ({want:.3f}, {got / want - 1:+.1%}, n {n})"
                      for phi, g, got, want, n, _ in cells), held(ok, complete))

    print()
    print("Reported without a prediction:")
    print("  Coverage per group and phi (pooled over unflagged fits): iid Wald | HAC at the rule [L = 1, 2 x rule] | "
          "fixed-b | profile (all fits)")
    for name in ("pd", "rho"):
        for g in "ABCD":
            for phi in PHIS:
                fits = [r for s in scen if grp[s] == g for r in by.get((phi, s), [])]
                prof = sum(int(r[f"prof_{name}"]) == 1 for r in fits) / len(fits) if fits else float("nan")
                print(f"    {name.upper()} {g} phi {phi}: {pool(phi, g, f'wald_{name}'):.3f} | "
                      f"{pool(phi, g, f'hac_{name}'):.3f} [{pool(phi, g, f'hac_{name}_l1'):.3f}, "
                      f"{pool(phi, g, f'hac_{name}_l2x'):.3f}] | {pool(phi, g, f'fixed_b_{name}'):.3f} | {prof:.3f}"
                      f"   (reference HAC {ref[(g, phi, f'{name}_hac')]:.3f} [{ref[(g, phi, f'{name}_hac_l1')]:.3f}, "
                      f"{ref[(g, phi, f'{name}_hac_l2x')]:.3f}])")
    print("  RMS HAC/iid SE ratio per group and phi: the rule [L = 1, 2 x rule] (reference)")
    for name in ("pd", "rho"):
        for g in "ABCD":
            for phi in PHIS:
                vals = [rms(phi, g, name, tag)[0] for tag in LAGS]
                refs = [ref[(g, phi, f"{name}_rms{tag}")] for tag in LAGS]
                print(f"    {name.upper()} {g} phi {phi}: {vals[0]:.3f} [{vals[1]:.3f}, {vals[2]:.3f}] "
                      f"({refs[0]:.3f} [{refs[1]:.3f}, {refs[2]:.3f}])")
    print("  Side of the HAC interval's misses (rule), per phi: truth above / below the interval, informative groups")
    for name in ("pd", "rho"):
        for phi in PHIS:
            above = below = 0
            for s in scen:
                if grp[s] == "A":
                    continue
                for r in unflagged(by.get((phi, s), [])):
                    if int(r[f"hac_{name}"]) == 0 and not int(r["hac_flags"]) & NOT_COMPUTED:
                        u_hat, u_true = logit(float(r[name])), logit(float(r[f"{name}_true"]))
                        above += u_true > u_hat
                        below += u_true < u_hat
            print(f"    {name.upper()} phi {phi}: above {above}, below {below}")

    summ = []
    for (phi, s), rs in sorted(by.items()):
        x = {"phi": phi, "scenario": s, "group": grp[s], "periods": rs[0]["periods"], "obligors": rs[0]["obligors"],
             "unflagged": len(unflagged(rs))}
        lo = hi = float("nan")
        if x["unflagged"] > 0:
            half = BAND_Z * math.sqrt(0.95 * 0.05 / x["unflagged"])
            lo, hi = 0.95 - half, 0.95 + half
        for name in ("pd", "rho"):
            x[f"wald_{name}"] = cov([int(r[f"wald_{name}"]) for r in rs])
            for tag in LAGS:
                x[f"hac_{name}{tag}"] = cov([int(r[f"hac_{name}{tag}"]) for r in rs])
            x[f"fixed_b_{name}"] = cov([int(r[f"fixed_b_{name}"]) for r in rs])
            c = x[f"hac_{name}"]
            x[f"hac_{name}_verdict"] = ("none" if math.isnan(c) else "below" if c < lo else "above" if c > hi
                                        else "inside")
        x["no_hac"] = sum(int(r["wald_pd"]) >= 0 and bool(int(r["hac_flags"]) & NOT_COMPUTED) for r in rs)
        summ.append(x)
    print("  HAC interval (rule) against the recovery band, per scenario: below / inside / above (no unflagged fit)")
    for name in ("pd", "rho"):
        for phi in PHIS:
            v = [x[f"hac_{name}_verdict"] for x in summ if x["phi"] == phi]
            print(f"    {name.upper()} phi {phi}: {v.count('below')} / {v.count('inside')} / {v.count('above')} "
                  f"({v.count('none')})")
    if a.summary_out:
        keys = list(summ[0])
        with open(a.summary_out, "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(keys)
            for x in summ:
                w.writerow([fmt6(x[k]) if isinstance(x[k], float) and k != "phi" else x[k] for k in keys])
    return 0


# ---- self-test -------------------------------------------------------------------------------------------------------
# A synthetic fixture: every scenario and phi, R replicates, no flagged fit unless tweaked. Per (group, phi) the HAC
# interval covers a fixed share of fits and the iid Wald another, chosen so each prediction's verdict is known.
FIX_HAC = {("B", 0.0): 0.86, ("B", 0.3): 0.82, ("B", 0.6): 0.71, ("C", 0.0): 0.90, ("C", 0.3): 0.87,
           ("C", 0.6): 0.80, ("D", 0.0): 0.93, ("D", 0.3): 0.91, ("D", 0.6): 0.86}
FIX_WALD_MINUS_HAC = {(0.0, "B"): 0.06, (0.0, "C"): 0.04, (0.0, "D"): 0.02,
                      (0.6, "B"): -0.03, (0.6, "C"): -0.10, (0.6, "D"): -0.14}


def fixture(path: Path, replicates: int, ref, tweak=None) -> None:
    cols = ["phi", "scenario", "replicate", "pd_true", "rho_true", "periods", "obligors", "pd", "rho", "flags",
            "se_pd_u", "se_rho_u", "wald_pd", "wald_rho", "prof_pd", "prof_rho", "prof_flags_pd", "prof_flags_rho"]
    for tag in LAGS:
        cols += [f"lag{tag}", f"hac_flags{tag}", f"hac_se_pd_u{tag}", f"hac_se_rho_u{tag}", f"hac_pd{tag}",
                 f"hac_rho{tag}"]
    cols += ["fixed_b_c", "fixed_b_pd", "fixed_b_rho"]
    with open(path, "w", newline="") as fh:
        fh.write(f"# S-48 (synthetic fixture); seed {SEED}; dgp_kernel none\n")
        w = csv.DictWriter(fh, fieldnames=cols, lineterminator="\n")
        w.writeheader()
        for phi in PHIS:
            for s in range(SCENARIOS):
                t = scenario_periods(s)
                g = group(s, t)
                h = FIX_HAC.get((g, phi), 0.90)
                wd = min(1.0, h + FIX_WALD_MINUS_HAC.get((phi, g), 0.0))
                for r in range(replicates):
                    hc = int(r < round(h * replicates))
                    v = {"phi": phi, "scenario": s, "replicate": r, "pd_true": 0.01, "rho_true": 0.1, "periods": t,
                         "obligors": 1000, "pd": 0.011 if r % 2 else 0.009, "rho": 0.1, "flags": 0,
                         "se_pd_u": 0.2, "se_rho_u": 0.5, "wald_pd": int(r < round(wd * replicates)),
                         "wald_rho": int(r < round(wd * replicates)), "prof_pd": int(r < 95 * replicates // 100),
                         "prof_rho": 1, "prof_flags_pd": 0, "prof_flags_rho": 0, "fixed_b_c": 2.26,
                         "fixed_b_pd": 1, "fixed_b_rho": 1}
                    for tag in LAGS:
                        v.update({f"lag{tag}": 2, f"hac_flags{tag}": 0,
                                  f"hac_se_pd_u{tag}": 0.2 * ref.get((g, phi, f"pd_rms{tag}"), 1.0),
                                  f"hac_se_rho_u{tag}": 0.5 * ref.get((g, phi, f"rho_rms{tag}"), 1.0),
                                  f"hac_pd{tag}": hc, f"hac_rho{tag}": hc})
                    v["hac_rho"] = int(r < round(ref.get((g, phi, "rho_hac"), 0.9) * replicates))
                    if tweak:
                        tweak(v)
                    w.writerow(v)


def fixture_summary(rows_path: Path, out: Path, tweak=None) -> None:
    by = {}
    for r in read_csv(rows_path):
        by.setdefault((float(r["phi"]), int(r["scenario"])), []).append(r)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["phi", "scenario", "flagged", "wald_pd", "prof_pd", "wald_rho", "prof_rho"])
        for (phi, s), rs in sorted(by.items()):
            row = [phi, s, sum(int(r["wald_pd"]) == -1 for r in rs)]
            for name in ("pd", "rho"):
                row += [fmt6(cov([int(r[f"wald_{name}"]) for r in rs])), fmt6(cov_all([int(r[f"prof_{name}"])
                                                                                        for r in rs]))]
            row = [row[0], row[1], row[2], row[3], row[4], row[5], row[6]]
            if tweak:
                tweak(row)
            w.writerow(row)


def self_test() -> int:
    global REPLICATES
    failures = []

    def expect(name, cond):
        if not cond:
            failures.append(name)

    tree = ast.parse((S35 / "compare.py").read_text(encoding="utf-8"))
    theirs = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign)
                  and any(isinstance(t, ast.Name) and t.id == "GROUP_A" for t in n.targets))
    expect("GROUP_A as compare.py's", theirs == GROUP_A)
    expect("summary.csv's hash", manifest_sha256_ok("summary.csv")[0])
    expect("summary_m6.csv's hash", manifest_sha256_ok("summary_m6.csv")[0])

    # The registered numbers are reference.csv's, rounded (PREDICTION.md's table).
    ref = reference_pools()
    table = {"pd_hac": {"B": (0.864, 0.815, 0.707), "C": (0.901, 0.871, 0.799), "D": (0.927, 0.909, 0.862)},
             "rho_hac": {"B": (0.881, 0.867, 0.804), "C": (0.918, 0.910, 0.876), "D": (0.935, 0.931, 0.911)}}
    for key, per in table.items():
        for g, vals in per.items():
            for phi, v in zip(PHIS, vals, strict=True):
                expect(f"reference {key} {g} {phi} rounds to {v}", round(ref[(g, phi, key)], 3) == v)
    for g, vals in {"B": (0.926, 1.170), "C": (0.950, 1.318), "D": (0.975, 1.425)}.items():
        for phi, v in zip((0.0, 0.6), vals, strict=True):
            expect(f"reference RMS {g} {phi} rounds to {v}", round(ref[(g, phi, "pd_rms")], 3) == v)
    m6 = m6_pools()
    for g, v in M6_REGISTERED.items():
        expect(f"M6 {g} rounds to {v}", round(m6[g], 3) == v)
    # The engine's flag bit, and the scenario layout, as the tool writes them.
    expect("kIntervalNotComputed is bit 2", NOT_COMPUTED == 4)
    expect("scenario 29 is T = 20", scenario_periods(29) == 20 and scenario_periods(32) == 40 and scenario_periods(71) == 100)

    def run(argv):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = main(argv)
        return rc, out.getvalue() + err.getvalue()

    saved = REPLICATES
    REPLICATES = 100
    try:
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            good, summ = tmp / "rows.csv", tmp / "s35.csv"
            fixture(good, REPLICATES, ref)
            fixture_summary(good, summ)
            out = tmp / "summary.csv"
            rc, text = run([str(good), "--s35-summary", str(summ), "--summary-out", str(out)])
            expect(f"clean fixture scores (rc {rc})", rc == 0)
            expect("full matrix", "NOT THE FULL MATRIX" not in text)
            expect("check 1 passes", "check 1 (the fits are S-35's): 0 failures" in text)
            expect("check 2 none", "check 2 (unflagged fits without a HAC SE, the rule): 0" in text)
            for k in range(1, 7):
                expect(f"H{k} held", re.search(rf"\| H{k} \|.*\| held \|", text) is not None)
            expect("H3 gains", "B +0.030; C +0.100; D +0.140" in text)
            expect("H4 shortfalls", "B +0.060; C +0.040; D +0.020" in text)
            expect("H6 ratio exact", "+0.0%, n 2" in text)
            expect("summary written", len(read_csv(out)) == len(PHIS) * SCENARIOS)
            expect("band verdicts printed", "PD phi 0.0: " in text and " / " in text)

            # H1: two cells out by more than 0.04 is not held; one is.
            off = tmp / "off.csv"
            fixture(off, REPLICATES, ref, tweak=lambda v: v.update(hac_pd=int(v["replicate"] < 80))
                    if v["phi"] == 0.0 and scenario_periods(v["scenario"]) == 100 and v["scenario"] not in GROUP_A
                    else None)
            fixture_summary(off, summ)
            rc, text = run([str(off), "--s35-summary", str(summ)])
            expect("one cell out still holds H1", re.search(r"\| H1 \| .* \| 8 of 9: .*D 0.0 0.800 \(0.927, OUT\)"
                                                            r".* \| held \|", text) is not None)
            expect("H4 then fails ordering", re.search(r"\| H4 \|.*\| NOT HELD \|", text) is not None)

            # Check 1 stops the scoring: a profile coverage that is not S-35's.
            fixture(good, REPLICATES, ref)
            fixture_summary(good, summ, tweak=lambda row: row.__setitem__(4, "0.9") if row[1] == 5 else None)
            rc, text = run([str(good), "--s35-summary", str(summ)])
            expect("a fit not S-35's stops", rc == 1 and "scenario 5: prof_pd 0.95 against S-35's 0.9" in text)

            # Check 2: an unflagged fit without a HAC SE is named, and misses.
            miss = tmp / "miss.csv"
            fixture(miss, REPLICATES, ref, tweak=lambda v: v.update(hac_flags=4, hac_pd=0, hac_rho=0)
                    if (v["phi"], v["scenario"], v["replicate"]) == (0.3, 40, 0) else None)
            fixture_summary(miss, summ)
            rc, text = run([str(miss), "--s35-summary", str(summ)])
            expect("missing HAC SE named", rc == 0 and "the rule): 1\n  phi 0.3 scenario 40 replicate 0" in text)

            # A flagged fit leaves the Wald and HAC denominators; check 1 compares its count.
            flag = tmp / "flag.csv"
            fixture(flag, REPLICATES, ref, tweak=lambda v: v.update(
                wald_pd=-1, wald_rho=-1, hac_pd=-1, hac_rho=-1, hac_pd_l1=-1, hac_rho_l1=-1, hac_pd_l2x=-1,
                hac_rho_l2x=-1, fixed_b_pd=-1, fixed_b_rho=-1, hac_flags=4, hac_flags_l1=4, hac_flags_l2x=4,
                flags=32) if v["replicate"] == 99 else None)
            fixture_summary(flag, summ)
            rc, text = run([str(flag), "--s35-summary", str(summ)])
            expect("flagged fits scored", rc == 0 and "the rule): 0" in text)
            fixture_summary(good, summ)
            rc, text = run([str(flag), "--s35-summary", str(summ)])
            expect("flagged count checked", rc == 1 and "flagged 1 against S-35's 0" in text)

            # Another seed is refused; a small matrix is marked.
            other = tmp / "other.csv"
            other.write_text(good.read_text().replace(SEED, "0x7A11DA7E"))
            rc, _ = run([str(other), "--s35-summary", str(summ)])
            expect("another seed refused", rc == 2)
            REPLICATES = 1000
            fixture_summary(good, summ)
            rc, text = run([str(good), "--s35-summary", str(summ)])
            expect("small matrix marked", "NOT THE FULL MATRIX" in text and "would hold" in text
                   and "| held |" not in text)
    finally:
        REPLICATES = saved
    for f in failures:
        print("self-test FAILED:", f)
    print(f"self-test: {'ok' if not failures else f'{len(failures)} failures'}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
