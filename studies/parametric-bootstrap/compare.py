# SPDX-License-Identifier: Apache-2.0
"""S-10 (K1-K5) and S-4 (H1-H6): score the registered predictions against the run.

Reads fits.parquet (export_fits.py, from study_parametric_bootstrap --out), S-9's rows
(studies/bayes-coverage/fits.parquet: the Jeffreys interval for the replicates the run did not
recompute) and the pinned recovery summary (D-136 groups; the iid percentile coverage for K4).
Committed before the subset run, so before any result of the run existed.

    uv run --no-project --with pyarrow==25.0.1 python studies/parametric-bootstrap/compare.py FITS_PARQUET \
        [--summary-out SUMMARY_CSV]

Checks first (the registrations require them; a failure stops the scoring):
  - Jeffreys: every recomputed replicate (0-49) equals S-9's row bit for bit (ends and flags);
  - profile: every replicate's recomputed ends and flags equal S-9's rows (the same panels and engine);
  - W0: the profile interval covers exactly when W0 <= c, replicate by replicate; mismatches are
    counted and reported (the two are solved independently to the solvers' tolerances).

Coverage is over all replicates; an interval not computed does not cover (D-131). The band is
0.95 +/- 3.29 sqrt(0.95 x 0.05 / R). S-4a and S-4b cover when W0 <= c k (c = 3.8415), k_a the
scenario's mean finite W0 (in-sample), k_b the replicate's own. Widths are median logit widths over
the replicates with an interval; a method's width ratio to the profile's is per scenario, the median
of the method's widths over the median of the profile's on the same replicates (for S-4, replicates
0-199, where the corrected intervals are solved).

Readings fixed before the subset run (the addenda of 2026-09-29 in both PREDICTION.md files):
  - K3 and H4 pool groups C and D together (38 scenarios), per method and parameter.
  - "Pooled" group coverage is over all replicates of the group's scenarios.
  - On a subset of scenarios every figure is computed on the scenarios present and printed with its
    denominators; "held" is marked "subset" rather than yes/no when a prediction's registered
    denominator is not complete.
"""

import argparse
import csv
import math
import sys
from pathlib import Path

import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[2]
NOT_COMPUTED = 1 << 2  # engine::kIntervalNotComputed
REFUSED, NUMERIC = 1 << 1, 1 << 3  # engine::kPosterior*
C = 3.841458820694124  # chi2_1 0.95 quantile = 2 x kProfileThreshold95
FINDINGS = [(29, "pd"), (29, "rho"), (55, "rho"), (68, "pd"), (72, "rho"), (74, "rho")]
PARAMS = ("pd", "rho")


def logit(v):
    return math.log(v) - math.log1p(-v)


def median(x):
    x = sorted(x)
    n = len(x)
    return float("nan") if n == 0 else x[n // 2] if n % 2 else 0.5 * (x[n // 2 - 1] + x[n // 2])


def summary_rows():
    with open(REPO / "tests/golden/recovery/summary.csv", newline="") as f:
        return {int(r["scenario"]): r for r in csv.DictReader(line for line in f if not line.startswith("#"))}


def group_of(r):
    sid = int(r["scenario"])
    return "A" if r["pd_verdict"] == "DEFERRED" else "BCD"[(sid // 3) % 3]


def check_and_join(rows):
    """The registration's checks, and the Jeffreys ends joined from S-9 for replicates 50 on."""
    wanted = {(r["scenario"], r["replicate"]) for r in rows}
    s9 = {}
    for r in pq.read_table(REPO / "studies/bayes-coverage/fits.parquet").to_pylist():
        if (r["scenario"], r["replicate"]) in wanted:
            s9[(r["scenario"], r["replicate"])] = r
    jef_bad = prof_bad = jef_checked = 0
    for r in rows:
        o = s9[(r["scenario"], r["replicate"])]
        for p in PARAMS:
            if (r[f"prof_lo_{p}"], r[f"prof_hi_{p}"], r[f"prof_flags_{p}"]) != (o[f"prof_lo_{p}"], o[f"prof_hi_{p}"],
                                                                                 o[f"prof_flags_{p}"]):
                if not (math.isnan(r[f"prof_lo_{p}"]) and math.isnan(o[f"prof_lo_{p}"])):
                    prof_bad += 1
        if r["replicate"] < 50:
            jef_checked += 1
            same = r["jef_flags"] == o["jeffreys_flags"] and all(
                r[f"jef_lo_{p}"] == o[f"jeffreys_et_lo_{p}"] and r[f"jef_hi_{p}"] == o[f"jeffreys_et_hi_{p}"] for p in PARAMS)
            jef_bad += not same
        else:
            r["jef_flags"] = o["jeffreys_flags"]
            for p in PARAMS:
                r[f"jef_lo_{p}"], r[f"jef_hi_{p}"] = o[f"jeffreys_et_lo_{p}"], o[f"jeffreys_et_hi_{p}"]
    print(f"check: Jeffreys recomputed on {jef_checked} replicates, {jef_bad} differ from S-9's rows")
    print(f"check: profile ends and flags, {prof_bad} of {2 * len(rows)} differ from S-9's rows")
    if jef_bad or prof_bad:
        sys.exit("a required check failed: the comparison stops until it is explained")


class Scenario:
    def __init__(self, rows, srow):
        r0 = rows[0]
        self.id, self.rows, self.R = r0["scenario"], rows, len(rows)
        self.truth = {"pd": r0["pd_true"], "rho": r0["rho_true"]}
        self.group = group_of(srow)
        self.iid = {p: float(srow[f"{p}_boot_coverage"]) for p in PARAMS}
        self.ka = {}
        for p in PARAMS:
            w = [r[f"w0_{p}"] for r in rows if math.isfinite(r[f"w0_{p}"])]
            self.ka[p] = sum(w) / len(w) if w else float("nan")
        half = 3.29 * math.sqrt(0.95 * 0.05 / self.R)
        self.band = (0.95 - half, 0.95 + half)

    def covers(self, r, m, p):
        t = self.truth[p]
        if m == "profile":
            return not r[f"prof_flags_{p}"] & NOT_COMPUTED and r[f"prof_lo_{p}"] <= t <= r[f"prof_hi_{p}"]
        if m == "jeffreys":
            return not r["jef_flags"] & (REFUSED | NUMERIC) and r[f"jef_lo_{p}"] <= t <= r[f"jef_hi_{p}"]
        if m in ("pct", "stud"):
            return r[f"{m}_lo_{p}"] <= t <= r[f"{m}_hi_{p}"]  # NaN (not computed) compares false
        w = r[f"w0_{p}"]
        k = self.ka[p] if m == "s4a" else r[f"kb_{p}"]
        return math.isfinite(w) and w <= C * k

    def coverage(self, m, p):
        return sum(self.covers(r, m, p) for r in self.rows) / self.R

    def in_band(self, c):
        return self.band[0] <= c <= self.band[1]

    def ends(self, r, m, p):
        key = {"profile": "prof", "jeffreys": "jef"}.get(m, m)
        return r[f"{key}_lo_{p}"], r[f"{key}_hi_{p}"]

    def width_ratio(self, m, p, ref="profile", replicates=None):
        num, den = [], []
        for r in self.rows:
            if replicates is not None and r["replicate"] >= replicates:
                continue
            a, b = self.ends(r, m, p), self.ends(r, ref, p)
            if all(math.isfinite(x) and 0 < x < 1 for x in a + b):
                num.append(logit(a[1]) - logit(a[0]))
                den.append(logit(b[1]) - logit(b[0]))
        return median(num) / median(den) if num and median(den) > 0 else float("nan")

    def paired(self, m, ref, p):
        n10 = n01 = 0
        for r in self.rows:
            x, y = self.covers(r, m, p), self.covers(r, ref, p)
            n10 += x and not y
            n01 += y and not x
        return n10, n01


def pooled(sc, m, p):
    tot = sum(s.R for s in sc)
    return sum(s.coverage(m, p) * s.R for s in sc) / tot if tot else float("nan")


def verdict(ok, complete):
    return ("yes" if ok else "no") if complete else f"subset ({'would hold' if ok else 'would not hold'})"


def score(sc):
    by = {g: [s for s in sc.values() if s.group == g] for g in "ABCD"}
    B, CD, BCD = by["B"], by["C"] + by["D"], by["B"] + by["C"] + by["D"]
    fullB, fullCD, fullBCD = len(B) == 12, len(CD) == 38, len(BCD) == 50
    print("\n| # | Prediction | Result | Held |\n|---|---|---|---|")

    prof_b = pooled(B, "profile", "rho")
    pct_b = pooled(B, "pct", "rho")
    below = sum(s.coverage("pct", "rho") < s.coverage("profile", "rho") for s in B)
    ok = 0.910 <= pct_b <= 0.935 and pct_b < prof_b and below >= 9
    print(f"| K1 | group B: pooled percentile rho coverage 0.910-0.935 and below the profile's; below it in >= 9 of 12 | "
          f"percentile {pct_b:.4f}, profile {prof_b:.4f}; below in {below} of {len(B)} | {verdict(ok, fullB)} |")

    stud_b, jef_b = pooled(B, "stud", "rho"), pooled(B, "jeffreys", "rho")
    wr = [s.width_ratio("stud", "rho") for s in B]
    k = sum(1.03 <= x <= 1.20 for x in wr)
    ok = jef_b <= stud_b <= 0.960 and k >= 10
    print(f"| K2 | group B: pooled studentised rho coverage >= Jeffreys' and <= 0.960; width ratio to profile 1.03-1.20 "
          f"in >= 10 of 12 | studentised {stud_b:.4f}, Jeffreys {jef_b:.4f}; width ratio in range in {k} of {len(B)} "
          f"({min(wr, default=float('nan')):.3f}-{max(wr, default=float('nan')):.3f}) | {verdict(ok, fullB)} |")

    parts, ok = [], True
    for m in ("pct", "stud"):
        for p in PARAMS:
            d = pooled(CD, m, p) - pooled(CD, "profile", p)
            wr = [s.width_ratio(m, p) for s in CD]
            k = sum(0.95 <= x <= 1.08 for x in wr)
            ok = ok and abs(d) <= 0.010 and k >= 32
            parts.append(f"{m} {p}: {d:+.4f}, widths in range {k} of {len(CD)}")
    print(f"| K3 | groups C-D: pooled coverage within 0.010 of the profile's; width ratio 0.95-1.08 in >= 32 of 38, per "
          f"method and parameter | {'; '.join(parts)} | {verdict(ok, fullCD)} |")

    k = sum(s.coverage("pct", "rho") > s.iid["rho"] for s in BCD)
    print(f"| K4 | groups B-D: parametric percentile rho coverage above the pinned iid percentile in >= 45 of 50 | "
          f"{k} of {len(BCD)} | {verdict(k >= 45, fullBCD)} |")

    parts, ok = [], True
    for m in ("pct", "stud"):
        n10 = sum(s.paired(m, "jeffreys", "rho")[0] for s in B)
        n01 = sum(s.paired(m, "jeffreys", "rho")[1] for s in B)
        z = (n10 - n01) / math.sqrt(n10 + n01) if n10 + n01 else 0.0
        eq = median([s.width_ratio(m, "rho", ref="jeffreys") for s in B])
        beats, equal = z > 2, 0.98 <= eq <= 1.02
        ok = ok and not (beats and equal)
        parts.append(f"{m}: paired z {z:+.2f} ({n10} vs {n01}), width vs Jeffreys {eq:.3f}")
    print(f"| K5 | neither bootstrap interval beats Jeffreys (paired z > 2, rho, group B) at equal width (0.98-1.02) | "
          f"{'; '.join(parts)} | {verdict(ok, fullB)} |")

    # S-4
    present = [f for f in FINDINGS if f[0] in sc]
    parts, ok = [], True
    for p in PARAMS:
        c = pooled(B, "s4a", p)
        ok = ok and 0.944 <= c <= 0.955
        parts.append(f"{p} {c:.4f}")
    inb = [f for f in present if sc[f[0]].in_band(sc[f[0]].coverage("s4a", f[1]))]
    ok = ok and len(inb) == len(present) == 6
    print(f"| H1 | oracle bound (in-sample; not a test): S-4a pooled group B 0.944-0.955, PD and rho; all six findings in "
          f"the band | {'; '.join(parts)}; findings in band {len(inb)} of {len(present)} | {verdict(ok, fullB and len(present) == 6)} |")

    parts, ok = [], True
    for p, (lo, hi) in (("rho", (1.05, 1.16)), ("pd", (1.03, 1.14))):
        m = median([s.ka[p] for s in B])
        ok = ok and lo <= m <= hi
        parts.append(f"B {p} {m:.3f}")
    for p in PARAMS:
        m = median([s.ka[p] for s in by["D"]])
        ok = ok and 0.98 <= m <= 1.05
        parts.append(f"D {p} {m:.3f}")
    print(f"| H2 | median k_a: group B rho 1.05-1.16, PD 1.03-1.14; group D 0.98-1.05 | {'; '.join(parts)} | "
          f"{verdict(ok, fullB and len(by['D']) == 21)} |")

    parts, ok = [], True
    for p in PARAMS:
        base = pooled(B, "profile", p)
        ga, gb = pooled(B, "s4a", p) - base, pooled(B, "s4b", p) - base
        ok = ok and gb >= 0.5 * ga
        parts.append(f"{p}: S-4a gain {ga:+.4f}, S-4b {gb:+.4f}")
    inb = [f for f in present if sc[f[0]].in_band(sc[f[0]].coverage("s4b", f[1]))]
    ok = ok and len(inb) >= 4
    print(f"| H3 | S-4b keeps at least half of S-4a's gain (group B, PD and rho); >= 4 of the six findings in the band | "
          f"{'; '.join(parts)}; findings in band {len(inb)} of {len(present)} | {verdict(ok, fullB and len(present) == 6)} |")

    parts, ok = [], True
    for p in PARAMS:
        d = pooled(CD, "s4b", p) - pooled(CD, "profile", p)
        ok = ok and abs(d) <= 0.005
        parts.append(f"{p} {d:+.4f}")
    left = sum(s.in_band(s.coverage("profile", p)) and not s.in_band(s.coverage("s4b", p)) for s in CD for p in PARAMS)
    ok = ok and left <= 2
    print(f"| H4 | groups C-D: S-4b pooled within 0.005 of the profile's; at most 2 in-band verdicts leave the band | "
          f"{'; '.join(parts)}; left the band {left} | {verdict(ok, fullCD)} |")

    s4b, jef = pooled(B, "s4b", "rho"), pooled(B, "jeffreys", "rho")
    print(f"| H5 | group B: S-4b pooled rho within +/-0.006 of Jeffreys' | S-4b {s4b:.4f}, Jeffreys {jef:.4f} "
          f"({s4b - jef:+.4f}) | {verdict(abs(s4b - jef) <= 0.006, fullB)} |")

    wr = [s.width_ratio("s4b", "rho", replicates=200) for s in B]
    k = sum(1.02 <= x <= 1.08 for x in wr)
    print(f"| H6 | group B: S-4b median rho width ratio to the profile 1.02-1.08 in >= 10 of 12 | {k} of {len(B)} "
          f"({min(wr, default=float('nan')):.3f}-{max(wr, default=float('nan')):.3f}) | {verdict(k >= 10, fullB)} |")


def write_summary(sc, path):
    methods = ("profile", "jeffreys", "pct", "stud", "s4a", "s4b")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["scenario", "group", "param", "method", "coverage", "in_band", "width_ratio_to_profile", "k_a"])
        for s in sc.values():
            for p in PARAMS:
                for m in methods:
                    c = s.coverage(m, p)
                    wr = s.width_ratio(m, p, replicates=200 if m.startswith("s4") else None) if m != "profile" else 1.0
                    w.writerow([s.id, s.group, p, m, f"{c:.4f}", int(s.in_band(c)), f"{wr:.4f}", f"{s.ka[p]:.4f}"])
    print(f"\nwrote {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fits")
    ap.add_argument("--summary-out")
    a = ap.parse_args()
    rows = pq.read_table(a.fits).to_pylist()
    check_and_join(rows)
    srows = summary_rows()
    by = {}
    for r in rows:
        by.setdefault(r["scenario"], []).append(r)
    sc = {sid: Scenario(sorted(v, key=lambda r: r["replicate"]), srows[sid]) for sid, v in sorted(by.items())}
    mism = sum((math.isfinite(r[f"w0_{p}"]) and r[f"w0_{p}"] <= C) != s.covers(r, "profile", p)
               for s in sc.values() for r in s.rows for p in PARAMS)
    print(f"check: W0 <= c against the profile interval's coverage, {mism} of {2 * len(rows)} replicates disagree")
    fb = sum(r[f"kb_fallback_{p}"] for r in rows for p in PARAMS)
    print(f"S-4b factor fell back to 1 in {fb} of {2 * len(rows)} (parameter, replicate) pairs")
    score(sc)
    if a.summary_out:
        write_summary(sc, a.summary_out)


if __name__ == "__main__":
    main()
