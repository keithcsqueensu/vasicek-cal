# SPDX-License-Identifier: Apache-2.0
"""Interval methods compared on the recovery panels: profile likelihood, Jeffreys equal-tailed, BCa, and the
parametric bootstrap's percentile and studentised intervals and the S-4b Bartlett-corrected profile interval.

A descriptive comparison, written after the studies it draws on (S-9, D-175; S-5, D-159); it scores
no prediction. For each method and parameter (PD, rho) it reports the 81 coverage verdicts against the
D-131 band (below / in / above), split into group A (near-uninformative, D-136) and groups B-D, and
the width of the intervals in logit units: per scenario the median over replicates, then across the
B-D scenarios the median and range of each method's width relative to the profile interval's.

Sources, all per replicate on the same pinned recovery panels (81 x 1,000):
  - profile likelihood and Jeffreys equal-tailed: studies/bayes-coverage/fits.parquet (S-9's rows; the
    profile ends there reproduce the pinned coverage exactly);
  - BCa (logit, B = 999 iid bootstrap and the jackknife acceleration): studies/jackknife-bias-rho/
    bca_ends.parquet (the S-5 method rerun with its ends written, see its MANIFEST entry);
  - parametric bootstrap percentile and studentised (B = 999, logit) and S-4b (the profile interval at
    the threshold c k_b, k_b the parametric bootstrap's mean W*): studies/parametric-bootstrap/
    fits.parquet (S-10 and S-4, D-181). S-4b covers when W0 <= c k_b, on every replicate; its ends, and
    so its widths, are solved on replicates 0-199 only. S-4a (the oracle factor) is not a usable
    interval and is not listed.
Coverage is over all replicates; an interval not computed does not cover (D-131). Widths are over the
replicates with an interval.

    uv run --no-project --with pyarrow==25.0.1 python studies/interval_comparison.py [--csv OUT]
"""

import argparse
import csv
import math
from pathlib import Path

import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
R = 1000
HALF = 3.29 * math.sqrt(0.95 * 0.05 / R)
NOT_COMPUTED = 1 << 2  # engine::kIntervalNotComputed
CHI2_95 = 3.841458820694124  # chi2_1 0.95 quantile, as studies/parametric-bootstrap/compare.py


def logit(v):
    return math.log(v) - math.log1p(-v)


def median(x):
    x = sorted(x)
    n = len(x)
    return float("nan") if n == 0 else x[n // 2] if n % 2 else 0.5 * (x[n // 2 - 1] + x[n // 2])


def groups():
    out = {}
    with open(REPO / "tests/golden/recovery/summary.csv", newline="") as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            sid = int(r["scenario"])
            out[sid] = "A" if r["pd_verdict"] == "DEFERRED" else "BCD"[(sid // 3) % 3]
    return out


def load():
    """(method, scenario, param) -> list of (covered, logit width or None); and (scenario, param) -> the
    profile interval's median logit width over replicates 0-199, S-4b's width reference."""
    res, ref200 = {}, {}
    for r in pq.read_table(ROOT / "bayes-coverage/fits.parquet").to_pylist():
        sid = r["scenario"]
        for p in ("pd", "rho"):
            t = r[f"{p}_true"]
            ok = not r[f"prof_flags_{p}"] & NOT_COMPUTED
            lo, hi = r[f"prof_lo_{p}"], r[f"prof_hi_{p}"]
            res.setdefault(("profile", sid, p), []).append(
                (ok and lo <= t <= hi, logit(hi) - logit(lo) if ok else None))
            lo, hi = r[f"jeffreys_et_lo_{p}"], r[f"jeffreys_et_hi_{p}"]
            fine = not r["jeffreys_flags"] & 0b1010  # refused or numeric: no interval
            res.setdefault(("jeffreys_et", sid, p), []).append(
                (fine and lo <= t <= hi, logit(hi) - logit(lo) if fine else None))
    bca = ROOT / "jackknife-bias-rho/bca_ends.parquet"
    if bca.exists():
        truth = {}
        for r in pq.read_table(ROOT / "bayes-coverage/fits.parquet", columns=["scenario", "pd_true", "rho_true"]).to_pylist():
            truth[r["scenario"]] = {"pd": r["pd_true"], "rho": r["rho_true"]}
        for r in pq.read_table(bca).to_pylist():
            sid = r["scenario"]
            for p in ("pd", "rho"):
                lo, hi = r[f"{p}_bca_lo_u"], r[f"{p}_bca_hi_u"]
                fine = math.isfinite(lo) and math.isfinite(hi)
                t = logit(truth[sid][p])
                res.setdefault(("bca", sid, p), []).append((fine and lo <= t <= hi, hi - lo if fine else None))
    pb = ROOT / "parametric-bootstrap/fits.parquet"
    if pb.exists():
        for r in pq.read_table(pb).to_pylist():
            sid = r["scenario"]
            for p in ("pd", "rho"):
                t = r[f"{p}_true"]
                for m in ("pct", "stud"):
                    lo, hi = r[f"{m}_lo_{p}"], r[f"{m}_hi_{p}"]
                    fine = math.isfinite(lo) and math.isfinite(hi)
                    res.setdefault((m, sid, p), []).append(
                        (fine and lo <= t <= hi, logit(hi) - logit(lo) if fine else None))
                w, k = r[f"w0_{p}"], r[f"kb_{p}"]
                lo, hi = r[f"s4b_lo_{p}"], r[f"s4b_hi_{p}"]
                fine = math.isfinite(lo) and math.isfinite(hi)
                res.setdefault(("s4b", sid, p), []).append(
                    (math.isfinite(w) and w <= CHI2_95 * k, logit(hi) - logit(lo) if fine else None))
                if r["replicate"] < 200 and not r[f"prof_flags_{p}"] & NOT_COMPUTED:  # S-4b's width reference
                    ref200.setdefault((sid, p), []).append(logit(r[f"prof_hi_{p}"]) - logit(r[f"prof_lo_{p}"]))
    else:
        print(f"note: {pb.relative_to(REPO)} not present (not committed, D-163): the S-10 and S-4b rows are "
              "left out; see studies/parametric-bootstrap/MANIFEST.json for the release and the regeneration")
    return res, {k: median(v) for k, v in ref200.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv")
    a = ap.parse_args()
    grp = groups()
    res, ref200 = load()
    methods = [m for m in ("profile", "jeffreys_et", "bca", "pct", "stud", "s4b") if any(k[0] == m for k in res)]
    rows = []
    for (m, sid, p), v in sorted(res.items()):
        c = sum(x[0] for x in v) / R
        cls = "below" if c < 0.95 - HALF else "above" if c > 0.95 + HALF else "in"
        w = median([x[1] for x in v if x[1] is not None])
        rows.append({"method": m, "scenario": sid, "group": grp[sid], "param": p, "coverage": c, "class": cls,
                     "width": w})
    by = {(r["method"], r["scenario"], r["param"]): r for r in rows}
    print("| Method | Parameter | Groups B-D (50): below / in / above | Group A (31): below / in / above | "
          "B-D width vs profile: median (range) | Below the band |")
    print("|---|---|---|---|---|---|")
    for m in methods:
        for p in ("pd", "rho"):
            sel = [r for r in rows if r["method"] == m and r["param"] == p]
            cnt = {g: [sum(r["class"] == k for r in sel if (r["group"] == "A") == (g == "A")) for k in ("below", "in", "above")]
                   for g in ("A", "BCD")}
            # S-4b's widths are solved on replicates 0-199: its reference is the profile's median there.
            ref = {r["scenario"]: ref200[(r["scenario"], p)] if m == "s4b" else by[("profile", r["scenario"], p)]["width"]
                   for r in sel}
            ratios = [r["width"] / ref[r["scenario"]] for r in sel if r["group"] != "A" and ref[r["scenario"]] > 0]
            below = ", ".join(f"{r['scenario']} ({r['group']})" for r in sel if r["class"] == "below") or "none"
            print(f"| {m} | {p} | {' / '.join(map(str, cnt['BCD']))} | {' / '.join(map(str, cnt['A']))} | "
                  f"{median(ratios):.3f} ({min(ratios):.3f}-{max(ratios):.3f}) | {below} |")
    if a.csv:
        with open(a.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
            w.writeheader()
            for r in rows:
                w.writerow({**r, "coverage": f"{r['coverage']:.3f}", "width": f"{r['width']:.6g}"})
        print(f"\nwrote {a.csv}")


if __name__ == "__main__":
    main()
