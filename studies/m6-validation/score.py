# SPDX-License-Identifier: Apache-2.0
"""Score the M6 validation run against PREDICTION.md (W1-W9; D-200).

    python studies/m6-validation/score.py ROWS.csv [--information FILE]

ROWS.csv is study_lrpd_serial_m6's output (--seed 0x7A11DA7E --replicates 100 --phi-intervals all). The expected
information (studies/lrpd-serial-correlation/information_m6.csv) marks the credible pairs and gives the SEs used to
standardise. W1 runs first and stops the scoring on any failure. Standard library only.
"""

import argparse
import csv
import math
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
INFO = ROOT / "studies" / "lrpd-serial-correlation" / "information_m6.csv"
SEED = "0x7A11DA7E"
PHIS = (0.0, 0.3, 0.6)
FLAG_QUAD = 1 << 2  # engine::kFlagQuadratureUnconverged (engine/refine.hpp): the doubled-rule check
FLAG_PLACEMENT = 1 << 8  # ar1::kFlagPlacementCheck (benchmark_models/ar1/search.hpp)
INTERVAL_NOT_COMPUTED = 1 << 2  # engine::kIntervalNotComputed


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def logit(p):
    return math.log(p) - math.log1p(-p)


def inv3(m):
    a, b, c, d, e, f, g, h, i = m
    det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    return (e * i - f * h) / det, (a * i - c * g) / det


def information(path):
    """(scenario, phi) -> credible flag, SEs, and the paired-difference variances at phi = 0."""
    out = {}
    for r in read_csv(path):
        key = (int(r["scenario"]), float(r["phi"]))
        ipp, irr, iff = float(r["i_pd_pd"]), float(r["i_rho_rho"]), float(r["i_phi_phi"])
        ipr, ipf, irf = float(r["i_pd_rho"]), float(r["i_pd_phi"]), float(r["i_rho_phi"])
        v3_pd, v3_rho = inv3((ipp, ipr, ipf, ipr, irr, irf, ipf, irf, iff))
        det2 = ipp * irr - ipr * ipr
        out[key] = {
            "credible": r["credible_pd"] == "1" and r["credible_rho"] == "1" and r["credible_phi"] == "1",
            "T": int(r["periods"]),
            "se": (float(r["se_pd"]), float(r["se_rho"]), float(r["se_phi"])),
            "se_logit": (float(r["se_logit_pd"]), float(r["se_logit_rho"])),
            "v3": (v3_pd, v3_rho),
            "v2": (irr / det2, ipp / det2),
            "se2_rho": math.sqrt(ipp / det2),
        }
    return out


def covers(r, a):
    name = ("pd", "rho", "phi")[a]
    if int(r[f"prof_flags_{name}"]) & INTERVAL_NOT_COMPUTED:
        return False
    lo, hi = float(r[f"prof_lo_{name}"]), float(r[f"prof_hi_{name}"])
    truth = (float(r["pd_true"]), float(r["rho_true"]), float(r["phi"]))[a]
    return lo <= truth <= hi


def verdict(ok):
    return "held" if ok else "NOT HELD"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rows")
    ap.add_argument("--information", default=str(INFO))
    a = ap.parse_args()
    with open(a.rows) as f:
        head = f.readline()
    if SEED.lower() not in head.lower():
        sys.exit(f"{a.rows}: not the validation run's rows (no seed {SEED} in its header)")
    rows = read_csv(a.rows)
    info = information(a.information)
    by_pair = {}
    for r in rows:
        by_pair.setdefault((int(r["scenario"]), float(r["phi"])), []).append(r)
    expected = {(s, p) for s in range(81) for p in PHIS}
    if set(by_pair) != expected or any(len(v) != 100 for v in by_pair.values()):
        sys.exit("the rows are not 81 scenarios x 3 phi x 100 replicates")

    # W1: checks before scoring.
    bad = []
    for r in rows:
        if r["fitted"] != "1" or not math.isfinite(float(r["loglik"])):
            bad.append((r["scenario"], r["phi"], r["replicate"], "not fitted"))
        elif float(r["phi0_loglik"]) < float(r["binom_loglik"]) - 1e-8:
            bad.append((r["scenario"], r["phi"], r["replicate"], "phi = 0 point below the binomial MLE"))
        elif float(r["loglik"]) < float(r["phi0_loglik"]) - 1e-10:
            bad.append((r["scenario"], r["phi"], r["replicate"], "maximum below the phi = 0 point"))
    print(f"W1 checks before scoring: {len(rows) - len(bad)} of {len(rows)} fits pass -> {verdict(not bad)}")
    if bad:
        for b in bad[:20]:
            print("   ", b)
        sys.exit("W1 failed: scoring stops")

    cred100 = sorted(k for k, v in info.items() if v["credible"] and v["T"] == 100)
    cred40 = sorted(k for k, v in info.items() if v["credible"] and v["T"] == 40)
    cred100_0 = [k for k in cred100 if k[1] == 0.0]
    held = {}

    # W2 and W3: phi = 0 against the binomial MLE, paired.
    worst = [0.0, 0.0]
    haus = [[], []]
    for k in cred100_0:
        inf = info[k]
        for ai, (est, ref) in enumerate((("pd", "phi0_pd"), ("rho", "phi0_rho"))):
            d_logit = sorted(abs(logit(float(r[est])) - logit(float(r[ref]))) / inf["se_logit"][ai] for r in by_pair[k])
            q95 = d_logit[94]  # the 95th of 100, type 1
            worst[ai] = max(worst[ai], q95)
            sd = math.sqrt(max(inf["v3"][ai] - inf["v2"][ai], 0.0))
            haus[ai] += [(float(r[est]) - float(r[ref])) / sd for r in by_pair[k]]
    held["W2"] = worst[0] <= 0.25 and worst[1] <= 0.25
    print(f"W2 phi = 0 paired, worst 95th percentile |logit diff| / SE: PD {worst[0]:.3f}, rho {worst[1]:.3f} "
          f"(<= 0.25) -> {verdict(held['W2'])}")
    hsd = [statistics.pstdev(h) for h in haus]
    held["W3"] = all(0.7 <= s <= 1.5 for s in hsd)
    print(f"W3 SD of (M6 - binomial) / sqrt(v3 - v2): PD {hsd[0]:.3f}, rho {hsd[1]:.3f} (0.7-1.5) -> {verdict(held['W3'])}")

    # W4: Kendall's bias at phi = 0.
    med_phi = statistics.median(float(r["phi_hat"]) for k in cred100_0 for r in by_pair[k])
    held["W4"] = -0.05 <= med_phi <= 0.01
    print(f"W4 median phi-hat at phi = 0 (15 cells): {med_phi:+.4f} (-0.05 to +0.01) -> {verdict(held['W4'])}")

    # W5 and W6: coverage.
    def pooled(keys, a):
        return statistics.fmean(1.0 if covers(r, a) else 0.0 for k in keys for r in by_pair[k])

    cov100 = [pooled(cred100, a) for a in range(3)]
    per_pair = {k: [pooled([k], a) for a in range(3)] for k in cred100}
    low = [(k, c) for k, c in per_pair.items() if min(c) < 0.86]
    held["W5"] = all(0.93 <= c <= 0.97 for c in cov100) and not low
    print(f"W5 coverage, 51 credible T = 100 pairs: PD {cov100[0]:.4f}, rho {cov100[1]:.4f}, phi {cov100[2]:.4f} "
          f"(0.93-0.97); pairs below 0.86: {len(low)} -> {verdict(held['W5'])}")
    cov40 = [pooled(cred40, a) for a in range(3)]
    in40 = all(0.90 <= c <= 0.97 for c in cov40)
    print(f"W6 coverage, 13 credible T = 40 pairs: PD {cov40[0]:.4f}, rho {cov40[1]:.4f}, phi {cov40[2]:.4f} "
          f"(expected range 0.90-0.97, not pass/fail): {'inside' if in40 else 'OUTSIDE'}")

    # W7: PD's median bias.
    pd_bias = {k: statistics.median((float(r["pd"]) - float(r["pd_true"])) / info[k]["se"][0] for r in by_pair[k])
               for k in cred100}
    held["W7"] = all(abs(b) <= 0.5 for b in pd_bias.values())
    print(f"W7 PD median bias in SE, worst over 51 pairs: {max(pd_bias.values(), key=abs):+.3f} (|.| <= 0.5) -> "
          f"{verdict(held['W7'])}")

    # W8: rho's bias against the iid MLE's at phi = 0 in the same scenario.
    diffs = {}
    for k in cred100:
        s = k[0]
        m6 = statistics.median((float(r["rho"]) - float(r["rho_true"])) / info[k]["se"][1] for r in by_pair[k])
        iid = statistics.median((float(r["phi0_rho"]) - float(r["rho_true"])) / info[(s, 0.0)]["se2_rho"]
                                for r in by_pair[(s, 0.0)])
        diffs[k] = m6 - iid
    share = sum(abs(d) <= 0.35 for d in diffs.values()) / len(diffs)
    held["W8"] = share >= 0.9
    print(f"W8 rho median bias minus the iid MLE's: {share:.1%} of 51 pairs within 0.35 SE (>= 90%) -> "
          f"{verdict(held['W8'])}")

    # W9: flags.
    def frac(rs, bit):
        return sum(1 for r in rs if int(r["flags"]) & bit) / len(rs)

    fp, fq = frac(rows, FLAG_PLACEMENT), frac(rows, FLAG_QUAD)
    worst_pair = max(max(frac(v, FLAG_PLACEMENT), frac(v, FLAG_QUAD)) for v in by_pair.values())
    held["W9"] = fp <= 0.01 and fq <= 0.01 and worst_pair <= 0.05
    print(f"W9 flags: placement {fp:.2%}, doubled rule {fq:.2%} (<= 1%); worst pair {worst_pair:.0%} (<= 5%) -> "
          f"{verdict(held['W9'])}")

    print("\nper credible pair (scenario, phi, T): coverage PD rho phi | PD median bias (SE) | rho bias minus iid (SE)")
    for k in cred100 + cred40:
        c = [pooled([k], a) for a in range(3)]
        extra = (f" | {pd_bias[k]:+.3f} | {diffs[k]:+.3f}" if k in pd_bias else "")
        print(f"  {k[0]:2d} {k[1]:.1f} {info[k]['T']:3d}: {c[0]:.2f} {c[1]:.2f} {c[2]:.2f}{extra}")
    strict = ("W2", "W3", "W4", "W5", "W7", "W8", "W9")
    n_held = sum(held[w] for w in strict) + 1  # W1 held to get here
    print(f"\n{n_held} of 8 strict predictions held (W1-W5, W7-W9); W6 reported against its range")
    return 0 if n_held == 8 else 1


if __name__ == "__main__":
    sys.exit(main())
