# SPDX-License-Identifier: Apache-2.0
"""S-35: score the run against PREDICTION.md (L1-L15), after the checks it requires.

Inputs: the rows of study_lrpd_serial (--out and --rate-out) and pinned files only:
  tests/golden/recovery/summary.csv               phi = 0 counts of the binomial MLE, the profile band
  tests/golden/vasicek_rate/recovery_summary.csv  phi = 0 counts of the rate MLE (n -> infinity)
  studies/mle-vs-mom/fits.parquet                 phi = 0 rows of the observed-rate arm (S-8's censor_*)
  reference.csv, reference_rate.csv, reference_observed.csv   the registered reference numbers

Checks (any failure stops the scoring, PREDICTION.md "Checks before scoring" and the addendum):
  - binomial, phi = 0: per scenario, the Wald, t, profile and bootstrap covered counts and the flagged count
    for PD and rho equal the pinned summary's;
  - rate MLE, phi = 0: per cell, the PD, rho and q covered counts equal the pinned rate summary's;
  - observed-rate arm, phi = 0: every row's estimates and coverage equal S-8's, bit for bit.

Pooling: "pooled" coverage is over all replicates of the scenarios pooled; with R equal across them it is the
mean of the per-scenario coverages, which is how the reference means are formed. Wald and t coverage is
among unflagged fits (indicator -1 excluded), as pinned; profile, bootstrap and every rate-MLE coverage is
over all replicates, a not-computed interval missing. A run on part of the matrix (the study subset) is
scored with the same rules on the scenarios present; every figure is printed with its denominator, and a
prediction whose registered denominator is incomplete reads "subset" instead of yes or no.

    uv run --no-project --with scipy==1.16.2 --with pyarrow==25.0.1 python \\
        studies/lrpd-serial-correlation/compare.py ROWS.csv RATE_ROWS.csv [--summary-out FILE]
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics as st
import sys
from pathlib import Path

import pyarrow.parquet as pq
from scipy.special import ndtr, ndtri

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
GROUP_A = {0, 1, 2, 3, 4, 6, 7, 9, 10, 12, 15, 18, 19, 20, 21, 22, 23, 24, 27, 28, 30, 33, 36, 39, 45, 46, 48, 54,
           57, 60, 72}
PHIS = (0.0, 0.3, 0.6)
INFORMATIVE = 50


def group(sid: int, periods: int) -> str:
    return "A" if sid in GROUP_A else {20: "B", 40: "C", 100: "D"}[periods]


def read_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(line for line in fh if not line.startswith("#")))


def logit(v: float) -> float:
    return math.log(v) - math.log1p(-v)


def mean(xs):
    xs = list(xs)
    return sum(xs) / len(xs) if xs else float("nan")


def cov(ind):
    """Coverage among indicators: 1 covers; 0 misses; -1 excluded (Wald/t flagged)."""
    kept = [x for x in ind if x >= 0]
    return sum(kept) / len(kept) if kept else float("nan")


def cov_all(ind):
    """Coverage over all replicates: -1 (not computed or not estimated) misses."""
    return sum(x == 1 for x in ind) / len(ind) if ind else float("nan")


def map_cov(c0: float, f: float) -> float:
    c0 = min(max(c0, 1e-6), 1 - 1e-6)
    z0 = float(ndtri((1 + c0) / 2))
    return float(ndtr(z0 / math.sqrt(f)) - ndtr(-z0 / math.sqrt(f)))


class Report:
    def __init__(self):
        self.lines, self.failed = [], False

    def row(self, label, prediction, result, held):
        self.lines.append((label, prediction, result, held))
        print(f"| {label} | {prediction} | {result} | {held} |")


def held(ok: bool, complete: bool) -> str:
    return ("yes" if ok else "**no**") if complete else ("subset: " + ("yes" if ok else "no"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("rows")
    ap.add_argument("rate_rows")
    ap.add_argument("--summary-out")
    a = ap.parse_args()

    rows = read_rows(a.rows)
    rate = read_rows(a.rate_rows)
    pinned = {int(r["scenario"]): r for r in read_rows(REPO / "tests/golden/recovery/summary.csv")}
    pinned_rate = {int(r["z_source_scenario"]): r for r in read_rows(REPO / "tests/golden/vasicek_rate/recovery_summary.csv")}
    ref = {(int(r["scenario"]), float(r["phi"])): r for r in read_rows(HERE / "reference.csv")}
    ref_rate = read_rows(HERE / "reference_rate.csv")
    ref_obs = {(int(r["scenario"]), float(r["phi"])): r for r in read_rows(HERE / "reference_observed.csv")}

    # Index the rows: by (phi, scenario) -> list of rows in replicate order.
    by = {}
    for r in rows:
        by.setdefault((float(r["phi"]), int(r["scenario"])), []).append(r)
    for v in by.values():
        v.sort(key=lambda r: int(r["replicate"]))
    rby = {}
    for r in rate:
        rby.setdefault((float(r["phi"]), int(r["z_source_scenario"])), []).append(r)
    for v in rby.values():
        v.sort(key=lambda r: int(r["replicate"]))
    scen = sorted({s for (_, s) in by})
    R = len(by[(0.0, scen[0])])
    meta = {s: by[(0.0, s)][0] for s in scen}
    per = {s: int(meta[s]["periods"]) for s in scen}
    grp = {s: group(s, per[s]) for s in scen}
    full = len(scen) == 81 and R == 1000
    cells = sorted({s for (_, s) in rby})
    print(f"rows: {len(scen)} scenarios x {R} replicates x {len(PHIS)} phi; rate cells: {len(cells)}; "
          f"{'full matrix' if full else 'PART OF THE MATRIX: predictions read subset'}")

    # ---- checks -------------------------------------------------------------------------------------
    bad = []
    if R == 1000:
        for s in scen:
            p, rs = pinned[s], by[(0.0, s)]
            for name, col in (("pd", "pd"), ("rho", "rho")):
                got = {
                    "covered": sum(int(r[f"wald_{col}"]) == 1 for r in rs),
                    "t_covered": sum(int(r[f"t_{col}"]) == 1 for r in rs),
                    "profile_covered": sum(int(r[f"prof_{col}"]) == 1 for r in rs),
                    "boot_covered": sum(int(r[f"boot_{col}"]) == 1 for r in rs),
                }
                for k, v in got.items():
                    if v != int(p[f"{name}_{k}"]):
                        bad.append(f"binomial scenario {s} {name}_{k}: {v} against pinned {p[f'{name}_{k}']}")
            flagged = sum(int(r["wald_pd"]) == -1 for r in rs)
            if flagged != int(p["flagged"]):
                bad.append(f"binomial scenario {s} flagged: {flagged} against pinned {p['flagged']}")
        for c in cells:
            p, rs = pinned_rate[c], rby[(0.0, c)]
            for k, col in (("pd_covered", "cover_pd"), ("rho_covered", "cover_rho"), ("q_covered", "cover_q")):
                v = sum(int(r[col]) == 1 for r in rs)
                if v != int(p[k]):
                    bad.append(f"rate cell {c} {k}: {v} against pinned {p[k]}")
    else:
        print("checks: the pinned counts are for R = 1,000; the count checks are skipped at this R")
    s8 = pq.read_table(REPO / "studies/mle-vs-mom/fits.parquet",
                       columns=["scenario", "replicate", "censor_pd", "censor_rho", "censor_cover_pd",
                                "censor_cover_rho"]).to_pylist()
    s8 = {(r["scenario"], r["replicate"]): r for r in s8 if r["scenario"] < 81}
    obs_checked = 0
    for s in scen:
        for r in by[(0.0, s)]:
            q = s8[(s, int(r["replicate"]))]
            obs_checked += 1
            if (float(r["obs_pd"]) != q["censor_pd"] or float(r["obs_rho"]) != q["censor_rho"]
                    or int(r["obs_cover_pd"]) != q["censor_cover_pd"] or int(r["obs_cover_rho"]) != q["censor_cover_rho"]):
                bad.append(f"observed-rate scenario {s} replicate {r['replicate']} differs from S-8's row")
    print(f"checks: {len(bad)} failures ({obs_checked} observed-rate rows compared with S-8 bit for bit)")
    for b in bad[:20]:
        print("  " + b)
    if bad:
        print("CHECKS FAILED: not scored", file=sys.stderr)
        return 1

    # ---- per (phi, scenario) summaries --------------------------------------------------------------
    summ = {}
    for (phi, s), rs in by.items():
        pd, rho, t = float(meta[s]["pd_true"]), float(meta[s]["rho_true"]), per[s]
        avg = [float(r["avg"]) for r in rs]
        x = {
            "avg_mean": mean(avg), "avg_sd": st.stdev(avg), "avg_below": mean(v < pd for v in avg),
            "avg_cov": mean(float(r["avg_lo"]) <= pd <= float(r["avg_hi"]) for r in rs),
            "avg_median_rel": st.median(avg) / pd - 1,
            "pd_median_rel": st.median(float(r["pd"]) for r in rs) / pd - 1,
            "rho_rel_bias": mean(float(r["rho"]) for r in rs) / rho - 1,
            "q_prof": cov_all([int(r["q_prof"]) for r in rs]),
            "flagged": sum(int(r["wald_pd"]) == -1 for r in rs),
            "obs_rho_u": mean(logit(float(r["obs_rho"])) for r in rs),
        }
        for col in ("pd", "rho"):
            x[f"wald_{col}"] = cov([int(r[f"wald_{col}"]) for r in rs])
            x[f"t_{col}"] = cov([int(r[f"t_{col}"]) for r in rs])
            x[f"prof_{col}"] = cov_all([int(r[f"prof_{col}"]) for r in rs])
            x[f"boot_{col}"] = cov_all([int(r[f"boot_{col}"]) for r in rs])
            x[f"obs_{col}"] = cov_all([int(r[f"obs_cover_{col}"]) for r in rs])
        x["periods"], x["obligors"], x["group"] = t, int(meta[s]["obligors"]), grp[s]
        summ[(phi, s)] = x
    rsumm = {}
    for (phi, c), rs in rby.items():
        rho, t = float(rs[0]["rho_true"]), int(rs[0]["periods"])
        s2 = rho / (1 - rho)
        ratios = [(float(r["rho"]) / (1 - float(r["rho"]))) / s2 for r in rs]
        rsumm[(phi, c)] = {
            "periods": t, "ratio_mean": mean(ratios), "ratio_mcse": st.stdev(ratios) / math.sqrt(len(ratios)),
            "pd": cov_all([int(r["cover_pd"]) for r in rs]), "rho": cov_all([int(r["cover_rho"]) for r in rs]),
            "flagged": sum(int(r["flags"]) != 0 for r in rs),
        }

    def pooled(phi, key, sel):
        return mean(summ[(phi, s)][key] for s in scen if sel(s))

    inf = [s for s in scen if grp[s] != "A"]
    n_inf = len(inf)
    print()
    print("| # | Prediction | Result | Held |")
    print("|---|---|---|---|")
    rep = Report()

    # L1
    ok = [abs(summ[(phi, s)]["avg_mean"] - float(meta[s]["pd_true"])) <= 3 * summ[(phi, s)]["avg_sd"] / math.sqrt(R)
          for phi in PHIS[1:] for s in scen]
    rep.row("L1", "the average within 3 MCSE of PD in >= 160 of 162", f"{sum(ok)} of {len(ok)}",
            held(sum(ok) >= 160, full))
    # L2
    up = sum(summ[(0.6, s)]["avg_below"] > summ[(0.0, s)]["avg_below"] for s in inf)
    gb = pooled(0.6, "avg_below", lambda s: grp[s] == "B")
    rep.row("L2", "P(avg < PD) higher at 0.6 than 0 in >= 40 of 50; group B mean at 0.6 in 0.56-0.70",
            f"{up} of {n_inf}; group B {gb:.3f}", held(up >= 40 and 0.56 <= gb <= 0.70, full))
    # L3
    f_mu = {(phi, t): 1 + 2 * sum((1 - k / t) * phi ** k for k in range(1, t)) for phi in PHIS for t in (20, 40, 100)}
    ok3 = [abs(rsumm[(phi, c)]["ratio_mean"] - (1 - f_mu[(phi, rsumm[(phi, c)]["periods"])] / rsumm[(phi, c)]["periods"]))
           <= 3 * rsumm[(phi, c)]["ratio_mcse"] for phi in PHIS[1:] for c in cells]
    boxed = sum(rsumm[(phi, c)]["flagged"] for phi in PHIS[1:] for c in cells)
    rep.row("L3", "rate MLE mean sigma2-hat/sigma2 within 3 MCSE of 1 - f_mu/T in >= 50 of 54",
            f"{sum(ok3)} of {len(ok3)} ({boxed} flagged fits)", held(sum(ok3) >= 50, len(cells) == 27 and R == 1000))
    # L4
    res4, ok4 = [], True
    for phi in PHIS[1:]:
        for t in (20, 40, 100):
            for key in ("pd", "rho"):
                got = mean(rsumm[(phi, c)][key] for c in cells if rsumm[(phi, c)]["periods"] == t)
                want = mean(float(r[f"profile_{key}_pred"]) for r in ref_rate
                            if float(r["phi"]) == phi and int(r["periods"]) == t)
                ok4 = ok4 and abs(got - want) <= 0.03
                res4.append(f"{key} T{t} phi {phi}: {got:.3f} ({want:.3f})")
    rep.row("L4", "rate MLE pooled profile coverage per T and phi within +/-0.03 of the reference",
            "; ".join(res4), held(ok4, len(cells) == 27 and R == 1000))
    # L5
    res5, ok5 = [], True
    for phi in PHIS[1:]:
        for g in "BCD":
            for key in ("pd", "rho"):
                got = pooled(phi, f"prof_{key}", lambda s, g=g: grp[s] == g)
                want = mean(float(ref[(s, phi)][f"profile_{key}_pred"]) for s in scen if grp[s] == g)
                ok5 = ok5 and abs(got - want) <= 0.03
                res5.append(f"{key} {g} {phi}: {got:.3f} ({want:.3f})")
    rep.row("L5", "binomial pooled profile coverage per group and phi within +/-0.03", "; ".join(res5), held(ok5, full))
    # L6
    fall = {(g, phi): pooled(phi, "prof_pd", lambda s, g=g: grp[s] == g) - pooled(0.0, "prof_pd", lambda s, g=g: grp[s] == g)
            for g in "ABCD" for phi in PHIS[1:]}
    ok6 = all(fall[("A", phi)] > fall[("D", phi)] for phi in PHIS[1:])
    rep.row("L6", "group A's PD profile fall smaller than group D's at both phi",
            "; ".join(f"phi {phi}: A {fall[('A', phi)]:+.3f}, D {fall[('D', phi)]:+.3f}" for phi in PHIS[1:]),
            held(ok6, full))
    # L7
    worst, ok7 = 0.0, True
    for phi in PHIS[1:]:
        for g in "BCD":
            for key in ("pd", "rho"):
                pf = pooled(phi, f"prof_{key}", lambda s, g=g: grp[s] == g) - pooled(0.0, f"prof_{key}", lambda s, g=g: grp[s] == g)
                for m in ("wald", "t", "boot"):
                    f_ = (pooled(phi, f"{m}_{key}", lambda s, g=g: grp[s] == g)
                          - pooled(0.0, f"{m}_{key}", lambda s, g=g: grp[s] == g))
                    worst = max(worst, abs(f_ - pf))
                    ok7 = ok7 and abs(f_ - pf) <= 0.04
    rep.row("L7", "Wald, t and bootstrap falls within +/-0.04 of the profile's", f"largest difference {worst:.3f}",
            held(ok7, full))
    # L8
    b6 = pooled(0.6, "prof_pd", lambda s: grp[s] == "B")
    d6 = pooled(0.6, "prof_pd", lambda s: grp[s] == "D")
    r20 = mean(rsumm[(0.6, c)]["pd"] for c in cells if rsumm[(0.6, c)]["periods"] == 20)
    r100 = mean(rsumm[(0.6, c)]["pd"] for c in cells if rsumm[(0.6, c)]["periods"] == 100)
    rep.row("L8", "phi 0.6: |B - D| < 0.04 (binomial PD profile) and |T20 - T100| < 0.04 (rate PD)",
            f"B {b6:.3f}, D {d6:.3f}; rate T20 {r20:.3f}, T100 {r100:.3f}",
            held(abs(b6 - d6) < 0.04 and abs(r20 - r100) < 0.04, full))
    # L9
    band_lo = {s: float.fromhex(pinned[s]["profile_band_lo_hex"]) for s in scen}
    below = {(phi, key): sum(summ[(phi, s)][f"prof_{key}"] < band_lo[s] for s in inf) for phi in PHIS[1:] for key in ("pd", "rho")}
    ok9 = (below[(0.3, "pd")] == n_inf and 15 <= below[(0.3, "rho")] <= 40 and below[(0.6, "pd")] == n_inf
           and below[(0.6, "rho")] >= 45)
    rep.row("L9", "below the band: phi 0.3 all 50 PD, 15-40 rho; phi 0.6 all 50 PD, >= 45 rho",
            f"phi 0.3: PD {below[(0.3, 'pd')]}, rho {below[(0.3, 'rho')]}; phi 0.6: PD {below[(0.6, 'pd')]}, "
            f"rho {below[(0.6, 'rho')]} (of {n_inf})", held(ok9, full))
    # L10
    ch = {g: mean(summ[(0.6, s)]["rho_rel_bias"] - summ[(0.0, s)]["rho_rel_bias"] for s in scen if grp[s] == g) for g in "BD"}
    rep.row("L10", "change in mean relative bias of rho-hat, 0 to 0.6: B in -0.17..-0.06, D in -0.05..-0.01",
            f"B {ch['B']:+.3f}, D {ch['D']:+.3f}", held(-0.17 <= ch["B"] <= -0.06 and -0.05 <= ch["D"] <= -0.01, full))
    # L11
    res11, ok11 = [], True
    for phi in PHIS[1:]:
        for g in "BCD":
            got = pooled(phi, "avg_cov", lambda s, g=g: grp[s] == g)
            want = mean(map_cov(summ[(0.0, s)]["avg_cov"], float(ref[(s, phi)]["f_avg"])) for s in scen if grp[s] == g)
            ok11 = ok11 and abs(got - want) <= 0.04
            res11.append(f"{g} {phi}: {got:.3f} ({want:.3f})")
    rep.row("L11", "average's t interval per group and phi within +/-0.04 of its mapping", "; ".join(res11),
            held(ok11, full))
    # L12
    res12, ok12 = [], True
    for phi in PHIS[1:]:
        for g in "ABCD":
            got = pooled(phi, "obs_pd", lambda s, g=g: grp[s] == g)
            want = mean(float(ref_obs[(s, phi)]["pd_cov_pred"]) for s in scen if grp[s] == g)
            ok12 = ok12 and abs(got - want) <= 0.04
            res12.append(f"{g} {phi}: {got:.3f} ({want:.3f})")
    rep.row("L12", "observed-rate PD coverage per group and phi within +/-0.04", "; ".join(res12), held(ok12, full))
    # L13
    res13, ok13 = [], True
    for phi in PHIS[1:]:
        for n in (100, 1000, 10000):
            got = pooled(phi, "obs_rho", lambda s, n=n: summ[(0.0, s)]["obligors"] == n)
            want = mean(float(ref_obs[(s, phi)]["rho_cov_pred"]) for s in scen if summ[(0.0, s)]["obligors"] == n)
            ok13 = ok13 and abs(got - want) <= 0.04
            res13.append(f"n {n} {phi}: {got:.3f} ({want:.3f})")
    rep.row("L13", "observed-rate rho coverage per n and phi within +/-0.04", "; ".join(res13), held(ok13, full))
    # L14
    d100 = (pooled(0.6, "obs_rho", lambda s: summ[(0.0, s)]["obligors"] == 100)
            - pooled(0.0, "obs_rho", lambda s: summ[(0.0, s)]["obligors"] == 100))
    d1e4 = (pooled(0.6, "obs_rho", lambda s: summ[(0.0, s)]["obligors"] == 10000)
            - pooled(0.0, "obs_rho", lambda s: summ[(0.0, s)]["obligors"] == 10000))
    rep.row("L14", "observed rho, 0 to 0.6: n = 100 within 0.03; n = 10,000 at least 0.06 lower",
            f"n 100 {d100:+.3f}; n 10,000 {d1e4:+.3f}", held(abs(d100) <= 0.03 and d1e4 <= -0.06, full))
    # L15
    down = sum(summ[(0.6, s)]["obs_rho_u"] < summ[(0.0, s)]["obs_rho_u"] for s in inf)
    rep.row("L15", "observed mean logit rho-hat lower at 0.6 than 0 in >= 45 of 50", f"{down} of {n_inf}",
            held(down >= 45, full))

    print()
    print("Reported without a prediction (median relative error of PD-hat and the average; q profile coverage; "
          "flagged fits), pooled over the informative groups:")
    for phi in PHIS:
        print(f"  phi {phi}: PD-hat median rel {pooled(phi, 'pd_median_rel', lambda s: grp[s] != 'A'):+.4f}; "
              f"average median rel {pooled(phi, 'avg_median_rel', lambda s: grp[s] != 'A'):+.4f}; "
              f"q profile {pooled(phi, 'q_prof', lambda s: grp[s] != 'A'):.3f}; "
              f"flagged {sum(summ[(phi, s)]['flagged'] for s in scen)}")

    if a.summary_out:
        keys = [k for k in next(iter(summ.values())) if k not in ("group",)]
        with open(a.summary_out, "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(["phi", "scenario", "group", *keys])
            for (phi, s), x in sorted(summ.items()):
                w.writerow([phi, s, x["group"], *[f"{x[k]:.6g}" if isinstance(x[k], float) else x[k] for k in keys]])
    return 0


if __name__ == "__main__":
    sys.exit(main())
