# SPDX-License-Identifier: Apache-2.0
"""S-10, K2's miss: a post-hoc diagnosis (written after the scoring; not a registered prediction).

Group B (T = 20), rho: for each interval, how often it was not computed, and on which side a
miss falls (the truth below the lower end, i.e. the interval too high, or above the upper end);
then, for the studentised interval, how the logit-scale analytic SE moves with the estimate
(a studentised pivot is only pivotal if the SE does not carry the estimate's error).

    uv run --no-project --with pyarrow==25.0.1 python studies/parametric-bootstrap/diagnose_k2.py FITS_PARQUET
"""

import math
import sys

import pyarrow.parquet as pq

sys.path.insert(0, __file__.rsplit("diagnose_k2.py", 1)[0])
import compare as cmp


def main():
    rows = pq.read_table(sys.argv[1]).to_pylist()
    summ = cmp.summary_rows()
    rows = [r for r in rows if cmp.group_of(summ[r["scenario"]]) == "B"]
    cmp.check_and_join(rows)
    sc = {}
    for r in rows:
        sc.setdefault(r["scenario"], []).append(r)
    scen = [cmp.Scenario(v, summ[k]) for k, v in sorted(sc.items())]
    n = sum(s.R for s in scen)
    print(f"group B: {len(scen)} scenarios, {n} replicates; rho\n")
    print("| interval | coverage | not computed | miss: truth below lower end | miss: truth above upper end |")
    print("|---|---|---|---|---|")
    for m in ("profile", "jeffreys", "pct", "stud"):
        cov = nc = lo = hi = 0
        for s in scen:
            t = s.truth["rho"]
            for r in s.rows:
                a, b = s.ends(r, m, "rho")
                if s.covers(r, m, "rho"):
                    cov += 1
                elif not (math.isfinite(a) and math.isfinite(b)):
                    nc += 1
                elif t < a:
                    lo += 1
                else:
                    hi += 1
        print(f"| {m} | {cov / n:.4f} | {nc / n:.4f} | {lo / n:.4f} | {hi / n:.4f} |")

    # The analytic SE against the estimate's error, on the logit scale, per scenario.
    print("\n| scenario | corr(logit rho_hat - logit rho, SE) | stud cov: rho_hat below truth | above truth |"
          " SE ratio, below / above |")
    print("|---|---|---|---|---|")
    for s in scen:
        t = cmp.logit(s.truth["rho"])
        e, se, cb, ca, sb, sa = [], [], [0, 0], [0, 0], [], []
        for r in s.rows:
            x, v = cmp.logit(r["rho"]) - t, r["se_analytic_u_rho"]
            if not math.isfinite(v):
                continue
            e.append(x)
            se.append(v)
            c = s.covers(r, "stud", "rho")
            if x < 0:
                cb[0] += c
                cb[1] += 1
                sb.append(v)
            else:
                ca[0] += c
                ca[1] += 1
                sa.append(v)
        me, ms = sum(e) / len(e), sum(se) / len(se)
        cov = sum((a - me) * (b - ms) for a, b in zip(e, se))
        corr = cov / math.sqrt(sum((a - me) ** 2 for a in e) * sum((b - ms) ** 2 for b in se))
        print(f"| {s.id} | {corr:+.3f} | {cb[0] / max(cb[1], 1):.3f} ({cb[1]}) | {ca[0] / max(ca[1], 1):.3f} ({ca[1]})"
              f" | {cmp.median(sb) / cmp.median(sa):.3f} |")


if __name__ == "__main__":
    main()
