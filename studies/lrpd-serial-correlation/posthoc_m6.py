# SPDX-License-Identifier: Apache-2.0
"""S-35's M6 arm, post-hoc diagnosis of L21 and L28 (D-208). Written after the run; not a prediction; changes no
verdict. The questions were committed in D-208 before this was run. Reads the phi = 0.6 rows only, the only ones
those two questions use.

L21: of the 50 informative scenarios, where did the AR(1) MLE's median relative error of PD-hat move closer to 0 than
the binomial MLE's? By T, n and PD: the count closer and the median shift (M6 minus binomial, relative error).

L28: per subset scenario (replicates 0-99), the percentile interval's PD coverage against the profile interval's, the
side of its misses, and its bias: the bootstrap PD-hat's median is not in the rows, so the percentile interval's
midpoint stands in for it, against the panel's PD-hat and against the true PD.

    python studies/lrpd-serial-correlation/posthoc_m6.py rows.csv
"""

from __future__ import annotations

import csv
import statistics as st
import sys
from collections import defaultdict

from score_m6 import (  # the scorer's definitions, unchanged
    BOOT_REPLICATES,
    SUBSET,
    group,
)


def main() -> int:
    cells = defaultdict(list)
    with open(sys.argv[1], newline="") as fh:
        for r in csv.DictReader(line for line in fh if not line.startswith("#")):
            if float(r["phi"]) != 0.6:
                continue
            cells[int(r["scenario"])].append(r)

    print("L21 (phi 0.6, groups B-D): median relative error of PD-hat, M6 against binomial")
    per = []
    for s, rs in sorted(cells.items()):
        pd, t, n = float(rs[0]["pd_true"]), int(rs[0]["periods"]), int(rs[0]["obligors"])
        if group(s, t) == "A":
            continue
        m6 = st.median(float(r["pd"]) for r in rs) / pd - 1
        bi = st.median(float(r["binom_pd"]) for r in rs) / pd - 1
        per.append({"s": s, "T": t, "n": n, "PD": pd, "m6": m6, "bi": bi, "closer": abs(m6) < abs(bi)})
    print(f"  {len(per)} scenarios, closer in {sum(p['closer'] for p in per)}")
    for key in ("T", "n", "PD"):
        print(f"  by {key}:")
        for v in sorted({p[key] for p in per}):
            sel = [p for p in per if p[key] == v]
            print(f"    {key} {v:g}: closer in {sum(p['closer'] for p in sel)} of {len(sel)}; median shift "
                  f"{st.median(p['m6'] - p['bi'] for p in sel):+.4f}; binomial's median relative error "
                  f"{st.median(p['bi'] for p in sel):+.4f}, M6's {st.median(p['m6'] for p in sel):+.4f}")
    print("  per scenario: s T n PD | binomial | M6 | shift | closer")
    for p in per:
        print(f"    {p['s']:2d} {p['T']:3d} {p['n']:5d} {p['PD']:.4f} | {p['bi']:+.4f} | {p['m6']:+.4f} | "
              f"{p['m6'] - p['bi']:+.4f} | {'yes' if p['closer'] else 'no'}")

    print("\nL28 (phi 0.6, subset, replicates 0-99): percentile against profile")
    print("  s T n PD | percentile cov | profile cov | gap | percentile misses: PD above / below | "
          "median midpoint/PD-hat - 1 | median midpoint/PD - 1 | median PD-hat/PD - 1")
    by_t = defaultdict(list)
    for s in SUBSET:
        rs = [r for r in cells[s] if int(r["replicate"]) < BOOT_REPLICATES]
        pd, t, n = float(rs[0]["pd_true"]), int(rs[0]["periods"]), int(rs[0]["obligors"])
        lo = [float(r["boot_lo_pd"]) for r in rs]
        hi = [float(r["boot_hi_pd"]) for r in rs]
        est = [float(r["pd"]) for r in rs]
        pct = sum(a <= pd <= b for a, b in zip(lo, hi)) / len(rs)
        prof = sum(float(r["prof_lo_pd"]) <= pd <= float(r["prof_hi_pd"]) for r in rs) / len(rs)
        above = sum(pd > b for b in hi)
        below = sum(pd < a for a in lo)
        mid = [(a + b) / 2 for a, b in zip(lo, hi)]
        by_t[t].append(pct - prof)
        print(f"    {s:2d} {t:3d} {n:5d} {pd:.4f} | {pct:.2f} | {prof:.2f} | {pct - prof:+.2f} | {above} / {below} | "
              f"{st.median(m / e - 1 for m, e in zip(mid, est)):+.4f} | {st.median(m / pd - 1 for m in mid):+.4f} | "
              f"{st.median(e / pd - 1 for e in est):+.4f}")
    print("  gap by T (mean over scenarios): " + "; ".join(f"T {t}: {st.fmean(v):+.3f} ({len(v)} scenarios)"
                                                        for t, v in sorted(by_t.items())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
