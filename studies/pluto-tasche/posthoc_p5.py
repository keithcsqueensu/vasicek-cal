# SPDX-License-Identifier: Apache-2.0
"""S-6 post-hoc: how far outside its band each P5 (i) miss lies. Written after the run; changes no verdict (D-217).

P5 (i) was scored NOT HELD by score.py (272 of 405 cells in band, against >= 403). The band ends ratio_lo and
ratio_hi are reference.py's bounds (scipy brentq) at the totals d_lo and d_hi; the simulated medians are the tool's
bounds (bisection to 1e-10 in logit, TOL_PT_SOLVE_LOGIT). score.py allows 1e-12 of slack, narrower than the two
solvers' agreement (TOL_SCIPY_PT_LOGIT = 2e-10). This script reports, for every cell out of band as scored:

  - its relative distance outside the band;
  - whether the totals behind its simulated median (the middle two of R = 1,000) lie in [d_lo, d_hi], the band in
    totals, where no floating point enters;

and the same two readings for P5 (iii) and (iv)'s medians over scenarios.

    python studies/pluto-tasche/posthoc_p5.py ROWS.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import score as S


def main(argv: list[str]) -> int:
    ref = S.read_reference()
    run = S.Run(Path(argv[1]))
    sc = S.scenarios(ref)
    med = {}
    out = []
    for s in sorted(sc):
        for g in S.GAMMAS:
            cell = run.cell(s, g)
            ds = sorted(r["d_total"] for r in cell)
            mid = (ds[len(ds) // 2 - 1], ds[len(ds) // 2])
            for m in S.METHODS:
                v = med[(s, m, g)] = S.median_ratio(cell, m)
                rr = ref[(s, m, g)]
                if not S.in_range(v, rr["ratio_lo"], rr["ratio_hi"]):
                    rel = (rr["ratio_lo"] - v) / rr["ratio_lo"] if v < rr["ratio_lo"] else (v - rr["ratio_hi"]) / rr["ratio_hi"]
                    in_d = rr["d_lo"] <= mid[0] and mid[1] <= rr["d_hi"]
                    out.append((s, m, g, rel, mid, int(rr["d_lo"]), int(rr["d_hi"]), in_d))
    print("## P5 (i): cells out of band as scored (score.py, slack 1e-12)")
    print()
    print(f"- out of band as scored: {len(out)} of 405")
    print(f"- relative distance outside the band: largest {max(abs(o[3]) for o in out):.2e}, "
          f"smallest {min(abs(o[3]) for o in out):.2e}")
    print(f"- of these, with the middle two totals outside [d_lo, d_hi]: {sum(not o[7] for o in out)}")
    print(f"- of these, with a band of zero width in totals (d_lo = d_hi): {sum(o[5] == o[6] for o in out)}")
    print()
    print("| scenario | method | gamma | relative distance | middle two D | [d_lo, d_hi] |")
    print("|---|---|---|---|---|---|")
    for s, m, g, rel, mid, lo, hi, _ in out:
        print(f"| {s} | {m} | {g} | {rel:+.2e} | {mid[0]}, {mid[1]} | [{lo}, {hi}] |")
    print()
    print("## P5 (iii) and (iv): the simulated median over scenarios at gamma 0.95 against its band")
    print()
    print("| n PD | method | rho | simulated | band | relative distance outside |")
    print("|---|---|---|---|---|---|")
    for tier, keep in (("1", S.tier1), ("0.1", lambda t: not S.tier1(t))):
        for m in ("oracle", "independent"):
            for rho in S.RHOS:
                ids = [s for s in sc if keep(sc[s]) and sc[s][1] == rho]
                sim = S.median(med[(s, m, 0.95)] for s in ids)
                lo = S.median(ref[(s, m, 0.95)]["ratio_lo"] for s in ids)
                hi = S.median(ref[(s, m, 0.95)]["ratio_hi"] for s in ids)
                dist = 0.0 if lo <= sim <= hi else ((lo - sim) / lo if sim < lo else (sim - hi) / hi)
                print(f"| {tier} | {m} | {rho} | {sim:.6f} | [{lo:.6f}, {hi:.6f}] | {dist:+.2e} |")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
