# SPDX-License-Identifier: Apache-2.0
"""M6 validation, post-hoc diagnosis of W3's miss (D-202). Written after the run; not a prediction; changes no verdict.

W3 standardised the paired difference theta3 - theta2 at phi = 0 by the Hausman SD sqrt(v3 - v2), taken from the
Monte Carlo expected information (500 panels per cell). Three checks, per credible T = 100 cell at phi = 0:

1. The first-order Hausman term: theta3 - theta2 ~= beta * phi_hat with beta = -(I_tt^-1 I_tphi). Compared: the
   regression slope of the difference on phi_hat against beta, its R^2, and phi_hat's SD against the information's.
2. The information's cross terms: corr(s_theta, s_phi) = i_theta_phi / sqrt(i_theta_theta i_phi_phi), scaled by
   sqrt(500), across cells, at phi = 0 (zero in theory, below) and at 0.3 and 0.6 for contrast.
3. Second order: corr(diff, phi_hat^2) and corr(|diff|, |phi_hat|).

Theory: at phi = 0 the periods are independent, so the phi score is a sum of products of adjacent periods' posterior
factor means, each mean zero; its covariance with a sum of per-period PD or rho scores, and the expected derivative,
vanish. So I_theta_phi = 0 and v3 - v2 = 0 at phi = 0; the difference is second order.

    python studies/m6-validation/posthoc_w3.py rows.csv
"""

from __future__ import annotations

import math
import statistics
import sys

import score  # the registered scorer's readers: the same cells, SEs and information

PANELS = 500


def main() -> int:
    rows = score.read_csv(sys.argv[1])
    info = score.information(str(score.INFO))
    raw = {(int(r["scenario"]), float(r["phi"])): r for r in score.read_csv(str(score.INFO))}
    by_pair: dict[tuple[int, float], list[dict[str, str]]] = {}
    for r in rows:
        by_pair.setdefault((int(r["scenario"]), float(r["phi"])), []).append(r)
    cells = sorted(k for k, v in info.items() if v["credible"] and v["T"] == 100 and k[1] == 0.0)

    print("1. first order, per credible phi = 0 cell (difference in units of the information's SE)")
    print("  cell par | slope / beta = ratio | SD(phi_hat)/se_phi | R^2 | corr(d, phi_hat^2) | corr(|d|, |phi_hat|)")
    second = {0: [], 1: []}
    for k in cells:
        r0 = raw[k]
        ipp, irr = float(r0["i_pd_pd"]), float(r0["i_rho_rho"])
        ipr, ipf, irf = float(r0["i_pd_rho"]), float(r0["i_pd_phi"]), float(r0["i_rho_phi"])
        det2 = ipp * irr - ipr * ipr
        beta = (-(irr * ipf - ipr * irf) / det2, -(ipp * irf - ipr * ipf) / det2)
        ph = [float(r["phi_hat"]) for r in by_pair[k]]
        sd_ph = statistics.pstdev(ph) / info[k]["se"][2]
        for ai, (est, ref, name) in enumerate((("pd", "phi0_pd", "PD "), ("rho", "phi0_rho", "rho"))):
            se = info[k]["se"][ai]
            d = [(float(r[est]) - float(r[ref])) / se for r in by_pair[k]]
            fit = statistics.linear_regression(ph, d)
            res = [y - fit.intercept - fit.slope * x for x, y in zip(ph, d)]
            r2 = 1 - statistics.pvariance(res) / statistics.pvariance(d)
            b = beta[ai] / se
            c2 = statistics.correlation([x * x for x in ph], d)
            ca = statistics.correlation([abs(x) for x in ph], [abs(y) for y in d])
            second[ai].append((c2, ca))
            print(f"  {k[0]:3d} {name} | {fit.slope:+.4f} / {b:+.4f} = {fit.slope / b:+6.2f} | {sd_ph:4.2f} | {r2:4.2f} |"
                  f" {c2:+.2f} | {ca:+.2f}")

    print(f"\n2. the information's cross-correlations, credible T = 100 cells; z = corr * sqrt({PANELS}), SD 1 if noise")
    for phi in (0.0, 0.3, 0.6):
        sel = [raw[k] for k, v in info.items() if v["credible"] and v["T"] == 100 and k[1] == phi]
        for a in ("pd", "rho"):
            c = [float(r[f"i_{a}_phi"]) / math.sqrt(float(r[f"i_{a}_{a}"]) * float(r["i_phi_phi"])) for r in sel]
            z = [x * math.sqrt(PANELS) for x in c]
            print(f"  phi {phi}  {a:3s}: {len(c):2d} cells, corr mean {statistics.fmean(c):+.4f}, z SD "
                  f"{statistics.pstdev(z):.2f}, |z| > 2 in {sum(abs(x) > 2 for x in z)}")

    print("\n3. second order, means over the phi = 0 cells")
    for ai, name in ((0, "PD "), (1, "rho")):
        print(f"  {name}: corr(d, phi_hat^2) {statistics.fmean(c for c, _ in second[ai]):+.2f}, "
              f"corr(|d|, |phi_hat|) {statistics.fmean(c for _, c in second[ai]):+.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
