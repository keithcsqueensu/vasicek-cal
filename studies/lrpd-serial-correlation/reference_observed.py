# SPDX-License-Identifier: Apache-2.0
"""S-35's observed-rate arm (appended registration, 2026-10-09): reference numbers from pinned results only.

The arm fits the Vasicek-rate MLE to the observed rates d_t / n with S-28's censored treatment of
zero-default periods (detection limit 1/(2n)), as the S-8 pass fitted it. At phi = 0 those fits are
pinned, replicate by replicate, in studies/mle-vs-mom/fits.parquet (the censor_* columns; on panels without
a zero- or full-default period the censored and refusing likelihoods are the same, and the pass fitted
those once). Unlike the n -> infinity arm, this estimator is biased at phi = 0 (S-28: binomial noise read
as factor variance, and censoring), so the mapping carries each scenario's own bias forward:

  - per scenario at phi = 0, in logit: the bias b0 = mean(u-hat) - u_true and spread s0 = sd(u-hat) of the
    censored fits, and the profile coverage c0;
  - the interval's half-width w is the one that gives coverage c0 for a normal estimate with bias b0 and
    sd s0: Phi((w - b0)/s0) - Phi((-w - b0)/s0) = c0;
  - at phi the spread inflates to s0 sqrt(f), with f the binomial MLE's inflation for the same parameter
    (reference.csv's f_pd, f_rho; the censored rate likelihood's per-period score depends on Z_t much as
    the count likelihood's does), the width w is unchanged (it comes from the iid curvature), and rho's
    bias moves down by (f_pd - 1)/T in logit, the extra downward bias of a variance about a mean of T
    dependent periods (reference.py's rate section; f_pd standing in for f_mu where binomial noise dilutes
    the dependence). PD's bias is unchanged: the marginal, and so the pseudo-true value, does not move.

    uv run --no-project --with scipy==1.16.2 --with pyarrow==25.0.1 python studies/lrpd-serial-correlation/reference_observed.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import pyarrow.parquet as pq
from scipy.optimize import brentq
from scipy.special import ndtr

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]


def logit(v: float) -> float:
    return math.log(v) - math.log1p(-v)


def coverage(w: float, b: float, s: float) -> float:
    return float(ndtr((w - b) / s) - ndtr((-w - b) / s))


def main() -> int:
    fits = pq.read_table(REPO / "studies/mle-vs-mom/fits.parquet",
                         columns=["scenario", "pd_true", "rho_true", "periods", "obligors", "censor_pd", "censor_rho",
                                  "censor_cover_pd", "censor_cover_rho"]).to_pylist()
    by = {}
    for r in fits:
        if r["scenario"] < 81:
            by.setdefault(r["scenario"], []).append(r)
    with open(HERE / "reference.csv", newline="") as fh:
        infl = {(int(r["scenario"]), float(r["phi"])): r for r in csv.DictReader(fh)}
    out = []
    for sid in sorted(by):
        rows = by[sid]
        t = rows[0]["periods"]
        stats = {}
        for name in ("pd", "rho"):
            truth = logit(rows[0][f"{name}_true"])
            u = [logit(r[f"censor_{name}"]) for r in rows]
            m = sum(u) / len(u)
            s0 = math.sqrt(sum((x - m) ** 2 for x in u) / (len(u) - 1))
            c0 = sum(r[f"censor_cover_{name}"] == 1 for r in rows) / len(rows)
            b0 = m - truth
            c = min(max(c0, 1e-6), 1 - 1e-6)
            w = brentq(lambda w_, b0=b0, s0=s0, c=c: coverage(w_, b0, s0) - c, 1e-9, 100 * s0 + abs(b0) + 10)
            stats[name] = (b0, s0, c0, w)
        for phi in (0.3, 0.6):
            ref = infl[(sid, phi)]
            f_pd, f_rho = float(ref["f_pd"]), float(ref["f_rho"])
            row = {"scenario": sid, "pd": rows[0]["pd_true"], "rho": rows[0]["rho_true"], "periods": t,
                   "obligors": rows[0]["obligors"], "phi": phi}
            for name, f, shift in (("pd", f_pd, 0.0), ("rho", f_rho, -(f_pd - 1) / t)):
                b0, s0, c0, w = stats[name]
                row[f"{name}_bias0_u"] = b0
                row[f"{name}_sd0_u"] = s0
                row[f"{name}_c0"] = c0
                row[f"{name}_bias_u_pred"] = b0 + shift
                row[f"{name}_cov_pred"] = coverage(w, b0 + shift, s0 * math.sqrt(f))
            out.append(row)
    path = HERE / "reference_observed.csv"
    with open(path, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(out[0]), lineterminator="\n")
        wr.writeheader()
        for row in out:
            wr.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in row.items()})
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
