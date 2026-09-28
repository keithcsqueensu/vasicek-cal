# SPDX-License-Identifier: Apache-2.0
"""S-13 targeted: score the pre-registered predictions (PREDICTION.md, Q1-Q7) against the run.

Reads the R = 10,000 summary (summary_r10000.csv, written by recovery_harness --summary-out; same
columns as tests/golden/recovery/summary.csv) and, for Q7 only, the committed per-replicate table
fits_r10000.parquet (export_fits.py): the summary does not carry the direction of PD's and rho's
profile misses. Prints one Markdown row per prediction: predicted, result, held or not.
Committed before the run's results were seen.

    uv run --no-project --with pyarrow==25.0.1 python studies/recovery-r10000/compare.py SUMMARY_CSV FITS_PARQUET
"""

import csv
import math
import sys

import pyarrow.parquet as pq

R = 10_000
BAND = (0.95 - 3.29 * math.sqrt(0.95 * 0.05 / R), 0.95 + 3.29 * math.sqrt(0.95 * 0.05 / R))
FINDINGS = [(29, "pd"), (29, "rho"), (55, "rho"), (68, "pd"), (72, "rho"), (74, "rho")]
GROUP_B = [(29, "pd"), (74, "rho"), (55, "rho"), (29, "rho"), (73, "pd"), (64, "rho"), (37, "rho"), (74, "pd"),
           (47, "pd"), (55, "pd"), (65, "pd")]
A_LOWER = [(2, "rho"), (45, "pd"), (45, "rho"), (22, "rho"), (60, "pd"), (72, "pd"), (2, "pd"), (28, "pd")]
A_CONSERVATIVE = [(6, "pd"), (21, "pd"), (27, "rho"), (9, "pd"), (12, "pd"), (18, "pd"), (30, "rho")]
A_UPPER_PASS = [(1, "pd"), (7, "rho"), (27, "pd"), (15, "pd"), (24, "pd"), (54, "rho")]
ALL38 = (GROUP_B + [(68, "pd"), (41, "pd"), (31, "rho"), (49, "rho"), (59, "pd"), (72, "rho")] + A_LOWER
         + A_CONSERVATIVE + A_UPPER_PASS)
PDS, RHOS = (0.001, 0.01, 0.05), (0.02, 0.12, 0.24)


def coverage(rows, sid, p):
    return int(rows[sid][f"{p}_profile_covered"]) / R


def rho_misses(parquet_path, scenarios):
    """Per scenario, [rho profile intervals entirely below the true rho, entirely above]."""
    t = pq.read_table(parquet_path, columns=["scenario", "rho_lo", "rho_hi", "rho_profile_flags"]).to_pydict()
    out = {sid: [0, 0] for sid in scenarios}
    for sid, lo, hi, flags in zip(t["scenario"], t["rho_lo"], t["rho_hi"], t["rho_profile_flags"]):
        if sid not in out or flags & 4:  # not computed: not covering, but in neither direction
            continue
        truth = RHOS[(sid // 9) % 3]
        out[sid][0] += hi < truth
        out[sid][1] += lo > truth
    return out



def main():
    summary, fits = sys.argv[1], sys.argv[2]
    with open(summary, newline="") as f:
        rows = {int(r["scenario"]): r for r in csv.DictReader(line for line in f if not line.startswith("#"))}
    assert all(int(r["replicates"]) == R for r in rows.values())
    cov = {k: coverage(rows, *k) for k in ALL38}

    def below(k):
        return cov[k] < BAND[0]

    def above(k):
        return cov[k] > BAND[1]

    def name(k):
        return f"{k[0]}/{'PD' if k[1] == 'pd' else 'ρ'}"

    def listed(keys):
        return ", ".join(f"{name(k)} {cov[k]:.4f}" for k in keys)

    out = []

    def row(n, predicted, result, held):
        out.append(f"| {n} | {predicted} | {result} | {'held' if held else 'not held'} |")

    row("Q1", "The six findings below the band, coverage 0.915-0.942", listed(FINDINGS),
        all(below(k) and 0.915 <= cov[k] <= 0.942 for k in FINDINGS))
    row("Q2", "Every targeted group B verdict (11) below the band, coverage 0.930-0.942", listed(GROUP_B),
        all(below(k) and 0.930 <= cov[k] <= 0.942 for k in GROUP_B))
    q3 = (below((68, "pd")) and below((41, "pd")) and not below((31, "rho")) and not above((31, "rho"))
          and not below((49, "rho")) and not above((49, "rho")) and not below((59, "pd")))
    row("Q3", "68/PD and 41/PD below; 31/ρ and 49/ρ inside; 59/PD inside or above",
        listed([(68, "pd"), (41, "pd"), (31, "rho"), (49, "rho"), (59, "pd")]), q3)
    n_a = sum(below(k) for k in A_LOWER)
    row("Q4", "72/ρ below (0.915-0.938); at least 6 of the other 8 group A lower-edge verdicts below",
        f"72/ρ {cov[(72, 'rho')]:.4f}; {n_a} of 8 below ({listed(A_LOWER)})",
        below((72, "rho")) and 0.915 <= cov[(72, "rho")] <= 0.938 and n_a >= 6)
    n_up = sum(above(k) for k in A_UPPER_PASS)
    row("Q5", "All seven CONSERVATIVE above; at least 4 of the six upper-edge PASS verdicts above",
        f"CONSERVATIVE above: {sum(above(k) for k in A_CONSERVATIVE)} of 7; PASS above: {n_up} of 6 "
        f"({listed(A_UPPER_PASS)})", all(above(k) for k in A_CONSERVATIVE) and n_up >= 4)
    nb, na = sum(below(k) for k in ALL38), sum(above(k) for k in ALL38)
    npass = len(ALL38) - nb - na
    row("Q6", "Totals over the 38: below 19-25; above 10-16; PASS 1-6", f"below {nb}; above {na}; PASS {npass}",
        19 <= nb <= 25 and 10 <= na <= 16 and 1 <= npass <= 6)
    b_rho_below = [k for k in GROUP_B if k[1] == "rho" and below(k)]
    misses = rho_misses(fits, sorted({k[0] for k in b_rho_below}))
    fr = {k: misses[k[0]][0] / (R - int(rows[k[0]]["rho_profile_covered"])) for k in b_rho_below}
    row("Q7", "Group B ρ verdicts below the band: more than 2/3 of non-covering intervals entirely below the truth",
        "; ".join(f"{name(k)} {v:.2f}" for k, v in fr.items()) or "none below", bool(fr) and all(v > 2 / 3 for v in fr.values()))
    print(f"Band at R = {R}: {BAND[0]:.4f}-{BAND[1]:.4f}.\n")
    print("| # | Prediction | Result | Held |\n|---|---|---|---|")
    print("\n".join(out))


if __name__ == "__main__":
    main()
