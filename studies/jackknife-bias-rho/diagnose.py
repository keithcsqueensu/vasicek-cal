# SPDX-License-Identifier: Apache-2.0
"""The shared jackknife run: post-hoc diagnostics of the q Wald sub-arms (J15-J18). Written after the
results, to diagnose them; it scores nothing (compare.py does, and was committed before the run).

For each of S-23's 17 below-band Wald scenarios, over the replicates without a Wald flag (the set the
sub-arms are scored on), from the committed per-replicate Parquet. e is the error of the corrected
centre, s-tilde - logit(q), and sd its standard deviation. For each SE, coverage of e within
+/- 1.96 x SE, adding one defect at a time:
  oracle   SE = sd in every replicate: steady, the right size, unrelated to e;
  size     SE = the SE's root mean square in every replicate: its average size, still steady;
  noise    the SE's own values, shuffled across replicates (mean of 200 shuffles): its size and its
           replicate-to-replicate variation, but unrelated to e;
  actual   the SE's own values: adds its link to e. This is the sub-arm's scored coverage.
Also the skewness of e, and the correlation of e with each SE. The oracle interval is symmetric,
so where it covers, skewness of e is not what the Wald intervals lack.

    uv run --no-project --with pyarrow==25.0.1 python studies/jackknife-bias-rho/diagnose.py
"""

import math
import random
from collections import defaultdict
from pathlib import Path

import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
WALD_BELOW = [11, 29, 32, 35, 37, 38, 47, 49, 55, 56, 59, 64, 65, 67, 68, 73, 74]  # S-23, D-156
Z = 1.96
SHUFFLES = 200


def mean(x):
    return sum(x) / len(x)


def corr(x, y):
    mx, my = mean(x), mean(y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    return sxy / math.sqrt(sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y))


def coverage(err, se):
    return mean([abs(e) <= Z * s for e, s in zip(err, se)])


def steps(err, se, sd):
    n = len(err)
    rms = math.sqrt(mean([s * s for s in se]))
    shuffled = []
    for k in range(SHUFFLES):
        perm = se[:]
        random.Random(k).shuffle(perm)
        shuffled.append(coverage(err, perm))
    return coverage(err, [sd] * n), coverage(err, [rms] * n), mean(shuffled), coverage(err, se)


def main():
    t = pq.read_table(HERE / "replicates.parquet",
                      columns=["scenario", "q_err_a", "se_delta", "se_jack", "qa_cover"]).to_pydict()
    rows = defaultdict(list)
    for i, s in enumerate(t["scenario"]):
        if t["qa_cover"][i] >= 0:
            rows[s].append(i)
    print("| Scenario | skew(e) | corr(e, SE_Δ) | corr(e, SE_J) | oracle | SE_Δ: size, noise, actual (a) "
          "| SE_J: size, noise, actual (b) |")
    print("|---|---|---|---|---|---|---|")
    for s in WALD_BELOW:
        err = [t["q_err_a"][i] for i in rows[s]]
        sed = [t["se_delta"][i] for i in rows[s]]
        sej = [t["se_jack"][i] for i in rows[s]]
        m = mean(err)
        sd = math.sqrt(mean([(e - m) ** 2 for e in err]))
        skew = mean([((e - m) / sd) ** 3 for e in err])
        oracle, *d = steps(err, sed, sd)
        _, *j = steps(err, sej, sd)
        print(f"| {s} | {skew:+.2f} | {corr(err, sed):.2f} | {corr(err, sej):.2f} | {oracle:.3f} | "
              f"{d[0]:.3f}, {d[1]:.3f}, {d[2]:.3f} | {j[0]:.3f}, {j[1]:.3f}, {j[2]:.3f} |")


if __name__ == "__main__":
    main()
