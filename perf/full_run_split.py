# SPDX-License-Identifier: Apache-2.0
"""Phase split of two full study runs, raw and normalised for machine state (perf/p24_p27_full_run.md, D-414).

Reads two --timing-out files of the same tool over the same scenarios and replicates (A, the baseline, then B)
and prints, for each phase, its thread-seconds and share in each run and the change from A to B.

Two runs taken one after the other on one machine can differ in speed by more than the change being measured:
on 2026-11-11 the phases whose code did not change ran 1-18% faster, cell by cell, in the second run. So it
also prints a normalised split. Within each (T, n) cell, the machine factor is B's thread-seconds over A's,
summed over the control phases (by default every phase whose code the change leaves alone), and B's cost in
that cell is divided by it before the cells are summed. This assumes the machine's speed moves every phase of
a cell in the same proportion; the control phases' own normalised changes (near zero, by construction, only
in total) are printed so a reader can see how far that holds phase by phase.

--log FILE --cut S restricts both runs to the scenarios B's log reports done within S seconds of its start
(lines "scenario N done (k left), S s"), to drop the part of B's run that shared the machine with other work.

    python perf/full_run_split.py A_timing.csv B_timing.csv [--control simulate,profile,...]
        [--log B_stderr.txt --cut 2981]
Standard library only.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict

CONTROL_DEFAULT = "simulate,profile,W0,Jeffreys,S-4 widths"  # unchanged by P-24 to P-27 (D-407 to D-410)
KEYS = ("scenario", "periods", "obligors", "replicate")


def read(path, keep):
    with open(path, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    phases = [k for k in rows[0] if k not in KEYS]
    cells = defaultdict(lambda: defaultdict(float))
    jobs = 0
    for r in rows:
        if keep is not None and int(r["scenario"]) not in keep:
            continue
        jobs += 1
        c = (int(r["periods"]), int(r["obligors"]))
        for p in phases:
            cells[c][p] += float(r[p])
    return phases, cells, jobs


def done_by(log, cut):
    keep = set()
    with open(log, encoding="utf-8", errors="replace") as f:
        for line in f:
            w = line.split()
            if len(w) >= 4 and w[0] == "scenario" and w[2] == "done" and int(w[-2]) <= cut:
                keep.add(int(w[1]))
    return keep


def table(title, phases, a, b):
    ta, tb = sum(a.values()), sum(b.values())
    print(title)
    for p in [*phases, "TOTAL"]:
        x, y = (ta, tb) if p == "TOTAL" else (a[p], b[p])
        change = f"{100 * (y / x - 1):+7.1f}%" if x else ""
        print(f"  {p:<18}{x:>11,.0f} ({100 * x / ta:4.1f}%) {y:>11,.0f} ({100 * y / tb:4.1f}%) {change}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--control", default=CONTROL_DEFAULT)
    ap.add_argument("--log")
    ap.add_argument("--cut", type=int)
    args = ap.parse_args()
    keep = done_by(args.log, args.cut) if args.log and args.cut is not None else None
    phases, ca, na = read(args.a, keep)
    phases_b, cb, nb = read(args.b, keep)
    if phases != phases_b or na != nb or set(ca) != set(cb):
        raise SystemExit("the two runs differ in phases, jobs or cells")
    control = args.control.split(",")
    for p in control:
        if p not in phases:
            raise SystemExit(f"unknown control phase {p!r}")
    print(f"{na} jobs per run in {len(ca)} cells" + (f" ({len(keep)} scenarios done by {args.cut} s)" if keep else ""))

    raw_a, raw_b = defaultdict(float), defaultdict(float)
    norm_b = defaultdict(float)
    print("machine factor by (T, n) cell (B / A over the control phases):")
    for c in sorted(ca):
        f = sum(cb[c][p] for p in control) / sum(ca[c][p] for p in control)
        print(f"  T={c[0]:<4} n={c[1]:<6} {f:.3f}")
        for p in phases:
            raw_a[p] += ca[c][p]
            raw_b[p] += cb[c][p]
            norm_b[p] += cb[c][p] / f
    table("raw thread-seconds: A, B, change", phases, raw_a, raw_b)
    table("normalised (B divided by its cell's machine factor): A, B, change", phases, raw_a, norm_b)


if __name__ == "__main__":
    main()
