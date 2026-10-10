# SPDX-License-Identifier: Apache-2.0
"""Project a full study run's runtime from a timing sample (perf/README.md, "Projecting a full run", D-187).

A study tool run with --replicates K --timing-out FILE records each (scenario, replicate) job's
thread-seconds by phase, with each phase's share of a full run's replicates. The projection:

  - groups the jobs by (T, n) cell, since runtime follows T x n far more than PD or rho;
  - a cell's cost per replicate is the sum over phases of share x the phase's mean seconds, the mean
    taken over the sampled jobs that run the phase in a full run (replicate < share x R_full); a phase
    with share "once" (a one-time cost, such as filling a surface-row cache) is counted once, in total;
  - the full run's thread-seconds are the sum over cells of that cost x the cell's scenarios x R_full,
    and its wall time that over the thread count (the study tools hand out jobs dynamically and have run
    at 99.7% utilisation, D-186);
  - the range is the 5th-95th percentile of 2,000 bootstrap draws, resampling jobs within each cell;
  - the subset-based projection (the subset's mean cost per job x every job), the method that missed S-10
    by 1.375x, is printed alongside for comparison.

Calibrated on S-10 (D-187): with 2 replicates per scenario (18 per cell) the projection was 6.11 h
(90% 5.68-6.54 h) against the actual 6.44 h; with 1 per scenario (9 per cell) 5.33 h (4.98-5.66 h), a
17% miss outside its range, because single-replicate costs are skewed (0.44-2.53x the scenario mean).
To be called validated after it lands in range on the next full run.

    python perf/project_runtime.py SAMPLE.csv [--replicates-full 1000] [--threads 24] [--subset 29,37,...]
    python perf/project_runtime.py --run "build/cpu-release/tests/study_parametric_bootstrap" \\
        [--start 2] [--target 0.10] [--max 32] [--threads 24] [--subset ...]

With --run the tool is run with --replicates K (K = --start, then doubled) and --timing-out, until the
range is within +/- target of the projection or K reaches --max; each run starts from replicate 0, so the
samples cost at most twice the last one. Standard library only.
"""

from __future__ import annotations

import argparse
import csv
import random
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

SUBSET_DEFAULT = "29,37,72,4,68,49,7,43,51"  # the study subset (D-150)


def read_sample(path):
    shares, tool, lines = {}, "", []
    with open(path, newline="") as f:
        for line in f:
            if line.startswith("# share,"):
                _, phase, share = line.rstrip("\n").split(",", 2)
                shares[phase] = share
            elif line.startswith("# tool,"):
                tool = line.rstrip("\n").split(",", 1)[1]
            elif not line.startswith("#"):
                lines.append(line)
    rows = list(csv.DictReader(lines))
    jobs = [{"scenario": int(r["scenario"]), "cell": (int(r["periods"]), int(r["obligors"])),
             "replicate": int(r["replicate"]), "sec": {p: float(r[p]) for p in shares}} for r in rows]
    return tool, shares, jobs


def cost_per_replicate(jobs, shares, r_full):
    """A cell's projected thread-seconds per replicate, from its sampled jobs."""
    total = 0.0
    for phase, share in shares.items():
        if share == "once":
            continue
        s = float(share)
        cut = s * r_full
        run = [j["sec"][phase] for j in jobs if j["replicate"] < cut]
        if run:
            total += s * sum(run) / len(run)
    return total


def project(jobs, shares, r_full, scenarios_per_cell, rng=None):
    once = sum(j["sec"][p] for j in jobs for p, s in shares.items() if s == "once")
    cells = {}
    for j in jobs:
        cells.setdefault(j["cell"], []).append(j)
    total = once
    for cell, cj in cells.items():
        if rng is not None:
            cj = [rng.choice(cj) for _ in cj]
        total += cost_per_replicate(cj, shares, r_full) * scenarios_per_cell[cell] * r_full
    return total


def summarise(jobs, shares, r_full, threads, subset, draws=2000):
    scen = {}
    for j in jobs:
        scen.setdefault(j["scenario"], j["cell"])
    per_cell = {}
    for cell in scen.values():
        per_cell[cell] = per_cell.get(cell, 0) + 1
    point = project(jobs, shares, r_full, per_cell)
    rng = random.Random(20261009)
    boot = sorted(project(jobs, shares, r_full, per_cell, rng) for _ in range(draws))
    lo, hi = boot[int(0.05 * draws)], boot[int(0.95 * draws)]
    sub_jobs = [j for j in jobs if j["scenario"] in subset]
    sub = None
    if sub_jobs:
        once = sum(j["sec"][p] for j in sub_jobs for p, s in shares.items() if s == "once")
        mean = cost_per_replicate(sub_jobs, shares, r_full)  # pooled over the subset's jobs
        sub = mean * len(scen) * r_full + once
    k = max(j["replicate"] for j in jobs) + 1
    return {"point": point, "lo": lo, "hi": hi, "half": max(point - lo, hi - point) / point, "subset": sub,
            "cells": len(per_cell), "scenarios": len(scen), "per_cell": len(jobs) / len(per_cell), "k": k,
            "threads": threads}


def report(tool, s):
    h = lambda x: x / s["threads"] / 3600
    print(f"{tool}: {s['scenarios']} scenarios in {s['cells']} (T, n) cells, {s['k']} replicates each "
          f"({s['per_cell']:.0f} jobs per cell)")
    print(f"  projection: {s['point']:,.0f} thread-s, {h(s['point']):.2f} h on {s['threads']} threads; "
          f"90% range {h(s['lo']):.2f}-{h(s['hi']):.2f} h (+/-{100 * s['half']:.0f}%)")
    if s["subset"] is not None:
        print(f"  subset-based, for comparison: {s['subset']:,.0f} thread-s, {h(s['subset']):.2f} h")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sample", nargs="?")
    ap.add_argument("--run", help="the study tool command; run with --replicates K --timing-out")
    ap.add_argument("--start", type=int, default=2)
    ap.add_argument("--max", type=int, default=32)
    ap.add_argument("--target", type=float, default=0.10)
    ap.add_argument("--replicates-full", type=int, default=1000)
    ap.add_argument("--threads", type=int, default=24)
    ap.add_argument("--subset", default=SUBSET_DEFAULT)
    a = ap.parse_args()
    subset = {int(x) for x in a.subset.split(",") if x}
    if a.sample:
        tool, shares, jobs = read_sample(a.sample)
        s = summarise(jobs, shares, a.replicates_full, a.threads, subset)
        report(tool, s)
        return 0 if s["half"] <= a.target else 3
    if not a.run:
        ap.error("give a SAMPLE file or --run")
    k = a.start
    with tempfile.TemporaryDirectory() as d:
        while True:
            out = Path(d) / f"timing_{k}.csv"
            cmd = shlex.split(a.run, posix=False) + ["--replicates", str(k), "--timing-out", str(out)]
            cmd[0] = str(Path(cmd[0]).resolve())  # Windows' CreateProcess does not resolve a relative a/b path
            print("running:", " ".join(cmd), file=sys.stderr)
            done = subprocess.run(cmd, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
            if done.returncode != 0:
                print(done.stderr, file=sys.stderr)
                return done.returncode
            tool, shares, jobs = read_sample(out)
            s = summarise(jobs, shares, a.replicates_full, a.threads, subset)
            report(tool, s)
            if s["half"] <= a.target:
                return 0
            if 2 * k > a.max:
                print(f"  range still wider than +/-{100 * a.target:.0f}% at {k} replicates per scenario: "
                      "report the projection with its range", file=sys.stderr)
                return 3
            k *= 2


if __name__ == "__main__":
    sys.exit(main())
