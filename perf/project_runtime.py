# SPDX-License-Identifier: Apache-2.0
"""Project a full study run's runtime from a timing sample (perf/README.md, "Projecting a full run", D-187).

A study tool run with --replicates K --timing-out FILE records each (scenario, replicate) job's
thread-seconds by phase, with each phase's share of a full run's replicates. The projection:

  - groups the jobs by (T, n) cell, since runtime follows T x n far more than PD or rho;
  - a cell's cost per replicate is the sum over phases of share x the phase's mean seconds, the mean
    taken over the sampled jobs that run the phase in a full run (replicate < share x R_full); a phase
    with share "once" (a one-time cost, such as filling a surface-row cache) is counted once, in total;
  - the full run's thread-seconds are the sum over cells of that cost x the cell's scenarios x R_full,
    and its wall time that over (threads x efficiency). The efficiency is measured from real runs, not
    assumed: thread-seconds / (threads x wall) of the last full run, which covers scenario overlap, cache
    fills and the tail where threads idle. S-10/S-4: 555,004 / (24 x 23,190) = 0.997, the default;
    each validation re-measures it;
  - the range is the 5th-95th percentile of 2,000 bootstrap draws, resampling jobs within each cell;
  - the subset-based projection (the subset's mean cost per job x every job), the method that missed S-10
    by 1.375x, is printed alongside for comparison.

Calibrated on S-10 (D-187): with 2 replicates per scenario (18 per cell) the projection was 6.13 h
(90% 5.70-6.56 h) against the actual 6.44 h; with 1 per scenario (9 per cell) 5.34 h (5.00-5.68 h), a
17% miss outside its range, because single-replicate costs are skewed (0.44-2.53x the scenario mean).
To be called validated after it lands in range on the next full run.

    python perf/project_runtime.py SAMPLE.csv [--replicates-full 1000] [--threads 24] [--efficiency 0.997]
        [--subset 29,37,...]
    python perf/project_runtime.py --run "build/cpu-release/tests/study_parametric_bootstrap" \\
        [--start 2] [--target 0.10] [--max 256] [--warmup 120] [--min-wall 120] [--threads 24] [--subset ...]

With --run the tool is run with --replicates K (K = --start, then doubled) and --timing-out. It stops when
all three hold, or when K reaches --max:
  - the range is within +/- target of the projection;
  - the projection agrees with the previous K's within the same target (D-192);
  - the sample itself ran at least --min-wall seconds (default 120; D-193).
Before sampling, --warmup seconds (default 120) of discarded runs bring the machine to sustained load. Each
run starts from replicate 0, so the samples cost at most twice the last one.

Two biases the range cannot see, both found on S-35 (0.293 h actual, 25,276 thread-s):
  - a one-time cost charged per replicate. The 2-replicate sample charged the surface-row cache's fills to
    its few jobs and projected 0.40 h; tools now record fills as a "once" phase;
  - boost clocks in a short sample. After that fix, cold samples of 2-32 replicates projected 0.20-0.25 h,
    a 64-replicate sample (67 s) 0.27 h, and a 128-replicate sample (145 s) 0.30 h.
Standard library only.
"""

from __future__ import annotations

import argparse
import csv
import random
import shlex
import subprocess
import sys
import tempfile
import time
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


def summarise(jobs, shares, r_full, threads, efficiency, subset, draws=2000):
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
            "threads": threads, "efficiency": efficiency}


def report(tool, s, allowance):
    h = lambda x: x / (s["threads"] * s["efficiency"]) / 3600
    print(f"{tool}: {s['scenarios']} scenarios in {s['cells']} (T, n) cells, {s['k']} replicates each "
          f"({s['per_cell']:.0f} jobs per cell)")
    print(f"  projection: {s['point']:,.0f} thread-s, {h(s['point']):.2f} h on {s['threads']} threads at efficiency {s['efficiency']}; "
          f"90% range {h(s['lo']):.2f}-{h(s['hi']):.2f} h (+/-{100 * s['half']:.0f}%)")
    print(f"  with the machine-state allowance (+/-{100 * allowance:.0f}%, D-193): "
          f"{h(s['lo'] * (1 - allowance)):.2f}-{h(s['hi'] * (1 + allowance)):.2f} h")
    if s["subset"] is not None:
        print(f"  subset-based, for comparison: {s['subset']:,.0f} thread-s, {h(s['subset']):.2f} h")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sample", nargs="?")
    ap.add_argument("--run", help="the study tool command; run with --replicates K --timing-out")
    ap.add_argument("--start", type=int, default=2)
    ap.add_argument("--max", type=int, default=256)
    ap.add_argument("--target", type=float, default=0.10)
    ap.add_argument("--replicates-full", type=int, default=1000)
    ap.add_argument("--threads", type=int, default=24)
    ap.add_argument("--efficiency", type=float, default=0.997,
                    help="thread-seconds / (threads x wall) measured on the last full run (S-10/S-4: 0.997)")
    ap.add_argument("--subset", default=SUBSET_DEFAULT)
    ap.add_argument("--warmup", type=float, default=120.0,
                    help="seconds of discarded runs before sampling, so the sample sees sustained clocks (D-193)")
    ap.add_argument("--allowance", type=float, default=0.05,
                    help="machine-state spread added to the sampling range (S-35: two 128-replicate samples on "
                         "different occasions, 24,059 and 25,990 thread-s, about +/-4%%; D-193)")
    ap.add_argument("--min-wall", type=float, default=120.0,
                    help="the accepted sample must itself run this long: short runs stay at boost clocks (D-193)")
    a = ap.parse_args()
    subset = {int(x) for x in a.subset.split(",") if x}
    if a.sample:
        tool, shares, jobs = read_sample(a.sample)
        s = summarise(jobs, shares, a.replicates_full, a.threads, a.efficiency, subset)
        report(tool, s, a.allowance)
        return 0 if s["half"] <= a.target else 3
    if not a.run:
        ap.error("give a SAMPLE file or --run")
    k = a.start
    previous = None  # the last sample's projection: a bias check the range cannot give (D-192)
    with tempfile.TemporaryDirectory() as d:
        # Warm-up (D-193): a cold machine runs a short sample at boost clocks. S-35's 32-replicate sample cost
        # 0.215 s per job cold and 0.279 s after sustained load, against 0.275 s in the full run.
        t_end = time.monotonic() + a.warmup
        while time.monotonic() < t_end:
            cmd = shlex.split(a.run, posix=False) + ["--replicates", str(a.start)]
            cmd[0] = str(Path(cmd[0]).resolve())
            print("warming up:", " ".join(cmd), file=sys.stderr)
            if subprocess.run(cmd, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode != 0:
                break  # the sampling run reports the error
        while True:
            out = Path(d) / f"timing_{k}.csv"
            cmd = shlex.split(a.run, posix=False) + ["--replicates", str(k), "--timing-out", str(out)]
            cmd[0] = str(Path(cmd[0]).resolve())  # Windows' CreateProcess does not resolve a relative a/b path
            print("running:", " ".join(cmd), file=sys.stderr)
            t0 = time.monotonic()
            done = subprocess.run(cmd, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
            wall = time.monotonic() - t0
            if done.returncode != 0:
                print(done.stderr, file=sys.stderr)
                return done.returncode
            tool, shares, jobs = read_sample(out)
            s = summarise(jobs, shares, a.replicates_full, a.threads, a.efficiency, subset)
            report(tool, s, a.allowance)
            drift = abs(s["point"] - previous) / s["point"] if previous is not None else None
            if drift is not None:
                print(f"  against {k // 2} replicates per scenario: {100 * drift:.0f}% apart")
            print(f"  sample wall {wall:.0f} s (minimum {a.min_wall:.0f} s)")
            if s["half"] <= a.target and drift is not None and drift <= a.target and wall >= a.min_wall:
                return 0
            if 2 * k > a.max:
                print(f"  range, agreement or sample length not met at {k} replicates per scenario: "
                      "report the projection with its range", file=sys.stderr)
                return 3
            previous = s["point"]
            k *= 2


if __name__ == "__main__":
    sys.exit(main())
