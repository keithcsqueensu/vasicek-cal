# SPDX-License-Identifier: Apache-2.0
"""P-13 (D-403): build study tools with profile-guided (and optionally link-time) optimisation, then check that
their output is byte-identical to a plain Release build's (perf/p13_lto_pgo.md).

The cycle, in one build directory:
  1. configure with -DVCAL_PGO=GENERATE (and -DVCAL_LTO=ON unless --no-lto), build the targets;
  2. run the training commands (each tool on scenarios the check does not use), which write the profiles;
  3. Clang only: merge the raw profiles into default.profdata (llvm-profdata merge);
  4. reconfigure with -DVCAL_PGO=USE and rebuild: the optimised binaries;
  5. with --check BASELINE_TESTS_DIR: run each check command with both builds and compare the outputs byte for
     byte; any difference is a failure (exit 1).

Training and check commands are fixed below (the recovery subset's scenarios, split so the check is not trained
on). Run from the repository root, in the shell you build from (on Windows: a developer prompt, so the
instrumented binaries find pgort140.dll):

    python perf/p13_pgo.py --build-dir build/pgo [--generator Ninja] [--no-lto] [--check build/cpu-release/tests]

Standard library only.
"""

from __future__ import annotations

import argparse
import filecmp
import glob
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

TARGETS = ["study_parametric_bootstrap", "study_lrpd_serial_m6"]
# (tool, arguments): scenarios 29, 72 and 49 of the subset for training; 37 and 4 for the check.
TRAIN = [
    ("study_parametric_bootstrap", ["--replicates", "1", "--scenarios", "29,72", "--threads", "1"]),
    ("study_lrpd_serial_m6", ["--replicates", "1", "--scenarios", "49", "--phis", "0.3", "--threads", "1"]),
]
CHECK = [
    ("study_parametric_bootstrap", ["--replicates", "1", "--scenarios", "37,4", "--threads", "1"]),
    ("study_lrpd_serial_m6", ["--replicates", "1", "--scenarios", "37", "--phis", "0.6", "--threads", "1"]),
]


def run(cmd: list[str], **kw) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, **kw)


def exe(tests_dir: Path, tool: str) -> str:
    for name in (tool, tool + ".exe"):
        p = tests_dir / name
        if p.exists():
            return str(p)
    sys.exit(f"{tool} not found in {tests_dir}")


def configure(build: Path, generator: str | None, pgo: str, lto: bool) -> None:
    cmd = ["cmake", "-S", ".", "-B", str(build), "-DCMAKE_BUILD_TYPE=Release", "-DVCAL_ENABLE_CUDA=OFF",
           f"-DVCAL_PGO={pgo}", f"-DVCAL_LTO={'ON' if lto else 'OFF'}", f"-DVCAL_PGO_DIR={build / 'pgo'}"]
    if generator:
        cmd += ["-G", generator]
    run(cmd)


def compiler_id(build: Path) -> str:
    for line in (build / "CMakeCache.txt").read_text(errors="replace").splitlines():
        if line.startswith("CMAKE_CXX_COMPILER_ID:"):
            return line.split("=", 1)[1]
    return ""


def timed(cmd: list[str], cwd: Path) -> float:
    t0 = time.perf_counter()
    run(cmd, cwd=cwd, stdout=subprocess.DEVNULL)
    return time.perf_counter() - t0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-dir", required=True)
    ap.add_argument("--generator", default=None)
    ap.add_argument("--no-lto", action="store_true")
    ap.add_argument("--check", metavar="BASELINE_TESTS_DIR", help="a plain Release build's tests directory")
    a = ap.parse_args()
    build = Path(a.build_dir).resolve()
    lto = not a.no_lto
    tests = build / "tests"
    work = build / "p13_runs"
    work.mkdir(parents=True, exist_ok=True)

    configure(build, a.generator, "GENERATE", lto)
    shutil.rmtree(build / "pgo", ignore_errors=True)
    (build / "pgo").mkdir()
    run(["cmake", "--build", str(build), "--target", *TARGETS])
    for tool, args in TRAIN:
        run([exe(tests, tool), *args, "--out", str(work / f"train_{tool}.csv")], cwd=work, stdout=subprocess.DEVNULL)
    if compiler_id(build) == "Clang":
        raw = glob.glob(str(build / "pgo" / "*.profraw"))
        if not raw:
            sys.exit("no .profraw files: did the training runs write their profiles?")
        run(["llvm-profdata", "merge", "-output", str(build / "pgo" / "default.profdata"), *raw])
    configure(build, a.generator, "USE", lto)
    run(["cmake", "--build", str(build), "--target", *TARGETS])

    if a.check:
        base = Path(a.check).resolve()
        bad = 0
        for tool, args in CHECK:
            outs = {}
            for tag, d in (("base", base), ("pgo", tests)):
                out = work / f"check_{tool}_{tag}.csv"
                s = timed([exe(d, tool), *args, "--out", str(out)], work)
                outs[tag] = out
                print(f"{tool} {tag}: {s:.1f} s", flush=True)
            same = filecmp.cmp(outs["base"], outs["pgo"], shallow=False)
            print(f"{tool}: output {'identical' if same else 'DIFFERS'}", flush=True)
            bad += 0 if same else 1
        sys.exit(1 if bad else 0)


if __name__ == "__main__":
    os.chdir(Path(__file__).resolve().parent.parent)
    main()
