# SPDX-License-Identifier: Apache-2.0
"""Fail if a large file is committed outside the named exceptions to D-163 (D-183).

Every file tracked at HEAD above LIMIT bytes must be listed in tools/committed_data_files.txt,
and every listed file must exist with exactly the listed size and SHA-256. Sizes and hashes are
of the committed blobs (git cat-file), so line-ending conversion in a working tree cannot affect
the result. Standard library only; needs git on PATH.

    python tools/check_data_files.py
"""

from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LIST = "tools/committed_data_files.txt"
LIMIT = 1 << 20  # 1 MiB: per-replicate study files are megabytes; summaries and goldens are below


def git(*args: str) -> bytes:
    return subprocess.run(["git", *args], cwd=ROOT, check=True, capture_output=True).stdout


def listed() -> dict[str, tuple[str, int]]:
    out = {}
    for line in (ROOT / LIST).read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            sha, size, path = line.split()
            out[path] = (sha, int(size))
    return out


def main() -> int:
    allowed = listed()
    tracked = {}
    for rec in git("ls-tree", "-r", "-l", "-z", "HEAD").split(b"\0"):
        if rec:
            meta, path = rec.decode("utf-8").split("\t", 1)
            kind, size = meta.split()[1], meta.split()[3]
            if kind == "blob":
                tracked[path] = int(size)
    bad = []
    for path, size in sorted(tracked.items()):
        if size > LIMIT and path not in allowed:
            bad.append(f"{path}: {size:,} bytes, above {LIMIT:,} and not in {LIST} (D-163: attach it to a release)")
    for path, (sha, size) in sorted(allowed.items()):
        if path not in tracked:
            bad.append(f"{path}: listed in {LIST} but not tracked (remove the line)")
            continue
        got = hashlib.sha256(git("cat-file", "blob", f"HEAD:{path}")).hexdigest()
        if tracked[path] != size or got != sha:
            bad.append(f"{path}: changed ({tracked[path]:,} bytes, sha256 {got}); a listed file is frozen")
    for b in bad:
        print(f"error: {b}")
    big = sum(1 for s in tracked.values() if s > LIMIT)
    print(f"{len(tracked)} tracked files, {big} above {LIMIT:,} bytes, {len(allowed)} listed exceptions: "
          + ("FAIL" if bad else "ok"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
