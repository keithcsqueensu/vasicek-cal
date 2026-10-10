# SPDX-License-Identifier: Apache-2.0
"""Fail if a MANIFEST.json records a malformed SHA-256, or one that does not match the committed file.

Every tracked MANIFEST.json is read and every SHA-256 in it found, in the three forms the manifests use:

- text "sha256 <hash>" anywhere in a string value (studies' "files" entries, including hashes in parentheses);
- a key ending in "sha256" (e.g. "executable_sha256"), whose value is the hash;
- a string value that is a bare run of 32 or more hex digits (the golden manifests' {path: hash} maps), unless its
  key names a commit ("tool_commit": a git object id, not a SHA-256).

Each must be exactly 64 lowercase hex characters (a 63-character hash once reached a manifest: S-35's rows.csv).
Where the hash belongs to a file in the repo, it must also match: the key names the file (a parenthetical such as
"(not committed)" is dropped; relative to the manifest's directory, then to the repo root), and the hash is the
bare value or the "sha256 <hash>" that opens the value. Hashes are of the committed blobs (git cat-file), so
line-ending conversion in a working tree cannot affect the result. Files not in the repo (release attachments,
executables) are checked for form only. Standard library only; needs git on PATH.

    python tools/check_manifests.py              check every tracked MANIFEST.json
    python tools/check_manifests.py --self-test  check that malformed entries are caught
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from collections.abc import Iterator
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parent.parent
VALID = re.compile(r"[0-9a-f]{64}")
IN_TEXT = re.compile(r"sha256\s+([^\s,;()]+)", re.IGNORECASE)
BARE = re.compile(r"[0-9A-Fa-f]{32,}")


def git(*args: str) -> bytes:
    return subprocess.run(["git", *args], cwd=ROOT, check=True, capture_output=True).stdout


def leaves(node: object, path: tuple[str, ...] = ()) -> Iterator[tuple[tuple[str, ...], str]]:
    if isinstance(node, dict):
        for k, v in node.items():
            yield from leaves(v, (*path, str(k)))
    elif isinstance(node, list):
        for i, v in enumerate(node):
            yield from leaves(v, (*path, str(i)))
    elif isinstance(node, str):
        yield path, node


def entries(doc: object) -> Iterator[tuple[tuple[str, ...], str, bool]]:
    """(key path, hash, whether the hash is the key's own file's) for every SHA-256 in a manifest."""
    for path, value in leaves(doc):
        key = path[-1] if path else ""
        if key.lower().endswith("sha256"):
            yield path, value.strip(), False
        elif BARE.fullmatch(value.strip()) and "commit" not in key.lower():
            yield path, value.strip(), True
        else:
            for m in IN_TEXT.finditer(value):
                yield path, m.group(1), m.start() == 0


def own_file(manifest: str, key: str, tracked: set[str]) -> str | None:
    """The tracked file a key names, or None: the key without a trailing parenthetical, beside the manifest or from the root."""
    name = re.sub(r"\s*\(.*\)\s*$", "", key).strip()
    return next((c for c in (str(PurePosixPath(manifest).parent / name), name) if c in tracked), None)


def check(manifest: str, doc: object, tracked: set[str]) -> tuple[list[str], int]:
    """Errors, and how many hashes were compared with a committed file."""
    bad, matched = [], 0
    for path, digest, own in entries(doc):
        where = f"{manifest}: {'/'.join(path)}"
        if not VALID.fullmatch(digest):
            bad.append(f"{where}: {digest!r} is not 64 lowercase hex characters ({len(digest)})")
            continue
        if not own:
            continue
        cand = own_file(manifest, path[-1], tracked)
        if cand is None:
            continue
        matched += 1
        got = hashlib.sha256(git("cat-file", "blob", f"HEAD:{cand}")).hexdigest()
        if got != digest:
            bad.append(f"{where}: recorded {digest}, but {cand} at HEAD is {got}")
    return bad, matched


def self_test() -> int:
    good = "855b5abcf7b0915e66f5b56d5334d3fed396b7f20feb981893a874c62b8b5e5b"
    cases = {
        "the S-35 63-character hash, in parentheses": {"f": {"rows.parquet": f"sha256 {good}; from rows.csv "
                                                              f"(sha256 {good[:-1]}, 54,884,780 bytes)"}},
        "uppercase hex": {"files": {"a.csv": f"sha256 {good.upper()}"}},
        "65 characters": {"run": {"executable_sha256": good + "0"}},
        "a bare 63-character value": {"outputs": {"x/y.csv": good[:-1]}},
        "a non-hex character": {"files": {"a.csv": f"sha256 {good[:-1]}g"}},
    }
    failed = [name for name, doc in cases.items() if len(check("t/MANIFEST.json", doc, set())[0]) != 1]
    clean = {"files": {"a.csv (not committed)": f"sha256 {good}; from b.csv (sha256 {good})"},
             "run": {"executable_sha256": good, "commit": "2b6aed3",
                     "tool_commit": "65e84208b434926a215dd253dcaa2339a4229a29"}}
    if check("t/MANIFEST.json", clean, set())[0]:
        failed.append("a well-formed manifest")
    for name in failed:
        print(f"self-test failed: {name}")
    print("self-test passed" if not failed else f"{len(failed)} self-test case(s) failed")
    return 1 if failed else 0


def main() -> int:
    ap = argparse.ArgumentParser(description="Check the SHA-256s recorded in every MANIFEST.json")
    ap.add_argument("--self-test", action="store_true")
    if ap.parse_args().self_test:
        return self_test()
    tracked = set(git("ls-files", "-z").decode("utf-8").split("\0")) - {""}
    manifests = sorted(p for p in tracked if PurePosixPath(p).name == "MANIFEST.json")
    bad, count, matched = [], 0, 0
    for m in manifests:
        doc = json.loads(git("cat-file", "blob", f"HEAD:{m}"))
        count += sum(1 for _ in entries(doc))
        errors, n = check(m, doc, tracked)
        bad += errors
        matched += n
    for b in bad:
        print(f"error: {b}")
    print(f"{len(manifests)} manifests, {count} SHA-256 entries, {matched} matched against committed files: "
          f"{'ok' if not bad else f'{len(bad)} error(s)'}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
