# SPDX-License-Identifier: Apache-2.0
"""Render tests/engine_flags.json from the engine's flag definitions (D-200).

Scoring scripts read the engine's flag bits from this file instead of writing them as literals, so a script cannot
disagree with the engine (a scorer once had kFlagQuadratureUnconverged on bit 3; the engine has it on bit 2).
Four groups: the estimate flags (kFlag*, in Estimate2/Estimate3's flags field), the interval flags (kInterval*,
in a profile interval's flags) and, from D-219, the grid posterior's flags (kPosterior*, in PosteriorResult's
flags) and the method of moments' flags (kMom*). A group may not use one bit twice. Standard library only.

    python tools/gen_engine_flags.py           write tests/engine_flags.json
    python tools/gen_engine_flags.py --check   fail if the committed file is stale
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TARGET = "tests/engine_flags.json"
SOURCES = (
    "engine/refine.hpp",
    "engine/profile.hpp",
    "engine/conditional_pd.hpp",
    "engine/posterior.hpp",
    "engine/moments.hpp",
    "benchmark_models/ar1/search.hpp",
)
PATTERN = re.compile(r"\b(k(Flag|Interval|Posterior|Mom)\w+)\s*=\s*1u\s*<<\s*(\d+)")
GROUP = {"Flag": "estimate", "Interval": "interval", "Posterior": "posterior", "Mom": "moments"}


def build() -> str:
    groups: dict[str, dict[str, dict[str, int | str]]] = {g: {} for g in GROUP.values()}
    for src in SOURCES:
        for line in (ROOT / src).read_text(encoding="utf-8").splitlines():
            m = PATTERN.search(line)
            if not m:
                continue
            name, kind, bit = m.group(1), m.group(2), int(m.group(3))
            group = groups[GROUP[kind]]
            if name in group:
                sys.exit(f"{name} is defined twice ({group[name]['source']} and {src})")
            clash = [n for n, v in group.items() if v["bit"] == bit]
            if clash:
                sys.exit(f"{name} ({src}) uses bit {bit}, as {clash[0]} does")
            group[name] = {"bit": bit, "value": 1 << bit, "source": src}
    doc = {
        "generated_by": "tools/gen_engine_flags.py (D-200); do not edit",
        **{g: dict(sorted(groups[g].items(), key=lambda kv: int(kv[1]["bit"]))) for g in GROUP.values()},
    }
    return json.dumps(doc, indent=2) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description="Render tests/engine_flags.json from the engine's flag definitions")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    text = build()
    path = ROOT / TARGET
    if a.check:
        if not path.exists() or path.read_text(encoding="utf-8") != text:
            print(f"{TARGET} is stale: run python tools/gen_engine_flags.py")
            return 1
        print("engine flags are up to date")
        return 0
    path.write_text(text, encoding="utf-8", newline="\n")
    print(f"wrote {TARGET}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
