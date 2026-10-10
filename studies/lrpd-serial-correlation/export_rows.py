# SPDX-License-Identifier: Apache-2.0
"""S-35: convert study_lrpd_serial's rows (--out, --rate-out) to Parquet for the release (D-163).

Every double is parsed from its 17-significant-digit text, an exact round trip; integers keep integer types.

    uv run --no-project --with pyarrow==25.0.1 python studies/lrpd-serial-correlation/export_rows.py \\
        ROWS.csv RATE_ROWS.csv OUT_DIR
"""

from __future__ import annotations

import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.csv as pc
import pyarrow.parquet as pq

INT8 = {"wald_pd", "wald_rho", "t_pd", "t_rho", "prof_pd", "prof_rho", "boot_pd", "boot_rho", "q_prof",
        "obs_cover_pd", "obs_cover_rho", "cover_pd", "cover_rho", "cover_q"}
UINT32 = {"scenario", "z_source_scenario", "replicate", "flags", "prof_flags_pd", "prof_flags_rho", "obs_flags"}
INT64 = {"periods", "obligors"}


def convert(src: Path, dst: Path) -> None:
    lines = src.read_text().splitlines()
    names = next(ln for ln in lines if not ln.startswith("#")).split(",")
    types = {n: pa.int8() if n in INT8 else pa.uint32() if n in UINT32 else pa.int64() if n in INT64 else pa.float64()
             for n in names}
    comments = sum(1 for ln in lines if ln.startswith("#"))
    table = pc.read_csv(src, read_options=pc.ReadOptions(skip_rows=comments),
                        convert_options=pc.ConvertOptions(column_types=types))
    pq.write_table(table, dst, compression="zstd", compression_level=19)
    print(f"wrote {dst}: {table.num_rows} rows")


def main() -> int:
    rows, rate, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    convert(rows, out / "rows.parquet")
    convert(rate, out / "rate_rows.parquet")
    return 0


if __name__ == "__main__":
    sys.exit(main())
