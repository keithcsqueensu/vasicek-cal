# SPDX-License-Identifier: Apache-2.0
"""S-5's BCa interval ends, for the interval comparison: CSV rows to a committed Parquet file.

The shared jackknife run (D-159) wrote each replicate's BCa coverage flags but not the interval ends.
To compare widths, the same run was repeated at its own commit (582895c) with only the ends added to
its --replicates-out rows (logit PD and rho, NaN where not computed); see MANIFEST.json. This keeps
the ends and the coverage flags, and checks that every flag equals the committed replicates.parquet.

    uv run --no-project --with pyarrow==25.0.1 python studies/jackknife-bias-rho/export_bca_ends.py ROWS_CSV OUT.parquet
"""

import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
KEEP = ["scenario", "replicate", "pd_bca_cover", "rho_bca_cover", "pd_bca_lo_u", "pd_bca_hi_u", "rho_bca_lo_u",
        "rho_bca_hi_u"]


def main():
    src, dst = sys.argv[1], sys.argv[2]
    types = {"scenario": pa.uint8(), "replicate": pa.uint16(), "pd_bca_cover": pa.int8(), "rho_bca_cover": pa.int8()}
    for c in KEEP[4:]:
        types[c] = pa.float64()
    opts = pacsv.ConvertOptions(column_types=types, include_columns=KEEP, null_values=[], strings_can_be_null=False)
    t = pacsv.read_csv(src, convert_options=opts).sort_by([("scenario", "ascending"), ("replicate", "ascending")])
    old = pq.read_table(HERE / "replicates.parquet", columns=["scenario", "replicate", "pd_bca_cover", "rho_bca_cover"])
    old = old.sort_by([("scenario", "ascending"), ("replicate", "ascending")])
    if t.num_rows != old.num_rows:
        sys.exit(f"{t.num_rows} rows, the committed run has {old.num_rows}")
    for c in ("scenario", "replicate", "pd_bca_cover", "rho_bca_cover"):
        diff = pc.sum(pc.not_equal(t.column(c).cast(pa.int64()), old.column(c).cast(pa.int64()))).as_py()
        if diff:
            sys.exit(f"{c}: {diff} rows differ from the committed replicates.parquet; not written")
    pq.write_table(t, dst, compression="zstd", compression_level=19)
    print(f"{t.num_rows} rows written to {dst}; every BCa coverage flag equals the committed run's")


if __name__ == "__main__":
    main()
