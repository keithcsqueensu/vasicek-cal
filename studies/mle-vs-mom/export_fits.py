# SPDX-License-Identifier: Apache-2.0
"""The estimator-comparison pass: convert study_estimator_pass's rows (CSV) to the committed Parquet file.

One row per (scenario, replicate), every column as written: the doubles are parsed from 17
significant digits, an exact round trip, so the Parquet values are the run's doubles bit for bit.

    uv run --no-project --with pyarrow==25.0.1 python studies/mle-vs-mom/export_fits.py ROWS_CSV OUT.parquet
"""

import sys

import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

INTS = {"scenario": pa.uint16(), "replicate": pa.uint16(), "periods": pa.uint16(), "obligors": pa.uint32(),
        "zero_periods": pa.uint16(), "full_periods": pa.uint16()}


def column_type(name):
    if name in INTS:
        return INTS[name]
    if name.endswith("_flags") or "_prof_flags_" in name:
        return pa.uint32()
    if name.endswith("_status") or "_cover_" in name:
        return pa.int8()
    return pa.float64()


def main():
    src, dst = sys.argv[1], sys.argv[2]
    header = open(src, encoding="ascii").readline().strip().split(",")
    # No null values: "nan" (a refused or unavailable value) must stay a NaN double, as the run wrote it,
    # not become a null (pyarrow's default null spellings include "nan").
    opts = pacsv.ConvertOptions(column_types={c: column_type(c) for c in header}, null_values=[],
                                strings_can_be_null=False)
    table = pacsv.read_csv(src, convert_options=opts)
    table = table.sort_by([("scenario", "ascending"), ("replicate", "ascending")])
    pq.write_table(table, dst, compression="zstd", compression_level=19)
    print(f"{table.num_rows} rows x {table.num_columns} columns written to {dst}")


if __name__ == "__main__":
    main()
