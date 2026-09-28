# SPDX-License-Identifier: Apache-2.0
"""S-34: convert study_severe_period --replicates-out CSV to committed Parquet.

One row per scenario, replicate and variant (0 the original panel; 1-4 as in PREDICTION.md): the
polished estimates, the three profile intervals, the fit's flags, coverage flags (1 covers, 0 misses,
-1 not computed) and the original panel's SEs. Doubles are kept exact (the CSV carries 17 significant
digits and is parsed back to the same doubles; "nan" stays a NaN, never a null).

    uv run --no-project --with pyarrow==25.0.1 python studies/severe-period-sensitivity/export_replicates.py IN.csv OUT.parquet
"""

import sys

import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

TYPES = {"scenario": pa.uint8(), "replicate": pa.uint16(), "variant": pa.uint8(), "flags": pa.uint32(),
         "pd_cover": pa.int8(), "rho_cover": pa.int8(), "q_cover": pa.int8()}


def main():
    names = pacsv.read_csv(sys.argv[1], read_options=pacsv.ReadOptions(block_size=1 << 26)).column_names
    types = {n: TYPES.get(n, pa.float64()) for n in names}
    table = pacsv.read_csv(sys.argv[1], convert_options=pacsv.ConvertOptions(column_types=types, null_values=[]))
    pq.write_table(table, sys.argv[2], compression="zstd", compression_level=19)
    print(f"{table.num_rows} rows, {table.num_columns} columns written to {sys.argv[2]}")


if __name__ == "__main__":
    main()
