# SPDX-License-Identifier: Apache-2.0
"""The shared jackknife run: convert study_jackknife_run --replicates-out CSV to committed Parquet.

Every per-replicate field of the run (estimates, corrected rho, coverage flags of each interval,
SEs, influence and refit gaps), with doubles kept exact (the CSV carries 17 significant digits and
is parsed back to the same doubles; "nan" stays a NaN, never a null). Coverage flags are int8:
1 covers, 0 misses, -1 not assessed.

    uv run --no-project --with pyarrow==25.0.1 python studies/jackknife-bias-rho/export_replicates.py IN.csv OUT.parquet
"""

import sys

import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

INTS = {"scenario": pa.uint8(), "replicate": pa.uint16(), "flags": pa.uint32()}
DOUBLES = {"pd", "rho", "rho_tilde", "q_err_a", "se_delta", "se_jack", "infl_rho", "infl_pd", "infl_pair_rho",
           "refit_gap", "refit_gap_clean", "rho_tilde_p"}


def main():
    names = pacsv.read_csv(sys.argv[1], read_options=pacsv.ReadOptions(block_size=1 << 26)).column_names
    types = {n: INTS.get(n, pa.float64() if n in DOUBLES else pa.int8()) for n in names}
    table = pacsv.read_csv(sys.argv[1], convert_options=pacsv.ConvertOptions(column_types=types, null_values=[]))
    pq.write_table(table, sys.argv[2], compression="zstd", compression_level=19)
    print(f"{table.num_rows} rows, {table.num_columns} columns written to {sys.argv[2]}")


if __name__ == "__main__":
    main()
