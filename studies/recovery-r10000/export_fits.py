# SPDX-License-Identifier: Apache-2.0
"""S-13 targeted: export the per-replicate facts its verdicts rest on to a committed Parquet file.

The harness's saved fits (recovery_harness --save-fits) are a raw 200-byte record per replicate and
scenario, about 160 MB at R = 10,000, too large to commit. This keeps, for the targeted scenarios
only, each replicate's PD and rho estimates and their profile-likelihood interval ends and flags:
enough to recompute every coverage and miss direction S-13 reports, with any tool that reads
Parquet. Values are the exact doubles (no rounding).

    uv run --no-project --with pyarrow==25.0.1 python studies/recovery-r10000/export_fits.py FITS OUT.parquet
"""

import struct
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

R = 10_000
SCENARIOS = [1, 2, 6, 7, 9, 12, 15, 18, 21, 22, 24, 27, 28, 29, 30, 31, 37, 41, 45, 47, 49, 54, 55, 59, 60, 64, 65,
             68, 72, 73, 74]
# The harness's Fit record (tests/recovery/recovery.hpp), 200 bytes, little-endian, x64 layout.
FIT = struct.Struct("<2d2dd I4x 2d2d 2I d 2d2d 2I d d d d I4x d d d")


def main():
    data = Path(sys.argv[1]).read_bytes()
    assert FIT.size == 200 and len(data) == R * 81 * FIT.size, "expected R = 10,000 x 81 fit records"
    cols = {k: [] for k in ("scenario", "replicate", "pd", "rho", "fit_flags", "pd_lo", "pd_hi", "pd_profile_flags",
                            "rho_lo", "rho_hi", "rho_profile_flags")}
    for sid in SCENARIOS:
        for r in range(R):
            f = FIT.unpack_from(data, (r * 81 + sid) * FIT.size)
            for k, v in zip(cols, (sid, r, f[0], f[1], f[5], f[6], f[8], f[10], f[7], f[9], f[11])):
                cols[k].append(v)
    types = {"scenario": pa.uint8(), "replicate": pa.uint16(), "fit_flags": pa.uint32(), "pd_profile_flags": pa.uint32(),
             "rho_profile_flags": pa.uint32()}
    table = pa.table({k: pa.array(v, type=types.get(k, pa.float64())) for k, v in cols.items()})
    pq.write_table(table, sys.argv[2], compression="zstd", compression_level=19)
    print(f"{table.num_rows} rows written to {sys.argv[2]}")


if __name__ == "__main__":
    main()
