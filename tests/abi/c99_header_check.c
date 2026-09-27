/* SPDX-License-Identifier: Apache-2.0
 *
 * vcal.h must be plain C99 (D-144): this translation unit is compiled as C99 with pedantic
 * errors (C11 on MSVC, which has no C99 mode), so a C++-only construct in the header fails the
 * build. Compile-only.
 */
#include "vcal/vcal.h"

int abi_c99_header_check(void);

int abi_c99_header_check(void) {
    vcal_estimate estimate;
    vcal_status status = VCAL_OK;
    estimate.struct_size = (uint32_t)sizeof estimate;
    estimate.flags = VCAL_FLAG_NEAR_BOUND | VCAL_FLAG_GRID_EDGE;
    return (int)estimate.struct_size + (int)estimate.flags + status + VCAL_ABI_VERSION;
}
