// SPDX-License-Identifier: Apache-2.0
// Must NOT compile. CTest builds it and asserts the static_assert message appears (D-063).
#include "core/precision.hpp"

struct FloatAccumulator {
    using eval_t = float;
    using accum_t = float;
};

static_assert(vcal::check_precision<FloatAccumulator>());
