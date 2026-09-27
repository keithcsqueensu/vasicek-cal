// SPDX-License-Identifier: Apache-2.0
#include "core/precision.hpp"

#include <type_traits>

#include "tests/harness/vcal_test.hpp"

namespace {
struct FloatAccumulator {  // violates D-039: accumulation must be FP64
    using eval_t = float;
    using accum_t = float;
};
struct IntegerEval {
    using eval_t = int;
    using accum_t = double;
};
struct MissingMembers {};
}  // namespace

static_assert(vcal::is_precision_v<vcal::PrecisionF64>);
static_assert(vcal::is_precision_v<vcal::PrecisionMixed>);
static_assert(!vcal::is_precision_v<FloatAccumulator>);
static_assert(!vcal::is_precision_v<IntegerEval>);
static_assert(!vcal::is_precision_v<MissingMembers>);

VCAL_TEST(precision_policies_types) {
    VCAL_CHECK((std::is_same_v<vcal::PrecisionF64::eval_t, double>));
    VCAL_CHECK((std::is_same_v<vcal::PrecisionMixed::eval_t, float>));
    VCAL_CHECK((std::is_same_v<vcal::PrecisionMixed::accum_t, double>));
}
