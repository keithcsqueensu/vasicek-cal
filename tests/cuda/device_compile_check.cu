// SPDX-License-Identifier: Apache-2.0
//
// Compile-only check that every device-visible function really compiles as device code
// (K-1, D-059, D-100). VCAL_HD inline functions only get device code when device code
// calls them, so compiling the headers alone proves little. This kernel calls each one:
// the special functions, all integrators (incl. the parity SplitRule) on the real integrand, the
// objective, the grid, weighted_sum and ArgMax. Nothing runs it (CI has no GPU). A host-only
// callee (e.g. a std:: function nvcc does not provide on the device) fails the build,
// because CUDA targets compile with --Werror cross-execution-space-call.
#include <cstdint>

#include "core/grid.hpp"
#include "core/model/binomial_mixture.hpp"
#include "core/model/vasicek.hpp"
#include "core/objectives/binomial_mixture.hpp"
#include "core/quadrature/composite_legendre.hpp"
#include "core/quadrature/gauss_hermite.hpp"
#include "core/reducers/argmax.hpp"
#include "core/special/erfcx.hpp"
#include "core/special/inverse_mills.hpp"
#include "core/special/lbinom.hpp"
#include "core/special/log_add_exp.hpp"
#include "core/special/log_phi.hpp"
#include "core/special/probit.hpp"
#include "engine/surface.hpp"

namespace {

using P = vcal::PrecisionF64;
using Objective = vcal::objectives::BinomialMixture<P>;

__global__ void exercise_device_code(const double* in, const double* node, const double* log_weight, int n_nodes,
                                     const double* gl_node, const double* gl_weight, int gl_n, double* out,
                                     std::int64_t* index_out) {
    const int i = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    const double x = in[i];
    namespace sp = vcal::special;

    double acc = sp::log_phi(x) + sp::erfcx(x) + sp::probit(0.3) + sp::probit_upper(0.3) + sp::inverse_mills(x) +
                 sp::log_add_exp(x, 1.0) + sp::lbinom(1000, 7) + sp::lbinom(40, 20) + sp::stirling_error(100) +
                 sp::stirling_error(5);

    const vcal::quadrature::GaussHermiteRule rule{node, log_weight, n_nodes};
    const vcal::quadrature::GaussHermiteAdaptive<P> adaptive{rule};
    const vcal::quadrature::GaussHermiteFixed<P> fixed{rule};
    const Objective objective;
    const Objective::Obs y{1000, 7};
    const double theta_values[2] = {0.01, 0.12};
    acc += objective.log_contrib(y, Objective::theta(theta_values), adaptive);
    acc += objective.log_contrib(y, Objective::theta(theta_values), fixed);
    const vcal::quadrature::CompositeLegendre<P> composite{{gl_node, gl_weight, gl_n}, 16};
    const vcal::quadrature::SplitRule<P> split{adaptive, composite};
    const Objective::Obs zero_defaults{1000, 0};
    acc += objective.log_contrib(zero_defaults, Objective::theta(theta_values), split);
    acc += Objective::rounding_scale(y, acc);

    const auto integrand =
        vcal::model::make_binomial_log_integrand(vcal::model::make_vasicek1f(0.02, 0.2), 100, 0);
    const auto hint = vcal::model::binomial_mixture_hint(integrand);
    acc += hint.mode + hint.scale + integrand(0.5);

    const vcal::Axis axis{1e-4, 0.5, 11, vcal::AxisScale::Logit};
    const vcal::Grid<2> grid{{axis, {1e-3, 0.5, 7, vcal::AxisScale::Probit}}};
    double values[2];
    grid.values(9, values);
    acc += values[0] + values[1] + axis.value_at(3) + axis.scaled_at(4) + axis.step() +
           vcal::grid::dvalue_dscaled(vcal::AxisScale::Probit, 0.2) +
           vcal::grid::to_scaled(vcal::AxisScale::Log, 0.3) + vcal::grid::from_scaled(vcal::AxisScale::Logit, -2.0);

    acc += vcal::engine::weighted_sum(in, 2, 1, in, 0);

    const vcal::reducers::ArgMax argmax;
    auto state = argmax.init();
    argmax.push(state, i, acc);
    auto merged = argmax.init();
    argmax.merge(merged, state);

    out[i] = acc;
    index_out[i] = merged.k;
}

}  // namespace

// A host entry point so the kernel is odr-used and fully compiled for every target architecture.
void vcal_device_compile_check_launch(const double* in, const double* node, const double* log_weight, int n_nodes,
                                      const double* gl_node, const double* gl_weight, int gl_n, double* out,
                                      std::int64_t* index_out) {
    exercise_device_code<<<1, 32>>>(in, node, log_weight, n_nodes, gl_node, gl_weight, gl_n, out, index_out);
}
