#include "inputs.hpp"

#include <benchmark/benchmark.h>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/domain/mbpt/biorthogonalization.hpp>

using namespace sequant;

static void biorthogonal_transform(benchmark::State &state,
                                   bool factor_out_nns_projector) {
  const ExprPtr &input =
      bench::spinfree_cc_residual(state.range(0), state.range(1));
  const auto external = external_indices(input->at(0));

  for (auto _ : state) {
    ExprPtr expression = input->clone();
    ExprPtr result = mbpt::biorthogonal_transform_pre_nnsproject(
        expression, external, factor_out_nns_projector);

    benchmark::DoNotOptimize(result);
  }
}

BENCHMARK_CAPTURE(biorthogonal_transform, factored_nns_projector, true)
    ->Apply(bench::cc_residual_args);
BENCHMARK_CAPTURE(biorthogonal_transform, full, false)
    ->Apply(bench::cc_residual_args);
