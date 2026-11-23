#include "inputs.hpp"

#include <benchmark/benchmark.h>

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index_basis_registry.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/rules/df.hpp>

#include <memory>

using namespace sequant;

// Density fitting of the integrals in spin-free CC residuals
static void density_fit_cc_residual(benchmark::State &state) {
  const ExprPtr &input =
      bench::spinfree_cc_residual(state.range(0), state.range(1));
  // N.B. the registry's "x" space is a union of orbital spaces, not an
  // auxiliary basis
  auto registry = std::make_shared<IndexBasisRegistry>(
      *get_default_context().index_basis_registry());
  mbpt::add_df_spaces(registry);
  const auto context_resetter =
      set_scoped_default_context(get_default_context_snapshot().set(registry));
  const IndexSpace aux = registry->retrieve_space(L"Κ");

  for (auto _ : state) {
    ExprPtr result = mbpt::density_fit(input->clone(), aux, L"g", L"B");

    benchmark::DoNotOptimize(result);
  }
}

BENCHMARK(density_fit_cc_residual)->Apply(bench::cc_residual_args);
