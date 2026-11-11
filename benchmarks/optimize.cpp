#include "inputs.hpp"

#include <benchmark/benchmark.h>

#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/optimize/common_subexpression_elimination.hpp>
#include <SeQuant/core/optimize/optimize.hpp>

#include <vector>

using namespace sequant;

// Contraction-order optimization of biorthogonal spin-free CC residuals
static void optimize_cc_result(benchmark::State &state) {
  const ResultExpr &input =
      bench::biorthogonal_cc_result(state.range(0), state.range(1));

  for (auto _ : state) {
    ResultExpr result = input.clone();
    optimize(result);

    benchmark::DoNotOptimize(result);
  }
}

BENCHMARK(optimize_cc_result)->Apply(bench::cc_residual_args);

// Common subexpression elimination across all residuals of a CC model
static void cse_cc_residuals(benchmark::State &state) {
  const std::size_t cc_rank = state.range(0);

  std::vector<ResultExpr> optimized;
  for (std::size_t p = 1; p <= cc_rank; ++p) {
    optimized.push_back(
        optimize(bench::biorthogonal_cc_result(cc_rank, p).clone()));
  }
  std::vector<ExportNode<>> trees;
  for (const ResultExpr &result : optimized) {
    trees.push_back(to_export_tree(result));
  }

  const auto to_tree = [](const auto &expr) { return to_export_tree(expr); };
  opt::CSEOptions<ExportNode<>> options;
  options.filter_predicate = opt::min_usage_filter<ExportNode<>>(2);

  for (auto _ : state) {
    std::vector<ExportNode<>> expressions = trees;
    auto positions =
        opt::eliminate_common_subexpressions(expressions, to_tree, options);

    benchmark::DoNotOptimize(positions);
    benchmark::DoNotOptimize(expressions);
  }
}

BENCHMARK(cse_cc_residuals)->DenseRange(2, bench::max_cc_rank);
