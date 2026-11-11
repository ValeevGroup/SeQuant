#include "inputs.hpp"

#include <benchmark/benchmark.h>

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/export/itf.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index_basis_registry.hpp>
#include <SeQuant/core/optimize/optimize.hpp>

#include <string>

using namespace sequant;

namespace {

// Optimized, biorthogonal spin-free CC residual
ResultExpr optimized_result(std::size_t cc_rank, std::size_t projection) {
  ResultExpr result =
      bench::biorthogonal_cc_result(cc_rank, projection).clone();
  optimize(result);
  return result;
}

}  // namespace

// Conversion of an optimized residual into an export tree
static void to_export_tree_cc_result(benchmark::State &state) {
  const ResultExpr input = optimized_result(state.range(0), state.range(1));

  for (auto _ : state) {
    ExportNode<> tree = to_export_tree(input);

    benchmark::DoNotOptimize(tree);
  }
}

// Preprocessing and code generation for an export tree
static void export_itf_cc_result(benchmark::State &state) {
  const ExportNode<> tree =
      to_export_tree(optimized_result(state.range(0), state.range(1)));

  const auto registry = get_default_context().index_basis_registry();
  const IndexSpace occupied = registry->retrieve_space(L"i");
  const IndexSpace unoccupied = registry->retrieve_space(L"a");

  for (auto _ : state) {
    ItfGenerator<ItfContext> generator;
    ItfContext context;
    context.set_tag(occupied, "c");
    context.set_tag(unoccupied, "e");
    context.set_name(occupied, "Closed");
    context.set_name(unoccupied, "External");

    export_expression(tree, generator, context);

    std::string code = generator.get_generated_code();
    benchmark::DoNotOptimize(code);
  }
}

BENCHMARK(to_export_tree_cc_result)->Apply(bench::cc_residual_args);
BENCHMARK(export_itf_cc_result)->Apply(bench::cc_residual_args);
