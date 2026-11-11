#include "inputs.hpp"

#include <benchmark/benchmark.h>

#include <SeQuant/core/expr.hpp>

#include <cstddef>

using namespace sequant;

// Multireference residuals differ structurally from the single-reference ones:
// tensors carry mixed index spaces in bra and ket and the projectors are
// partial or absent antisymmetrizers.

template <typename Operation>
static void nevpt2(benchmark::State &state, Operation operation) {
  const bench::Nevpt2 inputs;
  const std::size_t i = state.range(0);
  state.SetLabel(inputs.name(i));
  const ExprPtr &input = inputs.residual(i);

  for (auto _ : state) {
    ExprPtr expression = input->clone();
    operation(expression);

    benchmark::DoNotOptimize(expression);
  }
}

static void nevpt2_args(benchmark::Benchmark *benchmark) {
  benchmark->DenseRange(0, bench::Nevpt2::size() - 1);
}

BENCHMARK_CAPTURE(nevpt2, canonicalize, [](ExprPtr &expr) {
  canonicalize(expr);
})->Apply(nevpt2_args);
BENCHMARK_CAPTURE(nevpt2, simplify, [](ExprPtr &expr) {
  simplify(expr);
})->Apply(nevpt2_args);
