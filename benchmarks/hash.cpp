#include "inputs.hpp"

#include <benchmark/benchmark.h>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>

#include <cstddef>
#include <vector>

using namespace sequant;

namespace {

// Tensors of a spin-free CC residual, the way they are hashed in practice
std::vector<Tensor> collect_tensors(const ExprPtr &expr) {
  std::vector<Tensor> tensors;
  expr->visit(
      [&tensors](const ExprPtr &node) {
        if (node.is<Tensor>()) {
          tensors.push_back(node.as<Tensor>());
        }
      },
      /* atoms_only = */ true);
  return tensors;
}

}  // namespace

static void hash_index(benchmark::State &state) {
  const auto tensors = collect_tensors(
      bench::spinfree_cc_residual(state.range(0), state.range(1)));
  container::svector<Index> indices;
  for (const Tensor &tensor : tensors) {
    indices.insert(indices.end(), tensor.bra().begin(), tensor.bra().end());
    indices.insert(indices.end(), tensor.ket().begin(), tensor.ket().end());
    indices.insert(indices.end(), tensor.aux().begin(), tensor.aux().end());
  }

  for (auto _ : state) {
    for (const Index &index : indices) {
      auto hash = hash_value(index);
      benchmark::DoNotOptimize(hash);
    }
  }
  state.SetItemsProcessed(state.iterations() * indices.size());
}

// Hash values of tensors are memoized and invalidated on mutation. Re-setting
// the bra invalidates the hash such that it is computed from scratch (the
// bra copy is part of the timings)
static void hash_tensor_cold(benchmark::State &state) {
  auto tensors = collect_tensors(
      bench::spinfree_cc_residual(state.range(0), state.range(1)));

  for (auto _ : state) {
    for (Tensor &tensor : tensors) {
      auto bra = tensor.bra();
      tensor.set_bra(std::move(bra));
      auto hash = tensor.hash_value();
      benchmark::DoNotOptimize(hash);
    }
  }
  state.SetItemsProcessed(state.iterations() * tensors.size());
}

// Accumulation of the terms of a residual in a hash table, as done when
// summands are merged eagerly
static void hashing_accumulator(benchmark::State &state) {
  const ExprPtr &input =
      bench::unsimplified_cc_residual(state.range(0), state.range(1));

  for (auto _ : state) {
    // The accumulator merges duplicate terms in place
    const ExprPtr terms = input->clone();
    HashingAccumulator accumulator;
    for (const ExprPtr &term : *terms) {
      accumulator.append(term);
    }
    ExprPtr result = accumulator.make_expr();

    benchmark::DoNotOptimize(result);
  }
}

BENCHMARK(hash_index)->Apply(bench::cc_residual_args);
BENCHMARK(hash_tensor_cold)->Apply(bench::cc_residual_args);
BENCHMARK(hashing_accumulator)->Apply(bench::cc_residual_args);
