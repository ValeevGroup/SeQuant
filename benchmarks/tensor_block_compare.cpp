#include <benchmark/benchmark.h>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/utility/tensor.hpp>

#include <vector>

using namespace sequant;

std::vector<Tensor> get_tensors() {
  static std::vector<Tensor> tensors = {
      parse_expr(L"t{a1}")->as<Tensor>(),
      parse_expr(L"K{i1,a4,a2;i2,i3}")->as<Tensor>(),
      parse_expr(L"t{i1,a4,a2;i2,i3;i5,i6}")->as<Tensor>(),
      parse_expr(L"K{i2,a3,a2;i4,i1}")->as<Tensor>(),
      parse_expr(L"t{a1,a2;i1,i2}")->as<Tensor>(),
      parse_expr(L"t{a4,a5;i3,i7}")->as<Tensor>(),
      parse_expr(L"t{a4<i3,i7>,a5<i3,i7>;i3,i7}")->as<Tensor>(),
  };

  return tensors;
}

static void tensor_block_equal_comparison(benchmark::State &state) {
  TensorBlockEqualComparator cmp;
  for (auto _ : state) {
    for (const Tensor &lhs : get_tensors()) {
      for (const Tensor &rhs : get_tensors()) {
        bool result = cmp(lhs, rhs);

        benchmark::DoNotOptimize(result);
      }
    }
  }
}

static void tensor_block_less_than_comparison(benchmark::State &state) {
  TensorBlockLessThanComparator cmp;
  for (auto _ : state) {
    for (const Tensor &lhs : get_tensors()) {
      for (const Tensor &rhs : get_tensors()) {
        bool result = cmp(lhs, rhs);

        benchmark::DoNotOptimize(result);
      }
    }
  }
}

BENCHMARK(tensor_block_equal_comparison);
BENCHMARK(tensor_block_less_than_comparison);
