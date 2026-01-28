#include <benchmark/benchmark.h>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/parse.hpp>

using namespace sequant;

static constexpr std::size_t nInputs = 6;

static ExprPtr get_expression(std::size_t i) {
  switch (i) {
    case 1:
      return parse_expr(L"g{a1,i5;i6,p12}:A");
    case 2:
      return parse_expr(L"g{a1,i5;i6,p12}:A t{i6,i7;a1,a2}:A");
    case 3:
      return parse_expr(
          L"Â{i1,i2;a1,a2}:A g{i3,i4;a3,a4}:A t{a1,a3;i3,i4}:A "
          L"t{a2,a4;i1,i2}:A");
    case 4:
      return parse_expr(
          L"Â{i1,i2;a1,a2}:A DF{i3;a3;p1} DF{i4;a4;p1} t{a1;i3} t{a3;i4} "
          L"t{a2;i1} t{a4;i2}");
    case 5:
      return parse_expr(L"g{i3,i4;a3<i1,i4>,a4<i2>} s{a1<i1,i2>;a5<i3>}");
    case 6:
      return parse_expr(
          L"s{a2<i1,i2>;a6<i2,i4>} g{i3,i4;a3<i2,i4>,a4<i1,i3>} "
          L"t{a3<i2,i4>,a6<i2,i4>;i4,i2}");
  }

  throw "Invalid index";
}

static void bench_canonicalize(benchmark::State &state) {
  ExprPtr input = get_expression(state.range(0));

  for (auto _ : state) {
    ExprPtr canonicalized = canonicalize(input->clone());

    // Prevent canonicalization from being optimized away by the compiler
    benchmark::DoNotOptimize(canonicalized);
  }
}

BENCHMARK(bench_canonicalize)->Name("canonicalize")->DenseRange(1, nInputs);
