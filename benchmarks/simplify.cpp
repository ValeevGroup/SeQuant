#include <benchmark/benchmark.h>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/parse.hpp>

using namespace sequant;

static constexpr std::size_t nInputs = 2;

static ExprPtr get_expression(std::size_t i) {
  switch (i) {
    case 1:
      return parse_expr(
          L"1/2 Â{i1,i2;a1,a2}:A g{i3,i4;a3,a4}:A t{a1,a3;i3,i4}:A "
          L"t{a2,a4;i1,i2}:A"
          L"- 1/2 Â{i2,i1;a1,a2}:A g{i3,i4;a4,a3}:A t{a2,a3;i2,i1}:A "
          L"t{a4,a1;i3,i4}:A ");
    case 2:
      return parse_expr(
          L"Â{p4;p1;}:A 1/3 g{p1,p2;p3,p4}:A t{p3;p2} "
          L"+ Â{p1;p4;}:A 1/3 g{p4,p2;p3,p1}:A t{p3;p2} "
          L"+ Â{p3;p1;}:A 1/3 g{p2,p1;p3,p4}:A t{p4;p2} ");
  }

  throw "Invalid index";
}

template <bool rapid_only>
static void simplify(benchmark::State &state) {
  ExprPtr input = get_expression(state.range(0));

  for (auto _ : state) {
    ExprPtr expression = input->clone();
    if constexpr (rapid_only) {
      rapid_simplify(expression);
    } else {
      simplify(expression);
    }

    // Prevent simplification from being optimized away by the compiler
    benchmark::DoNotOptimize(expression);
  }
}

BENCHMARK(simplify<false>)->Name("simplify")->DenseRange(1, nInputs);

BENCHMARK(simplify<true>)->Name("rapid_simplify")->DenseRange(1, nInputs);
