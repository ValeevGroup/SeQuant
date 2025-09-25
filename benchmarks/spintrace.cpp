#include <benchmark/benchmark.h>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

using namespace sequant;

static constexpr std::size_t nInputs = 2;

static ExprPtr get_expression(std::size_t i) {
  switch (i) {
    case 1:
      return parse_expr(L"2 t{i1;a1<i1>} F{a1<i1>;i1}");
    case 2:
      return parse_expr(
          L"f{i1;a1} t{a1;i1} + 1/2 g{i1,i2;a1,a2}:A t{a1;i1} t{a2;i2} "
          L"+ 1/4 g{i1,i2;a1,a2}:A t{a1,a2;i1,i2}:A");
  }

  throw "Invalid index";
}

template <bool closed_shell>
static void spintrace(benchmark::State &state, bool assume_spinfree) {
  ExprPtr input = get_expression(state.range(0));

  for (auto _ : state) {
    ExprPtr result = [&]() {
      if constexpr (closed_shell) {
        return mbpt::closed_shell_spintrace(input->clone());
      } else {
        return mbpt::spintrace(input->clone(), {}, assume_spinfree);
      }
    }();

    // Prevent spintracing from being optimized away by the compiler
    benchmark::DoNotOptimize(result);
  }
}

BENCHMARK_CAPTURE(spintrace<false>, spintrace_remove_spin, true)
    ->DenseRange(1, nInputs);

BENCHMARK_CAPTURE(spintrace<false>, spintrace_keep_spin, false)
    ->DenseRange(1, nInputs);

BENCHMARK_CAPTURE(spintrace<true>, closed_shell_spintrace, true)
    ->DenseRange(1, nInputs);
