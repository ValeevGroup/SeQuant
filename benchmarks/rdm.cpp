#include <benchmark/benchmark.h>

#include <SeQuant/core/density.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/domain/mbpt/rdm.hpp>

#include <string>
#include <vector>

using namespace sequant;

static FNOperator make_nop(std::size_t rank) {
  std::vector<std::wstring> cre_labels, ann_labels;
  for (std::size_t i = 1; i <= rank; ++i) {
    cre_labels.push_back(L"i_" + std::to_wstring(i));
    ann_labels.push_back(L"i_" + std::to_wstring(rank + i));
  }
  return FNOperator(cre(cre_labels), ann(ann_labels));
}

// Expansion of a cumulant of the given rank in terms of densities
static void cumulants_to_densities(benchmark::State &state) {
  const ExprPtr cumulant = density::make_cumulant(make_nop(state.range(0)));

  for (auto _ : state) {
    ExprPtr result = mbpt::decompositions::cumulants_to_densities(cumulant);

    benchmark::DoNotOptimize(result);
  }
}

BENCHMARK(cumulants_to_densities)->DenseRange(2, 4);

// Decomposition of the 3-body density in terms of cumulants
static void three_body_decomp(benchmark::State &state) {
  const ExprPtr nop = ex<FNOperator>(make_nop(3));

  for (auto _ : state) {
    auto result = mbpt::decompositions::three_body_decomp(nop, false);

    benchmark::DoNotOptimize(result);
  }
}

BENCHMARK(three_body_decomp);
