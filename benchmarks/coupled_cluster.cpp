#include <benchmark/benchmark.h>

#include <SeQuant/domain/mbpt/models/cc.hpp>

static constexpr std::size_t maxRank = 15;

using namespace sequant;
using namespace sequant::mbpt;

static void cc_full_derivation(benchmark::State &state) {
  const std::size_t rank = state.range(0);

  for (auto _ : state) {
    CC cc(rank);
    auto equations = cc.t();

    benchmark::DoNotOptimize(equations);
  }
}

BENCHMARK(cc_full_derivation)->DenseRange(2, maxRank);
