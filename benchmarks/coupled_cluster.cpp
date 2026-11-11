#include <benchmark/benchmark.h>

#include "inputs.hpp"

#include <SeQuant/domain/mbpt/models/cc.hpp>

static constexpr std::size_t maxRank = sequant::bench::max_derivation_cc_rank;

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

static void cc_lambda_derivation(benchmark::State &state) {
  const std::size_t rank = state.range(0);

  for (auto _ : state) {
    CC cc(rank);
    auto equations = cc.λ();

    benchmark::DoNotOptimize(equations);
  }
}

BENCHMARK(cc_lambda_derivation)->DenseRange(2, bench::max_cc_rank);

// right-hand EOM-CC sigma equations for the CCSD model
static void cc_eom_derivation(benchmark::State &state) {
  const std::size_t np = state.range(0);
  const std::size_t nh = state.range(1);

  for (auto _ : state) {
    CC cc(2);
    auto equations = cc.eom_r(nₚ(np), nₕ(nh));

    benchmark::DoNotOptimize(equations);
  }
}

BENCHMARK(cc_eom_derivation)->Args({1, 1})->Args({2, 1})->Args({2, 2});

// unitary CC with a truncated BCH or Bernoulli expansion of H̄
static void ucc_derivation(benchmark::State &state,
                           CC::HbarExpansion expansion) {
  const std::size_t rank = state.range(0);
  const std::size_t comm_rank = state.range(1);

  for (auto _ : state) {
    CC cc(rank, {.ansatz = CC::Ansatz::U,
                 .hbar_comm_rank = comm_rank,
                 .hbar_expansion = expansion});
    auto equations = cc.t();

    benchmark::DoNotOptimize(equations);
  }
}

BENCHMARK_CAPTURE(ucc_derivation, bch, CC::HbarExpansion::BCH)
    ->ArgsProduct({{2}, {2, 3}});
BENCHMARK_CAPTURE(ucc_derivation, bernoulli, CC::HbarExpansion::Bernoulli)
    ->ArgsProduct({{2}, {2}});
