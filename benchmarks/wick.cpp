#include <benchmark/benchmark.h>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

#include <optional>
#include <vector>

using namespace sequant;
using namespace sequant::mbpt;

#ifdef SEQUANT_BENCH_MAX_WICK_INPUTS
static constexpr std::size_t nInputs = SEQUANT_BENCH_MAX_WICK_INPUTS;
#else
static constexpr std::size_t nInputs = 5;
#endif

template <Statistics stats>
ExprPtr get_op_sequence(std::size_t i) {
  using Op = NormalOperator<stats>;

  switch (i) {
    case 1:
      return deserialize<ExprPtr>("f{p1;p2} t{p3;p4}") *
             ex<Op>(cre({"p_1"}), ann({"p_2"})) *
             ex<Op>(cre({"p_3"}), ann({"p_4"}));
    case 2:
      return deserialize<ExprPtr>("Â{p1,p2;p3,p4}:A g{p5,p6;p7,p8}:A") *
             ex<Op>(cre({"p_1", "p_2"}), ann({"p_3", "p_4"})) *
             ex<Op>(cre({"p_5", "p_6"}), ann({"p_7", "p_8"}));
    case 3:
      return deserialize<ExprPtr>("Â{p1,p2,p3;p4,p5,p6}:A g{p7,p8;p9,p10}:A") *
             ex<Op>(cre({"p_1", "p_2", "p_3"}), ann({"p_4", "p_5", "p_6"})) *
             ex<Op>(cre({"p_7", "p_8"}), ann({"p_9", "p_10"}));
    case 4:
      return deserialize<ExprPtr>(
                 "Â{p1,p2,p3;p4,p5,p6}:A g{p7,p8;p9,p10}:A "
                 "t{p11,p12,p13;p14,p15,p16}:A") *
             ex<Op>(cre({"p_1", "p_2", "p_3"}), ann({"p_4", "p_5", "p_6"})) *
             ex<Op>(cre({"p_7", "p_8"}), ann({"p_9", "p_10"})) *
             ex<Op>(cre({"p_11", "p_12", "p_13"}),
                    ann({"p_14", "p_15", "p_16"}));
    case 5:
      return deserialize<ExprPtr>(
                 "Â{p1,p2;p3,p4}:A g{p5,p6;p7,p8}:A "
                 "t{p9,p10;p11,p12}:A t{p13,p14;p15,p16}:A") *
             ex<Op>(cre({"p_1", "p_2"}), ann({"p_3", "p_4"})) *
             ex<Op>(cre({"p_5", "p_6"}), ann({"p_7", "p_8"})) *
             ex<Op>(cre({"p_9", "p_10"}), ann({"p_11", "p_12"})) *
             ex<Op>(cre({"p_13", "p_14"}), ann({"p_15", "p_16"}));
  }

  throw sequant::Exception("Invalid index");
};

template <Statistics stats>
static void wick(benchmark::State &state, bool full_conractions_only,
                 bool use_topology) {
  const ExprPtr &input = get_op_sequence<stats>(state.range(0));

  std::size_t n_attempted = 0;
  std::size_t n_useful = 0;

  // Note: only what is contained in this loop will be part
  // of the benchmark timings
  for (auto _ : state) {
    WickTheorem<stats> wick(input->clone());
    wick.full_contractions(full_conractions_only);
    wick.use_topology(use_topology);

    // We only don't explicitly create the expressions for all contractions in
    // order to not blow up the memory consumption of the benchmark in the cases
    // where we compute all contractions.
    ExprPtr result = wick.compute(true);

    // Prevent the Wick computation from being optimized away by the compiler
    benchmark::DoNotOptimize(result);

    const auto &statistics = wick.stats();
    n_attempted = statistics.num_attempted_contractions;
    n_useful = statistics.num_useful_contractions;
  }

  state.counters["produced"] = n_useful;
  state.counters["attempted"] = n_attempted;
  state.counters["percentage"] =
      n_attempted > 0 ? static_cast<int>(100 * static_cast<float>(n_useful) /
                                         static_cast<float>(n_attempted))
                      : 0;
}

BENCHMARK_CAPTURE(wick<Statistics::FermiDirac>, full_only_with_topology, true,
                  true)
    ->DenseRange(1, nInputs);
BENCHMARK_CAPTURE(wick<Statistics::BoseEinstein>, full_only_with_topology, true,
                  true)
    ->DenseRange(1, nInputs);

BENCHMARK_CAPTURE(wick<Statistics::FermiDirac>, full_only_without_topology,
                  true, false)
    ->DenseRange(1, nInputs);
BENCHMARK_CAPTURE(wick<Statistics::BoseEinstein>, full_only_without_topology,
                  true, false)
    ->DenseRange(1, nInputs);

BENCHMARK_CAPTURE(wick<Statistics::FermiDirac>, all_with_topology, false, true)
    ->DenseRange(1, nInputs);
BENCHMARK_CAPTURE(wick<Statistics::BoseEinstein>, all_with_topology, false,
                  true)
    ->DenseRange(1, nInputs);

BENCHMARK_CAPTURE(wick<Statistics::FermiDirac>, all_without_topology, false,
                  false)
    ->DenseRange(1, nInputs);
BENCHMARK_CAPTURE(wick<Statistics::BoseEinstein>, all_without_topology, false,
                  false)
    ->DenseRange(1, nInputs);

struct VacAvPair {
  using Connections = std::vector<std::pair<int, int>>;
  ExprPtr expr;
  std::optional<Connections> connections;
};

static constexpr std::size_t nMbptInputs = 4;

VacAvPair get_mbpt_expr(std::size_t i) {
  namespace t = sequant::mbpt::tensor;

  switch (i) {
    case 1:
      return {t::A(nₚ(-2)) * t::H(2) * t ::T(2) * t::T(2),
              VacAvPair::Connections{{1, 2}, {1, 3}}};
    case 2:
      return {t::A(nₚ(-4)) * t::h(2) * t ::t(3) * t::t(3),
              VacAvPair::Connections{{1, 2}, {1, 3}}};
    case 3:
      return {t::A(-5) * t::h(2) * t::t(2) * t::t(2) * t::t(3),
              VacAvPair::Connections{{1, 2}, {1, 3}, {1, 4}}};
    case 4:
      return {t::l(nₚ(2), nₕ(1)) * t::h(2) * t::r(nₚ(1), nₕ(0)), std::nullopt};
  }

  throw sequant::Exception("Invalid index");
}

static void mbpt_vac_av(benchmark::State &state, bool csv) {
  auto base = mbpt::get_default_mbpt_context();
  base.set(csv ? CSV::Yes : CSV::No);
  auto ctx = sequant::mbpt::set_scoped_default_mbpt_context(base);

  VacAvPair input = get_mbpt_expr(state.range(0));

  // Note: only what is contained in this loop will be part
  // of the benchmark timings
  for (auto _ : state) {
    ExprPtr result = [&]() {
      if (input.connections) {
        return mbpt::tensor::vac_av(input.expr,
                                    {.connect = input.connections.value()});
      } else {
        return mbpt::tensor::vac_av(input.expr);
      }
    }();

    benchmark::DoNotOptimize(result);
  }
}

BENCHMARK_CAPTURE(mbpt_vac_av, mbpt_vac_av, false)->DenseRange(1, nMbptInputs);
BENCHMARK_CAPTURE(mbpt_vac_av, mbpt_vac_av_csv, true)
    ->DenseRange(1, nMbptInputs);
