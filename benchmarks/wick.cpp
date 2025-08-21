#include <benchmark/benchmark.h>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

#include <optional>
#include <vector>

using namespace sequant;
using namespace sequant::mbpt;

static constexpr std::size_t nInputs = 5;

template <Statistics stats>
ExprPtr get_op_sequence(std::size_t i) {
  using OpSeq = NormalOperatorSequence<stats>;
  using Op = NormalOperator<stats>;

  switch (i) {
    case 1:
      return parse_expr(L"f{p1;p2} t{p3;p4}") *
             ex<Op>(cre({L"p_1"}), ann({L"p_2"})) *
             ex<Op>(cre({L"p_3"}), ann({L"p_4"}));
    case 2:
      return parse_expr(L"A{p1,p2;p3,p4}:A g{p5,p6;p7,p8}:A") *
             ex<Op>(cre({L"p_1", L"p_2"}), ann({L"p_3", L"p_4"})) *
             ex<Op>(cre({L"p_5", L"p_6"}), ann({L"p_7", L"p_8"}));
    case 3:
      return parse_expr(L"A{p1,p2,p3;p4,p5,p6}:A g{p7,p8;p9,p10}:A") *
             ex<Op>(cre({L"p_1", L"p_2", L"p_3"}),
                    ann({L"p_4", L"p_5", L"p_6"})) *
             ex<Op>(cre({L"p_7", L"p_8"}), ann({L"p_9", L"p_10"}));
    case 4:
      return parse_expr(
                 L"A{p1,p2,p3;p4,p5,p6}:A g{p7,p8;p9,p10}:A "
                 L"t{p11,p12,p13;p14,p15,p16}:A") *
             ex<Op>(cre({L"p_1", L"p_2", L"p_3"}),
                    ann({L"p_4", L"p_5", L"p_6"})) *
             ex<Op>(cre({L"p_7", L"p_8"}), ann({L"p_9", L"p_10"})) *
             ex<Op>(cre({L"p_11", L"p_12", L"p_13"}),
                    ann({L"p_14", L"p_15", L"p_16"}));
    case 5:
      return parse_expr(
                 L"A{p1,p2;p3,p4}:A g{p5,p6;p7,p8}:A "
                 L"t{p9,p10;p11,p12}:A t{p13,p14;p15,p16}:A") *
             ex<Op>(cre({L"p_1", L"p_2"}), ann({L"p_3", L"p_4"})) *
             ex<Op>(cre({L"p_5", L"p_6"}), ann({L"p_7", L"p_8"})) *
             ex<Op>(cre({L"p_9", L"p_10"}), ann({L"p_11", L"p_12"})) *
             ex<Op>(cre({L"p_13", L"p_14"}), ann({L"p_15", L"p_16"}));
  }

  throw "Invalid index";
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
      return {t::A(nₚ(-4)) * t::H_(2) * t ::T_(3) * t::T_(3),
              VacAvPair::Connections{{1, 2}, {1, 3}}};
    case 3:
      return {t::A(-5) * t::H_(2) * t::T_(2) * t::T_(2) * t::T_(3),
              VacAvPair::Connections{{1, 2}, {1, 3}, {1, 4}}};
    case 4:
      return {t::L_(nₚ(2), nₕ(1)) * t::H_(2) * t::R_(nₚ(1), nₕ(0)),
              std::nullopt};
  }

  throw "Invalid index";
}

static void mbpt_vac_av(benchmark::State &state, bool csv) {
  auto ctx = sequant::mbpt::set_scoped_default_mbpt_context(
      mbpt::Context(csv ? mbpt::CSV::Yes : mbpt::CSV::No));

  VacAvPair input = get_mbpt_expr(state.range(0));

  // Note: only what is contained in this loop will be part
  // of the benchmark timings
  for (auto _ : state) {
    ExprPtr result = [&]() {
      if (input.connections) {
        return mbpt::tensor::vac_av(input.expr, input.connections.value());
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
