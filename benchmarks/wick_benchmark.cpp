#include <SeQuant/core/wick.hpp>
#include "SeQuant/core/op.hpp"
#include "SeQuant/core/runtime.hpp"
#include "SeQuant/core/space.hpp"
#include "SeQuant/core/tensor.hpp"
#include "SeQuant/domain/mbpt/convention.hpp"

#include <benchmark/benchmark.h>

#include <iostream>

using namespace sequant;

auto do_wick(bool use_nop_partitions, bool use_op_partitions) {
  constexpr Vacuum V = Vacuum::SingleProduct;

  auto opseq = FNOperatorSeq({FNOperator(IndexList{L"i_1", L"i_2"},
                                         {Index(L"a_1", {L"i_1", L"i_2"}),
                                          Index(L"a_2", {L"i_1", L"i_2"})},
                                         V),
                              FNOperator({L"p_1", L"p_2"}, {L"p_3", L"p_4"}),
                              FNOperator({Index(L"a_3", {L"i_3", L"i_4"}),
                                          Index(L"a_4", {L"i_3", L"i_4"})},
                                         IndexList{L"i_3", L"i_4"}),
                              FNOperator({Index(L"a_5", {L"i_5", L"i_6"}),
                                          Index(L"a_6", {L"i_5", L"i_6"})},
                                         IndexList{L"i_5", L"i_6"})});
  auto wick = FWickTheorem{opseq};
  wick.set_nop_connections({{1, 2}, {1, 3}}).use_topology(true);

  if (use_nop_partitions) wick.set_nop_partitions({{2, 3}});
  if (use_op_partitions)
    wick.set_op_partitions(
        {{0, 1}, {2, 3}, {4, 5}, {6, 7}, {8, 9}, {10, 11}, {12, 13}, {14, 15}});
  auto wick_result = wick.compute();

  // multiply tensor factors and expand
  auto wick_result_2 =
      ex<Constant>(rational{1, 256}) *
      ex<Tensor>(
          L"A", IndexList{L"i_1", L"i_2"},
          IndexList{{L"a_1", {L"i_1", L"i_2"}}, {L"a_2", {L"i_1", L"i_2"}}},
          IndexList{}, Symmetry::antisymm) *
      ex<Tensor>(L"g", WstrList{L"p_1", L"p_2"}, WstrList{L"p_3", L"p_4"},
                 WstrList{}, Symmetry::antisymm) *
      ex<Tensor>(
          L"t",
          IndexList{{L"a_3", {L"i_3", L"i_4"}}, {L"a_4", {L"i_3", L"i_4"}}},
          IndexList{L"i_3", L"i_4"}, IndexList{}, Symmetry::antisymm) *
      ex<Tensor>(
          L"t",
          IndexList{{L"a_5", {L"i_5", L"i_6"}}, {L"a_6", {L"i_5", L"i_6"}}},
          IndexList{L"i_5", L"i_6"}, IndexList{}, Symmetry::antisymm) *
      wick_result;
  expand(wick_result_2);
  wick.reduce(wick_result_2);
  rapid_simplify(wick_result_2);

  return std::make_pair(wick_result, wick_result_2);
}

static void BM_wick(benchmark::State &state) {
  for (auto _ : state) {
    benchmark::DoNotOptimize(do_wick(state.range(0), state.range(1)));
  }
}

BENCHMARK(BM_wick)->ArgsProduct({{0, 1}, {0, 1}});

int main(int argc, char **argv) {
  sequant::set_locale();
  sequant::detail::OpIdRegistrar op_id_registrar;
  sequant::set_default_context(
      Context(Vacuum::SingleProduct, IndexSpaceMetric::Unit,
              BraKetSymmetry::conjugate, SPBasis::spinorbital));
  mbpt::set_default_convention();
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
  set_num_threads(1);

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
}
