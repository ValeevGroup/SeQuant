#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/node_batch_annotation.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <unordered_map>

int main() {
  using namespace sequant;

  // start-snippet-1
  mbpt::load(mbpt::Convention::SR);

  // register a density-fitting auxiliary space ("Κ") and make it much larger
  // than the occupied/virtual spaces -- as in forming a density-fitted
  // 2-electron integral g[i,j;a,b] ≈ Σ_Κ B[i,a,Κ] B[j,b,Κ], where Κ ranges
  // over a large auxiliary fitting basis
  auto isr = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(isr);
  isr->retrieve_ptr(L"i")->approximate_size(4);
  isr->retrieve_ptr(L"a")->approximate_size(8);
  isr->retrieve_ptr(L"Κ")->approximate_size(100'000);

  auto expr = ex<Tensor>(L"B", bra{L"i_1"}, ket{L"a_1", L"Κ_1"}) *
              ex<Tensor>(L"B", bra{L"i_2"}, ket{L"a_2", L"Κ_1"}) *
              ex<Tensor>(L"S", bra{L"i_1"}, ket{L"i_2"});
  // end-snippet-1

  // start-snippet-2
  // a BatchPolicy declares which index spaces may be sliced (split by
  // contracted vs. external role), an upper bound on the per-slice block
  // size, and -- via peak_threshold -- the peak-memory budget that actually
  // turns batching on; the default peak_threshold of +infinity means "never
  // batch"
  BatchPolicy policy{
      .is_batchable_contracted_index =
          [](Index const& ix) { return ix.space().base_key() == L"Κ"; },
      .batch_target_size = [](Index const&) -> std::size_t { return 100; },
      .peak_threshold = 5'000'000.0};

  // term_batch_axes is an optional out-channel: when set, optimize() records
  // which indices it decided to slice at each contraction node
  auto batch_axes = std::make_shared<std::unordered_map<
      Expr const*, container::vector<NodeBatchAnnotation>>>();

  auto optimized = optimize(
      expr, {.objective_function = ObjectiveFunction::DenseSpaceTimeBatched,
             .batch_policy = policy,
             .term_batch_axes = batch_axes});

  bool batched_kappa = false;
  for (auto const& [node, annotations] : *batch_axes)
    for (auto const& annotation : annotations)
      for (auto const& [axis, mode_type] : annotation.axes)
        batched_kappa |= axis.space().base_key() == L"Κ" &&
                         mode_type == BatchModeType::Contracted;
  // end-snippet-2

  assert(batched_kappa);
  (void) batched_kappa;

  return 0;
}
