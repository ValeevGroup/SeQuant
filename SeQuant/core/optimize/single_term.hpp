#ifndef SEQUANT_CORE_OPTIMIZE_SINGLE_TERM_HPP
#define SEQUANT_CORE_OPTIMIZE_SINGLE_TERM_HPP

// The shared DP helpers (cost counters, OptRes/EvalSequence, init_results)
// live in single_term_detail.hpp; the generic CostModel driver and the
// objective models live in cost_model.hpp (which includes the detail header).
// Including cost_model.hpp here brings both, then the public
// single_term_opt<Metric> entry points below route through run_single_term_opt.
#include <SeQuant/core/optimize/cost_model.hpp>
#include <SeQuant/core/optimize/options.hpp>
#include <SeQuant/core/optimize/single_term_detail.hpp>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor_network.hpp>

#include <range/v3/view/filter.hpp>
#include <range/v3/view/reverse.hpp>

#include <functional>

namespace sequant::opt {
namespace detail {

///
/// \tparam Metric Objective function (ObjectiveFunction::DenseFLOPs or
///         ObjectiveFunction::DenseSize or ObjectiveFunction::DensePeakSize).
/// \tparam IdxToSz Invocable type mapping an Index to its extent.
/// \param network A TensorNetwork object.
/// \param idxsz An invocable on Index, that maps Index to its dimension.
/// \param subnet_cse Whether to recognize equivalent subnetworks to try
///        minimizing the ops counts.
/// \param cost Cost-model knobs (\ref CostParams): is_volatile_leaf,
///        volatile_weight, footprint_weight, peak_flops_tolerance, roofline.
///        is_volatile_leaf marks a leaf tensor as volatile (its value changes
///        on every replay); empty disables weighting. The predicate MUST be
///        invariant under slot/index canonicalization — key on tensor label or
///        structure, NOT on anonymous index identity — so that two subnetworks
///        deemed equivalent by the subnet-CSE canonicalization also agree on
///        volatility (the CSE path stores one cost per canonical subnet).
///        volatile_weight/footprint_weight apply to DenseFLOPs only;
///        roofline applies to the peak objectives only; peak_flops_tolerance
///        applies to DensePeakSize only (DensePeakSizeBatched's final
///        selection is instead threshold-gated by peak_threshold).
/// \param is_batchable_index Predicate marking an index as batchable (sliced);
///        ObjectiveFunction::DensePeakSizeBatched only.
/// \param batch_target_size Per-index per-batch slice size (an upper bound) for
///        batchable indices; ObjectiveFunction::DensePeakSizeBatched only.
/// \param inner_pow Optional k-aware CSV/PNO composite extent applied by every
///        cost counter; see \ref inner_aware_volume. Empty (default) sizes
///        composites by \p idxsz (k=1), which under-counts multi-composite
///        tensors.
/// \param batch_persistent_only When true, only persistent (volatile-leaf-free)
///        subnetworks are batched; ObjectiveFunction::DensePeakSizeBatched
///        only.
/// \param out_axes When non-null AND \p Metric ==
///        ObjectiveFunction::DensePeakSizeBatched, filled with the per-node
///        sliced-sets of the returned sequence's contraction (\c -1) nodes, in
///        the same left-first post-order the sequence itself was built in
///        (see \ref opt::detail::run_single_term_opt_axes). Cleared (left
///        empty) otherwise. Ignored (must stay null) for every other metric.
/// \return Optimal evaluation sequence under the chosen cost metric. If there
///         are equivalent optimal sequences then the result is the one that
///         keeps the order of tensors in the network as original as possible.
///
template <ObjectiveFunction Metric, has_index_extent IdxToSz>
EvalSequence single_term_opt(
    TensorNetwork const& network, IdxToSz&& idxsz, bool subnet_cse,
    CostParams const& cost = {},
    std::function<bool(Index const&)> const& is_batchable_index = {},
    std::function<std::size_t(Index const&)> batch_target_size = {},
    std::function<double(Index const&, std::size_t)> const& inner_pow = {},
    bool batch_persistent_only = false,
    container::vector<container::svector<std::pair<Index, AxisKind>>>*
        out_axes = nullptr) {
  if (out_axes) out_axes->clear();
  decltype(OptRes::indices) tidxs{};
  // Unpack the cost knobs into the names the recurrence arms below use, so the
  // arms are unchanged. (void)-cast all, since each if-constexpr arm uses only
  // a subset.
  auto const& is_volatile_leaf = cost.is_volatile_leaf;
  double const volatile_weight = cost.volatile_weight;
  double const footprint_weight = cost.footprint_weight;
  double const peak_flops_tolerance = cost.peak_flops_tolerance;
  double const accumulation_factor = cost.accumulation_factor;
  RooflineParams const& roofline = cost.roofline;
  bool const prune_outer_products = cost.prune_outer_products;
  (void)is_volatile_leaf;
  (void)volatile_weight;
  (void)footprint_weight;
  (void)peak_flops_tolerance;
  (void)accumulation_factor;
  (void)roofline;

  // Volatility weighting is a DenseFLOPs-only notion (persistent intermediates
  // cost MORE memory, not less, so it is wrong-signed for DenseSize). Build the
  // volatile-leaf bitmask in network.tensors() bit order so it aligns with the
  // DP's subset bits.
  size_t volatile_mask = 0;
  double nr = 1.0;
  if constexpr (Metric == ObjectiveFunction::DenseSpaceTime ||
                Metric == ObjectiveFunction::DenseTimeSpace) {
    SEQUANT_ASSERT(!subnet_cse &&
                   "subnet_cse not supported with DenseSpaceTime (Phase 1)");
    SEQUANT_ASSERT(!out_axes &&
                   "out_axes only supported with the batched peak objectives");
    (void)is_batchable_index;
    (void)batch_target_size;
    (void)batch_persistent_only;
    (void)footprint_weight;  // peak objectives use the roofline tie-break
    // is_volatile_leaf / volatile_weight / roofline feed only the secondary
    // tie-break among equal-peak schedules (peak itself ignores them).
    PeakModel model{idxsz,
                    inner_pow,
                    is_volatile_leaf,
                    volatile_weight,
                    roofline.machine_balance,
                    roofline.fast_mem_elems,
                    roofline.block_tiles,
                    roofline.block_prefactor,
                    peak_flops_tolerance};
    model.perf_first = (Metric == ObjectiveFunction::DenseTimeSpace);
    model.prune_outer_products = prune_outer_products;
    return run_single_term_opt(model, network, tidxs);
  } else if constexpr (Metric == ObjectiveFunction::DenseSpaceTimeBatched ||
                       Metric == ObjectiveFunction::DenseTimeSpaceBatched) {
    SEQUANT_ASSERT(
        !subnet_cse &&
        "subnet_cse not supported with DensePeakSizeBatched (Phase 2)");
    (void)footprint_weight;  // peak objectives use the roofline tie-break
    // is_volatile_leaf gates batching; volatile_weight / roofline feed the
    // secondary tie-break among equal-peak schedules.
    PeakBatchedModel model{idxsz,
                           is_batchable_index,
                           batch_target_size,
                           is_volatile_leaf,
                           inner_pow,
                           volatile_weight,
                           roofline.machine_balance,
                           roofline.fast_mem_elems,
                           roofline.block_tiles,
                           roofline.block_prefactor,
                           batch_persistent_only,
                           peak_flops_tolerance,
                           accumulation_factor,
                           cost.peak_threshold};
    model.perf_first = (Metric == ObjectiveFunction::DenseTimeSpaceBatched);
    model.prune_outer_products = prune_outer_products;
    model.batch_spectator_indices = cost.batch_spectator_indices;
    // charge_batch_recompute defaults to true (see PeakBatchedModel): the batch
    // recomputation cost is always reflected on the flops/exec axis.
    if (out_axes) {
      auto [seq, axes] = run_single_term_opt_axes(model, network, tidxs);
      *out_axes = std::move(axes);
      return seq;
    }
    return run_single_term_opt(model, network, tidxs);
  } else if constexpr (Metric == ObjectiveFunction::DenseFLOPs) {
    SEQUANT_ASSERT(!out_axes &&
                   "out_axes only supported with DensePeakSizeBatched");
    if (is_volatile_leaf && volatile_weight > 1.0) {
      size_t i = 0;
      for (auto&& t : network.tensors()) {
        auto tp = std::dynamic_pointer_cast<Tensor>(t);
        if (tp && is_volatile_leaf(*tp)) volatile_mask |= (size_t{1} << i);
        ++i;
      }
      nr = volatile_weight;
    }
    (void)is_batchable_index;
    (void)batch_target_size;
    (void)batch_persistent_only;
    (void)peak_flops_tolerance;
    AdditiveModel model{flops_counter(idxsz, inner_pow),
                        footprint_counter(idxsz, inner_pow),
                        volatile_mask,
                        nr,
                        footprint_weight,
                        subnet_cse};
    model.prune_outer_products = prune_outer_products;
    return run_single_term_opt(model, network, tidxs);
  } else {
    static_assert(Metric == ObjectiveFunction::DenseSize,
                  "Only DenseFLOPs, DenseSize, DenseSpaceTime, "
                  "DenseSpaceTimeBatched, DenseTimeSpace, and "
                  "DenseTimeSpaceBatched ObjectiveFunction supported.");
    SEQUANT_ASSERT(!out_axes &&
                   "out_axes only supported with DensePeakSizeBatched");
    (void)is_batchable_index;
    (void)batch_target_size;
    (void)batch_persistent_only;
    (void)peak_flops_tolerance;
    AdditiveModel model{memsize_counter(idxsz, inner_pow),
                        footprint_counter(idxsz, inner_pow),
                        volatile_mask,
                        nr,
                        footprint_weight,
                        subnet_cse};
    model.prune_outer_products = prune_outer_products;
    return run_single_term_opt(model, network, tidxs);
  }
}

}  // namespace detail

///
/// \tparam Metric Objective function (DenseFLOPs by default; DenseSize
///         minimizes total operand storage rather than flops; DensePeakSize
///         minimizes peak memory over the evaluation schedule -- see
///         ObjectiveFunction).
/// \param prod  Product to be optimized.
/// \param idxsz An invocable object that maps an Index object to size.
/// \return Parenthesized product expression.
///
/// @note @c prod is assumed to consist of only Tensor expressions
/// @note The remaining parameters (\c subnet_cse, \c cost,
///       \c is_batchable_index, \c batch_target_size, \c inner_pow,
///       \c batch_persistent_only) are forwarded verbatim to the detail
///       \ref single_term_opt overload; see it for their semantics.
/// \param out_axes When non-null AND \p Metric ==
///        ObjectiveFunction::DensePeakSizeBatched, filled with the per-node
///        sliced-sets of the returned Product tree's contraction nodes, in
///        the same left-first post-order the nested Product below is built
///        in (so \c (*out_axes)[j] annotates the j-th Product node formed by
///        the \c -1-handling arm of the loop below). Left empty if \p prod
///        has fewer than 3 factors (no factorization is performed) or if
///        \p Metric != DensePeakSizeBatched.
///
template <ObjectiveFunction Metric = ObjectiveFunction::DenseFLOPs,
          has_index_extent IdxToSz>
ExprPtr single_term_opt(
    Product const& prod, IdxToSz&& idxsz, bool subnet_cse = false,
    CostParams const& cost = {},
    std::function<bool(Index const&)> const& is_batchable_index = {},
    std::function<std::size_t(Index const&)> batch_target_size = {},
    std::function<double(Index const&, std::size_t)> const& inner_pow = {},
    bool batch_persistent_only = false,
    container::vector<container::svector<std::pair<Index, AxisKind>>>*
        out_axes = nullptr) {
  using ranges::views::filter;
  using ranges::views::reverse;

  if (out_axes) out_axes->clear();
  if (prod.factors().size() < 3)
    return ex<Product>(Product{prod.scalar(), prod.factors().begin(),
                               prod.factors().end(), Product::Flatten::No});
  auto const tensors =
      prod | filter(&ExprPtr::template is<Tensor>) | ranges::to_vector;
  auto seq = detail::single_term_opt<Metric>(
      TensorNetwork{tensors}, std::forward<IdxToSz>(idxsz), subnet_cse, cost,
      is_batchable_index, batch_target_size, inner_pow, batch_persistent_only,
      out_axes);
  auto result = container::svector<ExprPtr>{};
  for (auto i : seq)
    if (i == -1) {
      auto rexpr = *result.rbegin();
      result.pop_back();
      auto lexpr = *result.rbegin();
      result.pop_back();
      auto p = Product{1, ExprPtrList{lexpr, rexpr}, Product::Flatten::No};
      result.push_back(ex<Product>(Product{
          1, p.factors().begin(), p.factors().end(), Product::Flatten::No}));
    } else {
      result.push_back(tensors.at(i));
    }

  auto& p_ = (*result.rbegin()).as<Product>();
  for (auto&& v :
       prod | reverse | filter([](const ExprPtr& e) { return e->is_scalar(); }))
    p_.prepend(1, v, Product::Flatten::No);

  p_.scale(prod.scalar());
  return *result.rbegin();
}

}  // namespace sequant::opt

#endif  // SEQUANT_CORE_OPTIMIZE_SINGLE_TERM_HPP
