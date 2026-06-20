#ifndef SEQUANT_CORE_OPTIMIZE_SINGLE_TERM_HPP
#define SEQUANT_CORE_OPTIMIZE_SINGLE_TERM_HPP

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/optimize/options.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/external/bliss/graph.hh>

#include <range/v3/algorithm/find.hpp>
#include <range/v3/view/concat.hpp>

#include <cstdlib>
#include <iostream>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/indirect.hpp>
#include <range/v3/view/move.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/view/transform.hpp>

#include <algorithm>
#include <bit>
#include <functional>
#include <limits>
#include <type_traits>

namespace sequant::opt {

template <typename F>
concept has_index_extent = std::is_invocable_r_v<size_t, F, Index const&>;

namespace detail {

/// \brief Cost function returning the flop count of a binary tensor
/// contraction.
///
/// The returned callable computes the product of extents over the union of
/// indices on \c lhs, \c rhs, and \c result. A product of 1 (scalar
/// contraction) is reported as zero flops.
///
/// \param ixex Invocable mapping an Index to its extent.
/// \return A callable <tt>(lhs, rhs, result) -> double</tt> yielding the
/// flop count of the contraction.
auto flops_counter(has_index_extent auto&& ixex) {
  return [ixex = std::forward<decltype(ixex)>(ixex)](
             meta::range_of<Index> auto const& lhs,
             meta::range_of<Index> auto const& rhs,
             meta::range_of<Index> auto const& result) -> double {
    using ranges::views::concat;
    // <IndexSet> is required here: concatenating the three operands repeats
    // every contracted/shared index, so it must be deduplicated before taking
    // the extent product (cf. memsize_counter, which processes each operand
    // separately and so can use the default vector container).
    auto tot_idxs = tot_indices<IndexSet>(concat(lhs, rhs, result));
    double total_flops = ranges::accumulate(
        concat(tot_idxs.outer, tot_idxs.inner), 1., std::multiplies{}, ixex);
    // A product of exactly 1. means the index set was empty (the accumulation
    // init value), i.e. a scalar contraction => zero flops. Extents are
    // integer-valued, so this equality is exact.
    return total_flops == 1. ? 0. : total_flops;
  };
}

/// \brief Cost function returning the total memory footprint of a binary
/// tensor contraction.
///
/// The returned callable sums, over \c lhs, \c rhs, and \c result, the
/// product of extents of each operand's indices. Operands whose extent
/// product is 1 contribute zero.
///
/// \param ixex Invocable mapping an Index to its extent.
/// \return A callable <tt>(lhs, rhs, result) -> double</tt> yielding the
/// summed memory size of the three operands.
auto memsize_counter(has_index_extent auto&& ixex) {
  return [ixex = std::forward<decltype(ixex)>(ixex)](
             meta::range_of<Index> auto const& lhs,
             meta::range_of<Index> auto const& rhs,
             meta::range_of<Index> auto const& result) -> double {
    using ranges::views::concat;
    double total_mem{0};
    // Each operand is sized independently, so the default (vector) container of
    // tot_indices suffices -- a single operand's index list has no duplicates,
    // unlike the concatenated set flops_counter must dedup.
    for (auto&& tot_idxs :
         {tot_indices(lhs), tot_indices(rhs), tot_indices(result)}) {
      double mem = ranges::accumulate(concat(tot_idxs.outer, tot_idxs.inner),
                                      1., std::multiplies{}, ixex);
      // mem == 1. means this operand had no indices (the accumulation init
      // value), i.e. a scalar; it contributes no memory. Same exact-equality
      // convention as flops_counter above.
      if (mem != 1.) total_mem += mem;
    }
    return total_mem;
  };
}

/// \brief Cost function returning the storage footprint (element count) of a
/// single tensor: the product of the extents of its indices.
///
/// Used to apply a per-intermediate memory-footprint penalty in single-term
/// optimization (see OptimizeOptions::footprint_weight). A scalar (no indices)
/// contributes zero.
///
/// \param ixex Invocable mapping an Index to its extent.
/// \return A callable <tt>(result) -> double</tt> yielding the element count of
/// the result tensor.
auto footprint_counter(has_index_extent auto&& ixex) {
  return [ixex = std::forward<decltype(ixex)>(ixex)](
             meta::range_of<Index> auto const& result) -> double {
    using ranges::views::concat;
    auto tot_idxs = tot_indices(result);
    double mem = ranges::accumulate(concat(tot_idxs.outer, tot_idxs.inner), 1.,
                                    std::multiplies{}, ixex);
    return mem == 1. ? 0. : mem;
  };
}

///
/// Represents a result of optimization on a range of expressions
/// for a binary evaluation
///
struct OptRes {
  /// Free indices remaining upon evaluation
  IndexSet indices;

  /// The operations count of evaluation
  double ops;

  /// The evaluation sequence
  EvalSequence sequence;

  /// Bitmask splits that resulted into this OptRes
  size_t lp = 0;
  size_t rp = 0;

  /// unique canonical subnets in the optimal tree for this bitmask
  container::vector<size_t> subnets;
};

/// \brief Footprint (dense element count) of every subset's result tensor.
///
/// \c S[n] is the product of extents of the open indices of subset \c n
/// (those remaining after contracting the tensors in \c n, given \c tidxs as
/// the final target indices). \c S[0] (empty subset) and any scalar result are
/// 0. Shared by the peak DP and its tests so both agree on per-subset sizes.
template <typename TIdxs, typename IdxToSz>
container::vector<double> subset_footprints(TensorNetwork const& network,
                                            TIdxs const& tidxs,
                                            IdxToSz&& idxsz) {
  container::vector<OptRes> results((size_t{1} << network.tensors().size()));
  init_results(network, tidxs, results);
  auto fp = footprint_counter(std::forward<IdxToSz>(idxsz));
  container::vector<double> S(results.size(), 0.0);
  for (size_t n = 0; n < results.size(); ++n)
    S[n] = (n == 0) ? 0.0 : fp(results[n].indices);
  return S;
}

/// \brief Collects the distinct batchable indices (in appearance order) across
/// all tensors in \p network.
///
/// Iterates over every tensor slot (bra, ket, and aux) and appends an index to
/// the result the first time it is seen and \p is_batchable returns true for
/// it.  The returned list assigns each index a stable bit position: index at
/// position \c k is bit \c k of a sliced-set bitmask \c B.
///
/// \param network  The TensorNetwork to scan.
/// \param is_batchable  Predicate returning true for indices in a batchable
///        space (e.g. a DF/RI auxiliary space).
/// \return Ordered, deduplicated list of batchable indices.
inline container::vector<Index> batchable_index_list(
    TensorNetwork const& network,
    std::function<bool(Index const&)> const& is_batchable) {
  container::vector<Index> aux;
  if (!is_batchable) return aux;
  for (auto&& t : network.tensors()) {
    auto tp = std::dynamic_pointer_cast<Tensor>(t);
    for (auto&& ix : ranges::views::concat(tp->bra(), tp->ket(), tp->aux()))
      if (is_batchable(ix) && ranges::find(aux, ix) == ranges::end(aux))
        aux.push_back(ix);
  }
  return aux;
}

/// \brief Footprint tables for every sliced-set B of batchable indices.
///
/// Returns a vector of \c 2^m tables, where \c m = \c aux_list.size(). Table
/// \c B is the result of \ref subset_footprints evaluated with an extent
/// function that replaces the full extent of \c aux_list[k] with
/// \c min(full_extent, batch) whenever bit \c k is set in \c B.
///
/// \param network    The TensorNetwork.
/// \param tidxs      Target (open) indices of the network.
/// \param idxsz      Callable mapping an Index to its full extent.
/// \param is_batchable  Predicate identifying batchable indices.
/// \param batch_target_size  Per-index slice-size function: a sliced batchable
///        index \c ix contributes
///        min(full_extent, batch_target_size(aux_list[k])).
/// \param aux_list   Ordered list of distinct batchable indices (as returned
///        by \ref batchable_index_list).
/// \return \c tables[B][n] = footprint of subset \c n under sliced-set \c B.
template <typename TIdxs, typename IdxToSz>
container::vector<container::vector<double>> sliced_footprints(
    TensorNetwork const& network, TIdxs const& tidxs, IdxToSz&& idxsz,
    std::function<bool(Index const&)> const& is_batchable,
    std::function<std::size_t(Index const&)> const& batch_target_size,
    container::vector<Index> const& aux_list) {
  std::size_t const m = aux_list.size();
  container::vector<container::vector<double>> tables(std::size_t{1} << m);
  for (std::size_t B = 0; B < tables.size(); ++B) {
    auto extent = [&, B](Index const& ix) -> std::size_t {
      std::size_t e = idxsz(ix);
      if (is_batchable && is_batchable(ix)) {
        auto it = ranges::find(aux_list, ix);
        if (it != ranges::end(aux_list)) {
          std::size_t k =
              static_cast<std::size_t>(it - ranges::begin(aux_list));
          if (B & (std::size_t{1} << k))
            return std::min(e, batch_target_size(ix));
        }
      }
      return e;
    };
    tables[B] = subset_footprints(network, tidxs, extent);
  }
  return tables;
}

/// \brief Bitmask of volatile leaf tensors in \p network.
///
/// Bit \c i is set if tensor \c i (in \c network.tensors() order) satisfies
/// \p is_volatile_leaf.  Returns 0 when the predicate is empty (no tensor is
/// volatile, weighting disabled).
///
/// \param network           The TensorNetwork.
/// \param is_volatile_leaf  Predicate identifying volatile leaf tensors; may
///        be empty.
/// \return Bitmask with bit i set iff tensor i is a volatile leaf.
inline std::size_t leaf_volatile_mask(
    TensorNetwork const& network,
    std::function<bool(Tensor const&)> const& is_volatile_leaf) {
  std::size_t mask = 0;
  if (!is_volatile_leaf) return mask;
  std::size_t i = 0;
  for (auto&& t : network.tensors()) {
    auto tp = std::dynamic_pointer_cast<Tensor>(t);
    if (tp && is_volatile_leaf(*tp)) mask |= (std::size_t{1} << i);
    ++i;
  }
  return mask;
}

struct SubNetHash {
  size_t operator()(
      TensorNetwork::SlotCanonicalizationMetadata const& data) const noexcept {
    return data.hash_value();
  }
};

struct SubNetEqual {
  bool operator()(
      TensorNetwork::SlotCanonicalizationMetadata const& left,
      TensorNetwork::SlotCanonicalizationMetadata const& right) const {
    return bliss::ConstGraphCmp::cmp(*left.graph, *right.graph) == 0;
  }
};

/// \brief Seeds the DP table with per-subset open indices and singleton/empty
/// ops counts.
///
/// For each subset \c i of the input tensors, computes the indices that remain
/// open after evaluating that subset and stores them in \c results[i].indices.
/// Initializes \c results[i].ops to 0 for empty and singleton subsets, and
/// to \c max as a sentinel for subsets that will later be filled in by the DP.
/// Singleton subsets also get their one-element \c sequence pre-populated.
template <typename TIdxs>
void init_results(TensorNetwork const& network, TIdxs const& tidxs,
                  container::vector<OptRes>& results) {
  using IndexContainer = decltype(OptRes::indices);
  auto tensor_indices = network.tensors()          //
                        | ranges::views::indirect  //
                        | ranges::views::transform(slots);
  auto imed_indices = subset_target_indices(tensor_indices, tidxs);
  SEQUANT_ASSERT(ranges::distance(imed_indices) == ranges::distance(results));
  for (size_t i = 0; i < results.size(); ++i) {
    results[i].indices =
        imed_indices[i] | ranges::views::move | ranges::to<IndexContainer>;
    results[i].ops = std::popcount(i) > 1
                         ? std::numeric_limits<decltype(OptRes::ops)>::max()
                         : 0;
    if (std::popcount(i) == 1)
      results[i].sequence.emplace_back(std::countr_zero(i));
    // else results[i].sequence is left uninitialized
  }
}

struct SubnetMetadata {
  /// meta_ids[n] is the canonical-subnet id of subset n, or
  /// numeric_limits<size_t>::max() for subsets with popcount < 2.
  container::vector<size_t> meta_ids;
  /// Cost of evaluating one representative of each canonical subnet id,
  /// indexed by id. Populated lazily during the DP.
  container::vector<double> unique_meta_costs;
};

/// \brief Precomputes canonical-subnet identifiers for every subset of size
/// >= 2 so that structurally equivalent subnetworks share a CSE id.
///
/// Builds a `TensorNetwork` for each subset, canonicalizes it, and assigns a
/// dense integer id to each distinct canonical form. The returned
/// `unique_meta_costs` is sized to the number of distinct ids and zero-filled;
/// it is populated during the DP as each canonical subnet's optimal cost
/// becomes known.
///
/// Side effect: `results[n].indices` may be reordered by `canonicalize_slots`.
inline SubnetMetadata build_subnet_metadata(
    TensorNetwork const& network, container::vector<OptRes>& results) {
  SubnetMetadata out;
  // Use max as sentinel for entries with popcount < 2 (singletons/empty),
  // which are skipped below and never assigned a real meta ID.
  out.meta_ids.resize(results.size(), std::numeric_limits<size_t>::max());
  container::unordered_map<TensorNetwork::SlotCanonicalizationMetadata, size_t,
                           SubNetHash, SubNetEqual>
      meta_to_id;

  for (size_t n = 0; n < results.size(); ++n) {
    if (std::popcount(n) < 2) continue;
    auto ts = bits::on_bits_index(n) | bits::sieve(network.tensors());
    container::vector<ExprPtr> ts_expr;
    for (auto&& t : ts)
      ts_expr.emplace_back(std::dynamic_pointer_cast<Tensor>(t)->clone());

    auto tn = TensorNetwork{ts_expr};
    auto meta = tn.canonicalize_slots(
        TensorCanonicalizer::cardinal_tensor_labels(), &results[n].indices);

    auto [it, inserted] = meta_to_id.try_emplace(std::move(meta), 0);
    if (inserted) it->second = meta_to_id.size() - 1;

    out.meta_ids[n] = it->second;
  }
  out.unique_meta_costs.resize(meta_to_id.size(), 0.0);
  return out;
}

/// \brief Per-subset bitmask of batchable indices that are OPEN in that subset.
///
/// For each subset \c n of the input tensors, bit \c k of \c open_aux[n] is set
/// iff \c aux_list[k] is among the open (external) indices of subset \c n
/// (those that remain after contracting \c n's tensors, with \c tidxs as the
/// final targets). Used by the multi-mode batched DP and oracle to restrict the
/// sliced-set context to indices actually open in a sized subset, so that table
/// lookups (\ref sliced_footprints) are indexed consistently on both sides.
///
/// \param network   The TensorNetwork.
/// \param tidxs     Target (open) indices of the network.
/// \param aux_list  Ordered list of distinct batchable indices (as returned by
///        \ref batchable_index_list); index \c k maps to bit \c k.
/// \return \c open_aux[n] for every subset \c n.
template <typename TIdxs>
container::vector<std::size_t> subset_open_aux(
    TensorNetwork const& network, TIdxs const& tidxs,
    container::vector<Index> const& aux_list) {
  container::vector<OptRes> results(
      (std::size_t{1} << network.tensors().size()));
  init_results(network, tidxs, results);
  container::vector<std::size_t> open_aux(results.size(), 0);
  for (std::size_t n = 0; n < results.size(); ++n) {
    for (std::size_t k = 0; k < aux_list.size(); ++k) {
      if (ranges::find(results[n].indices, aux_list[k]) !=
          ranges::end(results[n].indices))
        open_aux[n] |= (std::size_t{1} << k);
    }
  }
  return open_aux;
}

/// Per-(subset, sliced-set) state for the multi-mode batched peak DP.
/// Indexed in the flat table by \c n*(2^m)+B (subset \c n, sliced-set \c B).
struct BatchedRes {
  double peak = std::numeric_limits<double>::max();  // min peak for (n, B)
  std::size_t lp = 0, rp = 0;  // winning bipartition (0 for singletons)
  bool lp_first = true;        // winning evaluation order
  std::size_t aprime = 0;  // winning subset of batchable indices sliced here
};

}  // namespace detail
}  // namespace sequant::opt

// The additive arms of single_term_opt<Metric> below are routed through the
// generic CostModel driver. cost_model.hpp includes this header (for the DP
// helpers + EvalSequence + OptRes); the include guards make the cycle a no-op,
// and every helper cost_model.hpp needs is defined above this point.
#include <SeQuant/core/optimize/cost_model.hpp>

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
/// \param is_volatile_leaf Predicate marking a leaf tensor as volatile (its
///        value changes on every replay); empty disables weighting. The
///        predicate MUST be invariant under slot/index canonicalization — key
///        on tensor label or structure, NOT on anonymous index identity — so
///        that two subnetworks deemed equivalent by the subnet-CSE
///        canonicalization also agree on volatility (the CSE path stores one
///        cost per canonical subnet). ObjectiveFunction::DenseFLOPs only.
/// \param volatile_weight Multiplier applied to the cost of each
///        volatile-result contraction (volatile contractions are re-evaluated
///        on every replay); persistent contractions are counted once. 1
///        (default) disables weighting.
/// \param footprint_weight Per-intermediate storage-footprint penalty added to
///        the cost (ObjectiveFunction::DenseFLOPs only); 0 (default) disables.
/// \return Optimal evaluation sequence under the chosen cost metric. If there
///         are equivalent optimal sequences then the result is the one that
///         keeps the order of tensors in the network as original as possible.
///
template <ObjectiveFunction Metric, has_index_extent IdxToSz>
EvalSequence single_term_opt(
    TensorNetwork const& network, IdxToSz&& idxsz, bool subnet_cse,
    std::function<bool(Tensor const&)> const& is_volatile_leaf = {},
    double volatile_weight = 1.0, double footprint_weight = 0.0,
    std::function<bool(Index const&)> const& is_batchable_index = {},
    std::function<std::size_t(Index const&)> batch_target_size = {}) {
  decltype(OptRes::indices) tidxs{};

  // Volatility weighting is a DenseFLOPs-only notion (persistent intermediates
  // cost MORE memory, not less, so it is wrong-signed for DenseSize). Build the
  // volatile-leaf bitmask in network.tensors() bit order so it aligns with the
  // DP's subset bits.
  size_t volatile_mask = 0;
  double nr = 1.0;
  if constexpr (Metric == ObjectiveFunction::DensePeakSize) {
    SEQUANT_ASSERT(!subnet_cse &&
                   "subnet_cse not supported with DensePeakSize (Phase 1)");
    (void)is_volatile_leaf;
    (void)volatile_weight;
    (void)footprint_weight;
    (void)is_batchable_index;
    (void)batch_target_size;
    return run_single_term_opt(PeakModel{idxsz}, network, tidxs);
  } else if constexpr (Metric == ObjectiveFunction::DensePeakSizeBatched) {
    SEQUANT_ASSERT(
        !subnet_cse &&
        "subnet_cse not supported with DensePeakSizeBatched (Phase 2)");
    (void)volatile_weight;
    (void)footprint_weight;
    if (std::getenv("SEQUANT_PEAK_DEBUG"))
      peak_batched_debug(network, tidxs, idxsz, is_batchable_index,
                         batch_target_size, is_volatile_leaf, std::cout);
    return run_single_term_opt(
        PeakBatchedModel{idxsz, is_batchable_index, batch_target_size,
                         is_volatile_leaf},
        network, tidxs);
  } else if constexpr (Metric == ObjectiveFunction::DenseFLOPs) {
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
    AdditiveModel model{flops_counter(idxsz), footprint_counter(idxsz),
                        volatile_mask,        nr,
                        footprint_weight,     subnet_cse};
    return run_single_term_opt(model, network, tidxs);
  } else {
    static_assert(Metric == ObjectiveFunction::DenseSize,
                  "Only DenseFLOPs, DenseSize, DensePeakSize, and "
                  "DensePeakSizeBatched ObjectiveFunction supported.");
    (void)is_batchable_index;
    (void)batch_target_size;
    AdditiveModel model{memsize_counter(idxsz), footprint_counter(idxsz),
                        volatile_mask,          nr,
                        footprint_weight,       subnet_cse};
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
///
template <ObjectiveFunction Metric = ObjectiveFunction::DenseFLOPs,
          has_index_extent IdxToSz>
ExprPtr single_term_opt(
    Product const& prod, IdxToSz&& idxsz, bool subnet_cse = false,
    std::function<bool(Tensor const&)> const& is_volatile_leaf = {},
    double volatile_weight = 1.0, double footprint_weight = 0.0,
    std::function<bool(Index const&)> const& is_batchable_index = {},
    std::function<std::size_t(Index const&)> batch_target_size = {}) {
  using ranges::views::filter;
  using ranges::views::reverse;

  if (prod.factors().size() < 3)
    return ex<Product>(Product{prod.scalar(), prod.factors().begin(),
                               prod.factors().end(), Product::Flatten::No});
  auto const tensors =
      prod | filter(&ExprPtr::template is<Tensor>) | ranges::to_vector;
  auto seq = detail::single_term_opt<Metric>(
      TensorNetwork{tensors}, std::forward<IdxToSz>(idxsz), subnet_cse,
      is_volatile_leaf, volatile_weight, footprint_weight, is_batchable_index,
      batch_target_size);
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
