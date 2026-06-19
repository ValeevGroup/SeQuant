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
/// \param batch      Shared slice size applied to every sliced batchable index.
/// \param aux_list   Ordered list of distinct batchable indices (as returned
///        by \ref batchable_index_list).
/// \return \c tables[B][n] = footprint of subset \c n under sliced-set \c B.
template <typename TIdxs, typename IdxToSz>
container::vector<container::vector<double>> sliced_footprints(
    TensorNetwork const& network, TIdxs const& tidxs, IdxToSz&& idxsz,
    std::function<bool(Index const&)> const& is_batchable, std::size_t batch,
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
          if (B & (std::size_t{1} << k)) return std::min(e, batch);
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
/// volatile, weighting disabled).  Mirrors the inline volatile-mask
/// construction inside \ref single_term_opt_impl.
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

/// \brief Finds the optimal evaluation sequence for a single-term tensor
/// contraction.
///
/// This function employs an exhaustive search using dynamic programming to
/// determine the contraction order that minimizes the total cost, as defined by
/// the provided cost function.
///
/// \tparam CostFn A function object type that computes the cost of a single
/// binary contraction.
/// Expected signature:
/// \code double(meta::range_of<Index> auto const& lhs,
///              meta::range_of<Index> auto const& rhs,
///             meta::range_of<Index> auto const& res)
/// \endcode
///
/// \param network The \ref TensorNetwork containing the tensors to be
/// contracted.
///  \param tidxs The set of indices that should remain open in the
/// final result.
///  \param cost_fn The cost model used to evaluate contractions
/// (e.g., flop count).
///  \param subnet_cse If true, enables Common Subexpression
/// Elimination (CSE) for
/// equivalent subnetworks. When enabled, the cost of
/// evaluating structurally identical subnetworks is counted
/// only once in the total cost of a contraction tree.
/// Equivalence is determined by canonicalizing the subnetwork
/// graph.
///
/// \return An \ref EvalSequence representing the optimal contraction order.
///
/// \details The optimization uses a bitmask-based dynamic programming approach
///  where each state represents a subnetwork (subset of tensors).
///  If \p subnet_cse is enabled, the algorithm precomputes canonical
///  metadata for every possible subnetwork to identify common
///  structures. This allows it to find trees that benefit from reusing
///  intermediate results, which is particularly effective for
///  expressions with repeating tensor patterns.
///
template <typename CostFn, typename FootprintFn>
  requires requires(CostFn&& fn, FootprintFn&& ffn,
                    decltype(OptRes::indices) const& ixs) {
    { fn(ixs, ixs, ixs) } -> std::floating_point;
    { ffn(ixs) } -> std::floating_point;
  }
EvalSequence single_term_opt_impl(TensorNetwork const& network,
                                  meta::range_of<Index> auto const& tidxs,
                                  CostFn&& cost_fn, bool subnet_cse,
                                  size_t volatile_mask, double volatile_weight,
                                  FootprintFn&& footprint_fn,
                                  double footprint_weight) {
  using ranges::views::concat;
  auto const nt = network.tensors().size();
  if (nt == 1) return EvalSequence{0};
  if (nt == 2) return EvalSequence{0, 1, -1};

  container::vector<OptRes> results((size_t{1} << nt));
  init_results(network, tidxs, results);

  // precompute all subnet_meta if subnet_cse is true
  // Note: the O(2^n) cost is bounded in practice — subset_target_indices above
  // asserts n <= 24, capping the number of subsets at ~16M.
  container::vector<size_t> meta_ids;
  container::vector<double> unique_meta_costs;
  if (subnet_cse) {
    auto md = build_subnet_metadata(network, results);
    meta_ids = std::move(md.meta_ids);
    unique_meta_costs = std::move(md.unique_meta_costs);
  }

  // find the optimal evaluation sequence
  for (size_t n = 0; n < results.size(); ++n) {
    if (std::popcount(n) < 2) continue;
    // A subset is volatile iff it contains any volatile leaf; the contraction
    // that forms it is then re-executed on every replay. volatile_mask == 0
    // (no predicate / DenseSize / volatile_weight<=1) makes w == 1 everywhere.
    double const w = (volatile_mask & n) ? volatile_weight : 1.0;
    for (auto& curr_cost = results[n].ops;
         auto&& [lp, rp] : bits::bipartitions(n)) {
      // do nothing with the trivial bipartition
      // i.e. one subset is the empty set and the other full
      if (lp == 0 || rp == 0) continue;

      // Per-intermediate memory-footprint penalty (storage of THIS result),
      // added once and NOT scaled by the replay weight w (peak footprint is a
      // one-time materialization cost). Zero when footprint_weight == 0.
      double const fp =
          footprint_weight != 0.0
              ? footprint_weight * footprint_fn(results[n].indices)
              : 0.0;

      double new_cost = 0;
      container::vector<size_t> combined_subnets;
      if (subnet_cse) {
        // subnets is always kept sorted; set_union requires sorted inputs and
        // produces sorted output — this invariant is maintained throughout.
        std::set_union(results[lp].subnets.begin(), results[lp].subnets.end(),
                       results[rp].subnets.begin(), results[rp].subnets.end(),
                       std::back_inserter(combined_subnets));
        new_cost = w * cost_fn(results[lp].indices,  //
                               results[rp].indices,  //
                               results[n].indices)   //
                   + fp;
        for (auto id : combined_subnets) {
          new_cost += unique_meta_costs[id];
        }
      } else {
        new_cost = w * cost_fn(results[lp].indices,  //
                               results[rp].indices,  //
                               results[n].indices)   //
                   + fp + results[lp].ops + results[rp].ops;
      }

      if (new_cost <= curr_cost) {
        curr_cost = new_cost;
        results[n].lp = lp;
        results[n].rp = rp;
        if (subnet_cse) {
          results[n].subnets = std::move(combined_subnets);
        }
      }
    }

    if (subnet_cse) {
      auto mid = meta_ids[n];
      // Canonically equivalent subnetworks share the same topology and index
      // sizes, so their cost is identical. Overwriting with a later bitmask's
      // cost is intentional and benign.
      unique_meta_costs[mid] =
          w * cost_fn(results[results[n].lp].indices,
                      results[results[n].rp].indices, results[n].indices) +
          (footprint_weight != 0.0
               ? footprint_weight * footprint_fn(results[n].indices)
               : 0.0);
      auto it = std::lower_bound(results[n].subnets.begin(),
                                 results[n].subnets.end(), mid);
      if (it == results[n].subnets.end() || *it != mid) {
        results[n].subnets.insert(it, mid);
      }
    }

    auto const& lseq = results[results[n].lp].sequence;
    auto const& rseq = results[results[n].rp].sequence;
    results[n].sequence =
        (lseq[0] < rseq[0] ? concat(lseq, rseq) : concat(rseq, lseq)) |
        ranges::to<EvalSequence>;
    results[n].sequence.push_back(-1);
  }

  return results.back().sequence;
}

/// Per-subset state for the peak DP.
struct PeakRes {
  double peak =
      std::numeric_limits<double>::max();  // min peak to build subtree
  size_t lp = 0, rp = 0;  // winning bipartition (0 for singletons)
  bool lp_first = true;   // winning evaluation order
};

/// Minimum peak memory to evaluate the whole network, plus the order.
/// All-co-resident model (realistic tensor peak): inputs are resident; a tensor
/// is live until consumed. Two build-independent tables: S[n] = result
/// footprint, L[n] = sum of leaf sizes in subset n. Fills `pr[n].peak` via:
///   peak[n] = min over (bipartition lp|rp, order) of
///     lp-first: max( L[rp] + peak[lp], S[lp] + peak[rp], S[lp]+S[rp]+S[n] )
///     rp-first: max( L[lp] + peak[rp], S[rp] + peak[lp], S[lp]+S[rp]+S[n] )
/// The L[other] term is the bystander cost: while one child evaluates, the
/// other child's inputs sit resident. S/L are build-independent and the
/// recurrence is monotone in child peaks, so the per-subset minimum is globally
/// optimal (optimal substructure).
template <typename TIdxs, typename IdxToSz>
container::vector<PeakRes> peak_dp(TensorNetwork const& network,
                                   TIdxs const& /*tidxs*/, IdxToSz&& /*idxsz*/,
                                   container::vector<double> const& S) {
  auto const nt = network.tensors().size();
  container::vector<PeakRes> pr(size_t{1} << nt);
  // L[n] = sum of leaf (singleton) sizes in subset n.
  container::vector<double> L(pr.size(), 0.0);
  for (size_t n = 0; n < pr.size(); ++n)
    for (size_t b = 0; b < nt; ++b)
      if (n & (size_t{1} << b)) L[n] += S[size_t{1} << b];
  for (size_t n = 0; n < pr.size(); ++n) {
    if (std::popcount(n) == 0) {
      pr[n].peak = 0.0;
      continue;
    }
    if (std::popcount(n) == 1) {
      pr[n].peak = S[n];  // a leaf, resident at its own size
      continue;
    }
    for (auto&& [lp, rp] : bits::bipartitions(n)) {
      if (lp == 0 || rp == 0) continue;
      double const both = S[lp] + S[rp] + S[n];
      double const lp_first =
          std::max({L[rp] + pr[lp].peak, S[lp] + pr[rp].peak, both});
      double const rp_first =
          std::max({L[lp] + pr[rp].peak, S[rp] + pr[lp].peak, both});
      double const cand = std::min(lp_first, rp_first);
      if (cand < pr[n].peak) {
        pr[n].peak = cand;
        pr[n].lp = lp;
        pr[n].rp = rp;
        pr[n].lp_first = (lp_first <= rp_first);
      }
    }
  }
  return pr;
}

/// Achieved minimum peak memory (the DensePeakSize objective value) for the
/// whole network under its optimal contraction order. Convenience wrapper that
/// runs \ref subset_footprints and \ref peak_dp and returns the full set's
/// (root subset's) peak. Used by tests to compare against the brute-force
/// oracle.
template <typename TIdxs, typename IdxToSz>
double peak_cost(TensorNetwork const& network, TIdxs const& tidxs,
                 IdxToSz&& idxsz) {
  auto S = subset_footprints(network, tidxs, idxsz);
  auto pr = peak_dp(network, tidxs, idxsz, S);
  return pr.back().peak;
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

/// \brief Multi-mode (per batchable index) batched peak-memory DP.
///
/// Implements the \c peak[n][B] recurrence of the batch-aware cost model design
/// (section 6.2, model A / all-co-resident). For a subset \c n contracted at
/// sliced-set context \c B, the batchable indices CONTRACTED at this node are
/// \c A = (open_aux[lp] | open_aux[rp]) & ~open_aux[n], gated by persistence
/// (only when \c (volatile_mask & n) == 0). The recurrence minimizes over
/// bipartition, evaluation order, and \c A' (a subset of \c A to batch here),
/// with the child context \c C = B | A'. Footprints come from
/// \ref sliced_footprints, indexed by a context restricted to what is open in
/// the sized subset (the \c & open_aux[s] mask), mirroring the oracle.
///
/// \param network        The TensorNetwork.
/// \param tidxs          Target (open) indices of the network.
/// \param idxsz          Callable mapping an Index to its full extent.
/// \param is_batchable   Predicate identifying batchable indices.
/// \param batch          Shared slice size for each sliced batchable index.
/// \param volatile_mask  Bitmask of volatile leaves (see \ref
///        leaf_volatile_mask); subset \c n is persistent iff
///        \c (volatile_mask & n) == 0.
/// \return Flat table \c pr of size \c (2^nt)*(2^m), indexed \c n*(2^m)+B, with
///         back-pointers (lp, rp, lp_first, aprime) per (n, B) for
///         reconstruction.
template <typename TIdxs, typename IdxToSz>
container::vector<BatchedRes> peak_dp_batched(
    TensorNetwork const& network, TIdxs const& tidxs, IdxToSz&& idxsz,
    std::function<bool(Index const&)> const& is_batchable, std::size_t batch,
    std::size_t volatile_mask) {
  auto const nt = network.tensors().size();
  auto aux = batchable_index_list(network, is_batchable);
  std::size_t const m = aux.size();
  std::size_t const nB = std::size_t{1} << m;
  auto tables =
      sliced_footprints(network, tidxs, idxsz, is_batchable, batch, aux);
  auto open_aux = subset_open_aux(network, tidxs, aux);

  // Flat (n, B) table.
  container::vector<BatchedRes> pr((std::size_t{1} << nt) * nB);
  auto idx = [nB](std::size_t n, std::size_t B) { return n * nB + B; };

  // Context-restricted size of subset s under sliced-set ctx: the table is
  // indexed by the part of ctx that is actually open in s (mirrors the oracle).
  auto sz = [&](std::size_t s, std::size_t ctx) {
    return tables[ctx & open_aux[s]][s];
  };
  // Per-context leaf-sum of subset s (sum of singleton sizes under ctx).
  auto Lof = [&](std::size_t s, std::size_t ctx) {
    double r = 0.0;
    for (std::size_t b = 0; b < nt; ++b)
      if (s & (std::size_t{1} << b)) r += sz(std::size_t{1} << b, ctx);
    return r;
  };

  for (std::size_t n = 0; n < (std::size_t{1} << nt); ++n) {
    if (std::popcount(n) == 0) {
      for (std::size_t B = 0; B < nB; ++B) pr[idx(n, B)].peak = 0.0;
      continue;
    }
    if (std::popcount(n) == 1) {
      for (std::size_t B = 0; B < nB; ++B) pr[idx(n, B)].peak = sz(n, B);
      continue;
    }
    for (std::size_t B = 0; B < nB; ++B) {
      auto& curr = pr[idx(n, B)];
      for (auto&& [lp, rp] : bits::bipartitions(n)) {
        if (lp == 0 || rp == 0) continue;
        // Batchable indices contracted at THIS node: open at children but not
        // at parent, gated by persistence.
        std::size_t const Acand =
            (volatile_mask & n)
                ? std::size_t{0}
                : ((open_aux[lp] | open_aux[rp]) & ~open_aux[n]);
        // Enumerate every subset A' of Acand (including the empty set).
        std::size_t Ap = Acand;
        while (true) {
          std::size_t const C = B | Ap;
          double const pl = pr[idx(lp, C)].peak;
          double const prr = pr[idx(rp, C)].peak;
          double const both = sz(lp, C) + sz(rp, C) + sz(n, B);
          double const lpf = std::max({Lof(rp, C) + pl, sz(lp, C) + prr, both});
          double const rpf = std::max({Lof(lp, C) + prr, sz(rp, C) + pl, both});
          double const cand = std::min(lpf, rpf);
          if (cand < curr.peak) {
            curr.peak = cand;
            curr.lp = lp;
            curr.rp = rp;
            curr.lp_first = (lpf <= rpf);
            curr.aprime = Ap;
          }
          if (Ap == 0) break;
          Ap = (Ap - 1) & Acand;
        }
      }
    }
  }
  return pr;
}

/// \brief Achieved minimum batched peak memory for the whole network: the
/// \c peak[root][B=0] objective of the multi-mode batched DP (section 6.2).
///
/// Builds the volatile-leaf bitmask from \p is_volatile_leaf, runs
/// \ref peak_dp_batched, and returns the root subset's peak at the empty
/// ancestor context (\c B = 0), per the design objective \c
/// peak[root][\emptyset].
///
/// \param network           The TensorNetwork.
/// \param tidxs             Target (open) indices of the network.
/// \param idxsz             Callable mapping an Index to its full extent.
/// \param is_batchable      Predicate identifying batchable indices.
/// \param batch             Shared slice size for each sliced batchable index.
/// \param is_volatile_leaf  Predicate marking a leaf as volatile; may be empty.
/// \return \c peak[root][0].
template <typename TIdxs, typename IdxToSz>
double peak_cost_batched(
    TensorNetwork const& network, TIdxs const& tidxs, IdxToSz&& idxsz,
    std::function<bool(Index const&)> const& is_batchable, std::size_t batch,
    std::function<bool(Tensor const&)> const& is_volatile_leaf) {
  std::size_t const vmask = leaf_volatile_mask(network, is_volatile_leaf);
  auto pr = peak_dp_batched(network, tidxs, idxsz, is_batchable, batch, vmask);
  auto aux = batchable_index_list(network, is_batchable);
  std::size_t const nB = std::size_t{1} << aux.size();
  std::size_t const root = (std::size_t{1} << network.tensors().size()) - 1;
  return pr[root * nB + 0].peak;
}

/// Reconstruct the EvalSequence from the peak DP back-pointers, honoring the
/// chosen evaluation order at each node (the lower-peak child is emitted last).
template <typename TIdxs, typename IdxToSz>
EvalSequence single_term_opt_peak_impl(TensorNetwork const& network,
                                       TIdxs const& tidxs, IdxToSz&& idxsz) {
  using ranges::views::concat;
  auto const nt = network.tensors().size();
  if (nt == 1) return EvalSequence{0};
  if (nt == 2) return EvalSequence{0, 1, -1};
  auto S = subset_footprints(network, tidxs, idxsz);
  auto pr = peak_dp(network, tidxs, idxsz, S);
  // bottom-up emit: a subset's sequence is (first-child seq)(second-child
  // seq)-1
  container::vector<EvalSequence> seq(pr.size());
  for (size_t n = 0; n < pr.size(); ++n) {
    if (std::popcount(n) == 1) {
      seq[n] = EvalSequence{static_cast<int>(std::countr_zero(n))};
    } else if (std::popcount(n) >= 2) {
      auto const& a = pr[n].lp_first ? seq[pr[n].lp] : seq[pr[n].rp];
      auto const& b = pr[n].lp_first ? seq[pr[n].rp] : seq[pr[n].lp];
      seq[n] = concat(a, b) | ranges::to<EvalSequence>;
      seq[n].push_back(-1);
    }
  }
  return seq.back();
}

/// \brief Reconstructs the optimal contraction \ref EvalSequence from the
/// multi-mode batched peak DP back-pointers (section 6.2, model A).
///
/// Runs \ref peak_dp_batched, then walks the chosen back-pointers from the
/// root subset at the empty ancestor context \c (root, B=0). At each node
/// \c (n, B) it reads the winning \c aprime, forms the child context
/// \c C = B | aprime, and recurses into the children at context \c C in
/// \c lp_first order (the lower-peak child is emitted last). A singleton emits
/// its tensor index; an internal node emits \c concat(first, second), -1.
///
/// \param network       The TensorNetwork.
/// \param tidxs         Target (open) indices of the network.
/// \param idxsz         Callable mapping an Index to its full extent.
/// \param is_batchable  Predicate identifying batchable indices.
/// \param batch         Shared slice size for each sliced batchable index.
/// \param is_volatile_leaf  Predicate marking a leaf as volatile; may be empty.
/// \return Optimal evaluation sequence under the batched peak objective.
template <typename TIdxs, typename IdxToSz>
EvalSequence single_term_opt_peak_batched_impl(
    TensorNetwork const& network, TIdxs const& tidxs, IdxToSz&& idxsz,
    std::function<bool(Index const&)> const& is_batchable, std::size_t batch,
    std::function<bool(Tensor const&)> const& is_volatile_leaf) {
  using ranges::views::concat;
  auto const nt = network.tensors().size();
  if (nt == 1) return EvalSequence{0};
  if (nt == 2) return EvalSequence{0, 1, -1};

  std::size_t const vmask = leaf_volatile_mask(network, is_volatile_leaf);
  auto pr = peak_dp_batched(network, tidxs, idxsz, is_batchable, batch, vmask);
  auto aux = batchable_index_list(network, is_batchable);
  std::size_t const nB = std::size_t{1} << aux.size();
  auto at = [&](std::size_t n, std::size_t B) -> BatchedRes const& {
    return pr[n * nB + B];
  };

  // Recursive back-pointer walk: emit the subsequence for subset n evaluated
  // at ancestor context B. The child context is C = B | aprime; children are
  // emitted in lp_first order (lower-peak child last).
  auto build = [&](auto&& self, std::size_t n, std::size_t B) -> EvalSequence {
    if (std::popcount(n) == 1)
      return EvalSequence{static_cast<int>(std::countr_zero(n))};
    auto const& r = at(n, B);
    std::size_t const C = B | r.aprime;
    auto first = self(self, r.lp_first ? r.lp : r.rp, C);
    auto second = self(self, r.lp_first ? r.rp : r.lp, C);
    auto seq = concat(first, second) | ranges::to<EvalSequence>;
    seq.push_back(-1);
    return seq;
  };

  std::size_t const root = (std::size_t{1} << nt) - 1;
  return build(build, root, 0);
}

/// \brief Independent memory-simulation recomputation of the chosen batched
/// reconstruction's model-A peak (the user-required full numeric reconstruction
/// check).
///
/// Runs \ref peak_dp_batched, walks the SAME chosen back-pointer tree as
/// \ref single_term_opt_peak_batched_impl, and recomputes the subtree peak by
/// DIRECT MEMORY SIMULATION rather than re-evaluating the DP's max/+ formula:
/// at a node \c (n, B) with child context \c C = B | aprime, the first child
/// (in \c lp_first order) is evaluated to completion, its result is held
/// resident while the second child is evaluated, and the all-co-resident peak
/// is the max over the schedule of the sum of live tensor sizes (leaf inputs
/// resident from the start, freed on consumption; each child's result held
/// until the parent contraction consumes it). All sizes come from
/// \ref sliced_footprints at the active context (\c tables[ctx & open_aux[s]],
/// matching the DP exactly). The returned value must EQUAL
/// \ref peak_cost_batched; a mismatch signals a bug in either the DP or the
/// reconstruction.
///
/// \param network       The TensorNetwork.
/// \param tidxs         Target (open) indices of the network.
/// \param idxsz         Callable mapping an Index to its full extent.
/// \param is_batchable  Predicate identifying batchable indices.
/// \param batch         Shared slice size for each sliced batchable index.
/// \param is_volatile_leaf  Predicate marking a leaf as volatile; may be empty.
/// \return The simulated model-A peak of the chosen reconstruction.
template <typename TIdxs, typename IdxToSz>
double reconstructed_batched_peak(
    TensorNetwork const& network, TIdxs const& tidxs, IdxToSz&& idxsz,
    std::function<bool(Index const&)> const& is_batchable, std::size_t batch,
    std::function<bool(Tensor const&)> const& is_volatile_leaf) {
  auto const nt = network.tensors().size();
  std::size_t const vmask = leaf_volatile_mask(network, is_volatile_leaf);
  auto pr = peak_dp_batched(network, tidxs, idxsz, is_batchable, batch, vmask);
  auto aux = batchable_index_list(network, is_batchable);
  std::size_t const nB = std::size_t{1} << aux.size();
  auto tables =
      sliced_footprints(network, tidxs, idxsz, is_batchable, batch, aux);
  auto open_aux = subset_open_aux(network, tidxs, aux);
  auto at = [&](std::size_t n, std::size_t B) -> BatchedRes const& {
    return pr[n * nB + B];
  };
  // Context-restricted size of subset s under sliced-set ctx (mirrors the DP).
  auto sz = [&](std::size_t s, std::size_t ctx) {
    return tables[ctx & open_aux[s]][s];
  };

  // Sum of leaf (singleton) sizes in subset s under ctx -- the inputs that are
  // resident throughout the evaluation of subtree s.
  auto Lof = [&](std::size_t s, std::size_t ctx) {
    double r = 0.0;
    for (std::size_t b = 0; b < nt; ++b)
      if (s & (std::size_t{1} << b)) r += sz(std::size_t{1} << b, ctx);
    return r;
  };

  // Simulate the peak of evaluating subtree n at ancestor context B by walking
  // the chosen back-pointers. A leaf is resident at its own size. For an
  // internal node, with child context C = B | aprime, evaluate the lp_first
  // child fully (peak: its own simulated peak), then hold its result (sized at
  // C) while evaluating the second child (whose inputs co-reside at Lof), then
  // both results co-reside while the parent result (sized at B) is formed.
  // This re-derives peak independent of the DP's max/+ formula.
  auto sim = [&](auto&& self, std::size_t n, std::size_t B) -> double {
    if (std::popcount(n) == 1) return sz(n, B);
    auto const& r = at(n, B);
    std::size_t const C = B | r.aprime;
    std::size_t const f = r.lp_first ? r.lp : r.rp;  // evaluated first
    std::size_t const s = r.lp_first ? r.rp : r.lp;  // evaluated second
    double const peak_f = self(self, f, C);
    double const peak_s = self(self, s, C);
    // While the first child evaluates, the second child's leaf inputs sit
    // resident (Lof(s, C)). While the second child evaluates, the first
    // child's result sits resident (sz(f, C)). When both results exist, the
    // parent result (sz(n, B)) is materialized alongside them.
    double const stage_first = Lof(s, C) + peak_f;
    double const stage_second = sz(f, C) + peak_s;
    double const stage_form = sz(f, C) + sz(s, C) + sz(n, B);
    return std::max({stage_first, stage_second, stage_form});
  };

  std::size_t const root = (std::size_t{1} << nt) - 1;
  if (nt == 0) return 0.0;
  return sim(sim, root, 0);
}

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
    std::size_t batch_target_size = 0) {
  decltype(OptRes::indices) tidxs{};

  // The per-intermediate footprint penalty needs idxsz too, so build a
  // footprint counter alongside the cost function (idxsz is copied into each).
  auto footprint_fn = footprint_counter(idxsz);

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
    return single_term_opt_peak_impl(network, tidxs, idxsz);
  } else if constexpr (Metric == ObjectiveFunction::DensePeakSizeBatched) {
    SEQUANT_ASSERT(
        !subnet_cse &&
        "subnet_cse not supported with DensePeakSizeBatched (Phase 2)");
    (void)volatile_weight;
    (void)footprint_weight;
    return single_term_opt_peak_batched_impl(
        network, tidxs, idxsz, is_batchable_index, batch_target_size,
        is_volatile_leaf);
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
    auto cost_fn = flops_counter(idxsz);
    return single_term_opt_impl(network, tidxs, cost_fn, subnet_cse,
                                volatile_mask, nr, footprint_fn,
                                footprint_weight);
  } else {
    static_assert(Metric == ObjectiveFunction::DenseSize,
                  "Only DenseFLOPs, DenseSize, DensePeakSize, and "
                  "DensePeakSizeBatched ObjectiveFunction supported.");
    (void)is_batchable_index;
    (void)batch_target_size;
    auto cost_fn = memsize_counter(idxsz);
    return single_term_opt_impl(network, tidxs, cost_fn, subnet_cse,
                                volatile_mask, nr, footprint_fn,
                                footprint_weight);
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
    std::size_t batch_target_size = 0) {
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
