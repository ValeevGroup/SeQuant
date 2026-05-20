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

#include <range/v3/view/concat.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/indirect.hpp>
#include <range/v3/view/move.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/view/transform.hpp>

#include <algorithm>
#include <bit>
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
    auto tot_idxs = tot_indices<IndexSet>(concat(lhs, rhs, result));
    double total_flops = ranges::accumulate(
        concat(tot_idxs.outer, tot_idxs.inner), 1., std::multiplies{}, ixex);
    // total_flops == 1. implies zero flops
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
    for (auto&& tot_idxs :
         {tot_indices(lhs), tot_indices(rhs), tot_indices(result)}) {
      double mem = ranges::accumulate(concat(tot_idxs.outer, tot_idxs.inner),
                                      1., std::multiplies{}, ixex);
      if (mem != 1.) total_mem += mem;
      // else 1. is assumed to be the accumulation init value;
      // skip adding it to the total.
    }
    return total_mem;
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
template <typename CostFn>
  requires requires(CostFn&& fn, decltype(OptRes::indices) const& ixs) {
    { fn(ixs, ixs, ixs) } -> std::floating_point;
  }
EvalSequence single_term_opt_impl(TensorNetwork const& network,
                                  meta::range_of<Index> auto const& tidxs,
                                  CostFn&& cost_fn, bool subnet_cse) {
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
    for (auto& curr_cost = results[n].ops;
         auto&& [lp, rp] : bits::bipartitions(n)) {
      // do nothing with the trivial bipartition
      // i.e. one subset is the empty set and the other full
      if (lp == 0 || rp == 0) continue;

      double new_cost = 0;
      container::vector<size_t> combined_subnets;
      if (subnet_cse) {
        // subnets is always kept sorted; set_union requires sorted inputs and
        // produces sorted output — this invariant is maintained throughout.
        std::set_union(results[lp].subnets.begin(), results[lp].subnets.end(),
                       results[rp].subnets.begin(), results[rp].subnets.end(),
                       std::back_inserter(combined_subnets));
        new_cost = cost_fn(results[lp].indices,  //
                           results[rp].indices,  //
                           results[n].indices);
        for (auto id : combined_subnets) {
          new_cost += unique_meta_costs[id];
        }
      } else {
        new_cost = cost_fn(results[lp].indices,  //
                           results[rp].indices,  //
                           results[n].indices)   //
                   + results[lp].ops + results[rp].ops;
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
          cost_fn(results[results[n].lp].indices,
                  results[results[n].rp].indices, results[n].indices);
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

///
/// \tparam OptFor Cost metric to optimize for (Flops or Memsize).
/// \tparam IdxToSz
/// \param network A TensorNetwork object.
/// \param idxsz An invocable on Index, that maps Index to its dimension.
/// \param subnet_cse Whether to recognize equivalent subnetworks to try
/// minimizing the ops counts.
/// \return Optimal evaluation sequence under the chosen cost metric. If there
///         are equivalent optimal sequences then the result is the one that
///         keeps the order of tensors in the network as original as possible.
///
template <OptFor Metric, has_index_extent IdxToSz>
EvalSequence single_term_opt(TensorNetwork const& network, IdxToSz&& idxsz,
                             bool subnet_cse) {
  decltype(OptRes::indices) tidxs{};
  if constexpr (Metric == OptFor::Flops) {
    auto cost_fn = flops_counter(std::forward<IdxToSz>(idxsz));
    return single_term_opt_impl(network, tidxs, cost_fn, subnet_cse);
  } else {
    static_assert(Metric == OptFor::Memsize,
                  "Only Flops and Memsize OptFor supported.");
    auto cost_fn = memsize_counter(std::forward<IdxToSz>(idxsz));
    return single_term_opt_impl(network, tidxs, cost_fn, subnet_cse);
  }
}

}  // namespace detail

///
/// \tparam Metric Cost metric to optimize for (Flops by default; Memsize
///         minimizes total operand memory rather than flops).
/// \param prod  Product to be optimized.
/// \param idxsz An invocable object that maps an Index object to size.
/// \return Parenthesized product expression.
///
/// @note @c prod is assumed to consist of only Tensor expressions
///
template <OptFor Metric = OptFor::Flops, has_index_extent IdxToSz>
ExprPtr single_term_opt(Product const& prod, IdxToSz&& idxsz,
                        bool subnet_cse = false) {
  using ranges::views::filter;
  using ranges::views::reverse;

  if (prod.factors().size() < 3)
    return ex<Product>(Product{prod.scalar(), prod.factors().begin(),
                               prod.factors().end(), Product::Flatten::No});
  auto const tensors =
      prod | filter(&ExprPtr::template is<Tensor>) | ranges::to_vector;
  auto seq = detail::single_term_opt<Metric>(
      TensorNetwork{tensors}, std::forward<IdxToSz>(idxsz), subnet_cse);
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
