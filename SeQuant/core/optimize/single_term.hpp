#ifndef SEQUANT_CORE_OPTIMIZE_SINGLE_TERM_HPP
#define SEQUANT_CORE_OPTIMIZE_SINGLE_TERM_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/external/bliss/graph.hh>

#include <range/v3/view.hpp>

#include <bit>
#include <type_traits>
#include <unordered_set>
#include "SeQuant/core/algorithm.hpp"
#include "SeQuant/core/tensor_canonicalizer.hpp"

namespace sequant::opt {
///
/// \param idxsz An invocable that returns size_t for Index argument.
/// \param idxs Index objects.
/// \return flops count
///
template <typename F>
concept has_index_extent = std::is_invocable_r_v<size_t, F, Index const&>;
namespace detail {

auto constexpr flops_counter(has_index_extent auto&& ixex) {
  return [ixex](meta::range_of<Index> auto const& lhs,
                meta::range_of<Index> auto const& rhs,
                meta::range_of<Index> auto const& result) {
    using ranges::views::concat;
    auto tot_idxs = tot_indices<container::set<Index, Index::LabelCompare>>(
        concat(lhs, rhs, result));
    double ops = 1.;
    for (auto&& idx : concat(tot_idxs.outer, tot_idxs.inner)) ops *= ixex(idx);
    // ops == 1 implies zero flops
    return ops == 1. ? 0. : ops;
  };
}

///
/// Represents a result of optimization on a range of expressions
/// for a binary evaluation
///
struct OptRes {
  /// Free indices remaining upon evaluation
  IndexSet indices;

  /// The flops count of evaluation
  double flops;

  /// The evaluation sequence
  EvalSequence sequence;
};

// constexpr auto cse_hasher = [](TNMeta const& data) -> size_t {
//   return data.hash_value();
// };

// constexpr auto cse_equal = [](TNMeta const& left,
//                               TNMeta const& right) -> bool {
//   return bliss::ConstGraphCmp::cmp(*left.graph, *right.graph) == 0;
// };

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

template <typename CostFn>
  requires requires(CostFn&& fn, decltype(OptRes::indices) const& ixs) {
    { std::forward<CostFn>(fn)(ixs, ixs, ixs) } -> std::floating_point;
  }
EvalSequence single_term_opt_impl(TensorNetwork const& network,
                                  meta::range_of<Index> auto const& tidxs,
                                  CostFn&& cost_fn, bool subnet_cse) {
  using ranges::views::concat;
  using ranges::views::indirect;
  using ranges::views::transform;
  using IndexContainer = decltype(OptRes::indices);
  auto const nt = network.tensors().size();
  if (nt == 1) return EvalSequence{0};
  if (nt == 2) return EvalSequence{0, 1, -1};

  using CanonSubnets =
      std::unordered_set<TensorNetwork::SlotCanonicalizationMetadata,
                         SubNetHash, SubNetEqual>;

  container::vector<OptRes> results((size_t{1} << nt));

  // initialize the intermediate results
  {
    auto tensor_indices = network.tensors()  //
                          | indirect         //
                          | transform(slots);
    auto imed_indices = subset_target_indices(tensor_indices, tidxs);
    SEQUANT_ASSERT(ranges::distance(imed_indices) == ranges::distance(results));
    for (size_t i = 0; i < results.size(); ++i) {
      results[i].indices =
          imed_indices[i] | ranges::views::move | ranges::to<IndexContainer>;
      results[i].flops =
          std::popcount(i) > 1
              ? std::numeric_limits<decltype(OptRes::flops)>::max()
              : 0;
      if (std::popcount(i) == 1)
        results[i].sequence.emplace_back(std::countr_zero(i));
      // else results[i].sequence is left uninitialized
    }
  }

  {
    container::vector<CanonSubnets> subnets_for_cse((size_t{1} << nt));
    for (size_t n = 0; n < subnets_for_cse.size(); ++n) {
      if (std::popcount(n) < 2) continue;
      auto ts = bits::on_bits_index(n) | bits::sieve(network.tensors());
      auto tn = TensorNetwork{ts};
      subnets_for_cse[n].emplace(tn.canonicalize_slots(
          TensorCanonicalizer::cardinal_tensor_labels(), &results[n].indices));
    }
  }

  // find the optimal evaluation sequence
  for (size_t n = 0; n < results.size(); ++n) {
    if (std::popcount(n) < 2) continue;
    std::pair<size_t, size_t> curr_parts{0, 0};
    for (auto& curr_cost = results[n].flops;
         auto&& [lp, rp] : bits::bipartitions(n)) {
      // do nothing with the trivial bipartition
      // i.e. one subset is the empty set and the other full
      if (lp == 0 || rp == 0) continue;
      auto new_cost = std::forward<CostFn>(cost_fn)(results[lp].indices,  //
                                                    results[rp].indices,  //
                                                    results[n].indices)   //
                      + results[lp].flops + results[rp].flops;
      if (new_cost <= curr_cost) {
        curr_cost = new_cost;
        curr_parts = decltype(curr_parts){lp, rp};
      }
    }
    auto const& lseq = results[curr_parts.first].sequence;
    auto const& rseq = results[curr_parts.second].sequence;
    results[n].sequence =
        (lseq[0] < rseq[0] ? concat(lseq, rseq) : concat(rseq, lseq)) |
        ranges::to<EvalSequence>;
    results[n].sequence.push_back(-1);
  }

  return results.back().sequence;
}

///
/// \tparam IdxToSz
/// \param network A TensorNetwork object.
/// \param idxsz An invocable on Index, that maps Index to its dimension.
/// \return Optimal evaluation sequence that minimizes flops. If there are
///         equivalent optimal sequences then the result is the one that keeps
///         the order of tensors in the network as original as possible.
///
template <has_index_extent IdxToSz>
EvalSequence single_term_opt(TensorNetwork const& network, IdxToSz&& idxsz,
                             bool subnet_cse) {
  auto cost_fn = flops_counter(std::forward<IdxToSz>(idxsz));
  decltype(OptRes::indices) tidxs{};
  return single_term_opt_impl(network, tidxs, cost_fn, subnet_cse);
}

}  // namespace detail

///
/// \param prod  Product to be optimized.
/// \param idxsz An invocable object that maps an Index object to size.
/// \return Parenthesized product expression.
///
/// @note @c prod is assumed to consist of only Tensor expressions
///
template <has_index_extent IdxToSz>
ExprPtr single_term_opt(Product const& prod, IdxToSz&& idxsz,
                        bool subnet_cse = false) {
  using ranges::views::filter;
  using ranges::views::reverse;

  if (prod.factors().size() < 3)
    return ex<Product>(Product{prod.scalar(), prod.factors().begin(),
                               prod.factors().end(), Product::Flatten::No});
  auto const tensors =
      prod | filter(&ExprPtr::template is<Tensor>) | ranges::to_vector;
  auto seq = detail::single_term_opt(TensorNetwork{tensors},
                                     std::forward<IdxToSz>(idxsz), subnet_cse);
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
  for (auto&& v : prod | reverse | filter(&Expr::template is<Variable>))
    p_.prepend(1, v, Product::Flatten::No);

  p_.scale(prod.scalar());
  return *result.rbegin();
}

}  // namespace sequant::opt

#endif  // SEQUANT_CORE_OPTIMIZE_SINGLE_TERM_HPP
