#ifndef SEQUANT_OPTIMIZE_OPTIMIZE_HPP
#define SEQUANT_OPTIMIZE_OPTIMIZE_HPP

#include <cmath>
#include <cstddef>
#include <functional>
#include <iterator>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <bit>

namespace sequant {

/// Optimize an expression assuming the number of virtual orbitals
/// greater than the number of occupied orbitals.

class Tensor;

/// \param expr Expression to be optimized.
/// \return EvalNode object.
// EvalNode optimize(ExprPtr const& expr);

namespace opt {

///
/// \param idxsz An invocable that returns size_t for Index argument.
/// \param idxs Index objects.
/// \return flops count
///
template <typename F>
concept has_index_extent = std::is_invocable_r_v<size_t, F, Index const&>;

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
  container::svector<sequant::Index> indices;

  /// The flops count of evaluation
  double flops;

  /// The evaluation sequence
  EvalSequence sequence;
};

template <typename CostFn>
  requires requires(CostFn&& fn, decltype(OptRes::indices) const& ixs) {
    { std::forward<CostFn>(fn)(ixs, ixs, ixs) } -> std::floating_point;
  }
EvalSequence single_term_opt_impl(TensorNetwork const& network,
                                  meta::range_of<Index> auto const& tidxs,
                                  CostFn&& cost_fn) {
  using ranges::views::concat;
  using ranges::views::indirect;
  using ranges::views::transform;
  using IndexContainer = decltype(OptRes::indices);
  auto const nt = network.tensors().size();
  if (nt == 1) return EvalSequence{0};
  if (nt == 2) return EvalSequence{0, 1, -1};

  container::vector<OptRes> results((1 << nt));

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
      // results[i].sequence is left uninitialized
    }
  }

  // find the optimal evaluation sequence
  for (size_t n = 0; n < results.size(); ++n) {
    if (std::popcount(n) < 2) continue;
    std::pair<size_t, size_t> curr_parts{0, 0};
    for (auto& curr_cost = results[n].flops;
         auto&& [lp, rp] : bits::bipartitions(n)) {
      auto new_cost = std::forward<CostFn>(cost_fn)(
          results[lp].indices, results[rp].indices, results[n].indices);
      if (new_cost < curr_cost) {
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
EvalSequence single_term_opt(TensorNetwork const& network, IdxToSz&& idxsz) {
  auto cost_fn = flops_counter(std::forward<IdxToSz>(idxsz));
  decltype(OptRes::indices) tidxs{};
  return single_term_opt_impl(network, tidxs, cost_fn);
}

///
/// Omit the first factor from the top level product from given expression.
/// Intended to drop "A" and "S" tensors from CC amplitudes as a preparatory
/// step for evaluation of the amplitudes.
///
ExprPtr tail_factor(ExprPtr const& expr) noexcept;

///
///
/// Pulls out scalar to the top level from a nested product.
/// If @c expr is not Product, does nothing.
void pull_scalar(sequant::ExprPtr expr) noexcept;

///
/// \param prod  Product to be optimized.
/// \param idxsz An invocable object that maps an Index object to size.
/// \return Parenthesized product expression.
///
/// @note @c prod is assumed to consist of only Tensor expressions
///
template <has_index_extent IdxToSz>
ExprPtr single_term_opt(Product const& prod, IdxToSz&& idxsz) {
  using ranges::views::filter;
  using ranges::views::reverse;

  if (prod.factors().size() < 3)
    return ex<Product>(Product{prod.scalar(), prod.factors().begin(),
                               prod.factors().end(), Product::Flatten::No});
  auto const tensors =
      prod | filter(&ExprPtr::template is<Tensor>) | ranges::to_vector;
  auto seq =
      single_term_opt(TensorNetwork{tensors}, std::forward<IdxToSz>(idxsz));
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

///
/// \brief Create clusters out of positions of terms in a sum that share common
///        intermediates.
///
/// \param expr A Sum to find clusters in.
/// \return A vector of clusters (vectors of position index of terms
///         in \c expr).
///
container::vector<container::vector<size_t>> clusters(Sum const& expr);

///
/// \brief Reorder summands so that terms having common intermediates appear
///        closer.
///
/// \param sum Expression to reorder.
/// \return Expression with summands re-ordered.
///
Sum reorder(Sum const& sum);

///
/// \param expr  Expression to be optimized.
/// \param idxsz An invocable object that maps an Index object to size.
/// \param reorder_sum If true, the summands are reordered so that terms with
///                    common sub-expressions appear closer to each other.
/// \return Optimized expression for lower evaluation cost.
template <typename IdxToSize, typename = std::enable_if_t<std::is_invocable_r_v<
                                  size_t, IdxToSize, const Index&>>>
ExprPtr optimize(ExprPtr const& expr, IdxToSize const& idx2size,
                 bool reorder_sum) {
  using ranges::views::transform;
  if (expr->is<Product>()) {
    if (ranges::all_of(*expr, [](auto&& x) {
          return x->template is<Tensor>() || x->template is<Variable>();
        }))
      return opt::single_term_opt(expr->as<Product>(), idx2size);
    else {
      auto const& prod = expr->as<Product>();

      container::svector<ExprPtr> non_tensors(prod.size());
      container::svector<ExprPtr> new_factors;

      for (auto i = 0; i < prod.size(); ++i) {
        auto&& f = prod.factor(i);
        if (f.is<Tensor>() || f.is<Variable>())
          new_factors.emplace_back(f);
        else {
          non_tensors[i] = f;
          auto target_idxs = get_unique_indices(f);
          new_factors.emplace_back(
              ex<Tensor>(L"I_" + std::to_wstring(i), bra(target_idxs.bra),
                         ket(target_idxs.ket), aux(target_idxs.aux)));
        }
      }

      auto result = opt::single_term_opt(
          Product(prod.scalar(), new_factors, Product::Flatten::No), idx2size);

      auto replacer = [&non_tensors](ExprPtr& out) {
        if (!out->is<Tensor>()) return;
        auto const& tnsr = out->as<Tensor>();
        auto&& label = tnsr.label();
        if (label.at(0) == L'I' && label.at(1) == L'_') {
          size_t suffix = std::stoi(std::wstring(label.data() + 2));
          out = non_tensors[suffix].clone();
        }
      };

      result->visit(replacer, /* atoms_only = */ true);
      return result;
    }
  } else if (expr->is<Sum>()) {
    auto smands = *expr | transform([&idx2size](auto&& s) {
      return optimize(s, idx2size, /* reorder_sum= */ false);
    }) | ranges::to_vector;
    auto sum = Sum{smands.begin(), smands.end()};
    return reorder_sum ? ex<Sum>(opt::reorder(sum)) : ex<Sum>(std::move(sum));
  } else
    return expr->clone();
}

}  // namespace opt

///
/// Optimize the expression using IndexSpace::aproximate_size() for reference
/// index extent.
///
/// \param expr  Expression to be optimized.
/// \param reorder_sum If true, the summands are reordered so that terms with
///                    common sub-expressions appear closer to each other.
///                    True by default.
/// \return Optimized expression for lower evaluation cost.
ExprPtr optimize(ExprPtr const& expr, bool reorder_sum = true);

/// Optimize the expression using IndexSpace::aproximate_size() for reference
/// index extent.
///
/// \param expr  Expression to be optimized.
/// \param reorder_sum If true, the summands are reordered so that terms with
///                    common sub-expressions appear closer to each other.
///                    True by default.
/// \return Optimized expression for lower evaluation cost.
ResultExpr& optimize(ResultExpr& expr, bool reorder_sum = true);

/// Optimize the expression using IndexSpace::aproximate_size() for reference
/// index extent.
///
/// \param expr  Expression to be optimized.
/// \param reorder_sum If true, the summands are reordered so that terms with
///                    common sub-expressions appear closer to each other.
///                    True by default.
/// \return Optimized expression for lower evaluation cost.
[[nodiscard]] ResultExpr& optimize(ResultExpr&& expr, bool reorder_sum = true);

}  // namespace sequant

#endif  // SEQUANT_OPTIMIZE_OPTIMIZE_HPP
