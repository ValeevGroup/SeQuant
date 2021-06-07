#ifndef SEQUANT_OPTIMIZE_OPTIMIZE_HPP
#define SEQUANT_OPTIMIZE_OPTIMIZE_HPP

#include <ios>
#include <utility>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/eval_seq.hpp>

namespace sequant::optimize {

/// Optimize an expression assuming the number of virtual orbitals
/// greater than the number of occupied orbitals.

/// \param expr Expression to be optimized.
/// \param canonize Whether to canonicalize the expression(s) during
///                 optimization. Leads to slightly(?) more common
///                 sub-expressions in the optimized result, at the cost of
///                 relabeling of the indices in the result making them
///                 different from the original.
/// \return EvalNode object.
EvalNode optimize(ExprPtr const& expr, bool canonize);

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
/// Result of the single term optimization of a term.
/// Holds operations count.
///
/// Iterable of one or more binary_expr<EvalExpr> root node pointers
/// that lead to the same operations count.
///
/// ie. degenerate evaluations leading to the minimal operations count are
/// stored as binary tree nodes.
///
struct STOResult {
  AsyCost cost;

  container::vector<EvalNode> optimal_seqs;
};

/// Perform single term optimization on a product.

/// @param canon whether to canonicalize each product before couting flops.
///              by canonicalizing before counting flops, we increase the
///              chance of encountering an intermediate whose hash value is
///              already present in @c imed_hash.
/// @return STOResult
template <typename F = std::function<bool(EvalNode const&)>,
          std::enable_if_t<std::is_invocable_r_v<bool, F, EvalNode const&>,
                           bool> = true>
STOResult single_term_opt(
    Product const& prod, bool canon,
    F&& pred = [](auto const&) { return true; }) {
  using ranges::to_vector;
  using ranges::views::iota;
  using ranges::views::take;
  using ranges::views::transform;
  using seq_t = EvalSeq<size_t>;

  struct {
    Product const& facs;

    ExprPtr operator()(size_t pos) { return facs.factor(pos)->clone(); }

    ExprPtr operator()(ExprPtr lf, ExprPtr rf) {
      auto p = Product{};
      p.append(std::move(lf));
      p.append(std::move(rf));

      return ex<Product>(std::move(p));
    }
  } fold_prod{prod};  // struct

  auto result = STOResult{AsyCost::max(), {}};

  auto finder = [&result, &fold_prod, &pred, prod, canon](auto const& seq) {
    auto expr = seq.evaluate(fold_prod);
    if (canon) {
      expr->canonicalize();
      pull_scalar(expr);
    }

    if (prod.scalar() != 1.) {
      if (!expr->template is<Product>())  // in case expr is non-product
        expr = ex<Product>(Product{std::move(expr)});
      *expr *= Constant{prod.scalar()};
    }

    auto node = to_eval_node(expr);

    auto cost = asy_cost(node, std::forward<F>(pred));

    if (cost == result.cost) {
      result.optimal_seqs.emplace_back(std::move(node));
    } else if (cost < result.cost) {
      result.optimal_seqs.clear();
      result.optimal_seqs.emplace_back(std::move(node));
      result.cost = std::move(cost);
    } else {
      // cost > optimal cost. do nothing.
    }
  };  // finder

  auto init_seq = iota(size_t{0}) | take(prod.size()) |
                  transform([](auto x) { return seq_t{x}; }) |
                  ranges::to_vector;

  enumerate_eval_seq<size_t>(init_seq, finder);
  return result;
}

}  // namespace sequant::optimize

#endif  // SEQUANT_OPTIMIZE_OPTIMIZE_HPP
