//
// Created by Bimal Gaudel on 5/24/21.
//

#ifndef SEQUANT_EVAL_NODE_HPP
#define SEQUANT_EVAL_NODE_HPP

#include "asy_cost.hpp"
#include "binary_node.hpp"
#include "eval_expr.hpp"

#include <boost/math/special_functions/factorials.hpp>

namespace sequant {

using EvalNode = FullBinaryNode<EvalExpr>;

EvalNode to_eval_node(ExprPtr const& expr);

[[deprecated("Unused function")]] EvalNode to_eval_node(EvalNode l, EvalNode r, EvalOp op);

[[deprecated]] EvalNode to_eval_node_antisymm(ExprPtr const& expr);

[[deprecated]] EvalNode to_eval_node_symm(ExprPtr const& expr);

template <typename Xpr,
          typename = std::enable_if_t<std::is_convertible_v<Xpr, EvalExpr>>>
ExprPtr to_expr(FullBinaryNode<Xpr> const& node) {
  auto const op = node->op();
  auto const& evxpr = *node;

  if (node.leaf()) {
    return evxpr.scalar() == Constant{1}
               ? evxpr.tensor().clone()
               : ex<Constant>(evxpr.scalar()) * evxpr.tensor().clone();
  }

  if (op == EvalOp::Prod || op == EvalOp::Symm || op == EvalOp::Antisymm) {
    auto prod = Product{};
    prod.scale(evxpr.scalar().value());

    ExprPtr lexpr = to_expr(node.left());
    ExprPtr rexpr = to_expr(node.right());

    if (lexpr->is<Tensor>())
      prod.append(1, lexpr);
    else
      prod.append(lexpr);

    if (rexpr->is<Tensor>())
      prod.append(1, rexpr);
    else
      prod.append(rexpr);

    return ex<Product>(std::move(prod));
  } else {
    assert(op == EvalOp::Sum && "unsupported operation type");
    return ex<Sum>(Sum{to_expr(node.left()), to_expr(node.right())});
  }
}

template <typename Xpr,
          typename = std::enable_if_t<std::is_convertible_v<Xpr, EvalExpr>>>
ExprPtr linearize_eval_node(FullBinaryNode<Xpr> const& node) {
  if (node.leaf()) return to_expr(node);

  ExprPtr lres = to_expr(node.left());
  ExprPtr rres = to_expr(node.right());
  if (node->op() == EvalOp::Sum) return ex<Sum>(ExprPtrList{lres, rres});

  if (node.left().leaf() && node.right().leaf()) {
    return ex<Product>(node->scalar().value(), ExprPtrList{lres, rres});
  } else if (!node.left().leaf() && !node.right().leaf()) {
    auto prod = Product(node->scalar().value(), ExprPtrList{});
    prod.append(lres);
    prod.append(rres);
    return ex<Product>(std::move(prod));
  } else if (node.left().leaf() && !node.right().leaf()) {
    auto prod = Product(node->scalar().value(), ExprPtrList{lres});
    prod.append(rres);
    return ex<Product>(std::move(prod));
  } else {  // (!node.left().leaf() && node.right().leaf())
    auto& res = lres->as<Product>();
    res.scale(node->scalar().value());
    res.append(rres);
    return lres;
  }
}

AsyCost asy_cost_single_node_symm_off(EvalNode const& node);

AsyCost asy_cost_single_node(EvalNode const& node);

namespace detail {
  ///
  ///  \tparam F function type that takes EvalNode const& argument and returns
  ///  bool.
  ///
  /// \param node Node to compute asymptotic cost on.
  ///
  /// \param exploit_symmetry Whether to use symmetry properties of an
  /// intermediate to get reduced cost.
  ///
  /// \param pred pred is called on every node and only those nodes that return
  /// true will be used to compute cost.
  ///
  /// \return Asymptotic cost of evaluation
  /// in terms of number of occupied and virtual orbitals.
  ///
  template <typename F = std::function<bool(EvalNode const&)>,
            std::enable_if_t<std::is_invocable_r_v<bool, F, EvalNode const&>,
                             bool> = true>
  AsyCost asy_cost_impl(EvalNode const& node, bool exploit_symmetry, F&& pred) {
    if (node.leaf() || !std::invoke(std::forward<F>(pred), node))
      return AsyCost::zero();

    return AsyCost{exploit_symmetry ? asy_cost_single_node(node)
                                    : asy_cost_single_node_symm_off(node)} +  //
           asy_cost_impl(node.left(), exploit_symmetry,
                         std::forward<F>(pred)) +  //
           asy_cost_impl(node.right(), exploit_symmetry, std::forward<F>(pred));
  }
}  // namespace detail

///
/// \param pred pred is called on every node and only those nodes that return
/// true will be used to compute cost. Default function: returns true.
///
template <typename F = std::function<bool(EvalNode const&)>,
          std::enable_if_t<std::is_invocable_r_v<bool, F, EvalNode const&>,
                           bool> = true>
AsyCost asy_cost_symm_off(
    EvalNode const& node, F&& pred = [](EvalNode const&) {
  return true; }) {
  return detail::asy_cost_impl(node, false, std::forward<F>(pred));
}

///
/// \param pred pred is called on every node and only those nodes that return
/// true will be used to compute cost. Default function: returns true.
///
template <typename F = std::function<bool(EvalNode const&)>,
          std::enable_if_t<std::is_invocable_r_v<bool, F, EvalNode const&>,
                           bool> = true>
AsyCost asy_cost(
    EvalNode const& node, F&& pred = [](EvalNode const&) {
  return true; }) {
  return detail::asy_cost_impl(node, true, std::forward<F>(pred));
}

}  // namespace sequant

#endif  // SEQUANT_EVAL_NODE_HPP
