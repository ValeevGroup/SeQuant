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

EvalNode to_eval_node(EvalNode l, EvalNode r, EvalOp op);

EvalNode to_eval_node_antisymm(ExprPtr const& expr);

EvalNode to_eval_node_symm(ExprPtr const& expr);

ExprPtr to_expr(EvalNode const& node);

ExprPtr linearize_eval_node(EvalNode const& node);

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
    EvalNode const& node, F&& pred = [](EvalNode const&) { return true; }) {
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
    EvalNode const& node, F&& pred = [](EvalNode const&) { return true; }) {
  return detail::asy_cost_impl(node, true, std::forward<F>(pred));
}

};  // namespace sequant

#endif  // SEQUANT_EVAL_NODE_HPP
