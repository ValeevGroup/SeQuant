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

EvalNode to_eval_node_antisymm(ExprPtr const& expr);

EvalNode to_eval_node_symm(ExprPtr const& expr);

ExprPtr to_expr(EvalNode const& node);

ExprPtr linearize_eval_node(EvalNode const& node);

AsyCostEntry asy_cost_single_node(EvalNode const& node);

AsyCostEntry asy_cost_single_node_symmetry(EvalNode const& node);

///
///  \tparam F function type that takes EvalNode const& argument and returns
///  bool.
///
/// \param node Node to compute asymptotic cost on.
///
/// \param exploit_symmetry Whether to use symmetry properties of an
/// intermediate to get reduced cost. Default: true.
///
/// \param pred pred is called
/// on every node and only those nodes that return true will be used to compute
/// cost. Default function: returns true.
///
/// \return Asymptotic cost of evaluation
/// in terms of number of occupied and virtual orbitals.
///
template <typename F = std::function<bool(EvalNode const&)>,
          std::enable_if_t<std::is_invocable_r_v<bool, F, EvalNode const&>,
                           bool> = true>
AsyCost asy_cost(
    EvalNode const& node, F&& pred = [](auto const&) { return true; },
    bool exploit_symmetry = true) {
  if (node.leaf() || !std::invoke(std::forward<F>(pred), node))
    return AsyCost::zero();

  return AsyCost{exploit_symmetry ? asy_cost_single_node_symmetry(node)
                                  : asy_cost_single_node(node)} +          //
         asy_cost(node.left(), std::forward<F>(pred), exploit_symmetry) +  //
         asy_cost(node.right(), std::forward<F>(pred), exploit_symmetry);
}

}  // namespace sequant

#endif  // SEQUANT_EVAL_NODE_HPP
