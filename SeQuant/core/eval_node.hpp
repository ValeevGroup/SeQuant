//
// Created by Bimal Gaudel on 5/24/21.
//

#ifndef SEQUANT_EVAL_NODE_HPP
#define SEQUANT_EVAL_NODE_HPP

#include "asy_cost.hpp"
#include "binary_node.hpp"
#include "eval_expr.hpp"

namespace sequant {

using EvalNode = BinaryNode<EvalExpr>;

EvalNode to_eval_node(ExprPtr const& expr);

EvalNode to_eval_node_antisymm(ExprPtr const& expr);

EvalNode to_eval_node_symm(ExprPtr const& expr);

ExprPtr to_expr(EvalNode const& node);

ExprPtr linearize_eval_node(EvalNode const& node);

AsyCost asy_cost_single_node(EvalNode const& node);

template <typename F = std::function<bool(EvalNode const&)>,
          std::enable_if_t<std::is_invocable_r_v<bool, F, EvalNode const&>,
                           bool> = true>
AsyCost asy_cost(
    EvalNode const& node, F&& pred = [](auto const&) { return true; }) {

  if (node.leaf() || !std::invoke(std::forward<F>(pred), node))
    return AsyCost::zero();

  return asy_cost_single_node(node) +                    //
         asy_cost(node.left(), std::forward<F>(pred)) +  //
         asy_cost(node.right(), std::forward<F>(pred));
}

}  // namespace sequant

#endif  // SEQUANT_EVAL_NODE_HPP
