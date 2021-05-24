#ifndef SEQUANT_UTILS_BINARIZE_EXPR_HPP
#define SEQUANT_UTILS_BINARIZE_EXPR_HPP

#include "binary_node.hpp"
#include "eval_expr.hpp"
#include "eval_seq.hpp"

namespace sequant::utils {

binary_node<eval_expr> binarize_expr(ExprPtr const& expr);

ExprPtr binarize_expr(binary_node<eval_expr> const& node);

ExprPtr linearize_expr(binary_node<eval_expr> const& node);

ExprPtr linearize_expr(ExprPtr const& expr);

}  // namespace sequant::utils

#endif  // SEQUANT_UTILS_BINARIZE_EXPR_HPP
