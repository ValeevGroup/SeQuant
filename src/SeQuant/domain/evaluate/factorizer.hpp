#ifndef SEQUANT_EVALUATE_FACTORIZER_HPP
#define SEQUANT_EVALUATE_FACTORIZER_HPP

///
/// Find common sub-networks between a pair of tensor networks.
///
/// @author Bimal Gaudel
/// @version March 2020
///

#include "eval_fwd.hpp"

#include <SeQuant/core/expr_fwd.hpp>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/evaluate/eval_tree.hpp>
#include <algorithm>
#include <tuple>

namespace sequant::evaluate {
///
/// Find the largest common sub network between a pair of tensor networks.
///
/// @param exprA Expression to find common subnetwork of.
/// @param exprB Expression to find common subnetwork of.
///
/// @return Tuple of positions of the common sub-networks.
///
std::tuple<container::svector<size_t>, container::svector<size_t>>
largest_common_subnet(const ExprPtr& exprA, const ExprPtr& exprB);

/// Generate a subexpression of an expr made of those present in indices
/// positions @c idx_vec.
/// @param expr The expression to extract subexpression from.
/// @param idx_vec Vector of indices of subexpressions in @c expr.
ExprPtr _getSubExpr(const ExprPtr& expr,
                    const container::svector<size_t>& idx_vec);

}  // namespace sequant::evaluate

#endif /* ifndef SEQUANT_EVALUATE_FACTORIZER_HPP */
