#ifndef SEQUANT_EXPR_UTILITIES_HPP
#define SEQUANT_EXPR_UTILITIES_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>

#include <string>

namespace sequant {

/// @returns A string describing (some of) the difference between the given
/// expressions. An empty diff means that they are equal. The produced diff is
/// meant to be (resonably) human-readable.
std::string diff(const Expr& lhs, const Expr& rhs);

/// @brief Applies index replacement rules to an ExprPtr
/// @param expr ExprPtr to transform
/// @param index_replacements index replacement map
/// @param scaling_factor to scale the result
/// @return a substituted and scaled expression pointer
ExprPtr transform_expr(const ExprPtr& expr,
                       const container::map<Index, Index>& index_replacements,
                       Constant::scalar_type scaling_factor = 1);

}  // namespace sequant

#endif
