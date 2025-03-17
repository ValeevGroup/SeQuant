#ifndef SEQUANT_TRANSFORM_EXPR_HPP
#define SEQUANT_TRANSFORM_EXPR_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>

namespace sequant {

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
