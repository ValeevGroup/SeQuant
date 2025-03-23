#ifndef SEQUANT_EXPR_UTILITIES_HPP
#define SEQUANT_EXPR_UTILITIES_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>

#include <optional>
#include <string>
#include <string_view>

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
[[nodiscard]] ExprPtr transform_expr(
    const ExprPtr& expr, const container::map<Index, Index>& index_replacements,
    Constant::scalar_type scaling_factor = 1);

/// @brief Searches for tensors with the given label and removes them from the
/// given expression Note: The function assumes that there don't exist multiple
/// tensors of that name that differ in their indexing.
///
/// @param expression The expression to modify
/// @param label The label of the tensor that shall be removed
/// @returns The removed tensor, if any occurrance has been found
std::optional<ExprPtr> pop_tensor(ExprPtr& expression, std::wstring_view label);

}  // namespace sequant

#endif
