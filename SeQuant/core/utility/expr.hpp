#ifndef SEQUANT_EXPR_UTILITIES_HPP
#define SEQUANT_EXPR_UTILITIES_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/expr_matcher.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/tensor.hpp>

#include <range/v3/algorithm/equal.hpp>
#include <range/v3/algorithm/find.hpp>
#include <range/v3/view/concat.hpp>
#include <range/v3/view/enumerate.hpp>

#include <algorithm>
#include <cassert>
#include <concepts>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

namespace sequant {

/// @returns A string describing (some of) the difference between the given
/// expressions. An empty diff means that they are equal. The produced diff is
/// meant to be (reasonably) human-readable.
std::string diff(const Expr &lhs, const Expr &rhs);

/// Checks whether the given expression is valid (i.e. uses
/// consistent indexing etc.)
/// @param expr The expression to validate
/// @param msg If given, the function will set the string to a message
///            describing why the provided expression is considered to
///            be invalid. If the expression is valid, the string will
///            be left unchanged.
/// @returns The validity of the expression
bool is_valid(const ExprPtr &expr, std::string *msg = nullptr);

/// Checks whether the given expression is valid (i.e. uses
/// consistent indexing etc.)
/// @param expr The expression to validate
/// @param msg If given, the function will set the string to a message
///            describing why the provided expression is considered to
///            be invalid. If the expression is valid, the string will
///            be left unchanged.
/// @returns The validity of the expression
bool is_valid(const Expr &expr, std::string *msg = nullptr);

/// Checks whether the given expression is valid (i.e. uses
/// consistent indexing etc.)
/// @param expr The expression to validate
/// @param msg If given, the function will set the string to a message
///            describing why the provided expression is considered to
///            be invalid. If the expression is valid, the string will
///            be left unchanged.
/// @returns The validity of the expression
bool is_valid(const ResultExpr &expr, std::string *msg = nullptr);

/// @brief Applies index replacement rules to an ExprPtr
/// @param expr ExprPtr to transform
/// @param index_replacements index replacement map
/// @param scaling_factor to scale the result
/// @return a substituted and scaled expression pointer
[[nodiscard]] ExprPtr transform_expr(
    const ExprPtr &expr, const container::map<Index, Index> &index_replacements,
    Constant::scalar_type scaling_factor = 1);
[[nodiscard]] ExprPtr transform_expr(
    const Expr &expr, const container::map<Index, Index> &index_replacements,
    Constant::scalar_type scaling_factor = 1);

/// @brief Searches for tensors with the given label and removes them from the
/// given expression Note: The function assumes that there don't exist multiple
/// tensors of that name that differ in their indexing.
///
/// @param expression The expression to modify
/// @param label The label of the tensor that shall be removed
/// @returns The removed tensor, if any occurrence has been found
std::optional<ExprPtr> pop_tensor(ExprPtr &expression, std::wstring_view label);

/// Replaces a given target expression by a given replacement
///
/// If target and replacement have common indices, the indices in the
/// replacement will be updated for each individual match of target. This is
/// relevant in cases the provided comparator doesn't account for index
/// equality. For instance, in t{a1;a2} -> r{a1;a5} applied to var * t{a3;a4}
/// the replacement shares the index a1 with the target. If the target gets
/// matched to t{a3;a4}, the replacement will be adapted via a1 -> a3, whereas
/// the index a5 will be left as is. Overall, this would lead to var * r{a3;a5}.
///
/// @param expr The expression to perform the replacements in
/// @param target The target expression to be replaced
/// @param replacement The expression to replace the target with
/// @returns A reference to expr, which has been modified in-place (useful for
/// chaining)
///
/// @note At this time, target must not be a composite expression
ExprPtr &replace(ExprPtr &expr, const ExprMatcher &target,
                 const Expr &replacement);

template <typename EqualityComparator = std::equal_to<>>
[[deprecated(
    "This is only a backwards-compat shim. Use overload using "
    "ExprMatcher instead")]] ExprPtr &
replace(ExprPtr &expr, const ExprPtr &target, const ExprPtr &replacement,
        EqualityComparator = {}) {
  ExprMatcherOptions options{.cross_comparisons = true};
  if constexpr (std::same_as<std::remove_cvref_t<EqualityComparator>,
                             std::equal_to<>>) {
    options.tensor_cmp = TensorComparison::Identity;
  } else if constexpr (std::same_as<std::remove_cvref_t<EqualityComparator>,
                                    TensorBlockEqualComparator>) {
    options.tensor_cmp = TensorComparison::Block;
  } else {
    static_assert(false,
                  "Compatibility shim can't deal with the provided comparator");
  }

  return replace(expr, ExprMatcher(*target, std::move(options)), *replacement);
}

/// Replaces a given target expression by a given replacement. Result indices
/// are adapted as needed.
///
/// If target and replacement have common indices, the indices in the
/// replacement will be updated for each individual match of target. This is
/// relevant in cases the provided comparator doesn't account for index
/// equality. For instance, in t{a1;a2} -> r{a1;a5} applied to var * t{a3;a4}
/// the replacement shares the index a1 with the target. If the target gets
/// matched to t{a3;a4}, the replacement will be adapted via a1 -> a3, whereas
/// the index a5 will be left as is. Overall, this would lead to var * r{a3;a5}.
///
/// @param expr The expression to perform the replacements in
/// @param target The target expression to be replaced
/// @param replacement The expression to replace the target with
/// @returns A reference to expr, which has been modified in-place (useful for
/// chaining)
///
/// @note At this time, target must not be a composite expression
ResultExpr &replace(ResultExpr &expr, const ExprMatcher &target,
                    const Expr &replacement);

template <typename EqualityComparator = std::equal_to<>>
[[deprecated(
    "This is only a backwards-compat shim. Use overload using "
    "ExprMatcher instead")]] ResultExpr &
replace(ResultExpr &expr, const ExprPtr &target, const ExprPtr &replacement,
        EqualityComparator = {}) {
  ExprMatcherOptions options{.cross_comparisons = true};
  if constexpr (std::same_as<std::remove_cvref_t<EqualityComparator>,
                             std::equal_to<>>) {
    options.tensor_cmp = TensorComparison::Identity;
  } else if constexpr (std::same_as<std::remove_cvref_t<EqualityComparator>,
                                    TensorBlockEqualComparator>) {
    options.tensor_cmp = TensorComparison::Block;
  } else {
    static_assert(false,
                  "Compatibility shim can't deal with the provided comparator");
  }

  return replace(expr, ExprMatcher(*target, std::move(options)), *replacement);
}

}  // namespace sequant

#endif
