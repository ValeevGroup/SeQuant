#ifndef SEQUANT_EXPR_UTILITIES_HPP
#define SEQUANT_EXPR_UTILITIES_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/indices.hpp>

#include <range/v3/view/concat.hpp>

#include <algorithm>
#include <cassert>
#include <optional>
#include <string>
#include <string_view>
#include <utility>

namespace sequant {

/// @returns A string describing (some of) the difference between the given
/// expressions. An empty diff means that they are equal. The produced diff is
/// meant to be (resonably) human-readable.
std::string diff(const Expr &lhs, const Expr &rhs);

/// @brief Applies index replacement rules to an ExprPtr
/// @param expr ExprPtr to transform
/// @param index_replacements index replacement map
/// @param scaling_factor to scale the result
/// @return a substituted and scaled expression pointer
[[nodiscard]] ExprPtr transform_expr(
    const ExprPtr &expr, const container::map<Index, Index> &index_replacements,
    Constant::scalar_type scaling_factor = 1);

/// @brief Searches for tensors with the given label and removes them from the
/// given expression Note: The function assumes that there don't exist multiple
/// tensors of that name that differ in their indexing.
///
/// @param expression The expression to modify
/// @param label The label of the tensor that shall be removed
/// @returns The removed tensor, if any occurrance has been found
std::optional<ExprPtr> pop_tensor(ExprPtr &expression, std::wstring_view label);

template <typename EqualityComparator = std::equal_to<>>
ExprPtr &replace(ExprPtr &expr, const ExprPtr &target,
                 const ExprPtr &replacement, EqualityComparator cmp = {}) {
  if (!target->is_atom()) {
    throw std::runtime_error(
        "Replacement of composite expressions is not yet implemented");
  }

  if (cmp(*expr, *target)) {
    expr = replacement->clone();
  } else {
    expr->visit(
        [&](ExprPtr &current) {
          if (cmp(*current, *target)) {
            current = replacement->clone();
          }
        },
        /*only_atoms*/ true);
  }

  return expr;
}

template <typename EqualityComparator = std::equal_to<>>
ResultExpr &replace(ResultExpr &expr, const ExprPtr &target,
                    const ExprPtr &replacement, EqualityComparator cmp = {}) {
  replace(expr.expression(), target, replacement, cmp);

  // We have to check whether the external indices have been modified by the
  // replacement and if they did, adapt the indices in the result
  IndexGroups<> externals = get_unique_indices(expr.expression());

  if (!std::ranges::equal(externals.bra, expr.bra()) ||
      !std::ranges::equal(externals.ket, expr.ket()) ||
      !std::ranges::equal(externals.aux, expr.aux())) {
    // Externals have changed -> update result
    // TODO: Is retaining result symmetry a reasonable thing to do? Generally
    // speaking, replacements could also change the result symmetry so in
    // principle we'd need a way to deduce result symmetry.
    expr =
        ResultExpr(bra(std::move(externals.bra)), ket(std::move(externals.ket)),
                   aux(std::move(externals.aux)), expr.symmetry(),
                   expr.braket_symmetry(), expr.column_symmetry(),
                   expr.has_label() ? std::optional<std::wstring>(expr.label())
                                    : std::nullopt,
                   std::move(expr.expression()));
  }

  return expr;
}

}  // namespace sequant

#endif
