#ifndef SEQUANT_EXPRESSIONS_LEXICOGRAPHIC_ORDER_HPP
#define SEQUANT_EXPRESSIONS_LEXICOGRAPHIC_ORDER_HPP

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/expressions/sum.hpp>

#include <algorithm>
#include <compare>
#include <concepts>
#include <utility>

namespace sequant {

/// Strict weak order of expressions based on Expr::lexicographic_compare,
/// i.e. on their contents rather than their hash values. Hence, unlike the
/// (faster) order defined by Expr::operator< , it does not change when the
/// hashing algorithm changes.
/// @note the comparison options (in particular the cardinal tensor labels)
/// are captured on construction
class LexicographicExprCompare {
 public:
  /// uses LexicographicCompareOptions::from_default_context()
  LexicographicExprCompare()
      : opts_(LexicographicCompareOptions::from_default_context()) {}

  explicit LexicographicExprCompare(LexicographicCompareOptions opts)
      : opts_(std::move(opts)) {}

  /// @return the order of @p a relative to @p b
  std::weak_ordering compare(const Expr &a, const Expr &b) const {
    return a.lexicographic_compare(b, opts_);
  }

  /// @return true if @p a precedes @p b
  bool operator()(const Expr &a, const Expr &b) const {
    return compare(a, b) < 0;
  }

  /// @return true if @p a precedes @p b
  bool operator()(const ExprPtr &a, const ExprPtr &b) const {
    return (*this)(*a, *b);
  }

  const LexicographicCompareOptions &options() const { return opts_; }

 private:
  LexicographicCompareOptions opts_;
};

/// strict weak order of ExprPtr objects
template <typename Compare>
concept expr_ptr_compare =
    std::predicate<const Compare &, const ExprPtr &, const ExprPtr &>;

/// Sorts @p expr recursively and in place according to @p cmp : the summands
/// of every Sum are sorted, as are the factors of every Product insofar as
/// they commute (i.e. the order of non-commuting factors is preserved). The
/// order of equivalent summands/factors is preserved.
/// @param expr the expression to sort
/// @param cmp the strict weak order of expressions to sort by; the default
/// orders expressions lexicographically (see LexicographicExprCompare), which,
/// unlike canonicalization, yields an order that does not depend on the hashing
/// algorithm
/// @return @p expr
template <expr_ptr_compare Compare = LexicographicExprCompare>
Expr &sort_expr(Expr &expr, const Compare &cmp = {}) {
  if (expr.is_atom()) return expr;

  // sort bottom-up so that the order of subexpressions is final before they
  // are compared to each other
  for (auto &subexpr : expr) sort_expr(*subexpr, cmp);

  if (expr.is<Sum>()) {
    std::stable_sort(expr.begin(), expr.end(), cmp);
  } else if (expr.is<Product>()) {
    // must do bubble sort to avoid swapping factors across a noncommuting
    // factor (cf. Product::canonicalize_impl)
    bubble_sort(expr.begin(), expr.end(),
                [&cmp](const ExprPtr &first, const ExprPtr &second) {
                  return first->commutes_with(*second) && cmp(first, second);
                });
  }

  return expr;
}

/// Sorts the expression of @p expr , see sort_expr(Expr&,const Compare&)
/// @return @p expr
template <expr_ptr_compare Compare = LexicographicExprCompare>
ResultExpr &sort_expr(ResultExpr &expr, const Compare &cmp = {}) {
  sort_expr(*expr.expression(), cmp);
  return expr;
}

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_LEXICOGRAPHIC_ORDER_HPP
