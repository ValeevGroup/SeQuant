#ifndef SEQUANT_TESTS_CATCH2_SEQUANT_H
#define SEQUANT_TESTS_CATCH2_SEQUANT_H

#include <catch2/catch_tostring.hpp>
#include <catch2/matchers/catch_matchers_templated.hpp>

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/latex.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/wstring.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

#include <cassert>
#include <string>
#include <type_traits>

namespace Catch {

// Make sure Catch uses proper string representation for SeQuant types

template <>
struct StringMaker<sequant::Expr> {
  static std::string convert(const sequant::Expr &expr) {
    try {
      return sequant::to_string(sequant::deparse(expr, true));
    } catch (const std::exception &) {
      // deparse doesn't support all kinds of expressions -> fall back to LaTeX
      // representation
      return sequant::to_string(sequant::to_latex(expr));
    }
  }
};
template <>
struct StringMaker<sequant::ExprPtr> {
  static std::string convert(const sequant::ExprPtr &expr) {
    return StringMaker<sequant::Expr>::convert(*expr);
  }
};

template <>
struct StringMaker<sequant::Index> {
  static std::string convert(const sequant::Index &idx) {
    return sequant::to_string(idx.full_label());
  }
};

}  // namespace Catch

/// Matches that the tested expression is equivalent to the given one. Two
/// expressions are considered equivalent if they both have the same canonical
/// form (i.e. they are the same expression after canonicalization)
class EquivalentToMatcher : public Catch::Matchers::MatcherGenericBase {
 public:
  /// Constructs the matcher with the expected expression. The constructor
  /// accepts either an actual expression object (as Expr & or ExprPtr) or
  /// a (w)string-like object which will then be parsed to yield the actual
  /// expression object.
  template <typename T>
  EquivalentToMatcher(T &&expression) {
    if constexpr (std::is_convertible_v<T, std::string>) {
      m_expr = sequant::parse_expr(
          sequant::to_wstring(std::string(std::forward<T>(expression))),
          sequant::Symmetry::nonsymm);
    } else if constexpr (std::is_convertible_v<T, std::wstring>) {
      m_expr = sequant::parse_expr(std::wstring(std::forward<T>(expression)),
                                   sequant::Symmetry::nonsymm);
    } else if constexpr (std::is_convertible_v<T, sequant::Expr>) {
      // Clone in order to not have to worry about later modification
      m_expr = expression.clone();
    } else {
      static_assert(std::is_convertible_v<T, sequant::ExprPtr>,
                    "Invalid type for expression");

      // Clone in order to not have to worry about later modification
      m_expr = expression->clone();
    }

    assert(m_expr);

    // Bring expression into canonical form
    m_expr->canonicalize();
  }

  bool match(const sequant::Expr &expr) const {
    // Never modify the expression that we are trying to check in order to avoid
    // side-effects
    sequant::ExprPtr clone = expr.clone();

    clone->canonicalize();

    return *clone == *m_expr;
  }

  std::string describe() const override {
    return "Equivalent to: " + Catch::Detail::stringify(m_expr);
  }

 private:
  sequant::ExprPtr m_expr;
};

template <typename T>
EquivalentToMatcher EquivalentTo(T &&expression) {
  return EquivalentToMatcher(std::forward<T>(expression));
}

#endif
