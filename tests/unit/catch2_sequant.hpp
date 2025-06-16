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
#include <SeQuant/core/result_expr.hpp>
#include <SeQuant/core/wstring.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

#include <range/v3/algorithm.hpp>

#include <algorithm>
#include <cassert>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <variant>

namespace Catch {

// Make sure Catch uses proper string representation for SeQuant types

template <>
struct StringMaker<sequant::Expr> {
  static std::string convert(const sequant::Expr &expr,
                             bool include_canonical = true) {
    std::string str;
    try {
      str = sequant::to_string(sequant::deparse(expr, true));
    } catch (const std::exception &) {
      // deparse doesn't support all kinds of expressions -> fall back to LaTeX
      // representation
      str = sequant::to_string(sequant::to_latex(expr));
    }

    if (include_canonical) {
      sequant::ExprPtr clone = expr.clone();
      canonicalize(clone);
      simplify(clone);
      std::string canon_str =
          StringMaker<sequant::Expr>::convert(*clone, false);

      if (canon_str != str) {
        str += " (canonicalized: " + canon_str + ")";
      }
    }

    return str;
  }
};
template <>
struct StringMaker<sequant::ExprPtr> {
  static std::string convert(const sequant::ExprPtr &expr) {
    return StringMaker<sequant::Expr>::convert(*expr);
  }
};
template <>
struct StringMaker<sequant::Tensor> {
  static std::string convert(const sequant::Tensor &tensor) {
    return StringMaker<sequant::Expr>::convert(tensor);
  }
};

template <>
struct StringMaker<sequant::Index> {
  static std::string convert(const sequant::Index &idx) {
    return sequant::to_string(idx.full_label());
  }
};

template <>
struct StringMaker<sequant::ResultExpr> {
  static std::string convert(const sequant::ResultExpr &res) {
    std::stringstream sstream;
    if (res.has_label()) {
      sstream << sequant::to_string(res.label());
    } else {
      sstream << "?";
    }

    using IndexString = StringMaker<sequant::Index>;
    using namespace std::literals;
    sstream << "{";
    sstream << (res.bra() | ranges::views::transform(IndexString::convert) |
                ranges::views::join(", "sv) | ranges::to<std::string>());
    sstream << ";";
    sstream << (res.ket() | ranges::views::transform(IndexString::convert) |
                ranges::views::join(", "sv) | ranges::to<std::string>());
    sstream << ";";
    sstream << (res.aux() | ranges::views::transform(IndexString::convert) |
                ranges::views::join(", "sv) | ranges::to<std::string>());
    sstream << "}";

    sstream << " = "
            << StringMaker<sequant::ExprPtr>::convert(res.expression());

    return sstream.str();
  }
};

}  // namespace Catch

namespace {

using ExprVar = std::variant<sequant::ExprPtr, sequant::ResultExpr>;

/// Converts the given expression-like object into an actual ExprPtr.
/// It accepts either an actual expression object (as Expr & or ExprPtr) or
/// a (w)string-like object which will then be parsed to yield the actual
/// expression object.
template <typename T>
ExprVar to_expression(T &&expression) {
  using std::begin;
  using std::end;

  if constexpr (std::is_convertible_v<T, std::string>) {
    std::wstring string = sequant::to_wstring(std::forward<T>(expression));

    if (std::find(begin(string), end(string), L'=') != end(string)) {
      return sequant::parse_result_expr(
          sequant::to_wstring(std::string(std::forward<T>(expression))),
          sequant::Symmetry::nonsymm);
    } else {
      return sequant::parse_expr(
          sequant::to_wstring(std::string(std::forward<T>(expression))),
          sequant::Symmetry::nonsymm);
    }
  } else if constexpr (std::is_convertible_v<T, std::wstring>) {
    if (std::find(begin(expression), end(expression), L'=') !=
        end(expression)) {
      return sequant::parse_result_expr(
          std::wstring(std::forward<T>(expression)),
          sequant::Symmetry::nonsymm);
    } else {
      return sequant::parse_expr(std::wstring(std::forward<T>(expression)),
                                 sequant::Symmetry::nonsymm);
    }
  } else if constexpr (std::is_convertible_v<T, sequant::ResultExpr>) {
    return expression;
  } else if constexpr (std::is_convertible_v<T, sequant::Expr>) {
    // Clone in order to not have to worry about later modification
    return expression.clone();
  } else {
    static_assert(std::is_convertible_v<T, sequant::ExprPtr>,
                  "Invalid type for expression");

    // Clone in order to not have to worry about later modification
    return expression->clone();
  }
}

template <typename Subclass>
class ExpressionMatcher : public Catch::Matchers::MatcherGenericBase {
 public:
  template <typename T>
  ExpressionMatcher(T &&expression)
      : m_expr(to_expression(std::forward<T>(expression))) {
    if (std::holds_alternative<sequant::ResultExpr>(m_expr)) {
      Subclass::pre_comparison(
          std::get<sequant::ResultExpr>(m_expr).expression());
    } else {
      Subclass::pre_comparison(std::get<sequant::ExprPtr>(m_expr));
    }
  }

  bool match(const sequant::ResultExpr &res) const {
    if (!std::holds_alternative<sequant::ResultExpr>(m_expr)) {
      return false;
    }

    const sequant::ResultExpr &self = std::get<sequant::ResultExpr>(m_expr);

    return self.has_label() == res.has_label() && self.label() == res.label() &&
           self.symmetry() == res.symmetry() &&
           self.braket_symmetry() == res.braket_symmetry() &&
           self.particle_symmetry() == res.particle_symmetry() &&
           do_match(*res.expression());
  }

  bool match(const sequant::ExprPtr &expr) const {
    return std::holds_alternative<sequant::ExprPtr>(m_expr) && do_match(*expr);
  }

  bool match(const sequant::Expr &expr) const {
    return std::holds_alternative<sequant::ExprPtr>(m_expr) && do_match(expr);
  }

  std::string describe() const override {
    return Subclass::comparison_requirement() + ": " + stringify(m_expr);
  }

 protected:
  ExprVar m_expr;

  bool do_match(const sequant::Expr &expr) const {
    // Never modify the expression that we are trying to check in order to avoid
    // side-effects
    sequant::ExprPtr clone = expr.clone();

    Subclass::pre_comparison(clone);

    const sequant::Expr &self = [&]() -> const sequant::Expr & {
      if (std::holds_alternative<sequant::ResultExpr>(m_expr)) {
        return *std::get<sequant::ResultExpr>(m_expr).expression();
      }

      return *std::get<sequant::ExprPtr>(m_expr);
    }();

    return *clone == self;
  }

  std::string stringify(const ExprVar &expr) const {
    if (std::holds_alternative<sequant::ResultExpr>(expr)) {
      return Catch::Detail::stringify(std::get<sequant::ResultExpr>(expr));
    }

    return Catch::Detail::stringify(std::get<sequant::ExprPtr>(expr));
  }
};

/// Matches that the tested expression is equivalent to the given one. Two
/// expressions are considered equivalent if they both have the same canonical
/// form (i.e. they are the same expression after canonicalization)
struct EquivalentToMatcher : ExpressionMatcher<EquivalentToMatcher> {
  using ExpressionMatcher::ExpressionMatcher;

  static void pre_comparison(sequant::ExprPtr &expr) {
    sequant::canonicalize(expr);
    sequant::simplify(expr);
  }

  static std::string comparison_requirement() { return "Equivalent to"; }
};

/// Matches that the tested expression simplifies (**without**
/// re-canonicalization!) to the same form as the given one
struct SimplifiesToMatcher : ExpressionMatcher<SimplifiesToMatcher> {
  using ExpressionMatcher::ExpressionMatcher;

  static void pre_comparison(sequant::ExprPtr &expr) {
    sequant::rapid_simplify(expr);
  }

  static std::string comparison_requirement() { return "Simplifies to"; }
};

}  // namespace

template <typename T>
EquivalentToMatcher EquivalentTo(T &&expression) {
  return EquivalentToMatcher(std::forward<T>(expression));
}

template <typename T>
SimplifiesToMatcher SimplifiesTo(T &&expression) {
  return SimplifiesToMatcher(std::forward<T>(expression));
}

#endif
