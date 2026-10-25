#ifndef SEQUANT_CORE_IO_FORMAT_HPP
#define SEQUANT_CORE_IO_FORMAT_HPP

#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <algorithm>
#include <concepts>
#include <format>
#include <memory>
#include <ostream>
#include <ranges>
#include <string>
#include <string_view>
#include <type_traits>

namespace sequant::io::detail {

// Deliberately not constexpr: when compile-time format-string checking meets
// a bad specifier, the compiler's diagnostic names this function.
inline void invalid_expression_format_specifier() {}

template <typename T>
  requires(std::same_as<T, Expr> || std::same_as<T, ResultExpr>)
std::string expression_string(const T* expr, bool serialize) {
  return expr ? (serialize ? serialization::to_string(*expr)
                           : latex::to_string(*expr))
              : std::string{};
}

template <bool SerializeDefault>
struct expression_formatter {
  constexpr auto parse(std::format_parse_context& ctx) {
    const auto begin = ctx.begin();
    auto end = begin;
    while (end != ctx.end() && *end != '}') ++end;
    const auto matches = [&](std::string_view name) {
      return std::equal(begin, end, name.begin(), name.end());
    };
    if (begin == end) {
      serialize_ = SerializeDefault;
    } else if (matches("l") || matches("latex")) {
      serialize_ = false;
    } else if (matches("s") || matches("serialize")) {
      serialize_ = true;
    } else {
      if (std::is_constant_evaluated()) invalid_expression_format_specifier();
      throw sequant::Exception("Invalid expression format specifier");
    }
    return end;
  }

  template <typename T, typename FormatContext>
  auto format(const T& expr, FormatContext& ctx) const {
    const auto text = [&] {
      if constexpr (std::same_as<T, sequant::ExprPtr>) {
        return expression_string(expr.get(), serialize_);
      } else if constexpr (std::same_as<T, sequant::ResultExpr>) {
        return expression_string(std::addressof(expr), serialize_);
      } else {
        return expression_string<sequant::Expr>(std::addressof(expr),
                                                serialize_);
      }
    }();
    return std::copy(text.begin(), text.end(), ctx.out());
  }

 private:
  bool serialize_ = SerializeDefault;
};

}  // namespace sequant::io::detail

namespace sequant {

/// Inserts the expression's serialized representation into a stream.
inline std::ostream& operator<<(std::ostream& stream, const Expr& expr) {
  return stream << io::detail::expression_string(&expr, true);
}

/// Inserts serialization, or an empty string for an empty expression pointer.
inline std::ostream& operator<<(std::ostream& stream, const ExprPtr& expr) {
  return stream << io::detail::expression_string(expr.get(), true);
}

/// Inserts the result expression's serialized representation.
inline std::ostream& operator<<(std::ostream& stream, const ResultExpr& expr) {
  return stream << io::detail::expression_string(&expr, true);
}

}  // namespace sequant

template <typename T>
  requires(std::derived_from<T, sequant::Expr> ||
           std::same_as<T, sequant::ExprPtr> ||
           std::same_as<T, sequant::ResultExpr>)
struct std::formatter<T, char>
    : sequant::io::detail::expression_formatter<true> {};

#if defined(__cpp_lib_format_ranges)
// Expressions iterate over their subexpressions; without this the standard
// range formatter and the formatter above are ambiguous. The library's own
// partial specialization is `template <ranges::input_range R> requires
// same_as<R, remove_cvref_t<R>>`; this one restates those constraints so that
// it subsumes them and is the more specialized of the two.
template <std::ranges::input_range T>
  requires std::same_as<T, std::remove_cvref_t<T>> &&
               std::derived_from<T, sequant::Expr>
inline constexpr std::range_format std::format_kind<T> =
    std::range_format::disabled;
#endif

#endif  // SEQUANT_CORE_IO_FORMAT_HPP
