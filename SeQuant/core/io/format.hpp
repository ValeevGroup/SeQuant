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
#include <string>
#include <string_view>

namespace sequant::io::detail {

template <typename Char, typename T>
  requires(std::same_as<T, Expr> || std::same_as<T, ResultExpr>)
std::basic_string<Char> expression_string(const T* expr, bool serialize) {
  auto text = expr ? (serialize ? serialization::to_string(*expr)
                                : latex::to_string(*expr))
                   : L"";
  if constexpr (std::same_as<Char, char>) {
    return toUtf8(text);
  } else {
    return text;
  }
}

template <typename T, typename Char, typename Traits>
std::basic_ostream<Char, Traits>& insert_expression(
    std::basic_ostream<Char, Traits>& stream, const T* expr, bool serialize) {
  const auto text = expression_string<Char>(expr, serialize);
  return stream << std::basic_string_view<Char, Traits>(text.data(),
                                                        text.size());
}

template <typename Char, bool SerializeDefault>
struct expression_formatter {
  constexpr auto parse(std::basic_format_parse_context<Char>& ctx) {
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
      throw sequant::Exception("Invalid expression format specifier");
    }
    return end;
  }

  template <typename T, typename FormatContext>
  auto format(const T& expr, FormatContext& ctx) const {
    const auto text = [&] {
      if constexpr (std::same_as<T, sequant::ExprPtr>) {
        return expression_string<Char>(expr.get(), serialize_);
      } else if constexpr (std::same_as<T, sequant::ResultExpr>) {
        return expression_string<Char>(std::addressof(expr), serialize_);
      } else {
        return expression_string<Char, sequant::Expr>(std::addressof(expr),
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

/// Inserts the expression's serialized representation into a narrow or wide
/// stream.
template <typename Char, typename Traits>
  requires(std::same_as<Char, char> || std::same_as<Char, wchar_t>)
std::basic_ostream<Char, Traits>& operator<<(
    std::basic_ostream<Char, Traits>& stream, const Expr& expr) {
  return io::detail::insert_expression(stream, &expr, true);
}

/// Inserts serialization, or an empty string for an empty expression pointer.
template <typename Char, typename Traits>
  requires(std::same_as<Char, char> || std::same_as<Char, wchar_t>)
std::basic_ostream<Char, Traits>& operator<<(
    std::basic_ostream<Char, Traits>& stream, const ExprPtr& expr) {
  return io::detail::insert_expression(stream, expr.get(), true);
}

/// Inserts the result expression's serialized representation.
template <typename Char, typename Traits>
  requires(std::same_as<Char, char> || std::same_as<Char, wchar_t>)
std::basic_ostream<Char, Traits>& operator<<(
    std::basic_ostream<Char, Traits>& stream, const ResultExpr& expr) {
  return io::detail::insert_expression(stream, &expr, true);
}

}  // namespace sequant

template <typename T, typename Char>
  requires((std::derived_from<T, sequant::Expr> ||
            std::same_as<T, sequant::ExprPtr> ||
            std::same_as<T, sequant::ResultExpr>) &&
           (std::same_as<Char, char> || std::same_as<Char, wchar_t>))
struct std::formatter<T, Char>
    : sequant::io::detail::expression_formatter<Char, true> {};

#if defined(__cpp_lib_format_ranges)
// Expressions iterate over their subexpressions; without this the standard
// range formatter and the formatter above are ambiguous.
template <typename T>
  requires std::derived_from<T, sequant::Expr>
inline constexpr std::range_format std::format_kind<T> =
    std::range_format::disabled;
#endif

#endif  // SEQUANT_CORE_IO_FORMAT_HPP
