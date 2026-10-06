#ifndef SEQUANT_CORE_IO_FORMAT_HPP
#define SEQUANT_CORE_IO_FORMAT_HPP

#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
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

template <typename Char>
std::basic_string<Char> expression_string(const Expr* expr, bool serialize) {
  auto text = expr ? (serialize ? serialization::to_string(*expr)
                                : latex::to_string(*expr))
                   : L"NULL";
  if constexpr (std::same_as<Char, char>) {
    return toUtf8(text);
  } else {
    return text;
  }
}

}  // namespace sequant::io::detail

namespace sequant {

/// Inserts the expression's LaTeX representation into a narrow or wide stream.
template <typename Char, typename Traits>
  requires(std::same_as<Char, char> || std::same_as<Char, wchar_t>)
std::basic_ostream<Char, Traits>& operator<<(
    std::basic_ostream<Char, Traits>& stream, const Expr& expr) {
  const auto text = io::detail::expression_string<Char>(&expr, false);
  return stream << std::basic_string_view<Char, Traits>(text.data(),
                                                        text.size());
}

/// Inserts LaTeX, or NULL for an empty expression pointer.
template <typename Char, typename Traits>
  requires(std::same_as<Char, char> || std::same_as<Char, wchar_t>)
std::basic_ostream<Char, Traits>& operator<<(
    std::basic_ostream<Char, Traits>& stream, const ExprPtr& expr) {
  const auto text = io::detail::expression_string<Char>(expr.get(), false);
  return stream << std::basic_string_view<Char, Traits>(text.data(),
                                                        text.size());
}

}  // namespace sequant

template <typename T, typename Char>
  requires((std::derived_from<T, sequant::Expr> ||
            std::same_as<T, sequant::ExprPtr>) &&
           (std::same_as<Char, char> || std::same_as<Char, wchar_t>))
struct std::formatter<T, Char> {
  constexpr auto parse(std::basic_format_parse_context<Char>& ctx) {
    const auto begin = ctx.begin();
    auto end = begin;
    while (end != ctx.end() && *end != '}') ++end;
    const auto matches = [&](std::string_view name) {
      return std::equal(begin, end, name.begin(), name.end());
    };
    if (begin == end || matches("l") || matches("latex")) {
      serialize_ = false;
    } else if (matches("s") || matches("serialize")) {
      serialize_ = true;
    } else {
      throw sequant::Exception("Invalid expression format specifier");
    }
    return end;
  }

  template <typename FormatContext>
  auto format(const T& expr, FormatContext& ctx) const {
    const sequant::Expr* ptr;
    if constexpr (std::same_as<T, sequant::ExprPtr>) {
      ptr = expr.get();
    } else {
      ptr = std::addressof(expr);
    }
    const auto text =
        sequant::io::detail::expression_string<Char>(ptr, serialize_);
    return std::copy(text.begin(), text.end(), ctx.out());
  }

 private:
  bool serialize_ = false;
};

#endif  // SEQUANT_CORE_IO_FORMAT_HPP
