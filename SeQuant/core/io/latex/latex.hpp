//
// Created by Eduard Valeyev on 3/30/18.
//

#ifndef SEQUANT_CORE_IO_LATEX_LATEX_HPP
#define SEQUANT_CORE_IO_LATEX_LATEX_HPP

#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <cmath>
#include <complex>
#include <limits>
#include <string>
#include <string_view>
#include <type_traits>

namespace sequant::io::latex {

template <typename T>
concept has_to_latex_member = requires(const T& t) { t.to_latex(); };

template <typename T>
concept pointer_can_call_to_latex = requires(const T& t) { t->to_latex(); };

template <typename T>
  requires(has_to_latex_member<T>)
std::wstring to_string(T&& t) {
  return t.to_latex();
}

template <typename T>
  requires(!has_to_latex_member<T> && pointer_can_call_to_latex<T>)
std::wstring to_string(T&& t) {
  return t->to_latex();
}

template <typename T>
  requires(!has_to_latex_member<T> && !pointer_can_call_to_latex<T> &&
           requires(const T& t) { t._to_latex(); })
std::wstring to_string(const T& t) {
  return t._to_latex();
}

template <typename T>
std::enable_if_t<std::is_arithmetic_v<std::decay_t<T>> &&
                     !std::is_floating_point_v<std::decay_t<T>>,
                 std::wstring>
to_string(T&& t) {
  std::wstring result = L"{";
  using ::sequant::to_wstring;
  result += to_wstring(t) + L"}";
  return result;
}

template <typename T>
std::enable_if_t<std::is_arithmetic_v<std::decay_t<T>> &&
                     std::is_floating_point_v<std::decay_t<T>> &&
                     !std::is_same_v<std::decay_t<T>, rational>,
                 std::wstring>
to_string(T&& t) {
  using Real = std::decay_t<T>;
  static const auto eps_sqrt = std::sqrt(std::numeric_limits<Real>::epsilon());

  std::wstring result = L"{";
  using ::sequant::to_wstring;

  const long round_t = std::lround(t);
  if (std::abs(round_t - t) < eps_sqrt)  // exact integer
    result += to_wstring(round_t) + L"}";
  else {
    const auto inv_t = Real(1) / t;
    const long round_inv_t = std::lround(inv_t);
    if (std::abs(round_inv_t - inv_t) < eps_sqrt) {  // exact inverse of an
                                                     // integer
      long denom = round_inv_t;
      using namespace std::literals;
      result += (std::signbit(t) ? L"-"s : L""s) + L"\\frac{1}{"s +
                std::to_wstring(std::abs(denom)) + L"}}"s;
    } else
      result += to_wstring(t) + L"}";
  }
  return result;
}

template <typename T>
std::wstring to_string(const std::complex<T>& t) {
  std::wstring result = L"{";
  result += to_string(t.real());
  if (t.imag() > 0) {
    result += L" + i " + to_string(t.imag());
  } else if (t.imag() < 0)
    result += L" - i " + to_string(-t.imag());
  result += L"}";
  return result;
}

std::wstring to_string(const rational& num);

namespace detail {

template <typename Char, typename Traits, typename Alloc>
std::basic_string<Char, Traits, Alloc> greek_characters_to_string_impl(
    std::basic_string_view<Char, Traits> str);

}  // namespace detail

/// replaces certain greek characters in a string with their (math-mode) LaTeX
/// equivalents
/// @tparam Char character type
/// @tparam Traits character traits type
/// @param str input string
/// @return string with greek characters contained within Unicode ranges
/// `[0x3B1,0x3C9]`
///         and `[0391,03A9]`replaced with their LaTeX equivalents
/// @warning if @p str contains non-ASCII characters `Char` must be `wchar_t`
/// @throw Exception if @p Char is narrow and @p str contains
/// non-ASCII characters
template <typename Char, typename Traits>
std::basic_string<Char, Traits> greek_characters_to_string(
    const std::basic_string_view<Char, Traits>& str) {
  return detail::greek_characters_to_string_impl<Char, Traits,
                                                 std::allocator<Char>>(str);
}

template <typename Char, typename Traits, typename Alloc>
std::basic_string<Char, Traits, Alloc> greek_characters_to_string(
    const std::basic_string<Char, Traits, Alloc>& str) {
  return detail::greek_characters_to_string_impl<Char, Traits, Alloc>(str);
}

namespace detail {

template <typename Char, typename Traits, typename Alloc>
std::basic_string<Char, Traits, Alloc> diactrics_to_string_impl(
    std::basic_string_view<Char, Traits> str);

}  // namespace detail

/// replaces certain diactric marks with their (math-mode) LaTeX equivalents
/// @tparam Char character type
/// @tparam Traits character traits type
/// @param str input string
/// @return with some diactrics replaced with their LaTeX equivalents
/// @warning if @p str contains non-ASCII characters `Char` must be `wchar_t`
/// @throw Exception if @p Char is narrow and @p str contains
/// non-ASCII characters
/// @note only some combined Unicode characters are currently supported
template <typename Char, typename Traits>
std::basic_string<Char, Traits> diactrics_to_string(
    const std::basic_string_view<Char, Traits>& str) {
  return detail::diactrics_to_string_impl<Char, Traits, std::allocator<Char>>(
      str);
}

template <typename Char, typename Traits, typename Alloc>
std::basic_string<Char, Traits, Alloc> diactrics_to_string(
    const std::basic_string<Char, Traits, Alloc>& str) {
  return detail::diactrics_to_string_impl<Char, Traits, Alloc>(str);
}

/// replaces certain Unicode characters with their (math-mode) LaTeX equivalents
/// @tparam Char character type
/// @tparam Traits character traits type
/// @param str input string
/// @return @p str with some Unicode characters replaced by their LaTeX
/// equivalents
/// @warning if @p str contains non-ASCII characters `Char` must be `wchar_t`
/// @throw Exception if @p Char is narrow and @p str contains
/// non-ASCII characters
/// @note only some combined Unicode characters are currently supported
/// @internal useful resources
/// - https://milde.users.sourceforge.net/LUCR/Math/unimathsymbols.pdf
/// - https://www.unicode.org/charts/PDF/U1D400.pdf
template <typename Char, typename Traits>
std::basic_string<Char, Traits> utf_to_string(
    const std::basic_string_view<Char, Traits>& str) {
  // replace diacritics first since it relies on wide character structure
  auto tmp = diactrics_to_string(str);
  return greek_characters_to_string(tmp);
}

template <typename Char, typename Traits, typename Alloc>
std::basic_string<Char, Traits, Alloc> utf_to_string(
    const std::basic_string<Char, Traits, Alloc>& str) {
  auto tmp = diactrics_to_string(str);
  return greek_characters_to_string(tmp);
}

}  // namespace sequant::io::latex

#endif  // SEQUANT_CORE_IO_LATEX_LATEX_HPP
