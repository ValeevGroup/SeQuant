//
// Created by Eduard Valeyev on 3/30/18.
//

#ifndef SEQUANT_CORE_LATEX_HPP
#define SEQUANT_CORE_LATEX_HPP

#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/wstring.hpp>
#include <type_traits>

namespace sequant {

template <typename T>
std::enable_if_t<meta::has_memfn_to_latex_v<std::decay_t<T>>, std::wstring>
to_latex(T&& t) {
  return t.to_latex();
}

template <typename T>
std::enable_if_t<std::is_arithmetic_v<std::decay_t<T>> &&
                     !std::is_floating_point_v<std::decay_t<T>>,
                 std::wstring>
to_latex(T&& t) {
  std::wstring result = L"{";
  using ::sequant::to_wstring;
  result += to_wstring(t) + L"}";
  return result;
}

template <typename T>
std::enable_if_t<std::is_arithmetic_v<std::decay_t<T>> &&
                     std::is_floating_point_v<std::decay_t<T>>,
                 std::wstring>
to_latex(T&& t) {
  using Real = std::decay_t<T>;
  static const auto eps_sqrt = std::sqrt(std::numeric_limits<Real>::epsilon());

  std::wstring result = L"{";
  using ::sequant::to_wstring;

  const long round_t = std::lround(t);
  if (std::abs(round_t - t) < eps_sqrt)  // exact integer
    result += to_wstring(round_t) + L"}";
  else {  // TODO detect rationals
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
std::wstring to_latex(const std::complex<T>& t) {
  std::wstring result = L"{";
  result += to_latex(t.real());
  if (t.imag() > 0) {
    result += L" + i " + to_latex(t.imag());
  } else if (t.imag() < 0)
    result += L" - i " + to_latex(-t.imag());
  result += L"}";
  return result;
}

namespace detail {

template <typename Char, typename Traits, typename Alloc>
std::basic_string<Char, Traits, Alloc> greek_characters_to_latex_impl(
    std::basic_string_view<Char, Traits> str);

}  // namespace detail

// chang-format off
/// replaces certain greek characters in a string with their latex equivalents
/// @tparam Char character type
/// @tparam Traits character traits type
/// @param str input string
/// @return string with greek characters contained within Unicode ranges
/// `[0x3B1,0x3C9]`
///         and `[0391,03A9]`replaced with their LaTeX equivalents
/// @warning if @p str contains non-ASCII characters `Char` must be `wchar_t`
/// @throw std::invalid_argument if @p Char is narrow and @p str contains
/// non-ASCII characters
// chang-format on
template <typename Char, typename Traits>
std::basic_string<Char, Traits> greek_characters_to_latex(
    const std::basic_string_view<Char, Traits>& str) {
  return detail::greek_characters_to_latex_impl<Char, Traits,
                                                std::allocator<Char>>(str);
}

template <typename Char, typename Traits, typename Alloc>
std::basic_string<Char, Traits, Alloc> greek_characters_to_latex(
    const std::basic_string<Char, Traits, Alloc>& str) {
  return detail::greek_characters_to_latex_impl<Char, Traits, Alloc>(str);
}

namespace detail {

template <typename Char, typename Traits, typename Alloc>
std::basic_string<Char, Traits, Alloc> diactrics_to_latex_impl(
    std::basic_string_view<Char, Traits> str);

}  // namespace detail

// chang-format off
/// replaces certain diactric marks with their latex equivalents
/// @tparam Char character type
/// @tparam Traits character traits type
/// @param str input string
/// @return with some diactrics replaced with their LaTeX equivalents
/// @warning if @p str contains non-ASCII characters `Char` must be `wchar_t`
/// @throw std::invalid_argument if @p Char is narrow and @p str contains
/// non-ASCII characters
/// @note only tilde is currently supported
// chang-format on
template <typename Char, typename Traits>
std::basic_string<Char, Traits> diactrics_to_latex(
    const std::basic_string_view<Char, Traits>& str) {
  return detail::diactrics_to_latex_impl<Char, Traits, std::allocator<Char>>(
      str);
}

template <typename Char, typename Traits, typename Alloc>
std::basic_string<Char, Traits, Alloc> diactrics_to_latex(
    const std::basic_string<Char, Traits, Alloc>& str) {
  return detail::diactrics_to_latex_impl<Char, Traits, Alloc>(str);
}

}  // namespace sequant

#endif  // SEQUANT_CORE_LATEX_HPP
