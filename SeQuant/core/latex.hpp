//
// Created by Eduard Valeyev on 3/30/18.
//

#ifndef SEQUANT_LATEX_HPP
#define SEQUANT_LATEX_HPP

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
    std::basic_string_view<Char, Traits> str) {
  // lower-case greek characters in the order of their appearance in Unicode
  // chart https://www.unicode.org/charts/PDF/U0370.pdf
  const std::vector<std::basic_string<Char, Traits, Alloc>> lc = {
      SQ_STRLIT(Char, "\\alpha"),   SQ_STRLIT(Char, "\\beta"),
      SQ_STRLIT(Char, "\\gamma"),   SQ_STRLIT(Char, "\\delta"),
      SQ_STRLIT(Char, "\\epsilon"), SQ_STRLIT(Char, "\\zeta"),
      SQ_STRLIT(Char, "\\eta"),     SQ_STRLIT(Char, "\\theta"),
      SQ_STRLIT(Char, "\\iota"),    SQ_STRLIT(Char, "\\kappa"),
      SQ_STRLIT(Char, "\\lambda"),  SQ_STRLIT(Char, "\\mu"),
      SQ_STRLIT(Char, "\\nu"),      SQ_STRLIT(Char, "\\xi"),
      SQ_STRLIT(Char, "o"),         SQ_STRLIT(Char, "\\pi"),
      SQ_STRLIT(Char, "\\rho"),     SQ_STRLIT(Char, "\\varsigma"),
      SQ_STRLIT(Char, "\\sigma"),   SQ_STRLIT(Char, "\\tau"),
      SQ_STRLIT(Char, "\\upsilon"), SQ_STRLIT(Char, "\\phi"),
      SQ_STRLIT(Char, "\\chi"),     SQ_STRLIT(Char, "\\psi"),
      SQ_STRLIT(Char, "\\omega")};
  const std::vector<std::basic_string<Char, Traits, Alloc>> uc = {
      SQ_STRLIT(Char, "A"),         SQ_STRLIT(Char, "B"),
      SQ_STRLIT(Char, "\\Gamma"),   SQ_STRLIT(Char, "\\Delta"),
      SQ_STRLIT(Char, "E"),         SQ_STRLIT(Char, "Z"),
      SQ_STRLIT(Char, "H"),         SQ_STRLIT(Char, "\\Theta"),
      SQ_STRLIT(Char, "I"),         SQ_STRLIT(Char, "K"),
      SQ_STRLIT(Char, "\\Lambda"),  SQ_STRLIT(Char, "M"),
      SQ_STRLIT(Char, "N"),         SQ_STRLIT(Char, "\\Xi"),
      SQ_STRLIT(Char, "O"),         SQ_STRLIT(Char, "\\Pi"),
      SQ_STRLIT(Char, "P"),         SQ_STRLIT(Char, ""),
      SQ_STRLIT(Char, "\\Sigma"),   SQ_STRLIT(Char, "T"),
      SQ_STRLIT(Char, "\\Upsilon"), SQ_STRLIT(Char, "\\Phi"),
      SQ_STRLIT(Char, "X"),         SQ_STRLIT(Char, "\\Psi"),
      SQ_STRLIT(Char, "\\Omega")};
  auto is_lc = [](Char ch) {
    return ch >= static_cast<Char>(0x3B1)     // alpha
           && ch <= static_cast<Char>(0x3C9)  // omega
        ;
  };
  auto is_uc = [](Char ch) {
    return ch >= static_cast<Char>(0x391)     // Alpha
           && ch <= static_cast<Char>(0x3A9)  // Omega
           && ch != static_cast<Char>(0x3A2)  // gap in the table
        ;
  };

  std::basic_string<Char, Traits, Alloc> result;
  const auto begin = cbegin(str);
  const auto end = cend(str);
  /// TODO iterate over characters (graphemes)
  for (auto it = begin; it != end; ++it) {
    const Char ch = *it;
    if (sizeof(Char) == 1 && static_cast<unsigned int>(ch) > 0x7F)
      throw std::invalid_argument(
          "greek_characters_to_latex<Char,...>(str): currently only supports "
          "non-ASCII characters in str if Char is a wide character (wchar_t, "
          "char16_t, or char32_t)");
    const bool ch_is_lc = is_lc(ch);
    const bool ch_is_uc = is_uc(ch);
    if (ch_is_lc) {
      const auto ch_addr = static_cast<long>(ch) - static_cast<long>(0x3B1);
      assert(ch_addr >= 0 && ch_addr < lc.size());
      const auto& lc_str = lc[static_cast<std::size_t>(ch_addr)];
      assert(lc_str.size() > 0);
      if (result.empty()) result = str.substr(0, it - begin);
      result += lc_str;
    } else if (ch_is_uc) {
      const auto ch_addr = static_cast<long>(ch) - static_cast<long>(0x391);
      assert(ch_addr >= 0 && ch_addr < uc.size());
      const auto& uc_str = uc[static_cast<std::size_t>(ch_addr)];
      assert(uc_str.size() > 0);
      if (result.empty()) result = str.substr(0, it - begin);
      result += uc_str;
    } else {
      if (!result.empty()) result.push_back(ch);
    }
  }

  if (!result.empty())
    return result;
  else
    return decltype(result)(str);
}
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

}  // namespace sequant

#endif  // SEQUANT_LATEX_HPP
