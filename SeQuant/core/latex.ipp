//
// Created by Eduard Valeyev on 3/30/18.
//

#ifndef SEQUANT_CORE_LATEX_IPP
#define SEQUANT_CORE_LATEX_IPP

#include <SeQuant/core/latex.hpp>
#include <SeQuant/core/wstring.hpp>
#include <map>
#include <optional>
#include <type_traits>
#include <vector>

namespace sequant {
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
    auto append = [&result,&str,&it,&begin](const auto& s) {
      if (result.empty()) result = str.substr(0, it - begin);
      result += s;
    };
    auto is_ascii = [](Char c) {
      return static_cast<unsigned int>(c) <= 0x7F;
    };

    const Char ch = *it;
    if (sizeof(Char) == 1 && !is_ascii(ch))
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
      append(lc_str);
    } else if (ch_is_uc) {
      const auto ch_addr = static_cast<long>(ch) - static_cast<long>(0x391);
      assert(ch_addr >= 0 && ch_addr < uc.size());
      const auto& uc_str = uc[static_cast<std::size_t>(ch_addr)];
      assert(uc_str.size() > 0);
      append(uc_str);
    } else {
      if (!result.empty()) result.push_back(ch);
    }
  }

  if (!result.empty())
    return result;
  else
    return decltype(result)(str);
}

template <typename Char, typename Traits, typename Alloc>
std::basic_string<Char, Traits, Alloc> diactrics_to_latex_impl(
    std::basic_string_view<Char, Traits> str) {
  using str_t = std::basic_string<Char, Traits, Alloc>;

  str_t result;
  const auto begin = cbegin(str);
  const auto end = cend(str);
  /// TODO iterate over characters (graphemes)
  std::optional<Char> next_ch;
  for (auto it = begin; it != end; ++it) {
    auto append = [&result, &str, &it, &begin](const auto& s) {
      if (result.empty()) result = str.substr(0, it - begin);
      result += s;
    };
    auto is_ascii = [](Char c) { return static_cast<unsigned int>(c) <= 0x7F; };

    const Char ch = *it;
    if (it + 1 != end) next_ch = *(it + 1);
    if (sizeof(Char) == 1 &&
        ((it == begin && !is_ascii(ch)) || (next_ch && !is_ascii(*next_ch)))) {
      throw std::invalid_argument(
          "diactrics_to_latex<Char,...>(str): currently only supports "
          "non-ASCII characters in str if Char is a wide character (wchar_t, "
          "char16_t, or char32_t)");
    }
    // Combining diacritics: https://www.ncbi.nlm.nih.gov/staff/beck/charents/accents.html
    if (next_ch) {
      // tilde
      if (*next_ch == static_cast<Char>(0x303)) {
        append(SQ_STRLIT(Char, "\\tilde{"));
      }
      // acute
      else if (*next_ch == static_cast<Char>(0x301)) {
        append(SQ_STRLIT(Char, "\\acute{"));
      }
      // grave
      else if (*next_ch == static_cast<Char>(0x300)) {
        append(SQ_STRLIT(Char, "\\grave{"));
      }
      // caron
      else if (*next_ch == static_cast<Char>(0x30C)) {
        append(SQ_STRLIT(Char, "\\check{"));
      }
      append(ch);
      append(SQ_STRLIT(Char, "}"));
      it += 1;
      continue;
    }

    {    // check for combined characters
      {  // tilde
         // lower-case characters with tilde
        const std::map<str_t, str_t> lc = {
            {SQ_STRLIT(Char, "ã"), SQ_STRLIT(Char, "\\tilde{a}")},
            {SQ_STRLIT(Char, "ẽ"), SQ_STRLIT(Char, "\\tilde{e}")},
            {SQ_STRLIT(Char, "ñ"), SQ_STRLIT(Char, "\\tilde{n}")},
            {SQ_STRLIT(Char, "õ"), SQ_STRLIT(Char, "\\tilde{o}")},
            {SQ_STRLIT(Char, "ũ"), SQ_STRLIT(Char, "\\tilde{u}")},
            {SQ_STRLIT(Char, "ṽ"), SQ_STRLIT(Char, "\\tilde{v}")}};
        auto lc_it = lc.find(str_t{ch});
        if (lc_it != lc.end()) {
          append(lc_it->second);
        } else {
          // upper-case characters with tilde
          const std::map<str_t, str_t> uc = {
              {SQ_STRLIT(Char, "Ã"), SQ_STRLIT(Char, "\\tilde{A}")},
              {SQ_STRLIT(Char, "Ẽ"), SQ_STRLIT(Char, "\\tilde{E}")},
              {SQ_STRLIT(Char, "Ñ"), SQ_STRLIT(Char, "\\tilde{N}")},
              {SQ_STRLIT(Char, "Õ"), SQ_STRLIT(Char, "\\tilde{O}")},
              {SQ_STRLIT(Char, "Ũ"), SQ_STRLIT(Char, "\\tilde{U}")},
              {SQ_STRLIT(Char, "Ṽ"), SQ_STRLIT(Char, "\\tilde{V}")}};
          auto uc_it = uc.find(str_t{ch});
          if (uc_it != uc.end()) {
            append(uc_it->second);
          }
        }
      }  // tilde
    }
  }

  if (!result.empty())
    return result;
  else
    return decltype(result)(str);
}

}  // namespace detail

}  // namespace sequant

#endif  // SEQUANT_CORE_LATEX_IPP
