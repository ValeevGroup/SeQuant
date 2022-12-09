//
// Created by Eduard Valeyev on 4/30/20.
//

#include "../../SeQuant/core/index.hpp"

#include "../../SeQuant/core/sequant.hpp"

namespace sequant {

const std::size_t Index::min_tmp_index() {
  return get_default_context().first_dummy_index_ordinal();
}

void Index::reset_tmp_index() { tmp_index_accessor() = min_tmp_index() - 1; }

std::wstring Index::to_latex() const {
  auto protect_subscript = [](const std::wstring_view str) {
    auto subsc_pos = str.find(L'_');
    if (subsc_pos == std::wstring_view::npos)
      return std::wstring(str);
    else {
      assert(subsc_pos + 1 < str.size());
      if (subsc_pos + 2 == str.size())  // don't protect single character
        return std::wstring(str);
      std::wstring result = std::wstring(str.substr(0, subsc_pos + 1)) + L"{" +
                            std::wstring(str.substr(subsc_pos + 1)) + L"}";
      return result;
    }
  };

  // replaces greek chars
  auto greek_characters_to_latex = [](std::wstring& str) {
    auto lc = {L"alpha",   L"beta", L"gamma",    L"delta", L"epsilon",
               L"zeta",    L"eta",  L"theta",    L"iota",  L"kappa",
               L"lambda",  L"mu",   L"nu",       L"xi",    L"",
               L"pi",      L"rho",  L"varsigma", L"sigma", L"tau",
               L"upsilon", L"phi",  L"chi",      L"psi",   L"omega"};
    auto is_lc = [](wchar_t ch) {
      return ch >= static_cast<wchar_t>(0x3B1)     // alpha
             && ch <= static_cast<wchar_t>(0x3C9)  // omega
          ;
    };
    auto is_uc = [](wchar_t ch) { return false; };

    std::wstring result;
    const auto begin = cbegin(str);
    const auto end = cend(str);
    for (auto it = begin; it != end; ++it) {
      const wchar_t ch = *it;
      const bool ch_is_lc = is_lc(ch);
      const bool ch_is_uc = !ch_is_lc ? false : is_uc(ch);
      if (ch_is_lc || ch_is_uc) {
        auto* lc_c_str =
            cbegin(lc) + (static_cast<long>(ch) - static_cast<long>(0x3B1));
        if (wcslen(*lc_c_str) > 0) {
          if (result.empty()) result = str.substr(0, it - begin);
          result.push_back(L'\\');
          result += *lc_c_str;
        }
      } else {
        if (!result.empty()) result.push_back(ch);
      }
    }

    if (!result.empty()) str = std::move(result);
  };

  std::wstring result;
  result = L"{";
  result += protect_subscript(this->label());
  if (this->has_proto_indices()) {
    result += L"^{";
    for (const auto& pi : this->proto_indices()) {
      result += pi.to_latex();
    }
    result += L"}";
  }
  result += L"}";
  greek_characters_to_latex(result);
  return result;
}

}  // namespace sequant
