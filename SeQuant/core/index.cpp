//
// Created by Eduard Valeyev on 4/30/20.
//

#include "../../SeQuant/core/index.hpp"

#include "../../SeQuant/core/context.hpp"

#include "../../SeQuant/core/wstring.hpp"

namespace sequant {

const std::size_t Index::min_tmp_index() {
  return get_default_context().first_dummy_index_ordinal();
}

void Index::reset_tmp_index() { tmp_index_accessor() = min_tmp_index() - 1; }

std::wstring Index::to_latex() const {
  auto protect_subscript = [](const std::wstring_view str) {
    auto subsc_pos = str.rfind(L'_');
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

std::string Index::ascii_label() const {
  static const std::unordered_map<wchar_t, std::string> greek_to_english_name =
      {{L'Α', "ALPHA"},   {L'Β', "BETA"},  {L'Γ', "GAMMA"},   {L'Δ', "DELTA"},
       {L'Ε', "EPSILON"}, {L'Ζ', "ZETA"},  {L'Η', "ETA"},     {L'Θ', "THETA"},
       {L'Ι', "IOTA"},    {L'Κ', "KAPPA"}, {L'Λ', "LAMBDA"},  {L'Μ', "MU"},
       {L'Ν', "NU"},      {L'Ξ', "XI"},    {L'Ο', "OMICRON"}, {L'Π', "PI"},
       {L'Ρ', "RHO"},     {L'Σ', "SIGMA"}, {L'Τ', "TAU"},     {L'Υ', "UPSILON"},
       {L'Φ', "PHI"},     {L'Χ', "CHI"},   {L'Ψ', "PSI"},     {L'Ω', "OMEGA"},
       {L'α', "alpha"},   {L'β', "beta"},  {L'γ', "gamma"},   {L'δ', "delta"},
       {L'ε', "epsilon"}, {L'ζ', "zeta"},  {L'η', "eta"},     {L'θ', "theta"},
       {L'ι', "iota"},    {L'κ', "kappa"}, {L'λ', "lambda"},  {L'μ', "mu"},
       {L'ν', "nu"},      {L'ξ', "xi"},    {L'ο', "omicron"}, {L'π', "pi"},
       {L'ρ', "rho"},     {L'σ', "sigma"}, {L'τ', "tau"},     {L'υ', "upsilon"},
       {L'φ', "phi"},     {L'χ', "chi"},   {L'ψ', "psi"},     {L'ω', "omega"}};

  std::wstring label(label_);

  std::replace(label.begin(), label.end(), L'↑', L'a');
  std::replace(label.begin(), label.end(), L'↓', L'b');
  std::string label_ascii;
  for (auto it = label.begin(); it != label.end(); ++it) {
    auto pos = greek_to_english_name.find(*it);
    if (pos != greek_to_english_name.end()) {
      label_ascii.append(pos->second);
    } else {
      label_ascii.push_back(*it);
    }
  }
  return label_ascii;
}

std::string Index::to_string() const {
  return sequant::to_string(this->label());
}

}  // namespace sequant
