//
// Created by Eduard Valeyev on 4/30/20.
//

#include "../../SeQuant/core/index.hpp"

#include "../../SeQuant/core/context.hpp"

#include "../../SeQuant/core/wstring.hpp"

#include "../../SeQuant/core/latex.hpp"

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
  return utf_to_latex(result);
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
