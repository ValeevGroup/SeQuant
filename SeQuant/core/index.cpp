//
// Created by Eduard Valeyev on 4/30/20.
//

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/latex.hpp>
#include <SeQuant/core/wstring.hpp>

#include <format>
#include <unordered_map>

namespace sequant {

const std::size_t Index::min_tmp_index() {
  return get_default_context().first_dummy_index_ordinal();
}

void Index::reset_tmp_index() { tmp_index_accessor() = min_tmp_index() - 1; }

std::wstring Index::to_latex() const {
  std::wstring protos{};
  if (has_proto_indices()) {
    protos += L"^{";
    for (auto&& pidx : proto_indices()) protos += pidx.to_latex();
    protos += L"}";
  }
  auto [lbl, sfx_] = split_label();

  std::wstring sfx{};
  if (!sfx_.empty())
    sfx = std::format(L"_{}",
                      sfx_.size() == 1 ? sfx_ : std::format(L"{{{}}}", sfx_));

  return std::format(L"{{{}{}{}}}", utf_to_latex(lbl), sfx, protos);
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

  std::wstring label(this->label());

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
