//
// Created by Eduard Valeyev on 4/30/20.
//

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <format>
#include <unordered_map>

namespace sequant {

std::size_t Index::min_tmp_index() noexcept {
  return get_default_context().first_dummy_index_ordinal();
}

void Index::reset_tmp_index() noexcept {
  tmp_index_accessor() = min_tmp_index() - 1;
}

std::string Index::to_latex() const noexcept {
  std::string protos{};
  const auto instance = basis_.unnamed_instance();
  if (has_proto_indices() || instance) {
    protos += "^{";
    for (auto&& pidx : proto_indices()) protos += pidx.to_latex();
    if (instance) protos += basis_.instance_suffix();
    protos += "}";
  }
  std::string sfx;
  if (ordinal_)
    sfx = std::format("_{}", *ordinal_ < 10 ? std::to_string(*ordinal_)
                                            : std::format("{{{}}}", *ordinal_));
  return std::format("{{{}{}{}}}", io::latex::utf_to_string(basis_.base_key()),
                     sfx, protos);
}

std::string Index::ascii_label() const {
  static const std::unordered_map<char32_t, std::string> greek_to_english_name =
      {{U'Α', "ALPHA"},   {U'Β', "BETA"},  {U'Γ', "GAMMA"},   {U'Δ', "DELTA"},
       {U'Ε', "EPSILON"}, {U'Ζ', "ZETA"},  {U'Η', "ETA"},     {U'Θ', "THETA"},
       {U'Ι', "IOTA"},    {U'Κ', "KAPPA"}, {U'Λ', "LAMBDA"},  {U'Μ', "MU"},
       {U'Ν', "NU"},      {U'Ξ', "XI"},    {U'Ο', "OMICRON"}, {U'Π', "PI"},
       {U'Ρ', "RHO"},     {U'Σ', "SIGMA"}, {U'Τ', "TAU"},     {U'Υ', "UPSILON"},
       {U'Φ', "PHI"},     {U'Χ', "CHI"},   {U'Ψ', "PSI"},     {U'Ω', "OMEGA"},
       {U'α', "alpha"},   {U'β', "beta"},  {U'γ', "gamma"},   {U'δ', "delta"},
       {U'ε', "epsilon"}, {U'ζ', "zeta"},  {U'η', "eta"},     {U'θ', "theta"},
       {U'ι', "iota"},    {U'κ', "kappa"}, {U'λ', "lambda"},  {U'μ', "mu"},
       {U'ν', "nu"},      {U'ξ', "xi"},    {U'ο', "omicron"}, {U'π', "pi"},
       {U'ρ', "rho"},     {U'σ', "sigma"}, {U'τ', "tau"},     {U'υ', "upsilon"},
       {U'φ', "phi"},     {U'χ', "chi"},   {U'ψ', "psi"},     {U'ω', "omega"}};

  std::string label_ascii;
  for (const char32_t ch : to_utf32(this->label())) {
    if (ch == U'↑') {
      label_ascii.push_back('a');
    } else if (ch == U'↓') {
      label_ascii.push_back('b');
    } else if (auto pos = greek_to_english_name.find(ch);
               pos != greek_to_english_name.end()) {
      label_ascii.append(pos->second);
    } else {
      label_ascii.append(to_utf8(ch));
    }
  }
  return label_ascii;
}

std::string Index::to_string() const { return std::string(this->label()); }

std::shared_ptr<const IndexBasisRegistry>
Index::obtain_default_index_registry() {
  return get_default_context().index_basis_registry();
}

Index Index::replace_basis_instance(
    IndexBasis::optional_instance basis_instance) const {
  Index result(*this);
  result.basis_ =
      default_registry_resolved(IndexBasis(space(), basis_instance));
  result.label_.reset();
  result.full_label_.reset();
  return result;
}

void Index::rebase(IndexSpace space) {
  // a name and the metadata belong to the instance of the old space; in
  // another space the instance is whatever the default registry knows it as
  if (space != basis_.space())
    basis_ = default_registry_resolved(
        IndexBasis(std::move(space), basis_.basis_instance()));
  label_.reset();
  full_label_.reset();
}

}  // namespace sequant
