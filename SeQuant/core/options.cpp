//
// Created by Eduard Valeyev on 8/25/25.
//

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/options.hpp>

namespace sequant {

CanonicalizationMethod operator&(CanonicalizationMethod m1,
                                 CanonicalizationMethod m2) {
  return static_cast<CanonicalizationMethod>(static_cast<int>(m1) &
                                             static_cast<int>(m2));
}

CanonicalizationMethod operator|(CanonicalizationMethod m1,
                                 CanonicalizationMethod m2) {
  return static_cast<CanonicalizationMethod>(static_cast<int>(m1) |
                                             static_cast<int>(m2));
}

std::wstring to_wstring(CanonicalizationMethod m) {
  switch (m) {
    case CanonicalizationMethod::Topological:
      return L"topological";
    case CanonicalizationMethod::Lexicographic:
      return L"lexicographic";
    case CanonicalizationMethod::Complete:
      return L"complete";
  }

  SEQUANT_UNREACHABLE;
}

CanonicalizeOptions CanonicalizeOptions::default_options() {
  return sequant::get_default_context().canonicalization_options().value_or(
      CanonicalizeOptions{});
}

CanonicalizeOptions CanonicalizeOptions::copy_and_set(
    CanonicalizationMethod arg) const {
  auto result = *this;
  result.method = arg;
  return result;
}

CanonicalizeOptions CanonicalizeOptions::copy_and_set(
    std::optional<std::initializer_list<Index>> arg) const {
  auto result = *this;
  result.named_indices = arg;
  return result;
}

CanonicalizeOptions CanonicalizeOptions::copy_and_set(
    IgnoreNamedIndexLabel arg) const {
  auto result = *this;
  result.ignore_named_index_labels = arg;
  return result;
}

SimplifyOptions SimplifyOptions::default_options() {
  auto result =
      sequant::get_default_context().canonicalization_options().value_or(
          CanonicalizeOptions{});
  return {result};
}
SimplifyOptions::SimplifyOptions(CanonicalizeOptions opts)
    : CanonicalizeOptions(opts) {}

}  // namespace sequant
