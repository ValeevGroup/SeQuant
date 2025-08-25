//
// Created by Eduard Valeyev on 8/25/25.
//

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/options.hpp>

namespace sequant {

bool logical_and(CanonicalizationMethod m1, CanonicalizationMethod m2) {
  return (static_cast<int>(m1) & static_cast<int>(m2)) != 0;
}

std::wstring to_wstring(CanonicalizationMethod m) {
  switch (m) {
    case CanonicalizationMethod::Topological:
      return L"topological";
    case CanonicalizationMethod::Lexicographic:
      return L"lexicographic";
    case CanonicalizationMethod::Complete:
      return L"complete";
    default:
      abort();
  }
}

CanonicalizeOptions CanonicalizeOptions::default_options() {
  return sequant::get_default_context().canonicalization_options().value_or(
      CanonicalizeOptions{});
}

SimplifyOptions SimplifyOptions::default_options() {
  auto result =
      sequant::get_default_context().canonicalization_options().value_or(
          CanonicalizeOptions{});
  return {result};
}

}  // namespace sequant
