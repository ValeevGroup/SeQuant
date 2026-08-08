#ifndef SEQUANT_CORE_OPTIMIZE_GA_INTERMEDIATE_LABEL_HPP
#define SEQUANT_CORE_OPTIMIZE_GA_INTERMEDIATE_LABEL_HPP

// The naming convention for the GA optimizer's shared intermediates, and
// nothing else. Deliberately standalone (only <string_view>): consumers that
// merely need to RECOGNIZE an emitted intermediate -- backends deciding whether
// a tensor is a user array or a factorizer temporary -- must not have to
// include emit.hpp, which drags in the whole evaluator (fitness.hpp -> genome
// codec -> key_table.hpp).

#include <string_view>

namespace sequant::opt::ga {

/// Label prefix of named GA intermediates ("IGA1", "IGA2", ...).
inline constexpr std::wstring_view named_intermediate_prefix = L"IGA";

/// Whether \p label names a GA intermediate (i.e. starts with
/// named_intermediate_prefix followed by digits).
inline bool is_named_intermediate(std::wstring_view label) {
  auto const& p = named_intermediate_prefix;
  if (label.size() <= p.size() || label.substr(0, p.size()) != p) return false;
  for (auto c : label.substr(p.size()))
    if (c < L'0' || c > L'9') return false;
  return true;
}

}  // namespace sequant::opt::ga

#endif  // SEQUANT_CORE_OPTIMIZE_GA_INTERMEDIATE_LABEL_HPP
