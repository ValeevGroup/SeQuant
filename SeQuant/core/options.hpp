//
// Created by Eduard Valeyev on 8/25/25.
//

#ifndef SEQUANT_CORE_OPTIONS_HPP
#define SEQUANT_CORE_OPTIONS_HPP

#include <string>

namespace sequant {
/// canonicalization methods
enum class CanonicalizationMethod {
  /// Enables use of expression topology in the canonicalization, may be
  /// expensive.
  /// @note Canonicalization of tensor networks must take into account the
  ///       topology of the TN; this enables the use of canonical graph
  ///       sort of all TN elements to produces complete
  ///       canonicalization. The result may be aesthetically poor.
  Topological = 0b01,
  /// Enables the use of rapid canonicalization based on lexicographic sort.
  /// @note for TN this performs lexicographic sort of tensors, slots and
  /// slot bundles, and indices. Produces aesthetically pleasing result,
  /// but incomplete if some tensors are identical
  Lexicographic = 0b10,
  /// Enables use of Topological, then Lexicographic. Guaranteed to work even
  /// if some tensors are identical, and the result is aesthetically pleasing.
  Complete = Topological | Lexicographic,
  /// Lexicographic = quick-and-dirty
  Rapid = Lexicographic
};

bool operator&(CanonicalizationMethod m1, CanonicalizationMethod m2);
std::wstring to_wstring(CanonicalizationMethod m);

struct CanonicalizeOptions {
  CanonicalizationMethod method = CanonicalizationMethod::Topological;
  static CanonicalizeOptions default_options();

  friend constexpr bool operator==(const CanonicalizeOptions& a,
                                   const CanonicalizeOptions& b) {
    return a.method == b.method;
  }
};

struct SimplifyOptions : public CanonicalizeOptions {
  static SimplifyOptions default_options();

  friend constexpr bool operator==(const SimplifyOptions& a,
                                   const SimplifyOptions& b) {
    return static_cast<const CanonicalizeOptions&>(a) ==
           static_cast<const CanonicalizeOptions&>(b);
  }
};

}  // namespace sequant

#endif  // SEQUANT_CORE_OPTIONS_HPP
