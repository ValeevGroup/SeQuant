//
// Created by Eduard Valeyev on 8/25/25.
//

#ifndef SEQUANT_CORE_OPTIONS_HPP
#define SEQUANT_CORE_OPTIONS_HPP

#include <string>

namespace sequant {

class Index;

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

/// @brief options that control behavior of `canonicalize()`
struct CanonicalizeOptions {
  /// TN canonicalization method
  CanonicalizationMethod method = CanonicalizationMethod::Topological;
  /// specifies named indices; by default all indices that appear only once are
  /// deduced to be named, but this may be misleading if e.g. single
  /// summed-over dummy index appears in an expression
  std::optional<std::initializer_list<Index>> named_indices = std::nullopt;
  /// whether to ignore the labels of named indices. Setting
  /// to false will cause named indices to be treated as equivalent slots, which
  /// the result to be independent of their labels. This does not make sense in
  /// contexts where labels are meaningful, e.g. when canonicalizing sum of
  /// tensor networks.
  bool ignore_named_index_labels = true;

  static CanonicalizeOptions default_options();
  CanonicalizeOptions copy_and_set_method(CanonicalizationMethod method) const;

  friend constexpr bool operator==(const CanonicalizeOptions& a,
                                   const CanonicalizeOptions& b) {
    return a.method == b.method;
  }
};

/// @brief options that control behavior of `simplify()`
/// @note this is a superset of CanonicalizeOptions
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
