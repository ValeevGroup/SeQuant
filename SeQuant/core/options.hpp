//
// Created by Eduard Valeyev on 8/25/25.
//

#ifndef SEQUANT_CORE_OPTIONS_HPP
#define SEQUANT_CORE_OPTIONS_HPP

#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/aggregate.hpp>

#include <optional>
#include <string>
#include <vector>

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

CanonicalizationMethod operator&(CanonicalizationMethod m1,
                                 CanonicalizationMethod m2);
CanonicalizationMethod operator|(CanonicalizationMethod m1,
                                 CanonicalizationMethod m2);
std::wstring to_wstring(CanonicalizationMethod m);

/// @brief options that control behavior of `canonicalize()`
struct CanonicalizeOptions {
  SEQUANT_DESIGNATED_INIT_ONLY;
  enum class IgnoreNamedIndexLabel : bool { Yes = true, No = false };
  enum class FoldKramers : bool { Yes = true, No = false };
  enum class FoldKramersEvalLeaves : bool { Yes = true, No = false };

  /// TN canonicalization method
  /// @internal
  /// TODO revert to CanonicalizationMethod::Topological once issues with
  /// external index handling exemplified by
  /// https://github.com/ValeevGroup/SeQuant/pull/406
  /// https://github.com/ValeevGroup/SeQuant/issues/426
  /// and https://github.com/ValeevGroup/SeQuant/issues/287
  /// is fixed
  /// @endinternal
  CanonicalizationMethod method = CanonicalizationMethod::Complete;
  /// specifies named indices; by default all indices that appear only once are
  /// deduced to be named, but this may be misleading if e.g. single
  /// summed-over dummy index appears in an expression
  std::optional<container::set<Index>> named_indices = std::nullopt;
  /// whether to ignore the labels of named indices. Setting
  /// to false will cause named indices to be treated as equivalent slots, which
  /// the result to be independent of their labels. This does not make sense in
  /// contexts where labels are meaningful, e.g. when canonicalizing sum of
  /// tensor networks and will be therefore ignored.
  IgnoreNamedIndexLabel ignore_named_index_labels = IgnoreNamedIndexLabel::Yes;
  /// whether to apply the Kramers (time-reversal) network fold: every
  /// connected component of KramersSymmetry::TimeReversal tensors joined by
  /// shared flavored dummy indices is given one orientation (the one with
  /// fewer down-first tensors), flipped tensors acquiring the conjugation
  /// marker and the phase; components touching a named index are pinned.
  /// Engaged only when the default context's registry has Kramers partner
  /// spaces (IndexSpaceRegistry::add_kramers_partners).
  FoldKramers fold_kramers = FoldKramers::No;
  /// whether an EvalExpr LEAF applies the single-tensor Kramers fold (a
  /// down-first leaf shares the up-first partner's hash/cache slot with a
  /// {conj, phase} transform). Off by default: the leaf provider serves the
  /// leaf's own spelling and the transform cannot relabel flavors, so this
  /// is only correct when the provider aliases the Kramers partner block
  /// (serving-level aliasing).
  FoldKramersEvalLeaves fold_kramers_eval_leaves = FoldKramersEvalLeaves::No;

  static CanonicalizeOptions default_options();
  CanonicalizeOptions copy_and_set(CanonicalizationMethod) const;
  CanonicalizeOptions copy_and_set(std::optional<container::set<Index>>) const;
  CanonicalizeOptions copy_and_set(IgnoreNamedIndexLabel) const;
  CanonicalizeOptions copy_and_set(FoldKramers) const;
  CanonicalizeOptions copy_and_set(FoldKramersEvalLeaves) const;

  friend constexpr bool operator==(const CanonicalizeOptions& a,
                                   const CanonicalizeOptions& b) {
    return a.method == b.method;
  }
};

/// @brief options that control behavior of `simplify()`
/// @note this is a superset of CanonicalizeOptions
struct SimplifyOptions : public CanonicalizeOptions {
  enum class FoldConjugatePairs : bool { Yes = true, No = false };

  /// whether simplify() folds complex-conjugate-related summand pairs
  /// (A + A* -> 2 Re(A), A - A* -> 2i Im(A)) after canonicalization;
  /// engaged only when the default context's registry contains a
  /// complex-field base space (in a real field conjugation is trivial and
  /// plain canonicalization already merges such pairs).
  /// @note default is No for now: consumers that predate RealPart/ImagPart
  ///       (tensor-network construction, evaluation) ingests the folded
  ///       nodes via EvalOp::RealPart/ImagPart, so the fold is on by
  ///       default.
  FoldConjugatePairs fold_conjugate_pairs = FoldConjugatePairs::Yes;

  // the base overloads are hidden by the FoldConjugatePairs overload below
  using CanonicalizeOptions::copy_and_set;
  SimplifyOptions copy_and_set(FoldConjugatePairs) const;

  static SimplifyOptions default_options();
  SimplifyOptions(CanonicalizeOptions opts);

  friend constexpr bool operator==(const SimplifyOptions& a,
                                   const SimplifyOptions& b) {
    return static_cast<const CanonicalizeOptions&>(a) ==
           static_cast<const CanonicalizeOptions&>(b);
  }
};

}  // namespace sequant

#endif  // SEQUANT_CORE_OPTIONS_HPP
