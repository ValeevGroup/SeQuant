#ifndef SEQUANT_CORE_EVAL_CANON_TRANSFORM_HPP
#define SEQUANT_CORE_EVAL_CANON_TRANSFORM_HPP

#include <cstddef>
#include <cstdint>

namespace sequant {

///
/// \brief Canonicalization byproduct mapping a node's _cached_ canonical result
///        to the value the node denotes. Applied on retrieval
///        (see apply_canon_transform in eval.hpp) and excluded from the node's
///        own (slot) hash, so that every transform of one canonical value
///        shares its slot. conj/braket_swap do enter the parent's structural
///        hash via structural_salt(): a uniform conjugation hoists, a mixed
///        one salts (see doc/developer/conjugation.rst, "The eval boundary").
struct CanonTransform {
  std::int8_t phase = 1;     ///< +/-1 linear byproduct (antisymmetric reorder)
  bool conj = false;         ///< elementwise complex conjugation
  bool braket_swap = false;  ///< bra<->ket transposition of the canonical slots

  [[nodiscard]] constexpr bool trivial() const noexcept {
    return phase == 1 && !conj && !braket_swap;
  }
  /// salt for the _parent_'s hash combination: conj/swap only -- phase is
  /// multiplicatively hoistable and never enters structural identity
  /// (a product folds its children's phases into its own transform; a sum
  /// hoists a uniform phase and salts a mixed one with phase_salt)
  [[nodiscard]] constexpr std::size_t structural_salt() const noexcept {
    return (conj ? 1u : 0u) | (braket_swap ? 2u : 0u);
  }
  /// salt a sum combines into a negated summand's hash when its summands'
  /// phases are mixed (disjoint from the structural_salt() bits)
  static constexpr std::size_t phase_salt = 4u;
  friend constexpr bool operator==(CanonTransform, CanonTransform) = default;
};

/// @return whether @p tr hoists out of a product or a sum node: elementwise
/// conjugation distributes over contraction and addition
/// (`(A·B)꙳ = A꙳·B꙳`, `(Σ T)꙳ = Σ T꙳`), while a bra<->ket exchange respells
/// the node's own result -- the partition its placeholder is built from --
/// so a transform carrying one salts the parent's hash instead
[[nodiscard]] constexpr bool hoistable(CanonTransform tr) noexcept {
  return tr.conj && !tr.braket_swap;
}

/// composition of two transforms: phases multiply, conj/swap compose as Z2
[[nodiscard]] constexpr CanonTransform compose(CanonTransform a,
                                               CanonTransform b) noexcept {
  return {static_cast<std::int8_t>(a.phase * b.phase), a.conj != b.conj,
          a.braket_swap != b.braket_swap};
}

}  // namespace sequant

#endif  // SEQUANT_CORE_EVAL_CANON_TRANSFORM_HPP
