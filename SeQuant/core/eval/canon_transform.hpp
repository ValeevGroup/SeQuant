#ifndef SEQUANT_CORE_EVAL_CANON_TRANSFORM_HPP
#define SEQUANT_CORE_EVAL_CANON_TRANSFORM_HPP

#include <cstddef>
#include <cstdint>

namespace sequant {

///
/// \brief Canonicalization byproduct mapping a node's CACHED canonical result
///        to the value the node denotes. Applied on retrieval
///        (see apply_canon_transform in eval.hpp); excluded from the node's
///        own (slot) hash, exactly as the former standalone canon_phase was.
///        conj/braket_swap DO enter the parent's structural hash via
///        structural_salt() -- see the design spec's uniform-conj-hoists /
///        mixed-conj-salts rule
///        (doc/dev/specs/2026-09-01-lazy-conj-eval-design.md).
struct CanonTransform {
  std::int8_t phase = 1;     ///< +/-1 linear byproduct (antisymmetric reorder)
  bool conj = false;         ///< elementwise complex conjugation
  bool braket_swap = false;  ///< bra<->ket transposition of the canonical slots

  [[nodiscard]] constexpr bool trivial() const noexcept {
    return phase == 1 && !conj && !braket_swap;
  }
  /// salt for the PARENT's hash combination: conj/swap only -- phase is
  /// multiplicatively hoistable and never enters structural identity
  [[nodiscard]] constexpr std::size_t structural_salt() const noexcept {
    return (conj ? 1u : 0u) | (braket_swap ? 2u : 0u);
  }
  friend constexpr bool operator==(CanonTransform, CanonTransform) = default;
};

/// composition of two transforms: phases multiply, conj/swap compose as Z2
[[nodiscard]] constexpr CanonTransform compose(CanonTransform a,
                                               CanonTransform b) noexcept {
  return {static_cast<std::int8_t>(a.phase * b.phase), a.conj != b.conj,
          a.braket_swap != b.braket_swap};
}

}  // namespace sequant

#endif  // SEQUANT_CORE_EVAL_CANON_TRANSFORM_HPP
