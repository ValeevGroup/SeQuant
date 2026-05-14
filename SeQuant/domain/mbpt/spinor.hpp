//
// Spinor decomposition passes for relativistic 2-component theories.
//
// This module provides two symbolic transforms over tensor expressions:
//
//   1. kramers_trace(expr): rewrite a closed-shell expression with
//      all-spinor (Kramers::any) indices into a sum over canonical
//      Kramers-block representatives (one per orbit under
//      {column-swap, Hermitian, time-reversal-with-(-1)^k}). Folds
//      the 2^n raw Kramers blocks of a rank-n tensor into ~5 unique
//      classes for rank-4, with conjugation/sign tracking baked into
//      the per-tensor lookup.
//
//   2. quaternion_decompose(expr): (Phase 4 — not yet implemented)
//      further decompose Kramers-up indices into the four real
//      (s, x, y, z) quaternion components (Helmich-Paris, Repisky,
//      Visscher 2016 Eq. 21), aggressively folding via Pauli/quaternion
//      algebra to yield ≤ ~10 surviving terms for an MP2-style input.
//
// Both passes operate on the canonical SeQuant Tensor representation
// and consume its `Kramers` index quantum number (defined in
// SeQuant/domain/mbpt/space_qns.hpp).
//

#ifndef SEQUANT_DOMAIN_MBPT_SPINOR_HPP
#define SEQUANT_DOMAIN_MBPT_SPINOR_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/space_qns.hpp>

#include <string>
#include <string_view>
#include <type_traits>

namespace sequant::mbpt {

// =====================================================================
// Kramers index annotations.
//
// Convention (mirrors the spin annotations using ↑/↓):
//   - L'⇑' (U+21D1) = Kramers-up
//   - L'⇓' (U+21D3) = Kramers-down
// Distinct from spin's ↑/↓ so the two QN families compose without
// label collision.
// =====================================================================

/// Extract the Kramers QN from a QuantumNumbersAttr.
inline Kramers to_kramers(const QuantumNumbersAttr& t) {
  SEQUANT_ASSERT((t.to_int32() & mask_v<Kramers>) != 0);
  return static_cast<Kramers>(t.to_int32() & mask_v<Kramers>);
}

/// Strip Kramers QN bits from a QuantumNumbersAttr (mirrors
/// `spinannotation_remove(QuantumNumbersAttr)` in `spin.hpp`).
inline QuantumNumbersAttr kramers_annotation_remove(
    const QuantumNumbersAttr& t) {
  static_assert((~(~mask_v<Kramers> & ~bitset::reserved) & ~bitset::reserved) ==
                    mask_v<Kramers>,
                "Kramers bitmask uses reserved bits");
  return t.intersection(QuantumNumbersAttr(~mask_v<Kramers> & ~bitset::reserved));
}

/// Strip a trailing ⇑/⇓ Kramers annotation from a label.
template <typename WS, typename = std::enable_if_t<
                          meta::is_wstring_or_view_v<std::decay_t<WS>>>>
std::wstring kramers_annotation_remove(WS&& label) {
  auto view = to_basic_string_view(label);
  const auto has_annotation = !view.empty() &&
                              (view.back() == L'⇑' || view.back() == L'⇓');
  return std::wstring{view.data(),
                      view.data() + view.size() - (has_annotation ? 1 : 0)};
}

/// Add a Kramers annotation (⇑ or ⇓) to a label (must not already
/// carry one). Kramers::any returns the label unchanged.
template <typename WS, typename = std::enable_if_t<
                          meta::is_wstring_or_view_v<std::decay_t<WS>>>>
std::wstring kramers_annotation_add(WS&& label, Kramers k) {
  [[maybe_unused]] auto view = to_basic_string_view(label);
  SEQUANT_ASSERT(view.empty() ||
                 (view.back() != L'⇑' && view.back() != L'⇓'));
  switch (k) {
    case Kramers::any:
      return std::wstring(std::forward<WS>(label));
    case Kramers::up:
      return std::wstring(std::forward<WS>(label)) + L'⇑';
    case Kramers::down:
      return std::wstring(std::forward<WS>(label)) + L'⇓';
    case Kramers::null:
      SEQUANT_ABORT("Invalid Kramers quantum number");
  }
  SEQUANT_UNREACHABLE;
}

/// Replace any existing ⇑/⇓ annotation with @p k (or strip if Kramers::any).
template <typename WS, typename = std::enable_if_t<
                          meta::is_wstring_or_view_v<std::decay_t<WS>>>>
std::wstring kramers_annotation_replace(WS&& label, Kramers k) {
  auto label_kf = kramers_annotation_remove(std::forward<WS>(label));
  return kramers_annotation_add(label_kf, k);
}

/// Construct a Kramers-up Index from @p idx (preserves type, ordinal,
/// proto-indices; sets the Kramers QN bits and decorates the label).
[[nodiscard]] Index make_kramers_up(const Index& idx);

/// Construct a Kramers-down Index from @p idx.
[[nodiscard]] Index make_kramers_dn(const Index& idx);

/// Construct a Kramers-free (Kramers::any) Index from @p idx.
[[nodiscard]] Index make_kramers_free(const Index& idx);

/// Apply an Index→Index replacement map to all tensors in @p expr.
/// Mirrors `append_spin` in `spin.hpp`; provided as a separate symbol
/// so callers don't drag in the spin tracer.
ExprPtr append_kramers(const ExprPtr& expr,
                       const container::map<Index, Index>& index_replacements);

// =====================================================================
// Complex-arithmetic operators on tensor expressions.
//
// Groundwork for the Kramers tracer (which uses TRS to rewrite barred
// configurations as conjugates of unbarred ones) and the upcoming
// quaternion decomposer (which separates real/imaginary parts).
//
// Convention: complex conjugation of a tensor is encoded as a `*`
// suffix on the tensor's label (e.g., `g` → `g*`, `t` → `t*`).
// Evaluator-side dispatch (e.g., MPQC's SQ_MP2 yielder) translates the
// suffix into a `.conj()` call on the underlying numeric tensor.
// =====================================================================

/// Test whether a tensor's label encodes complex conjugation
/// (trailing `*` suffix per the convention above).
inline bool has_conj_suffix(std::wstring_view label) {
  return !label.empty() && label.back() == L'*';
}

/// Append `*` to a tensor's label (or strip if already present, since
/// (z*)* = z). Pure label rewrite — does NOT affect Symmetry tags or
/// IndexSpace QNs. Caller's responsibility to ensure the underlying
/// numeric tensor's complex conjugate is what's wanted.
[[nodiscard]] std::wstring toggle_conj_suffix(std::wstring_view label);

/// Symbolic complex conjugation of an expression.
///
/// Recursively walks @p expr and applies:
///   - Tensor: toggle its `*` suffix (z → z*, z* → z).
///   - Constant: numeric complex conjugate.
///   - Variable: not yet supported (throws).
///   - Product: conjugate the scalar AND each factor. (Π zᵢ)* = Π zᵢ*.
///   - Sum: conjugate each summand. (Σ zᵢ)* = Σ zᵢ*.
///
/// `conjugate(conjugate(expr))` recovers @p expr.
[[nodiscard]] ExprPtr conjugate(const ExprPtr& expr);

/// Symbolic real part: returns `(expr + conjugate(expr)) / 2`.
/// Currently materialized as that two-term Sum (no native `Re` node);
/// downstream callers can apply numerical optimizations after evaluation.
[[nodiscard]] ExprPtr real_part(const ExprPtr& expr);

/// Symbolic imaginary part: returns `(expr - conjugate(expr)) / (2i)`.
[[nodiscard]] ExprPtr imaginary_part(const ExprPtr& expr);

// =====================================================================
// Kramers tracer.
// =====================================================================

/// Rewrite an all-spinor tensor expression into a sum over canonical
/// Kramers-block representatives. Implements the closed-shell relativistic
/// MP2/CC algebra: enumerate per-tensor Kramers blocks, apply per-tensor
/// TRS canonicalization (with (-1)^k sign), and fold via the standard
/// canonicalize+rapid_simplify pipeline. Single entry handles MP2 energy,
/// CC residuals, and arbitrary operator expressions; expands any
/// antisymmetrizer (@c A) tensors internally before tracing.
///
/// @param expr expression with indices in `IndexSpace`s carrying
///             `Kramers::any` (i.e., not yet specialized to up/down)
/// @return     a Sum<Product> with canonical-orbit-block tensors;
///             tensor labels with conjugation get a `*` suffix
ExprPtr kramers_trace(const ExprPtr& expr);

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_SPINOR_HPP
