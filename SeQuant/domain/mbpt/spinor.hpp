//
// Spinor decomposition passes for relativistic 2-component theories.
//
// This module provides symbolic transforms over tensor expressions:
//
//   1. kramers_trace(expr): rewrite a closed-shell expression with
//      all-spinor (Kramers::any) indices into a sum over canonical
//      Kramers-block representatives (one per orbit under
//      {column-swap, Hermitian, time-reversal-with-(-1)^k}), with
//      conjugation/sign tracking baked into the per-tensor lookup.
//
//   2. quaternion_decompose(expr): (not yet implemented) further
//      decompose Kramers-up indices into the four real (s, x, y, z)
//      quaternion components (Helmich-Paris, Repisky, Visscher 2016
//      Eq. 21), aggressively folding via Pauli/quaternion algebra to a
//      compact surviving-term set.
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

/// @brief Extracts the Kramers quantum number from a QuantumNumbersAttr.
/// @param t the quantum-numbers attribute to read
/// @return the Kramers QN encoded in @p t
/// @pre @p t has at least one Kramers bit set
inline Kramers to_kramers(const QuantumNumbersAttr& t) {
  SEQUANT_ASSERT((t.to_int32() & mask_v<Kramers>) != 0);
  return static_cast<Kramers>(t.to_int32() & mask_v<Kramers>);
}

/// @brief Strips the Kramers QN bits from a QuantumNumbersAttr.
/// @param t the quantum-numbers attribute to clear
/// @return @p t with all Kramers bits unset
/// @sa spinannotation_remove(const QuantumNumbersAttr&) in spin.hpp
inline QuantumNumbersAttr kramers_annotation_remove(
    const QuantumNumbersAttr& t) {
  static_assert((~(~mask_v<Kramers> & ~bitset::reserved) & ~bitset::reserved) ==
                    mask_v<Kramers>,
                "Kramers bitmask uses reserved bits");
  return t.intersection(
      QuantumNumbersAttr(~mask_v<Kramers> & ~bitset::reserved));
}

/// @brief Strips a trailing ⇑/⇓ Kramers annotation from a label.
/// @tparam WS a wide-string or wide-string-view type
/// @param label the (possibly Kramers-annotated) label
/// @return @p label with any trailing ⇑/⇓ glyph removed
template <typename WS, typename = std::enable_if_t<
                           meta::is_wstring_or_view_v<std::decay_t<WS>>>>
std::wstring kramers_annotation_remove(WS&& label) {
  auto view = to_basic_string_view(label);
  const auto has_annotation =
      !view.empty() && (view.back() == L'⇑' || view.back() == L'⇓');
  return std::wstring{view.data(),
                      view.data() + view.size() - (has_annotation ? 1 : 0)};
}

/// @brief Adds a Kramers annotation (⇑ or ⇓) to a label.
/// @tparam WS a wide-string or wide-string-view type
/// @param label the label to annotate
/// @param k the Kramers state; `Kramers::any` returns @p label unchanged
/// @return @p label decorated with the matching ⇑/⇓ glyph
/// @pre @p label does not already carry a ⇑/⇓ annotation
template <typename WS, typename = std::enable_if_t<
                           meta::is_wstring_or_view_v<std::decay_t<WS>>>>
std::wstring kramers_annotation_add(WS&& label, Kramers k) {
  [[maybe_unused]] auto view = to_basic_string_view(label);
  SEQUANT_ASSERT(view.empty() || (view.back() != L'⇑' && view.back() != L'⇓'));
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

/// @brief Replaces any existing ⇑/⇓ annotation on a label with @p k.
/// @tparam WS a wide-string or wide-string-view type
/// @param label the (possibly Kramers-annotated) label
/// @param k the Kramers state to set; `Kramers::any` strips the annotation
/// @return @p label re-annotated with @p k
template <typename WS, typename = std::enable_if_t<
                           meta::is_wstring_or_view_v<std::decay_t<WS>>>>
std::wstring kramers_annotation_replace(WS&& label, Kramers k) {
  auto label_kf = kramers_annotation_remove(std::forward<WS>(label));
  return kramers_annotation_add(label_kf, k);
}

/// @brief Constructs a Kramers-up Index from @p idx.
/// @param idx the source index (its type, ordinal and proto-indices are
///            preserved)
/// @return a copy of @p idx with the Kramers-up QN bit set and the label
///         decorated with ⇑
[[nodiscard]] Index make_kramers_up(const Index& idx);

/// @brief Constructs a Kramers-down Index from @p idx.
/// @param idx the source index
/// @return a copy of @p idx with the Kramers-down QN bit set and the
///         label decorated with ⇓
[[nodiscard]] Index make_kramers_dn(const Index& idx);

/// @brief Constructs a Kramers-free (`Kramers::any`) Index from @p idx.
/// @param idx the source index
/// @return a copy of @p idx carrying the `Kramers::any` QN
[[nodiscard]] Index make_kramers_free(const Index& idx);

/// @brief Applies an Index→Index replacement map to all tensors in an
///        expression.
/// @param expr the expression to rewrite
/// @param index_replacements the Index→Index substitution map
/// @return a copy of @p expr with every tensor index substituted
/// @sa append_spin in spin.hpp — this is the Kramers analogue, provided
///     as a separate symbol so callers need not pull in the spin tracer
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
// Evaluator-side dispatch translates the suffix into a `.conj()` call
// on the underlying numeric tensor.
// =====================================================================

/// @brief Tests whether a tensor label encodes complex conjugation.
/// @param label the tensor label to inspect
/// @return true iff @p label carries a trailing `*` suffix
inline bool has_conj_suffix(std::wstring_view label) {
  return !label.empty() && label.back() == L'*';
}

/// @brief Toggles the `*` conjugation suffix on a tensor label.
/// @param label the tensor label
/// @return @p label with a `*` appended, or stripped if already present
///         (since `(z*)* = z`)
/// @note pure label rewrite — does not touch Symmetry tags or IndexSpace
///       QNs; the caller must ensure the underlying numeric tensor's
///       complex conjugate is what is wanted
[[nodiscard]] std::wstring toggle_conj_suffix(std::wstring_view label);

/// @brief Symbolic complex conjugation of an expression.
///
/// Recursively walks @p expr and applies:
///   - Tensor: toggles its `*` suffix (`z → z*`, `z* → z`);
///   - Constant: numeric complex conjugate;
///   - Variable: not yet supported (throws);
///   - Product: conjugates the scalar and each factor, `(Π zᵢ)* = Π zᵢ*`;
///   - Sum: conjugates each summand, `(Σ zᵢ)* = Σ zᵢ*`.
///
/// @param expr the expression to conjugate
/// @return the symbolic complex conjugate of @p expr
/// @note `conjugate(conjugate(expr))` recovers @p expr
[[nodiscard]] ExprPtr conjugate(const ExprPtr& expr);

/// @brief Symbolic real-part wrapper: `RealPart(E)` represents `Re(E)`
///        for a scalar-valued Expr `E`.
/// @note real-valued itself and idempotent under conjugation; the
///       evaluator extracts @c inner, evaluates it, and takes the real
///       part of the resulting scalar
class RealPart : public sequant::Expr {
 public:
  RealPart() = delete;
  RealPart(const RealPart&) = default;
  RealPart(RealPart&&) = default;
  ~RealPart() override = default;

  /// @brief Wraps a scalar-valued expression in a real-part marker.
  /// @param inner the scalar-valued expression
  /// @pre @p inner is non-null and `inner->is_scalar()`
  explicit RealPart(ExprPtr inner) : inner_{std::move(inner)} {
    SEQUANT_ASSERT(inner_);
    SEQUANT_ASSERT(inner_->is_scalar());
  }

  /// @return the wrapped expression
  const ExprPtr& inner() const { return inner_; }

  bool is_scalar() const override { return true; }

  type_id_type type_id() const override { return get_type_id<RealPart>(); }

  ExprPtr clone() const override { return ex<RealPart>(inner_->clone()); }

  /// @brief Adjoint of `Re(E)` — a no-op, since `Re(E)` is real.
  void adjoint() override {}

  std::wstring to_latex() const override {
    return L"\\Re\\left[" + inner_->to_latex() + L"\\right]";
  }

 private:
  ExprPtr inner_;

  hash_type memoizing_hash() const override {
    auto compute = [this]() {
      auto v = hash::value(*inner_);
      hash::combine(v, std::size_t{0xC0FFEE01ull});  // RealPart tag
      return v;
    };
    if (!hash_value_) hash_value_ = compute();
    return *hash_value_;
  }

  bool static_equal(const Expr& that) const override {
    return *inner_ == *static_cast<const RealPart&>(that).inner_;
  }

  bool static_less_than(const Expr& that) const override {
    return *inner_ < *static_cast<const RealPart&>(that).inner_;
  }
};

/// @brief Symbolic imaginary-part wrapper: `ImagPart(E)` represents
///        `Im(E)` for a scalar-valued Expr `E`.
/// @note real-valued and idempotent under conjugation
class ImagPart : public sequant::Expr {
 public:
  ImagPart() = delete;
  ImagPart(const ImagPart&) = default;
  ImagPart(ImagPart&&) = default;
  ~ImagPart() override = default;

  /// @brief Wraps a scalar-valued expression in an imaginary-part marker.
  /// @param inner the scalar-valued expression
  /// @pre @p inner is non-null and `inner->is_scalar()`
  explicit ImagPart(ExprPtr inner) : inner_{std::move(inner)} {
    SEQUANT_ASSERT(inner_);
    SEQUANT_ASSERT(inner_->is_scalar());
  }

  /// @return the wrapped expression
  const ExprPtr& inner() const { return inner_; }

  bool is_scalar() const override { return true; }

  type_id_type type_id() const override { return get_type_id<ImagPart>(); }

  ExprPtr clone() const override { return ex<ImagPart>(inner_->clone()); }

  /// @brief Adjoint of `Im(E)` — a no-op, since `Im(E)` is real.
  void adjoint() override {}

  std::wstring to_latex() const override {
    return L"\\Im\\left[" + inner_->to_latex() + L"\\right]";
  }

 private:
  ExprPtr inner_;

  hash_type memoizing_hash() const override {
    auto compute = [this]() {
      auto v = hash::value(*inner_);
      hash::combine(v, std::size_t{0xC0FFEE02ull});  // ImagPart tag
      return v;
    };
    if (!hash_value_) hash_value_ = compute();
    return *hash_value_;
  }

  bool static_equal(const Expr& that) const override {
    return *inner_ == *static_cast<const ImagPart&>(that).inner_;
  }

  bool static_less_than(const Expr& that) const override {
    return *inner_ < *static_cast<const ImagPart&>(that).inner_;
  }
};

/// @brief Constructs `Re(expr)` as a RealPart wrapper.
/// @param expr the scalar-valued expression to wrap
/// @return a RealPart wrapping @p expr
[[nodiscard]] inline ExprPtr re(ExprPtr expr) {
  return ex<RealPart>(std::move(expr));
}

/// @brief Constructs `Im(expr)` as an ImagPart wrapper.
/// @param expr the scalar-valued expression to wrap
/// @return an ImagPart wrapping @p expr
[[nodiscard]] inline ExprPtr im(ExprPtr expr) {
  return ex<ImagPart>(std::move(expr));
}

/// @brief Symbolic real part of an expression.
/// @param expr the scalar-valued expression
/// @return `Re(expr)` as a RealPart wrapper Expr
/// @note replaces the older `(expr + conj(expr))/2` Sum form so callers
///       can detect RealPart structurally and dispatch real-part
///       evaluation
[[nodiscard]] inline ExprPtr real_part(const ExprPtr& expr) {
  return ex<RealPart>(expr);
}

/// @brief Symbolic imaginary part of an expression.
/// @param expr the scalar-valued expression
/// @return `Im(expr)` as an ImagPart wrapper Expr
[[nodiscard]] inline ExprPtr imaginary_part(const ExprPtr& expr) {
  return ex<ImagPart>(expr);
}

// =====================================================================
// Kramers tracer.
// =====================================================================

/// @brief Rewrites an all-spinor tensor expression into a sum over
///        canonical Kramers-block representatives.
///
/// Implements the closed-shell relativistic trace algebra: enumerate the
/// per-tensor Kramers blocks, apply per-tensor TRS canonicalization (with
/// the `(-1)^k` sign), and fold via the standard
/// `canonicalize` + `rapid_simplify` pipeline. A single entry point
/// handles energy expressions, amplitude residuals, and arbitrary
/// operator expressions; antisymmetric tensors are expanded internally
/// before tracing.
///
/// @param expr expression whose indices carry `Kramers::any` (i.e. not
///             yet specialized to up/down)
/// @return a `Sum<Product>` of canonical-orbit-block tensors; tensor
///         labels carrying conjugation get a `*` suffix
ExprPtr kramers_trace(const ExprPtr& expr);

/// @brief Flips every Kramers state (`up` ↔ `down`) in an expression.
///
/// Walks every tensor in @p expr, collects each Kramers-tagged index,
/// builds an Index→Index replacement that swaps the Kramers QN bit (and
/// the ⇑/⇓ label glyph), then applies the swap. Indices carrying
/// `Kramers::any` are left untouched.
///
/// @param expr the expression to flip
/// @return a copy of @p expr with every specialized Kramers index flipped
/// @sa fold_conj_pairs — uses this to find TRS-conjugate-pair Products
[[nodiscard]] ExprPtr flip_kramers(const ExprPtr& expr);

/// @brief Folds TRS conjugate-pair Products in a Sum.
///
/// For each summand `A`, builds `A_flipped` via flip_kramers() (the
/// full-bar TRS partner), canonicalizes it, and checks structural
/// equality against the remaining summands:
///   - `A == A_flipped` (self-conjugate after canonicalize): emits
///     `Re(A)`, since the term is necessarily real-valued;
///   - `A_flipped == B` for some unused `B`: emits `2·Re(A)` and marks
///     both `A` and `B` consumed (term count drops by one per fold);
///   - no match: keeps `A` unchanged.
///
/// Imaginary-pair detection (`A - A_flipped → 2i·Im(A)`) is supported
/// when the partner `B` carries the opposite scalar prefactor.
///
/// @param expr a `Sum` to fold (typically the output of kramers_trace())
/// @return a new `Sum`, typically with fewer summands, whose scalar
///         value equals that of @p expr
[[nodiscard]] ExprPtr fold_conj_pairs(const ExprPtr& expr);

// =====================================================================
// Cycle-based Kramers enumeration.
//
// The closed-shell spin-trace algorithm in SeQuant decomposes a Product's
// contraction graph into cycles and exploits cycle structure to count
// independent traces. We adapt the same idea for Kramers: each cycle in
// the contraction graph carries its own Kramers labelling, and the
// number of *canonical* Kramers patterns per cycle (under cyclic rotation
// + reflection) is typically << 2^L for length-L cycles — a long cycle
// visiting d distinct indices twice each has 2^d raw labellings that
// collapse to a much smaller canonical-pattern set.
// =====================================================================

/// @brief One contraction cycle in a Product's index graph.
///
/// The cycle alternates intra-tensor edges (bra[k] ↔ ket[k] of the same
/// tensor) and inter-tensor edges (the same Index on two tensors).
/// `nodes` walks the cycle in visiting order.
struct ContractionCycle {
  /// @brief One step of the cycle walk.
  struct Node {
    std::size_t tensor_idx;  ///< which factor of the Product
    std::size_t braket_pos;  ///< position in the tensor's const_braket()
    Index idx;               ///< index occupying that slot
  };
  container::svector<Node> nodes;  ///< the cycle walk, in visiting order
};

/// @brief Decomposes a Product's index-contraction structure into cycles.
/// @param product a fully-contracted scalar Product
/// @return the contraction cycles of @p product
/// @pre every Index in @p product appears in exactly two `(tensor, slot)`
///      positions; each tensor's bra[k]/ket[k] positions are paired as
///      intra-tensor edges
[[nodiscard]] container::svector<ContractionCycle> kramers_cycles(
    const sequant::Product& product);

/// @brief Diagnostic dump of the cycle structure of an expression.
///
/// Walks @p expr (a `Sum` of Products), decomposes each Product into
/// cycles, and prints each cycle's length, distinct indices visited,
/// per-cycle Kramers labelling, and the number of canonical Kramers
/// patterns under cyclic-rotation symmetry.
///
/// @param expr the expression to inspect
/// @param os the wide-character stream to write to
void kramers_cycle_dump(const ExprPtr& expr, std::wostream& os);

/// @brief Canonical Kramers label of a cycle.
/// @param c the contraction cycle
/// @return the lex-smallest rotation of @p c 's raw label string (e.g.
///         `BUUB` and `UBBU` both canonicalize to `BBUU`)
[[nodiscard]] std::wstring cycle_canonical_label(const ContractionCycle& c);

/// @brief Cycle-canonical signature of a Product.
/// @param product the Product to fingerprint
/// @return the sorted multiset of per-cycle canonical Kramers labels,
///         joined into one key string (e.g. `BBUU|UUUU`)
/// @note Products with the same signature are equivalent under per-cycle
///       cyclic-rotation symmetry
[[nodiscard]] std::wstring cycle_canonical_signature(
    const sequant::Product& product);

/// @brief Buckets Products by cycle-canonical signature and sums each
///        bucket into one representative.
///
/// Within a bucket the scalar prefactors are summed and one
/// representative is emitted. Acts independently of fold_conj_pairs()
/// and on a raw `Sum`-of-Products without any Re/Im wrappers.
///
/// @param expr a `Sum` to fold
/// @return a `Sum` with one summand per cycle-canonical signature
/// @warning NOT value-preserving in general — the cyclic-rotation
///          signature is not a complete contraction invariant, so two
///          genuinely different contractions can share it. Superseded by
///          kramers_trace_cycles(); kept only for diagnostic comparison.
[[nodiscard]] ExprPtr cycle_canonical_fold(const ExprPtr& expr);

// =====================================================================
// Standalone cycle-driven Kramers tracer.
//
// Unlike `kramers_trace` (which enumerates all 2^n per-index Kramers
// configurations then folds post-hoc), this tracer is *driven* by the
// contraction-graph cycle structure — analogous to how
// `closed_shell_spintrace` is driven by `count_cycles`.
//
// Algorithm (cycles decomposed on the ANTISYMMETRIZED form, before
// expand_antisymm):
//   1. Decompose the input Product's contraction graph into cycles.
//   2. Group cycles into structural-equivalence classes: cycles that map
//      to one another under the antisymmetrizer's index permutations
//      (structurally-identical cycles are one class — the antisymmetric
//      bra/ket index permutations swap them).
//   3. Enumerate Kramers labelings as MULTISETS within each class (not
//      ordered assignments) — this is where the reduction over naive
//      per-index 2^n enumeration comes from: equivalent cycles with
//      swapped labels are the same contraction.
//   4. For each multiset configuration: substitute Kramers labels,
//      expand_antisymm, canonicalize; attach the multiset multiplicity.
//   5. Collect into a Sum; the caller may apply `fold_conj_pairs` for
//      additional TRS conjugate-pair folding.
//
// Correct-by-construction: works from contraction topology, never guesses
// equivalences. Coexists with `kramers_trace` so the two can be
// cross-validated.
// =====================================================================

/// @brief Cycle-driven Kramers tracer.
///
/// See the block comment above for the algorithm.
///
/// @param expr a closed-shell all-spinor expression with `Kramers::any`
///             indices and (optionally) antisymmetric tensors
/// @return a `Sum<Product>` of canonical Kramers-block contractions whose
///         scalar value equals the full 2^n per-index trace
[[nodiscard]] ExprPtr kramers_trace_cycles(const ExprPtr& expr);

/// @brief Recombines direct/exchange Product pairs into antisymmetric
///        tensors — the inverse of `expand_antisymm`.
///
/// After `kramers_trace` + `fold_conj_pairs` the expression is a `Sum` of
/// fully-expanded NonSymm contractions. Many of these are the `direct`
/// and `exchange` halves of a single antisymmetrized contraction: two
/// summands recombine into one when
///   - one tensor is structurally identical between them, and
///   - the other tensor differs by a column permutation, with the
///     summand scalars in the ratio `sign(permutation)`.
/// Then `c·X·Y_direct + (sign·c)·X·Y_exchange  →  c·X·Ȳ`, where `Ȳ` is
/// `Y` re-tagged `Symmetry::Antisymm` over the Kramers-restricted index
/// space (NOT full spinor).
///
/// Applied iteratively to a fixed point, so a second pass over the
/// once-recombined terms catches doubly-antisymmetrized blocks. Re/Im
/// wrappers from fold_conj_pairs() are preserved — recombination acts on
/// the wrapped inner Product.
///
/// @param expr a `Sum` of (optionally Re/Im-wrapped) NonSymm contractions
/// @return a `Sum` in minimal mixed antisymm/non-antisymm form, with the
///         same scalar value as @p expr
/// @note value-preserving by construction — it only re-groups terms that
///       sum to the same scalar; the result is the analogue of what the
///       open-shell spin tracer yields, but driven by index wiring
///       rather than spin-label uniformity
[[nodiscard]] ExprPtr antisymm_recombine(const ExprPtr& expr);

/// @brief Burnside-enumeration Kramers tracer.
///
/// Like `kramers_trace`, but instead of walking all `2^n` per-index
/// Kramers configurations it enumerates only the orbit representatives
/// under the index-permutation symmetry group induced by the
/// antisymmetrizer structure of the input (the group generated by
/// transpositions within each tensor's antisymmetric bra/ket groups).
/// Each representative is emitted once, scaled by its orbit size; group
/// elements act as bit-permutations on the n-bit config and orbits are
/// found by bit-ops.
///
/// A transposition `(p,q)` is admitted as a group generator iff in every
/// tensor containing both indices they sit in the same antisymmetric
/// bra/ket group and the accumulated permutation sign is `+1`.
///
/// @param expr a closed-shell all-spinor expression with `Kramers::any`
///             indices
/// @return a `Sum<Product>` value-equivalent (after `canonicalize`) to
///         the output of kramers_trace(); feed it the same
///         `fold_conj_pairs → antisymm_recombine` pipeline
[[nodiscard]] ExprPtr kramers_trace_burnside(const ExprPtr& expr);

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_SPINOR_HPP
