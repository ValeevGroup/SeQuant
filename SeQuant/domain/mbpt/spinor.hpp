//
// Spinor decomposition passes for relativistic 2-component theories.
//
// Provides three symbolic transforms over closed-shell relativistic
// tensor expressions:
//
//   1. kramers_trace / kramers_trace_cycles / kramers_trace_burnside
//      Rewrite an all-spinor (Kramers::any) expression into a sum over
//      canonical Kramers-block representatives.
//   2. fold_conj_pairs + antisymm_recombine
//      Post-trace folds: collapse TRS-conjugate Product pairs to
//      `2·Re(A)`/`2i·Im(A)`, recombine direct/exchange pairs into
//      Kramers-restricted antisymmetric tensors.
//   3. complex_split
//      Split each complex Kramers tensor `g = g~r + i·g~i` so the
//      result is over real tensors only (`Re(g·t) → g~r·t~r −
//      g~i·t~i`).
//
// All passes operate on the canonical SeQuant Tensor representation
// and consume the `Kramers` index quantum number defined in
// SeQuant/domain/mbpt/space_qns.hpp.
//

#ifndef SEQUANT_DOMAIN_MBPT_SPINOR_HPP
#define SEQUANT_DOMAIN_MBPT_SPINOR_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/space_qns.hpp>

#include <optional>
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
inline Kramers to_kramers(const QuantumNumbersAttr& t) {
  SEQUANT_ASSERT((t.to_int32() & mask_v<Kramers>) != 0);
  return static_cast<Kramers>(t.to_int32() & mask_v<Kramers>);
}

/// @brief Strips the Kramers QN bits from a QuantumNumbersAttr.
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
template <typename WS, typename = std::enable_if_t<
                           meta::is_wstring_or_view_v<std::decay_t<WS>>>>
std::wstring kramers_annotation_replace(WS&& label, Kramers k) {
  auto label_kf = kramers_annotation_remove(std::forward<WS>(label));
  return kramers_annotation_add(label_kf, k);
}

/// @brief Constructs a Kramers-up Index from @p idx (label decorated with ⇑).
[[nodiscard]] Index make_kramers_up(const Index& idx);

/// @brief Constructs a Kramers-down Index from @p idx (label decorated with ⇓).
[[nodiscard]] Index make_kramers_dn(const Index& idx);

/// @brief Constructs a Kramers-free (`Kramers::any`) Index from @p idx.
[[nodiscard]] Index make_kramers_free(const Index& idx);

/// @brief Applies an Index→Index replacement map to all tensors in @p expr.
///
/// @sa append_spin in spin.hpp — append_kramers is a synonym; both are
///     thin wrappers around `Tensor::transform_indices`.
ExprPtr append_kramers(const ExprPtr& expr,
                       const container::map<Index, Index>& index_replacements);

// =====================================================================
// Re/Im wrappers + complex-conjugation label convention.
//
// Convention: complex conjugation of a tensor is encoded as a `*`
// suffix on the tensor's label (e.g., `g` → `g*`). Evaluator-side
// dispatch translates the suffix into a `.conj()` call on the
// underlying numeric tensor.
//
// `RealPart`/`ImagPart` are scalar Expr nodes that wrap a contracted
// (closed-form) Product; they are opaque to `simplify`/`canonicalize`,
// so the wrapped Antisymm sub-expressions are not re-permuted.
// =====================================================================

/// @brief Tests whether a tensor label encodes complex conjugation.
inline bool has_conj_suffix(std::wstring_view label) {
  return !label.empty() && label.back() == L'*';
}

/// @brief Toggles the `*` conjugation suffix on a tensor label.
[[nodiscard]] std::wstring toggle_conj_suffix(std::wstring_view label);

/// @brief Symbolic real-part wrapper: `RealPart(E)` = `Re(E)` for a
///        scalar-valued Expr `E`.
class RealPart : public sequant::Expr {
 public:
  RealPart() = delete;
  RealPart(const RealPart&) = default;
  RealPart(RealPart&&) = default;
  ~RealPart() override = default;

  /// @pre @p inner is non-null and scalar-valued
  ///
  /// We do not enforce `Expr::is_scalar()` because Tensor (an atom)
  /// returns false unconditionally and a closed-contraction Product of
  /// two Tensors inherits the same answer — the typical inner here.
  explicit RealPart(ExprPtr inner) : inner_{std::move(inner)} {
    SEQUANT_ASSERT(inner_);
  }

  const ExprPtr& inner() const { return inner_; }
  bool is_scalar() const override { return true; }
  type_id_type type_id() const override { return get_type_id<RealPart>(); }
  ExprPtr clone() const override { return ex<RealPart>(inner_->clone()); }
  void adjoint() override {}  // Re(E) is real, self-adjoint
  std::wstring to_latex() const override {
    return L"\\Re\\left[" + inner_->to_latex() + L"\\right]";
  }

 private:
  ExprPtr inner_;

  hash_type memoizing_hash() const override {
    auto compute = [this]() {
      auto v = hash::value(*inner_);
      hash::combine(v, std::size_t{0xC0FFEE01ull});
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

/// @brief Symbolic imaginary-part wrapper: `ImagPart(E)` = `Im(E)`.
class ImagPart : public sequant::Expr {
 public:
  ImagPart() = delete;
  ImagPart(const ImagPart&) = default;
  ImagPart(ImagPart&&) = default;
  ~ImagPart() override = default;

  /// @pre @p inner is non-null and scalar-valued (see RealPart for why
  ///      we don't use `Expr::is_scalar()` as the precondition)
  explicit ImagPart(ExprPtr inner) : inner_{std::move(inner)} {
    SEQUANT_ASSERT(inner_);
  }

  const ExprPtr& inner() const { return inner_; }
  bool is_scalar() const override { return true; }
  type_id_type type_id() const override { return get_type_id<ImagPart>(); }
  ExprPtr clone() const override { return ex<ImagPart>(inner_->clone()); }
  void adjoint() override {}  // Im(E) is real, self-adjoint
  std::wstring to_latex() const override {
    return L"\\Im\\left[" + inner_->to_latex() + L"\\right]";
  }

 private:
  ExprPtr inner_;

  hash_type memoizing_hash() const override {
    auto compute = [this]() {
      auto v = hash::value(*inner_);
      hash::combine(v, std::size_t{0xC0FFEE02ull});
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

/// @brief Wraps @p expr as `Re(expr)`.
[[nodiscard]] inline ExprPtr real_part(ExprPtr expr) {
  return ex<RealPart>(std::move(expr));
}

/// @brief Wraps @p expr as `Im(expr)`.
[[nodiscard]] inline ExprPtr imaginary_part(ExprPtr expr) {
  return ex<ImagPart>(std::move(expr));
}

// =====================================================================
// Kramers tracers.
//
// Three coexisting implementations, all returning expressions whose
// scalar value equals the full 2^n per-Kramers-index trace of @p expr.
// Provided in parallel for cross-validation; pick whichever is fastest
// in production.
//
// All assume @p expr is fully contracted (every barred index appears
// an even number of times). CC residuals carry external bars and need
// `(-1)^(external bars)` sign bookkeeping that is not yet implemented.
// =====================================================================

/// @brief 2^n exhaustive enumeration over Kramers configurations of
///        every Kramers::any index in @p expr.
ExprPtr kramers_trace(const ExprPtr& expr);

/// @brief ResultExpr overload: kramers-trace a full equation
///        `D^{...}_{...} = RHS`.
///
/// Enumerates 2^N external Kramers configurations; per configuration,
/// substitutes the externals on the RHS, then runs `kramers_trace` +
/// `antisymm_recombine` on the result. Two folding levels are
/// supported via @p fold_klein:
///   - **true** (default): emit only the canonical representative per
///     orbit under the Klein 4-group {e, α (bra-flip), β (ket-flip),
///     αβ (full TRS)}. For an N-external Kramers-restricted result
///     this reduces 2^N raw blocks by a factor of 4 (typically).
///   - **false**: emit only canonical representatives under full TRS
///     alone (cfg ≤ ~cfg). For partial-Kramers-flip equivalence
///     verification — emits the 2^(N-1) TRS-canonical configurations
///     so a caller can pair them into Klein orbits and numerically
///     check the partial-flip equivalence.
///
/// @return one ResultExpr per non-vanishing canonical configuration.
[[nodiscard]] container::svector<ResultExpr> kramers_trace(
    const ResultExpr& expr, bool fold_klein = true);

/// @brief Cycle-driven multiset enumeration: decompose the
///        contraction graph into cycles, group cycles by shape, and
///        enumerate Kramers-labeling multisets per group with
///        multinomial multiplicity.
[[nodiscard]] ExprPtr kramers_trace_cycles(const ExprPtr& expr);

/// @brief Burnside-orbit enumeration over the index-permutation
///        symmetry group induced by the antisymmetrizer structure.
[[nodiscard]] ExprPtr kramers_trace_burnside(const ExprPtr& expr);

/// @brief Flips every specialized Kramers state (`up` ↔ `down`) in @p expr.
///
/// `Kramers::any` indices are left untouched.
[[nodiscard]] ExprPtr flip_kramers(const ExprPtr& expr);

/// @brief Folds TRS conjugate-pair Products in a Sum.
///
/// For each summand `A`, builds `A_flipped = flip_kramers(A)`,
/// canonicalizes it, and looks for a structural match among the
/// remaining summands:
///   - `A == A_flipped` (self-conjugate): emits `Re(A)`;
///   - matching `B` with same scalar: emits `2·Re(A)`, drops B;
///   - matching `B` with opposite scalar: emits `2i·Im(A)`, drops B;
///   - no match: keeps `A` unchanged.
///
/// O(n) hash-bucketed; equivalent to pairwise search but much cheaper.
[[nodiscard]] ExprPtr fold_conj_pairs(const ExprPtr& expr);

/// @brief Recombines direct/exchange Product pairs into antisymmetric
///        tensors — the inverse of `expand_antisymm`.
///
/// Iterates to a fixed point. Recombines only when one tensor matches
/// structurally between two summands, the other differs by a single
/// column transposition, and the summand scalars are negatives. Re/Im
/// wrappers from `fold_conj_pairs` are preserved.
///
/// @todo CC support requires generalizing `detect_single_column_swap`
///       past rank-(2,2) for CCSDT triples.
[[nodiscard]] ExprPtr antisymm_recombine(const ExprPtr& expr);

// =====================================================================
// Complex (Re/Im) split.
//
// For closed-shell 1eX2C the Coulomb operator spin-traces each
// per-electron transition density to a *complex* scalar; no quaternion
// (4-component) structure survives. The minimal real-arithmetic form
// is the 2-component split: `g → g~r + i·g~i`, both real, propagated
// symbolically (`Re(g·t) = g~r·t~r − g~i·t~i`). Components are encoded
// as a `~r`/`~i` label suffix.
// =====================================================================

/// @brief Encodes a real/imaginary-part marker onto a tensor label.
[[nodiscard]] inline std::wstring complex_label_add(std::wstring_view core,
                                                    bool imag) {
  return std::wstring{core} + (imag ? L"~i" : L"~r");
}

/// @brief Parses the real/imaginary-part marker off a tensor label.
/// @return true for `~i`, false for `~r`, std::nullopt if no marker
///         (a trailing `*` conjugation marker is tolerated).
[[nodiscard]] inline std::optional<bool> complex_label_parse(
    std::wstring_view label) {
  const auto pos = label.find(L'~');
  if (pos == std::wstring_view::npos) return std::nullopt;
  auto tail = label.substr(pos + 1);
  if (!tail.empty() && tail.back() == L'*')
    tail = tail.substr(0, tail.size() - 1);
  if (tail == std::wstring_view{L"r"}) return false;
  if (tail == std::wstring_view{L"i"}) return true;
  return std::nullopt;
}

/// @brief Returns the marker-free core of a (possibly complex-split)
///        tensor label.
[[nodiscard]] inline std::wstring complex_label_core(std::wstring_view label) {
  return std::wstring{label.substr(0, label.find(L'~'))};
}

/// @brief Splits each complex tensor `g = g~r + i·g~i` and propagates
///        complex arithmetic symbolically.
///
/// Conjugated tensors (`g*`) split as `g~r − i·g~i`. Real
/// sub-expressions are distributed and constant-folded but
/// deliberately NOT canonicalized — the split tensors inherit the
/// input's index wiring, and canonicalizing `Antisymm` `~r`/`~i`
/// tensors would emit signs the evaluator cannot see.
///
/// @pre every summand of @p expr is wrapped in `RealPart`/`ImagPart`
///      (i.e. run `fold_conj_pairs` first)
[[nodiscard]] ExprPtr complex_split(const ExprPtr& expr);

// =====================================================================
// V1-based open-shell spin tracer.
//
// SeQuant's default `mbpt::open_shell_spintrace` runs the
// canonicalization pipeline through `TensorNetworkV3`, which assumes
// a CC-style index wiring (every contracted dummy appears once in a
// bra and once in a ket). Expressions in which a dummy appears across
// two bras (or two kets) — e.g. the PNO pair density
// `D^{ij}_{ab} = t̄^{ac}_{ij} t̄^{bc}_{ij}`, where `c` is summed across
// two bras — trip V3's graph-build assertion on a worker thread,
// bypassing user-side `try`/`catch`.
//
// `open_shell_spintrace_v1` is a drop-in V1-canonicalization variant
// of `open_shell_spintrace_impl` (see `spin.cpp:1367`). The algorithm
// is identical; only the canonicalization back-end is swapped, which
// V1's permissive graph builder accepts.
//
// `canonicalize_v1` and `simplify_v1` are the underlying primitives.
// =====================================================================

/// @brief Like `sequant::canonicalize`, but forces every all-tensor
///        Product subterm through `TensorNetworkV1`.
///
/// Use when an expression's contraction wiring trips
/// `TensorNetworkV3`'s graph-build assumptions (typically: a dummy
/// summed across two bras or two kets).
ExprPtr& canonicalize_v1(ExprPtr& expr,
                         CanonicalizeOptions opts = CanonicalizeOptions{
                             .method = CanonicalizationMethod::Lexicographic});

/// @brief Like `sequant::simplify`, but uses `canonicalize_v1`.
ExprPtr& simplify_v1(ExprPtr& expr);

/// @brief V1-canonicalization variant of `mbpt::open_shell_spintrace`.
///
/// Same algorithm as `open_shell_spintrace_impl` in `spin.cpp:1367`;
/// only difference is that `simplify`/`canonicalize` are routed
/// through `simplify_v1`/`canonicalize_v1` so the pipeline tolerates
/// non-CC index wiring (open dummies across two bras / two kets).
///
/// @param expr the spinorbital expression to trace
/// @param ext_index_groups groups of external (free) indices; one
///        per particle (e.g. `{{i,a}, {j,b}}` for a rank-(2,2) RHS)
/// @param target_spin_case if set, returns only the αα...β...β
///        configuration with @p target_spin_case beta-tagged
///        externals; otherwise returns all such configurations
[[nodiscard]] std::vector<ExprPtr> open_shell_spintrace_v1(
    const ExprPtr& expr,
    container::svector<container::svector<Index>> ext_index_groups,
    std::optional<int> target_spin_case = std::nullopt);

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_SPINOR_HPP
