//
// Kramers tracing for closed-shell relativistic (2-/4-component) theories.
//
// Rewrites an all-spinor (Kramers-free) expression — e.g. the MP2 energy
// `1/4 g-bar t-bar` over full spinor indices — into a sum over the
// time-reversal-symmetry (TRS) canonical Kramers-block representatives.
//
// Design (round 2): a single function `closed_shell_kramers_trace` that
// mirrors the naive spin tracer `spintrace_impl` but, instead of the
// Ms-conservation filter (`can_expand`, which is invalid relativistically),
// folds whole Kramers configurations related by global time reversal
// (T = flip every Kramers label + complex-conjugate) into a single
// `2 Re[...]` representative. Particle-interchange (sigma) folding is left
// to `sequant::canonicalize`; external-index antisymmetry folding (CC,
// rank-general) and the per-column internal-T-reach are deferred.
//
// Kramers labels reuse the `Spin` quantum number (space_qns.hpp):
//   Spin::alpha (== Spin::up)   = Kramers-up   (label suffix U+2191)
//   Spin::beta  (== Spin::down) = Kramers-down (label suffix U+2193)
//
// `RealPart`/`ImagPart` are symbolic scalar Expr markers: they wrap a
// closed (scalar-valued) contraction so the evaluator takes its real/
// imaginary part. They are opaque to `simplify`/`canonicalize`, so the
// wrapped sub-expression is never re-permuted. SeQuant's own evaluation
// engine does not (yet) evaluate `Re()`/`Im()` of a tensor network — that
// is a TODO; for now the marker is consumed by the MPQC evaluator.
//

#ifndef SEQUANT_DOMAIN_MBPT_SPINOR_HPP
#define SEQUANT_DOMAIN_MBPT_SPINOR_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <SeQuant/domain/mbpt/fwd.hpp>

#include <string>
#include <utility>

namespace sequant::mbpt {

// =====================================================================
// Re/Im evaluation markers + complex-conjugation label convention.
//
// Convention: complex conjugation of a tensor is encoded as a `*` suffix
// on the tensor label (e.g. `g` -> `g*`). Evaluator-side dispatch
// translates the suffix into a `.conj()` call on the numeric tensor.
// =====================================================================

/// @brief Tests whether a tensor label encodes complex conjugation.
inline bool has_conj_suffix(std::wstring_view label) {
  return !label.empty() && label.back() == L'*';
}

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
// Kramers tracer.
// =====================================================================

// clang-format off
/// @brief Traces an all-spinor (Kramers-free) closed-shell expression into a
///        sum over time-reversal-canonical Kramers-block representatives.
/// @details Mirrors the naive spin tracer `spintrace` in structure: it
/// enumerates the 2^n Kramers configurations of the term's indices and assigns
/// Kramers-up/down labels (reusing `Spin::alpha`/`Spin::beta`). Unlike
/// `spintrace`, it applies NO Ms-conservation filter (Kramers is not conserved
/// relativistically). Configurations related by global time reversal
/// (T: flip every Kramers label, complex-conjugate) are folded pairwise into a
/// single `RealPart`-wrapped representative (`A + A* = 2 Re A`). For the
/// integral `g` the bra antisymmetry is expanded (to the NonSymm form);
/// the amplitude `t` is kept antisymmetric (the "level-1" form). The caller is
/// expected to run `canonicalize` + `rapid_simplify` on the result to perform
/// the particle-interchange (sigma) merge.
/// @param expr an all-spinor expression (indices have no Kramers/Spin label)
/// @param ext_index_groups groups of external indices (empty for a fully
///        contracted energy; external-antisymmetry folding is not yet
///        implemented and these are treated like internal groups for now)
/// @return the Kramers-traced expression (unsimplified)
// clang-format on
ExprPtr closed_shell_kramers_trace(
    const ExprPtr& expr,
    const container::svector<container::svector<Index>>& ext_index_groups = {});

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_SPINOR_HPP
