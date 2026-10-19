//
// Created by Kshitij Surjuse on 2026-08-31.
//

#ifndef SEQUANT_CORE_EXPRESSIONS_COMPLEX_HPP
#define SEQUANT_CORE_EXPRESSIONS_COMPLEX_HPP

#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/product.hpp>

#include <utility>

/// @file complex.hpp
/// Symbolic real/imaginary-part wrappers for scalar-valued expressions.
///
/// `RealPart(E)` = `Re(E)` and `ImagPart(E)` = `Im(E)` for a scalar-valued
/// Expr `E`; both wrappers are _real_-valued by convention (`E = Re(E) +
/// i*Im(E)`), hence self-adjoint and conjugation-invariant. They are general
/// expression nodes (not domain-specific): any pipeline that folds a sum of
/// conjugate pairs (`A + A* = 2 Re(A)`, `A - A* = 2i Im(A)`) produces them.
///
/// Construction goes through the smart builders real_part()/imaginary_part(),
/// which apply the eager composition rules
///   `Re(Re x) = Re x`, `Re(Im x) = Im x`, `Im(Re x) = 0`, `Im(Im x) = 0`
/// (both wrappers are real-valued, so a second projection is the identity on
/// `Re`-of and annihilates on `Im`-of) and evaluate complex `Constant`s via
/// their `Complex` ring value.
///
/// @note the evaluation engine ingests both nodes: `EvalOp::RealPart` and
///       `EvalOp::ImagPart` evaluate the wrapped expression and take the
///       real or the imaginary part of the resulting scalar, so an inner
///       that is not scalar-valued has no evaluation.

namespace sequant {

namespace detail {

/// canonicalizes @p inner in place (rapidly if @p rapid), folds a Constant
/// byproduct back into it and hoists a real scalar out of it
/// @return the hoisted real scalar as a Constant, or nullptr
[[nodiscard]] ExprPtr canonicalize_projected(ExprPtr& inner,
                                             CanonicalizeOptions opts,
                                             bool rapid);

/// @brief the common implementation of RealPart and ImagPart: a real-valued
///        projection of a scalar-valued Expr, hence self-adjoint and
///        conjugation-invariant. `Derived` supplies the LaTeX name and the
///        hash salt that tell the two apart.
template <typename Derived>
class ProjectionExpr : public Expr {
 public:
  ProjectionExpr() = delete;
  ProjectionExpr(const ProjectionExpr&) = default;
  ProjectionExpr(ProjectionExpr&&) = default;
  ~ProjectionExpr() override = default;

  /// @pre @p inner is non-null and scalar-valued
  ///
  /// We do not enforce `Expr::is_scalar()` because Tensor (an atom)
  /// returns false unconditionally and a closed-contraction Product of
  /// two Tensors inherits the same answer -- the typical inner here.
  explicit ProjectionExpr(ExprPtr inner) : inner_{std::move(inner)} {
    SEQUANT_ASSERT(inner_);
  }

  const ExprPtr& inner() const { return inner_; }
  bool is_scalar() const override { return true; }
  type_id_type type_id() const override { return get_type_id<Derived>(); }
  ExprPtr clone() const override { return ex<Derived>(inner_->clone()); }
  std::int8_t adjoint() override { return 1; }     // real, self-adjoint
  std::int8_t kconjugate() override { return 1; }  // real
  std::wstring to_latex() const override {
    return std::wstring(Derived::latex_name) + L"\\left[" + inner_->to_latex() +
           L"\\right]";
  }

  /// canonicalizes the wrapped expression in place, then re-applies the eager
  /// hoist, so that no canonical wrapper wraps a Product with a real scalar
  /// @return the hoisted real scalar as a Constant, or nullptr when there is
  ///         none. A byproduct of the inner canonicalization is folded back
  ///         into the wrapped expression first: the wrapper is not linear
  ///         over a complex scalar, so only a real one can leave it
  ExprPtr canonicalize(CanonicalizeOptions opts =
                           CanonicalizeOptions::default_options()) override {
    auto hoisted = canonicalize_projected(inner_, opts, /*rapid=*/false);
    reset_hash_value();
    return hoisted;
  }

  /// @copydoc canonicalize()
  ExprPtr rapid_canonicalize(
      CanonicalizeOptions opts =
          CanonicalizeOptions::default_options().copy_and_set(
              CanonicalizationMethod::Rapid)) override {
    SEQUANT_ASSERT(opts.method == CanonicalizationMethod::Rapid);
    auto hoisted = canonicalize_projected(inner_, opts, /*rapid=*/true);
    reset_hash_value();
    return hoisted;
  }

  /// the wrapped expression is this node's only subexpression, so that
  /// Expr::visit(), index transforms and relabeling reach it
  ExprIterator begin_subexpr() override {
    // N.B. a mutable iterator into inner_ invalidates the memoized hash
    reset_hash_value();
    return ExprIterator{&inner_};
  }
  ExprIterator end_subexpr() override {
    reset_hash_value();
    return ExprIterator{&inner_ + 1};
  }
  ConstExprIterator begin_subexpr() const override {
    return ConstExprIterator{&inner_};
  }
  ConstExprIterator end_subexpr() const override {
    return ConstExprIterator{&inner_ + 1};
  }

 private:
  ExprPtr inner_;

  hash_type memoizing_hash() const override {
    if (!hash_value_) {
      auto v = hash::value(*inner_);
      hash::combine(v, Derived::hash_salt);
      hash_value_ = v;
    }
    return *hash_value_;
  }
  bool static_equal(const Expr& that) const override {
    return *inner_ == *static_cast<const ProjectionExpr&>(that).inner_;
  }
  bool static_less_than(const Expr& that) const override {
    return *inner_ < *static_cast<const ProjectionExpr&>(that).inner_;
  }
};

}  // namespace detail

/// @brief Symbolic real-part wrapper: `RealPart(E)` = `Re(E)` for a
///        scalar-valued Expr `E`.
class RealPart : public detail::ProjectionExpr<RealPart> {
 public:
  using ProjectionExpr::ProjectionExpr;
  static constexpr const wchar_t* latex_name = L"\\Re";
  static constexpr std::size_t hash_salt = 0xC0FFEE01ull;
};

/// @brief Symbolic imaginary-part wrapper: `ImagPart(E)` = `Im(E)`.
class ImagPart : public detail::ProjectionExpr<ImagPart> {
 public:
  using ProjectionExpr::ProjectionExpr;
  static constexpr const wchar_t* latex_name = L"\\Im";
  static constexpr std::size_t hash_salt = 0xC0FFEE02ull;
};

namespace detail {
/// strips a Product's scalar: returns the same factors with scalar 1
[[nodiscard]] ExprPtr strip_scalar(const Product& prod);

/// Takes a real scalar out of @p inner in place, as `Re(c X) = c Re(X)` and
/// `Im(c X) = c Im(X)` allow for a real @c c. A complex scalar stays put:
/// neither wrapper is linear over it.
/// @return the scalar taken out, or nullptr if there was none to take
[[nodiscard]] ExprPtr hoist_real_scalar(ExprPtr& inner);
}  // namespace detail

/// @brief Wraps @p expr as `Re(expr)`, applying the eager rules.
///
/// `Re(Constant c)` evaluates to `Constant(c.real())`; `Re(Re x) = Re x` and
/// `Re(Im x) = Im x` (both wrappers are real-valued). Scalar prefactors:
/// a real scalar hoists (`Re(c X) = c Re(X)`), a purely imaginary scalar
/// rotates (`Re(i b X) = -b Im(X)`); a general complex scalar stays wrapped
/// (`Re(c X) = Re(c) Re(X) - Im(c) Im(X)` is recognized, not auto-expanded).
[[nodiscard]] ExprPtr real_part(ExprPtr expr);

/// @brief Wraps @p expr as `Im(expr)`, applying the eager rules.
///
/// `Im(Constant c)` evaluates to `Constant(c.imag())`; `Im(Re x) = 0` and
/// `Im(Im x) = 0` (both wrappers are real-valued). Scalar prefactors:
/// a real scalar hoists (`Im(c X) = c Im(X)`), a purely imaginary scalar
/// rotates (`Im(i b X) = b Re(X)`); a general complex scalar stays wrapped.
[[nodiscard]] ExprPtr imaginary_part(ExprPtr expr);

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPRESSIONS_COMPLEX_HPP
