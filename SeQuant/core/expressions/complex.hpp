//
// Created by Kshitij Surjuse on 2026-08-31.
//

#ifndef SEQUANT_CORE_EXPRESSIONS_COMPLEX_HPP
#define SEQUANT_CORE_EXPRESSIONS_COMPLEX_HPP

#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/product.hpp>

#include <utility>

/// @file complex.hpp
/// Symbolic real/imaginary-part wrappers for scalar-valued expressions.
///
/// `RealPart(E)` = `Re(E)` and `ImagPart(E)` = `Im(E)` for a scalar-valued
/// Expr `E`; both wrappers are REAL-valued by convention (`E = Re(E) +
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
/// @note evaluation of `Re`/`Im` of a tensor network is not implemented yet;
///       consumers evaluate the inner expression and take the real/imaginary
///       part of the resulting scalar.

namespace sequant {

/// @brief Symbolic real-part wrapper: `RealPart(E)` = `Re(E)` for a
///        scalar-valued Expr `E`.
class RealPart : public Expr {
 public:
  RealPart() = delete;
  RealPart(const RealPart&) = default;
  RealPart(RealPart&&) = default;
  ~RealPart() override = default;

  /// @pre @p inner is non-null and scalar-valued
  ///
  /// We do not enforce `Expr::is_scalar()` because Tensor (an atom)
  /// returns false unconditionally and a closed-contraction Product of
  /// two Tensors inherits the same answer -- the typical inner here.
  explicit RealPart(ExprPtr inner) : inner_{std::move(inner)} {
    SEQUANT_ASSERT(inner_);
  }

  const ExprPtr& inner() const { return inner_; }
  bool is_scalar() const override { return true; }
  type_id_type type_id() const override { return get_type_id<RealPart>(); }
  ExprPtr clone() const override;
  void adjoint() override {}  // Re(E) is real, self-adjoint
  std::wstring to_latex() const override;

 private:
  ExprPtr inner_;

  hash_type memoizing_hash() const override;
  bool static_equal(const Expr& that) const override;
  bool static_less_than(const Expr& that) const override;
};

/// @brief Symbolic imaginary-part wrapper: `ImagPart(E)` = `Im(E)`.
class ImagPart : public Expr {
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
  ExprPtr clone() const override;
  void adjoint() override {}  // Im(E) is real, self-adjoint
  std::wstring to_latex() const override;

 private:
  ExprPtr inner_;

  hash_type memoizing_hash() const override;
  bool static_equal(const Expr& that) const override;
  bool static_less_than(const Expr& that) const override;
};

namespace detail {
/// strips a Product's scalar: returns the same factors with scalar 1
[[nodiscard]] ExprPtr strip_scalar(const Product& prod);
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
