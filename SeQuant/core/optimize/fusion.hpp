#ifndef SEQUANT_OPT_FUSION_HPP
#define SEQUANT_OPT_FUSION_HPP

#include <SeQuant/core/expr.hpp>

namespace sequant::opt {

///
/// Use this class to test the fusibility of two products.
///
/// The pattern of factors should match in both products for them
/// to have a common factor.
/// Fusion is either possible from the left hand side or the right hand side.
///
/// eg. abcd + abef => ab(cd + ef) from left. nullptr from right.
///     abef + cdef => nullptr from left. (ab + cd)ef from right.
///
/// Only common scalars are factored out as of now.
/// eg. (1/2)abcd + (1/2)abef => (1/2)ab(cd + ef)
///     (1/2)abcd + (1/4)abef => ab((1/2)cd + (1/4)ef)
///     (1/2)abcd - (1/2)abef => ab((1/2)cd - (1/2)ef)

class Fusion {
 public:
  Fusion(Product const& lhs, Product const& rhs);

  /// the result of fusion from left hand side.
  /// returns nullptr if no fusion possible.
  ExprPtr left() const;

  /// the result of fusion from right hand side.
  /// returns nullptr if no fusion possible.
  ExprPtr right() const;

  static ExprPtr fuse_left(Product const& lhs, Product const& rhs);

  static ExprPtr fuse_right(Product const& lhs, Product const& rhs);

  ///
  /// Get the greatest common divisor of two rational numbers.
  ///
  static rational gcd_rational(rational const& left, rational const& right);

  ///
  /// Fuse scalars @param left and @param right and the return the result
  /// as an array of three elements: first is the greatest common factor,
  /// second the fused sub-factor of @param left and the third is that
  /// of @param right.
  ///
  static std::array<rational, 3> fuse_scalar(rational const& left,
                                             rational const& right);

 private:
  ExprPtr left_;

  ExprPtr right_;
};

}  // namespace sequant::opt

#endif  // SEQUANT_OPT_FUSION_HPP
