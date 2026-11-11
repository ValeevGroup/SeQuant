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

///
/// Factors out a factor shared by the top-level contractions of summands:
///
///   c_1 F X_1 + c_2 F' X_2 + ...  ->  F (c_1 X_1 + c_2 X_2' + ...)
///
/// Each summand must be a single-term optimized product (as produced by
/// \ref single_term_opt), i.e. a Tensor, or a Product of two factors that are
/// Tensors or Products of Tensors; any other summand is kept as is.
/// \c F' must equal \c F up to a renaming of its dummy (non-external) indices,
/// and \c X_1 and \c X_2 must have the same external indices once that renaming
/// is applied to \c X_2 (which gives \c X_2'). Dummy indices internal to
/// \c X_2 are renamed when they would clash. Each fused group saves all but one
/// of its top-level contractions. The summed \c X's are fused recursively.
///
/// Each summand is fused through the side (left or right factor) that matches
/// the most summands; ties go to the left factor.
///
/// \pre no summand has an index that appears in both factors of its top
///      contraction without being summed over (a Hadamard-style index), and
///      the summands do not rely on broadcasting: a renamed dummy may take the
///      label of a result index that a group's summands do not carry
/// \param summands the summands of a Sum
/// \return the summands with each group of two or more fused summands replaced
///         by one summand, placed at the position of the group's first member;
///         a fused summand is a Product of the shared factor and a Sum
///
Sum::summands_type fuse_top_factors(Sum::summands_type const& summands);

}  // namespace sequant::opt

#endif  // SEQUANT_OPT_FUSION_HPP
