#ifndef SEQUANT_CORE_DERIV_DERIV_HPP
#define SEQUANT_CORE_DERIV_DERIV_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>

namespace sequant::deriv {

/// @return the derivative of the explicit function @p E by each of
/// @p variables in turn (any number, so any order): d^n E / dx_1 ... dx_n,
/// taken left to right, so `diff(E, {x, y})` is `diff(diff(E, {x}), {y})`
///
/// A scalar `Variable` x is complex: x and its adjoint x^* are independent.
///
/// For a tensor leaf x, let G be the positional gradient of @p E by x: the
/// leaves of @p E with x's label and ranks and, where x^+ holds the same
/// elements as x but the exchange of bra and ket is no operation of x (see
/// below), those with the label and ranks of x^+, at the exchanged position,
/// are read, as they are written, as elements of one unconstrained array,
/// and G(J) is the derivative by the element at
/// position J. The operations of x are the permutations its declared symmetry
/// allows: for Symmetry::Antisymm the permutations of the bra slots times
/// those of the ket slots, each with its sign; for Symmetry::Nonsymm with
/// ColumnSymmetry::Symm the permutations of the columns; each alone and,
/// where that is an operation (see below), followed by the exchange of bra
/// and ket; and, a protoindex list being unordered, each of these with every
/// order of each distinct protoindex list of x whose labels are not slots of
/// x. At a position J of x's block (the list of the bases of its slots)
/// the result is the sum of G over all operations g of x:
/// `diff(E, x)(J) = sum_g sign(g) G(g J)`. This is the derivative with
/// respect to the independent elements of x: at J it is the derivative along
/// J and its symmetry partners, times the order of J's stabilizer.
/// The sum is taken over the spellings of x, its indices in each slot order
/// an operation makes, as the operation's sign times the derivative by that
/// spelling; the operations that permute like slots are summed by one Â or Ŝ
/// factor.
///
/// If x^+ holds other numbers than x (no bra-ket symmetry, or
/// BraKetSymmetry::Conjugate, over a complex field), the derivative is
/// Wirtinger's: leaves spelled x^+ are constants. If x^+ holds the same
/// elements as x (BraKetSymmetry::Symm, or a real base field), the exchange of
/// bra and ket is an operation of x where x has bra-ket symmetry and bra and
/// ket have one nonzero rank; otherwise the result is the derivative by x plus
/// the derivative by x^+.
///
/// The free indices of the result are those of @p E plus those of each
/// variable, under the names the variable is spelled with. The result is a
/// flat sum of products with its scalars combined; it is not canonicalized.
/// With no variables, the result is @p E as it is.
/// Symmetrizers are never expanded: those of @p E stay factors, and that of
/// x is put on the result as a factor.
/// @throw Exception for an input that is rejected
[[nodiscard]] ExprPtr diff(const ExprPtr& E,
                           const container::svector<ExprPtr>& variables);

/// @return @p E where every one of @p variables is zero: the terms that hold
/// one of them, in one of its spellings or those of its adjoint, are dropped.
/// A leaf with a tensor variable's label vanishes if it lies inside that one
/// variable: a leaf that meets a spelling vanishes if each of its slots lies
/// within the index space of that spelling's slot
/// @throw Exception if a variable is not a Variable or a tensor leaf (such as a
/// sum, a product, a symmetrizer or a Kronecker delta), repeats an index, has
/// an index whose protoindices mix its slots with other labels, or declares
/// Symmetry::Symm; if a variable stands under a non-positive power; or if a
/// leaf holds elements of a variable and lies inside none of the variables
/// (also where it lies inside their union), or declares other symmetries than
/// the variable whose elements it may hold, or its protoindices do not follow
/// its slots as the variable's do
[[nodiscard]] ExprPtr at_zero(const ExprPtr& E,
                              const container::svector<ExprPtr>& variables);

}  // namespace sequant::deriv

#endif  // SEQUANT_CORE_DERIV_DERIV_HPP
