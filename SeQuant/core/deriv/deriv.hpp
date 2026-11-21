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
/// leaves of @p E with x's label and ranks
/// are read, as they are written, as elements of one unconstrained array,
/// and G(J) is the derivative by the element at
/// position J. The operations of x are the permutations its declared symmetry
/// allows: for Symmetry::Antisymm the permutations of the bra slots times
/// those of the ket slots, each with its sign; for Symmetry::Nonsymm with
/// ColumnSymmetry::Symm the permutations of the columns. At a position J of
/// x's block (the list of the bases of its slots) the result is the sum of G
/// over the operations g of x that permute like slots (slots of one basis and
/// protoindex list): `diff(E, x)(J) = sum_g sign(g) G(g J)`, summed by one Â
/// or Ŝ factor.
///
/// The free indices of the result are those of @p E plus those of each
/// variable, under the names the variable is spelled with. The result is a
/// flat sum of products with its scalars combined; it is not canonicalized.
/// With no variables, the result is @p E as it is.
/// Symmetrizers are never expanded: those of @p E stay factors, and that of
/// x is put on the result as a factor.
[[nodiscard]] ExprPtr diff(const ExprPtr& E,
                           const container::svector<ExprPtr>& variables);

}  // namespace sequant::deriv

#endif  // SEQUANT_CORE_DERIV_DERIV_HPP
