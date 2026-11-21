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
/// position J. At a position J of x's block (the list of the bases of its
/// slots) the result is G(J).
///
/// The free indices of the result are those of @p E plus those of each
/// variable, under the names the variable is spelled with. The result is a
/// flat sum of products with its scalars combined; it is not canonicalized.
/// With no variables, the result is @p E as it is.
[[nodiscard]] ExprPtr diff(const ExprPtr& E,
                           const container::svector<ExprPtr>& variables);

}  // namespace sequant::deriv

#endif  // SEQUANT_CORE_DERIV_DERIV_HPP
