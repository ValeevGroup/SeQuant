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
/// The result is a flat sum of products with its scalars combined; it is not
/// canonicalized. With no variables, the result is @p E as it is.
[[nodiscard]] ExprPtr diff(const ExprPtr& E,
                           const container::svector<ExprPtr>& variables);

}  // namespace sequant::deriv

#endif  // SEQUANT_CORE_DERIV_DERIV_HPP
