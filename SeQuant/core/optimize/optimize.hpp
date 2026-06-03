#ifndef SEQUANT_OPTIMIZE_OPTIMIZE_HPP
#define SEQUANT_OPTIMIZE_OPTIMIZE_HPP

#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/optimize/options.hpp>

namespace sequant {

/// Optimize the expression for lower evaluation cost.
///
/// \param expr  Expression to be optimized.
/// \param opts  Optimization parameters; see \c OptimizeOptions. By default:
///              the cost metric is flop count, index extents are taken from
///              \c IndexSpace::approximate_size(), and the summands of a sum
///              are reordered to cluster terms that share intermediates.
/// \return Optimized expression.
ExprPtr optimize(ExprPtr const& expr, OptimizeOptions opts = {});

/// \copydoc optimize(ExprPtr const&, OptimizeOptions)
ResultExpr& optimize(ResultExpr& expr, OptimizeOptions opts = {});

/// \copydoc optimize(ExprPtr const&, OptimizeOptions)
[[nodiscard]] ResultExpr& optimize(ResultExpr&& expr,
                                   OptimizeOptions opts = {});

// Overloads for backwards compatibility

/// \deprecated Use the \c OptimizeOptions overload instead.
///
/// Equivalent to calling the primary overload with default \c OptimizeOptions
/// and \c OptimizeOptions::reorder_sum set to \p reorder_sum.
///
/// \param expr         Expression to be optimized.
/// \param reorder_sum  If true, reorder the summands of a sum to cluster terms
///                     that share intermediates.
/// \return Optimized expression.
[[deprecated(
    "use the OptimizeOptions"
    " overload of optimize() instead")]] ExprPtr
optimize(ExprPtr const& expr, bool reorder_sum);

/// \copydoc optimize(ExprPtr const&, bool)
[[deprecated(
    "use the OptimizeOptions"
    " overload of optimize() instead")]] ResultExpr&
optimize(ResultExpr& expr, bool reorder_sum);

/// \copydoc optimize(ExprPtr const&, bool)
[[nodiscard, deprecated("use the OptimizeOptions"
                        " overload of optimize() instead")]] ResultExpr&
optimize(ResultExpr&& expr, bool reorder_sum);

}  // namespace sequant

#endif  // SEQUANT_OPTIMIZE_OPTIMIZE_HPP
