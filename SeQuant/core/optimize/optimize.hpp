#ifndef SEQUANT_OPTIMIZE_OPTIMIZE_HPP
#define SEQUANT_OPTIMIZE_OPTIMIZE_HPP

#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/optimize/options.hpp>

namespace sequant {

/// Optimize the expression for lower evaluation cost.
///
/// \param expr  Expression to be optimized.
/// \param opts  Options controlling the cost metric, summand reordering, and
///              the index → extent provider. An empty
///              \c OptimizeOptions::idx_to_extent falls back to
///              \c IndexSpace::approximate_size().
/// \return Optimized expression.
ExprPtr optimize(ExprPtr const& expr, OptimizeOptions opts = {});

ResultExpr& optimize(ResultExpr& expr, OptimizeOptions opts = {});

[[nodiscard]] ResultExpr& optimize(ResultExpr&& expr,
                                   OptimizeOptions opts = {});

}  // namespace sequant

#endif  // SEQUANT_OPTIMIZE_OPTIMIZE_HPP
