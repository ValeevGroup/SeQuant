#ifndef SEQUANT_OPTIMIZE_OPTIMIZE_HPP
#define SEQUANT_OPTIMIZE_OPTIMIZE_HPP

#include <SeQuant/core/expr_fwd.hpp>

namespace sequant {

///
/// Optimize the expression using IndexSpace::approximate_size() for reference
/// index extent.
///
/// \param expr  Expression to be optimized.
/// \param reorder_sum If true, the summands are reordered so that terms with
///                    common sub-expressions appear closer to each other.
///                    True by default.
/// \return Optimized expression for lower evaluation cost.
ExprPtr optimize(ExprPtr const& expr, bool reorder_sum = true);

/// Optimize the expression using IndexSpace::approximate_size() for reference
/// index extent.
///
/// \param expr  Expression to be optimized.
/// \param reorder_sum If true, the summands are reordered so that terms with
///                    common sub-expressions appear closer to each other.
///                    True by default.
/// \return Optimized expression for lower evaluation cost.
ResultExpr& optimize(ResultExpr& expr, bool reorder_sum = true);

/// Optimize the expression using IndexSpace::approximate_size() for reference
/// index extent.
///
/// \param expr  Expression to be optimized.
/// \param reorder_sum If true, the summands are reordered so that terms with
///                    common sub-expressions appear closer to each other.
///                    True by default.
/// \return Optimized expression for lower evaluation cost.
[[nodiscard]] ResultExpr& optimize(ResultExpr&& expr, bool reorder_sum = true);

}  // namespace sequant

#endif  // SEQUANT_OPTIMIZE_OPTIMIZE_HPP
