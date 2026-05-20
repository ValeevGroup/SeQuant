#ifndef SEQUANT_OPTIMIZE_OPTIMIZE_HPP
#define SEQUANT_OPTIMIZE_OPTIMIZE_HPP

#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/optimize/flags.hpp>
#include <SeQuant/core/optimize/single_term.hpp>

namespace sequant {

///
/// Optimize the expression using IndexSpace::approximate_size() as the index
/// extent provider.
///
/// \param expr  Expression to be optimized.
/// \param opt_for  Cost metric to optimize for (Flops by default).
/// \param reorder  Whether to reorder summands; on by default.
/// \return Optimized expression for lower evaluation cost.
ExprPtr optimize(ExprPtr const& expr, OptFor opt_for = OptFor::Flops,
                 ReorderSum reorder = ReorderSum::Reorder);

/// Optimize using a caller-supplied index extent provider. Pass an empty
/// \c idx2size to fall back to IndexSpace::approximate_size().
ExprPtr optimize(ExprPtr const& expr, index_to_extent_t idx2size,
                 OptFor opt_for = OptFor::Flops,
                 ReorderSum reorder = ReorderSum::Reorder);

ResultExpr& optimize(ResultExpr& expr, OptFor opt_for = OptFor::Flops,
                     ReorderSum reorder = ReorderSum::Reorder);

ResultExpr& optimize(ResultExpr& expr, index_to_extent_t idx2size,
                     OptFor opt_for = OptFor::Flops,
                     ReorderSum reorder = ReorderSum::Reorder);

[[nodiscard]] ResultExpr& optimize(ResultExpr&& expr,
                                   OptFor opt_for = OptFor::Flops,
                                   ReorderSum reorder = ReorderSum::Reorder);

[[nodiscard]] ResultExpr& optimize(ResultExpr&& expr,
                                   index_to_extent_t idx2size,
                                   OptFor opt_for = OptFor::Flops,
                                   ReorderSum reorder = ReorderSum::Reorder);

}  // namespace sequant

#endif  // SEQUANT_OPTIMIZE_OPTIMIZE_HPP
