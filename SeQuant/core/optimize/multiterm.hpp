#ifndef SEQUANT_CORE_OPTIMIZE_MULTITERM_HPP
#define SEQUANT_CORE_OPTIMIZE_MULTITERM_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/optimize/options.hpp>

namespace sequant {
class Sum;
}

namespace sequant::opt {

/// \brief Multi-term factorization over single-term-optimized summands.
///
/// Pulls shared factors across the summands of \p sum (e.g.
/// \c A*B + A*C -> A*(B + C)), generalized to N terms at once and to the
/// two-sided biclique case \c (A+B)*(X+Y), and emits the result in place as a
/// nested \c ExprPtr (a \c Product whose factor is a \c Sum). Factorizations
/// are applied only when they lower the evaluation cost (cost-driven biclique
/// selection); structurally-shareable factors whose saving is non-positive are
/// left untouched.
///
/// \param sum   The single-term-optimized sum to factor.
/// \param nodes Per-summand binarized eval nodes, positionally aligned with
///              \c sum.summands() (\c nodes[i] is the binary tree of
///              \c sum.summand(i)). Read-only inputs to the matcher.
/// \param opts  Optimization options. \c opts.idx_to_extent must be populated
///              (asserted); the same cost model that drives single-term
///              optimization drives the biclique \c saving here.
/// \return The factored sum as an \c ExprPtr.
ExprPtr factorize_multiterm(
    Sum const& sum, container::vector<FullBinaryNode<EvalExpr>> const& nodes,
    OptimizeOptions const& opts);

}  // namespace sequant::opt

#endif  // SEQUANT_CORE_OPTIMIZE_MULTITERM_HPP
