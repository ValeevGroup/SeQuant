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
/// \par Scope (what is and isn't pulled out)
/// The matcher operates on the two factors of each summand's *top* binary
/// contraction, as fixed by single-term optimization. Two consequences:
///  - A shared *composite* factor is pulled out only when single-term
///    optimization already exposes it as one of those two top-level factors
///    (e.g. \c A*B*C + A*B*D folds to \c A*B*(C + D) because \c A*B is the
///    shared top factor). A composite that single-term brackets *away* into a
///    deeper subtree is not recovered: \c A*B*C*D + A*B*E*F is left as two
///    products when the chosen evaluation order buries \c A*B (e.g.
///    \c ((A*B)*C)*D), since the two top factors then differ.
///  - Only two-sided bicliques are emitted, and partner sums are not
///    recursively re-factored. An input that mathematically factors into three
///    or more groups, e.g. \c (A+B)(C+D)(E+F) expanded to eight terms, yields a
///    single two-sided fold such as \c (AC+AD+BC+BD)*(E+F), not the fully
///    nested form.
///
/// \param sum   The single-term-optimized sum to factor.
/// \param nodes Per-summand binarized eval nodes, positionally aligned with
///              \c sum.summands() (\c nodes[i] is the binary tree of
///              \c sum.summand(i)). Read-only inputs to the matcher.
/// \param opts  Optimization options. \c opts.idx_to_extent must be populated
///              (asserted). The biclique \c saving is scored with the same base
///              per-contraction cost counter single-term optimization uses --
///              \c flops_counter (DenseFLOPs) or \c memsize_counter
///              (DenseSize), selected by \c opts.objective_function. It does
///              not, however, apply single-term's \c volatile_weight (with \c
///              is_volatile_leaf) or \c footprint_weight adjustments: those
///              fields do not currently influence multi-term factorization.
/// \return The factored sum as an \c ExprPtr.
ExprPtr factorize_multiterm(
    Sum const& sum, container::vector<FullBinaryNode<EvalExpr>> const& nodes,
    OptimizeOptions const& opts);

}  // namespace sequant::opt

#endif  // SEQUANT_CORE_OPTIMIZE_MULTITERM_HPP
