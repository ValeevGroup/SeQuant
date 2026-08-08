#ifndef SEQUANT_CORE_OPTIMIZE_GA_OPTIMIZE_GA_HPP
#define SEQUANT_CORE_OPTIMIZE_GA_OPTIMIZE_GA_HPP

// Public entry of the genetic factorizer: jointly optimizes the contraction
// order of every summand (L1) and the distributive factoring across summands
// (L2) of one or more targets, scoring candidates over the whole cross-target
// DAG so common-subexpression sharing is inside the objective.
//
// The objective is FLOPs and nothing else -- no storage, footprint or peak
// term enters it. Its one piece of iteration awareness is
// OptimizeOptions::volatile_weight: a merge is charged `volatile_weight` times
// if its cluster contains a volatile leaf and once otherwise, which with no
// volatile-leaf predicate reduces exactly to a single-shot flop count. All
// flop figures are upper bounds: producer resolution is greedy by default,
// exact resolution being exponential in the key-fibre product.

#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/optimize/ga/emit.hpp>
#include <SeQuant/core/optimize/ga/ga.hpp>
#include <SeQuant/core/optimize/options.hpp>

namespace sequant::opt::ga {

struct GAResult {
  /// factorized expression per target; shared arrays appear as named leaf
  /// tensors (see emit.hpp)
  container::svector<ExprPtr> exprs;
  /// named-intermediate definitions, in dependency order; each must be
  /// evaluated before any tree that references its name
  container::svector<ResultExpr> definitions;
  /// total flops of the schedule (shared arrays counted once). With a
  /// volatile-leaf predicate this is the REPLAY-WEIGHTED total,
  /// `persistent_flops + volatile_weight * volatile_flops`.
  double flops = 0;
};

/// Search over the given targets. \p opts supplies idx_to_extent (defaults
/// to IndexSpace::approximate_size).
GAResult optimize_ga(container::svector<TargetInput> const& targets,
                     OptimizeOptions const& opts = {},
                     GAOptions const& ga_opts = {},
                     CostModel const& cost = CostModel::native(),
                     ProducerResolution resolution =
                         ProducerResolution::Greedy);

}  // namespace sequant::opt::ga

namespace sequant {

/// Result of the ResultExpr-based joint GA optimization: the factorized
/// targets (input order) plus the named-intermediate definitions they
/// reference (dependency order).
struct GAOptimized {
  container::vector<ResultExpr> targets;
  container::vector<ResultExpr> definitions;
};

/// Joint GA optimization of several targets; each ResultExpr's expression must
/// be a Sum of Products of Tensors (or a single Product). On return the
/// targets' expressions are their factorized forms, which reference the
/// returned definitions by name.
GAOptimized optimize_ga(container::vector<ResultExpr> exprs,
                        OptimizeOptions const& opts = {},
                        opt::ga::GAOptions const& ga_opts = {});

}  // namespace sequant

#endif  // SEQUANT_CORE_OPTIMIZE_GA_OPTIMIZE_GA_HPP
