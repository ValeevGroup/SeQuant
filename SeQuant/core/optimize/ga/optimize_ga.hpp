#ifndef SEQUANT_CORE_OPTIMIZE_GA_OPTIMIZE_GA_HPP
#define SEQUANT_CORE_OPTIMIZE_GA_OPTIMIZE_GA_HPP

// Public entry of the genetic factorizer: jointly optimizes the contraction
// order of every summand (L1) and the distributive factoring across summands
// (L2) of one or more targets, scoring candidates over the whole cross-target
// DAG so common-subexpression sharing is inside the objective. The staged
// pipeline (optimize()) remains the per-term seed and the baseline.
//
// The objective is FLOPs and nothing else -- no storage, footprint or peak term
// enters it (peak is a batching concern, handled downstream). The one piece of
// iteration awareness it has is OptimizeOptions::volatile_weight: for an
// iterative solver the DAG is replayed once per iteration but only its
// amplitude-dependent part is rebuilt, so a merge is charged
// `volatile_weight` times if its cluster contains a volatile leaf and once
// otherwise (see CostModel::volatile_weight). With no volatile-leaf predicate
// this reduces exactly to the previous single-shot flop count.
//
// All flop figures are upper bounds: producer resolution (which cluster
// builds a shared array) is greedy by default; exact resolution is
// exponential in the key-fibre product and only feasible for small universes.

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
  /// evaluated (once per set of leaf values) before any tree that references
  /// its name
  container::svector<ResultExpr> definitions;
  /// the winning genome and its decoded schedule
  Genome genome;
  Schedule schedule;
  /// total flops of the schedule (shared arrays counted once). When a
  /// volatile-leaf predicate is supplied this is the REPLAY-WEIGHTED total,
  /// `persistent_flops + volatile_weight * volatile_flops` -- still a pure flop
  /// count, just with amplitude-dependent work charged once per replay.
  double flops = 0;
  /// total flops with no cross-term sharing (each use counted); what the
  /// emitted forest would cost with every definition inlined at every use.
  /// Weighted on the same convention as \ref flops, except that its L2
  /// volatility is read off the cluster mask rather than propagated through the
  /// Val tree -- it is a diagnostic, not the objective.
  double flops_no_sharing = 0;
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
/// targets (in input order) plus the named-intermediate definitions they
/// reference (in dependency order; see opt::ga::emit_named).
struct GAOptimized {
  container::vector<ResultExpr> targets;
  container::vector<ResultExpr> definitions;
};

/// Joint GA optimization of several targets; each ResultExpr's expression
/// must be a Sum of Products of Tensors (or a single Product). On return the
/// targets' expressions are replaced by their factorized forms, which
/// reference the returned definitions by name.
GAOptimized optimize_ga(container::vector<ResultExpr> exprs,
                        OptimizeOptions const& opts = {},
                        opt::ga::GAOptions const& ga_opts = {});

}  // namespace sequant

#endif  // SEQUANT_CORE_OPTIMIZE_GA_OPTIMIZE_GA_HPP
