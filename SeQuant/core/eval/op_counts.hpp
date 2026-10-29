//
// Op-count and cost probes read by the evaluation tests and the dry-run cost
// model: thread-local diagnostic counters stashed by the dry-run ops and by the
// ordered executor.
//

#ifndef SEQUANT_CORE_EVAL_OP_COUNTS_HPP
#define SEQUANT_CORE_EVAL_OP_COUNTS_HPP

#include <cstddef>

namespace sequant::eval {

namespace detail {

/// Diagnostic: the flops the last product op computed (cm->flops), stashed by
/// DryRunOps::record_op_cost and read by \c CacheManager::tally_build (called
/// from eval.hpp) and the unit tests, so the per-build tally carries the same
/// cost the replay computed -- no second flop path. Thread-local;
/// single-threaded dry-run replay.
inline double& last_op_flops() noexcept {
  static thread_local double f = 0.0;
  return f;
}

/// Diagnostic: the roofline exec-cost estimate (cm->exec_cost) for the last
/// product op, stashed by DryRunOps::record_op_cost alongside
/// \c last_op_flops() and read by \c CacheManager::tally_build (called from
/// eval.hpp) so each per-build tally entry carries the same time estimate the
/// replay computed. Thread-local; single-threaded dry-run replay.
inline double& last_op_exec() noexcept {
  static thread_local double f = 0.0;
  return f;
}

/// Test-facing (and cheap enough to stay live): what the most recent ordered
/// evaluation did, by op kind.
///
///  - \c builds     -- \c compute_cell calls (ordered_executor.hpp): one per
///                     Build production the walk actually performed.
///  - \c assembles  -- Assemble steps performed (a per-batch partial folded
///                     into its output's running result).
///  - \c probes     -- entries into the generic recursive evaluator
///                     \c evaluate_impl (eval.hpp) while the evaluation ran.
///                     The ordered executor computes every cell itself from
///                     table reads, so a completed ordered run leaves this at
///                     zero; a nonzero count means some production fell back
///                     to the tree-walking engine.
///
/// Reset at \c run_ordered_schedule_pre_results entry; thread-local, so a
/// reader must run on the thread that drove the evaluation.
struct OrderedOpCounts {
  std::size_t builds = 0;
  std::size_t assembles = 0;
  std::size_t probes = 0;
};

/// The per-thread counters behind \c ordered_op_counts; the executor tallies
/// into this slot directly.
inline OrderedOpCounts& ordered_op_counts_slot() noexcept {
  static thread_local OrderedOpCounts c;
  return c;
}

/// \return a snapshot of \c ordered_op_counts_slot (see \c OrderedOpCounts).
[[nodiscard]] inline OrderedOpCounts ordered_op_counts() noexcept {
  return ordered_op_counts_slot();
}

}  // namespace detail

}  // namespace sequant::eval

#endif  // SEQUANT_CORE_EVAL_OP_COUNTS_HPP
