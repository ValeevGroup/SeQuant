#ifndef SEQUANT_CORE_BATCH_POLICY_HPP
#define SEQUANT_CORE_BATCH_POLICY_HPP

#include <SeQuant/core/utility/aggregate.hpp>

#include <cstddef>
#include <functional>
#include <limits>

namespace sequant {

class Index;
class Tensor;

/// The two runtime execution models for batched evaluation (see
/// `doc/dev/specs/2026-09-12-batched-array-dag-eval-as-built.md`,
/// sections 8 and 12.1):
///   - \c forest_descent (default): one tree at a time,
///     `sequant::evaluate(Nodes const&, ...)`, unchanged.
///   - \c ordered: one fused, table-driven walk over the whole forest,
///     driven by the `eval::OrderedSchedule` IR
///     (`sequant::eval::evaluate_ordered_schedule`), so a value shared across
///     trees is built once per home block and reused, rather than rebuilt
///     per tree.
enum class BatchScheduler { forest_descent, ordered };

/// One batchability policy shared by the single-term optimizer and the runtime
/// batched evaluator (make_evaluator, Task A3). All predicates default empty.
struct BatchPolicy {
  SEQUANT_DESIGNATED_INIT_ONLY;
  /// Spaces batchable in the contracted role: a mode of such a space is
  /// batchable where it is summed. Companion to \ref
  /// is_batchable_external_index (the external role). Splitting batchability by
  /// role lets a caller admit a space only where batching it is meaningful --
  /// e.g. a space batchable only as an external spectator contributes none of
  /// its contracted occurrences to the optimizer's 2^m search. Building block;
  /// the derived "batchable in any role" query is \ref is_batchable_index().
  /// Defaults to decline every index; a caller opts spaces in explicitly.
  std::function<bool(Index const&)> is_batchable_contracted_index =
      [](Index const&) { return false; };
  /// Spaces batchable in the external role: a mode of such a space is batchable
  /// where it is open on the term root (a spectator carried to the result), not
  /// where it is contracted. Building block; declared adjacent to its
  /// contracted companion. Defaults to decline every index; a caller that wants
  /// external batching sets this predicate explicitly (there is no fallback to
  /// the contracted role).
  std::function<bool(Index const&)> is_batchable_external_index =
      [](Index const&) { return false; };

  /// Derived "batchable in any role": the union of the two building-block
  /// predicates. This is never a settable field -- it is computed from
  /// \ref is_batchable_contracted_index and \ref is_batchable_external_index.
  /// The runtime batched evaluator's accept predicate is this union (a mode is
  /// accepted at runtime if it is batchable in either role); the factorizer's
  /// role filters instead consume the individual building blocks. The building
  /// blocks default-decline, so both are always callable here.
  std::function<bool(Index const&)> is_batchable_index() const {
    auto contracted = is_batchable_contracted_index;
    auto external = is_batchable_external_index;
    return [contracted, external](Index const& ix) {
      return contracted(ix) || external(ix);
    };
  }
  /// Per-index per-batch slice size (in elements) for a batchable index -- an
  /// UPPER BOUND, not a goal. Both the single-term optimizer and the runtime
  /// batched evaluator treat it as a ceiling: the realized whole-tile batch is
  /// rounded *down* to a tile multiple and never exceeds this value, except the
  /// one-tile floor (a lone tile larger than the target forms its own batch).
  std::function<std::size_t(Index const&)> batch_target_size = {};
  std::function<bool(Tensor const&)> is_volatile_leaf = {};

  /// If true, an external/spectator index -- open on the whole network's result
  /// yet contracted at no node -- is eligible for batching; its per-slice size
  /// comes from \c batch_target_size(ix) like any batchable index. Default
  /// false = no spectator batching.
  /// Necessary but not sufficient: the DP opens externals per node, and only
  /// where \c peak_threshold is finite -- the gate inside \c
  /// PeakBatchedModel::relax is exactly
  /// `batch_spectator_indices && std::isfinite(peak_threshold)`
  /// (\c optimize/cost_model.hpp), with no objective condition and no
  /// post-DP placement pass. Both batched objectives therefore admit
  /// spectator axes; an infinite budget admits none.
  bool batch_spectator_indices = false;

  /// If true, restrict batching to persistent (amplitude-independent) subtrees,
  /// declining to batch any subtree that contains a volatile leaf. If false
  /// (the default), batch ACROSS THE BOARD: slicing the batch axis shrinks any
  /// intermediate carrying it regardless of volatility (footprint objective)
  /// and leaves flops unchanged, so the persistence gate would only ever raise
  /// the modelled/realized peak. Set true to recover the persistent-only
  /// behavior (amortizes the per-replay partition + relaxed-screening cost over
  /// many reuses, at the price of a higher peak for volatile intermediates).
  /// Read identically by the single-term optimizer and the runtime evaluator.
  bool persistent_only = false;

  /// Footprint multiplier for the in-flight batch contribution that co-resides
  /// with a batch-accumulated intermediate (K += contribution). 0 = ignore
  /// (default); ~1 = full contribution materialized; backend-specific (TA's
  /// eager tile accumulation lowers it, multiple in-flight Summa steps raise it
  /// ~30%). Read by the single-term optimizer's PeakBatchedModel to price the
  /// accumulator + contribution co-residency of a node that contracts a
  /// batchable index.
  double accumulation_factor = 0.0;

  /// Selects between the two runtime execution models (\ref BatchScheduler
  /// above). Consulted by the
  /// `sequant::evaluate(Nodes const&, BatchPolicy const&, ...)` driver
  /// overload (`ordered_executor.hpp`) to select the driver. Default
  /// \c forest_descent selects the forest-descent evaluator.
  BatchScheduler scheduler = BatchScheduler::forest_descent;

  /// Peak-memory budget in bytes. It is a feasibility ceiling under both
  /// batched objectives, and it is the single knob that turns batching on:
  /// \c PeakBatchedModel::relax opens neither a contracted nor an external
  /// loop unless `std::isfinite(peak_threshold)`, so the default +infinity
  /// means no batching at all.
  ///
  /// - space-first (\c DenseSpaceTimeBatched): among the frontier points whose
  ///   modeled byte peak is <= peak_threshold, minimize flops, ties broken by
  ///   lower peak; fall back to global min-peak (best effort) when none fit.
  ///
  /// - time-first (\c DenseTimeSpaceBatched): among the frontier points whose
  ///   modeled byte peak is <= peak_threshold, minimize flops, ties broken
  ///   toward the least-sliced realization (\c nsl) and then lower peak --
  ///   so a schedule is not sliced for free below the ceiling. When nothing
  ///   fits, the fallback keeps the perf-first character: global min flops,
  ///   ties by min peak (accepting the overage).
  ///
  /// See \c PeakBatchedModel::select_root (\c optimize/cost_model.hpp) and
  /// the as-built design, \c
  /// doc/dev/specs/2026-09-12-batched-array-dag-eval-as-built.md section 4.4.
  double peak_threshold = std::numeric_limits<double>::infinity();
};

}  // namespace sequant

#endif  // SEQUANT_CORE_BATCH_POLICY_HPP
