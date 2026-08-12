#ifndef SEQUANT_CORE_BATCH_POLICY_HPP
#define SEQUANT_CORE_BATCH_POLICY_HPP

#include <cstddef>
#include <functional>
#include <limits>

namespace sequant {

class Index;
class Tensor;

/// One batchability policy shared by the single-term optimizer and the runtime
/// batched evaluator (make_evaluator, Task A3). All predicates default empty.
struct BatchPolicy {
  /// Spaces batchable in the CONTRACTED role: a mode of such a space is
  /// batchable where it is summed. Companion to \ref
  /// is_batchable_external_index (the EXTERNAL role). Splitting batchability by
  /// role lets a caller admit a space only where batching it is meaningful --
  /// e.g. a space batchable only as an external spectator contributes none of
  /// its contracted occurrences to the optimizer's 2^m search. Building block;
  /// the derived "batchable in any role" query is \ref is_batchable_index().
  /// Defaults to decline every index; a caller opts spaces in explicitly.
  std::function<bool(Index const&)> is_batchable_contracted_index =
      [](Index const&) { return false; };
  /// Spaces batchable in the EXTERNAL role: a mode of such a space is batchable
  /// where it is open on the term root (a spectator carried to the result), not
  /// where it is contracted. Building block; declared adjacent to its
  /// contracted companion. Defaults to decline every index; a caller that wants
  /// external batching sets this predicate explicitly (there is no fallback to
  /// the contracted role).
  std::function<bool(Index const&)> is_batchable_external_index =
      [](Index const&) { return false; };

  /// Derived "batchable in ANY role": the union of the two building-block
  /// predicates. This is NEVER a settable field -- it is computed from
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
  /// false = no spectator batching (byte-identical to non-spectator behavior).
  /// Necessary but not sufficient: spectator axes are emitted only under a
  /// TIME-FIRST objective (DenseTimeSpaceBatched) and only when the selected
  /// root's modeled peak exceeds \c peak_threshold. Spectator batching is
  /// therefore currently unavailable under the space-first objectives.
  bool batch_spectator_indices = false;

  /// Enable the order-aware multilevel recompute cost model (resident-scan peak
  /// + ordered-key flops recompute). SELECTION knob ONLY: it makes the DP
  /// charge recompute realistically and thus pick a different (better-batching)
  /// factorization. It does NOT control external-mode EMISSION -- that is the
  /// independent \ref node_level_placement. Consulted only by the batched
  /// objectives (threaded via CostParams). Default TRUE: the recompute-aware
  /// model is the more realistic cost for selection. This is SAFE precisely
  /// because it is now selection-only -- the node-level emission it used to
  /// force is separately gated by \ref node_level_placement (default off), so
  /// the emission stays the correct, cheap root-level forest seed. (Before the
  /// decouple, defaulting this true forced the node-level runtime regression.)
  bool order_aware_recompute = true;

  /// Emission-placement knob for external (spectator) modes, INDEPENDENT of the
  /// order-aware cost model. Only meaningful with \ref batch_spectator_indices.
  /// If true, the emit uses node-level placement (per-node External stamps); if
  /// false (default) it uses the root-level forest seed (one global spectator
  /// loop). Node-level placement is currently a net runtime REGRESSION -- ~6x
  /// wall time and ~8x batch scopes on water-8, and it produces a wrong
  /// residual on water-20 -- because it nests a batch scope at every carrying
  /// node and the batched evaluator replays each. It stays OFF by default until
  /// that is fixed; the root-seed emission is correct and cheap regardless of
  /// order_aware_recompute.
  bool node_level_placement = false;

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

  /// Coexistence switch (Task 6 of the whole-scope batched DAG execution
  /// design, `doc/dev/specs/2026-08-10-whole-scope-batched-dag-execution-
  /// design.md`) between the two RUNTIME EXECUTION MODELS: forest descent
  /// (default false -- one tree at a time, `sequant::evaluate(Nodes const&,
  /// ...)`, unchanged) and whole-scope descent (true -- one fused scope-tree
  /// walk over the whole forest, `sequant::eval::evaluate_whole_scope`, so a
  /// value shared across trees is built once per home block and reused,
  /// rather than rebuilt per tree). Consulted by the `sequant::evaluate(
  /// Nodes const&, BatchPolicy const&, ...)` driver overload
  /// (`scope_executor.hpp`) to select the driver, and by
  /// `sequant::eval::dryrun::cost_profile()` to select the matching peak
  /// model: the co-residency oracle (`peak_profile_sweep` over `home_modes`)
  /// when true, since that model is what predicts the whole-scope realized
  /// peak, vs the batched-scratch replay high-watermark (models forest
  /// descent) when false. Purely additive: false reproduces today's behavior
  /// on both call sites byte-for-byte.
  bool whole_scope_execution = false;

  /// SP3 gating switch (`doc/dev/specs/2026-08-05-dryrun-wetrun-schedule-
  /// equivalence-design.md` follow-on, the ordered-scope batched-eval
  /// design): between forest descent / whole-scope descent (both selected
  /// above via \ref whole_scope_execution) and the new ORDERED executor
  /// (true -- `sequant::eval::evaluate_ordered_schedule`, driven by the SP2
  /// `eval::OrderedSchedule` IR rather than the narrow `ScopeSchedule` scope
  /// tree). Consulted by the `sequant::evaluate(Nodes const&, BatchPolicy
  /// const&, ...)` driver overload (`scope_executor.hpp`) BEFORE \ref
  /// whole_scope_execution, so it takes priority when both are set (setting
  /// both is well-defined: the ordered executor wins) and the two pre-
  /// existing dispatch arms (forest descent / whole-scope descent) are
  /// reached, byte-identically, only when this flag is false. Default false
  /// reproduces today's dispatch byte-for-byte on every existing caller.
  bool ordered_schedule_execution = false;

  /// Peak-memory budget in BYTES for the batched objectives. Its meaning
  /// DIFFERS between them:
  ///
  /// - SPACE-FIRST (DenseSpaceTimeBatched): a hard feasibility gate. The
  ///   single-term optimizer minimizes flops among schedules whose modeled peak
  ///   is <= peak_threshold, falling back to min-peak (best effort) when none
  ///   fit. Default +infinity => every schedule feasible => min flops => no
  ///   batching, i.e. here a finite value is the *enable* trigger for batching.
  ///
  /// - TIME-FIRST (DenseTimeSpaceBatched): NOT a feasibility gate. Root
  ///   selection ignores it entirely (peak breaks exact flop ties only), so it
  ///   can neither constrain the schedule's peak nor enable CONTRACTED-axis
  ///   batching (which is emitted regardless). Its ONLY effect is to trigger
  ///   EXTERNAL (spectator) axis emission, together with
  ///   \c batch_spectator_indices: axes are emitted iff the selected root's
  ///   modeled peak exceeds this threshold.
  double peak_threshold = std::numeric_limits<double>::infinity();
};

}  // namespace sequant

#endif  // SEQUANT_CORE_BATCH_POLICY_HPP
