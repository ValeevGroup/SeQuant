#ifndef SEQUANT_CORE_BATCH_POLICY_HPP
#define SEQUANT_CORE_BATCH_POLICY_HPP

#include <SeQuant/core/utility/aggregate.hpp>

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
  /// + ordered-key flops recompute). false (default) => byte-identical
  /// set-keyed DP. Consulted only by the batched objectives (threaded via
  /// CostParams).
  bool order_aware_recompute = false;

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
