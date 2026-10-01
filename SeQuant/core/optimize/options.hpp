#ifndef SEQUANT_CORE_OPTIMIZE_OPTIONS_HPP
#define SEQUANT_CORE_OPTIMIZE_OPTIONS_HPP

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/fwd.hpp>
#include <SeQuant/core/utility/aggregate.hpp>

#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <unordered_map>
#include <utility>

namespace sequant {

class Index;
class Tensor;
class Expr;

/// Objective function to minimize in single-term and top-level optimize
/// routines. The `Dense*` models assume dense tensors:
/// - `DenseFLOPs` counts floating-point operations.
/// - `DenseSize` counts result-tensor storage elements (summed over
///   intermediates) -- a gross-traffic proxy, not a peak.
/// - `DenseSpaceTime` (also spelled `DensePeakSize`, a deprecated alias) is
///   peak-first, perf-second: it minimizes peak memory --
///   the maximum over the evaluation schedule of the combined size of all
///   simultaneously-live tensors (intermediates and resident input leaves, the
///   all-co-resident model) -- and breaks ties by the roofline perf cost.
///   Unlike the order-independent `DenseFLOPs`/`DenseSize`, the contraction
///   order is a real lever here. Note: does not yet support
///   common-subexpression elimination (`CSEOptions::subnet` must be false).
/// - `DenseSpaceTimeBatched` (deprecated alias `DensePeakSizeBatched`) extends
///   `DenseSpaceTime` with a per-index batchability model: each index
///   satisfying `OptimizeOptions::batch_policy.is_batchable_index()` (the
///   derived union of the contracted- and external-role building blocks) is
///   treated as independently sliced to `min(extent,
///   batch_policy.batch_target_size(ix))` elements per index --
///   `batch_target_size` is an upper bound, so this is a conservative
///   (over-)estimate of the realized whole-tile batch, which the backend rounds
///   *down* to a tile multiple (never above the target; see
///   `mode_batches_of_trange1`). The
///   DP minimises peak over the worst-case sliced configuration. Only consulted
///   by the batched oracle and DP; requires
///   a batchability role predicate
///   (`batch_policy.is_batchable_contracted_index` and/or
///   `is_batchable_external_index`) and `batch_policy.batch_target_size` to be
///   set. Final selection is ceiling-gated by `peak_threshold`, which must be
///   finite for any loop to be opened at all.
/// - `DenseTimeSpace` / `DenseTimeSpaceBatched` are the perf-first, peak-second
///   duals of `DenseSpaceTime` / `DenseSpaceTimeBatched`: they select a
///   factorization by roofline perf first and peak second (same Pareto-frontier
///   and roofline machinery, opposite lexicographic order). Because slicing is
///   perf-neutral, a perf-first primary never prefers a flops-catastrophic
///   factorization merely for its sliceability. Naming:
///   `Dense{Primary}{Secondary}`, `Space` = peak/size, `Time` = perf.
///
///   `peak_threshold` and these objectives: it is a feasibility ceiling on the
///   perf-first branch of `select_root` too -- among the frontier points that
///   fit the budget it takes the fewest flops, ties broken toward the
///   least-sliced realization (`nsl`) and then lower peak, so nothing is
///   sliced for free below the ceiling. When no point fits, the fallback stays
///   perf-first: global min flops, ties by min peak. It cannot force a
///   flops-catastrophic factorization the way the space-first branch can.
///
///   `peak_threshold` is also the switch that turns batching on at all:
///   `PeakBatchedModel::relax` opens neither a contracted nor an external loop
///   unless `std::isfinite(peak_threshold)`. External opens are emitted per
///   node, gated by exactly
///   `batch_policy.batch_spectator_indices && std::isfinite(peak_threshold)`
///   -- under both batched objectives, with no objective condition and no
///   post-DP external placement pass.
///
/// Leaves room for `Sparse*` models later.
enum class ObjectiveFunction {
  DenseFLOPs,
  DenseSize,
  /// Peak-first, perf-second. Deprecated alias: `DensePeakSize`.
  DenseSpaceTime,
  /// Batched variant of `DenseSpaceTime`. Deprecated alias:
  /// `DensePeakSizeBatched`.
  DenseSpaceTimeBatched,
  /// Perf-first, peak-second: never prefers a flops-catastrophic factorization
  /// for its sliceability.
  DenseTimeSpace,
  /// Batched variant of `DenseTimeSpace`.
  DenseTimeSpaceBatched,
  /// Deprecated aliases, sharing the underlying values of the constants above
  /// so that code and JSON inputs spelling them ("dense_peak_size") keep
  /// working: a `== DensePeakSize` guard catches `DenseSpaceTime`.
  DensePeakSize = DenseSpaceTime,
  DensePeakSizeBatched = DenseSpaceTimeBatched
};

/// Whether to reorder summands so terms with shared intermediates appear
/// closer to each other.
enum class ReorderSum { Reorder, NoReorder };

/// Common-subexpression-elimination (CSE) options for single-term
/// optimization. `subnet` recognizes equivalent subnetworks while searching for
/// an evaluation order, trading extra search time for potentially lower op
/// counts. (Room to grow, e.g. a maximum subnet size to consider.)
struct CSEOptions {
  SEQUANT_DESIGNATED_INIT_ONLY;
  bool subnet = false;
};

/// Roofline parameters for the peak objectives' secondary (tie-break) cost.
/// When \c machine_balance > 0, the per-contraction tie-break cost becomes the
/// roofline wall-time proxy \c max(flops, machine_balance * Q), where the data
/// movement \c Q = max(operand+result footprint, block_prefactor * flops /
/// sqrt(fast_mem_elems / block_tiles)) combines the compulsory single-pass
/// traffic with the finite-cache (Hong-Kung) re-read bound. This charges
/// bandwidth-bound contractions (e.g. single-PNO-index ones) their true memory
/// traffic while leaving compute-bound (dense) contractions at \c flops, so it
/// is inert in the dense case. \c machine_balance == 0 (default) recovers the
/// pure-flop tie-break. See doc/dev/specs/2026-06-23-roofline-tiebreak-cost.md.
struct RooflineParams {
  SEQUANT_DESIGNATED_INIT_ONLY;
  /// Machine balance beta = 8*F/B in FLOPs per element of traffic. 0 = off.
  /// Calibrated in 8-byte (real double) elements, as is \c fast_mem_elems: a
  /// complex network is priced with twice the traffic width and four real
  /// flops per multiply-add (see opt::detail::FieldCostFactors), so one
  /// calibration serves both fields.
  double machine_balance = 0.0;
  /// Capacity M of the binding fast memory level, in elements (e.g. LLC/8).
  double fast_mem_elems = 0.0;
  /// Resident-tile count c0 in the blocking bound (b ~ sqrt(M/c0)); ~3 for the
  /// A,B,C tiles of a blocked GEMM. Calibratable.
  double block_tiles = 3.0;
  /// Prefactor kappa folding FMA/packing/BLAS constants into the re-read term.
  double block_prefactor = 1.0;
};

/// Cost-model tuning knobs consumed by \ref opt::single_term_opt and forwarded
/// to the objective models. Bundled so callers (and \ref sequant::optimize via
/// \ref OptimizeOptions) pass one object rather than five positional arguments.
/// All have neutral defaults: an empty \c is_volatile_leaf disables replay
/// weighting, and \c roofline.machine_balance == 0 keeps the pure-flop
/// tie-break.
struct CostParams {
  SEQUANT_DESIGNATED_INIT_ONLY;
  /// Marks a LEAF tensor as volatile (amplitude-dependent), so the contraction
  /// forming any subset that contains it is replayed every iteration. Empty =>
  /// nothing volatile (replay weighting off).
  std::function<bool(Tensor const&)> is_volatile_leaf = {};
  /// Replay weight on volatile contractions (conceptually the replay count).
  double volatile_weight = 1.0;
  /// Per-intermediate storage-footprint penalty (DenseFLOPs/DenseSize only; see
  /// OptimizeOptions::footprint_weight). Not used by the peak objectives.
  double footprint_weight = 0.0;
  /// Relative peak tolerance for DenseSpaceTime's final selection; see
  /// OptimizeOptions::peak_flops_tolerance. Unused by DenseSpaceTimeBatched,
  /// whose final selection is instead threshold-gated by \c peak_threshold.
  double peak_flops_tolerance = 0.10;
  /// Roofline parameters for the peak objectives' secondary cost; see
  /// \ref RooflineParams. machine_balance == 0 => pure-flop tie-break.
  RooflineParams roofline = {};
  /// In-flight batch-contribution footprint multiplier for the batched
  /// objectives; see BatchPolicy::accumulation_factor. 0 (default) = no
  /// penalty.
  double accumulation_factor = 0.0;
  /// Peak-memory budget in bytes for the batched objectives; see
  /// BatchPolicy::peak_threshold. A feasibility ceiling on root selection under
  /// both batched objectives (space-first minimizes flops among the points that
  /// fit, ties by lower peak; time-first minimizes flops, ties by fewer slices
  /// then lower peak), and the enable switch for batching itself: +infinity
  /// (the default) opens no loop at all.
  double peak_threshold = std::numeric_limits<double>::infinity();
  /// Prune disconnected (outer-product) subsets from the single-term DP; see
  /// OptimizeOptions::prune_outer_products. true (default) = prune.
  bool prune_outer_products = true;
  /// Gate for external-index batching; see
  /// BatchPolicy::batch_spectator_indices. false (default) => the DP opens no
  /// external loop, so \ref opt::detail::PeakBatchedModel emits no \c
  /// BatchModeType::External entries, so a non-external-aware caller sees
  /// none. Note this is necessary but not sufficient:
  /// the per-node gate in \c PeakBatchedModel::relax is
  /// `batch_spectator_indices && std::isfinite(peak_threshold)`, so a finite
  /// \c peak_threshold is required too. It is not conditioned on the
  /// objective, and there is no post-DP external placement pass.
  bool batch_spectator_indices = false;

  /// Spaces batchable in the contracted role, threaded from
  /// BatchPolicy::is_batchable_contracted_index. Building block consumed by the
  /// single-term optimizer's DP contracted-role filter (a mode of such a space
  /// is sliced where it is summed). Declared adjacent to its external
  /// companion. Defaults to decline every index.
  std::function<bool(Index const&)> is_batchable_contracted_index =
      [](Index const&) { return false; };
  /// Spaces batchable in the external role, threaded from
  /// BatchPolicy::is_batchable_external_index. Defaults to decline every index;
  /// a caller that wants external batching sets it explicitly (there is no
  /// fallback to the contracted-role predicate).
  std::function<bool(Index const&)> is_batchable_external_index =
      [](Index const&) { return false; };
  /// Per-index per-batch slice size (an upper bound) for a batchable index;
  /// threaded from BatchPolicy::batch_target_size. Consulted only by the
  /// batched objectives. Empty (default) => no target => no slicing.
  std::function<std::size_t(Index const&)> batch_target_size = {};
  /// k-aware inner (CSV/PNO composite) extent applied by every cost counter;
  /// threaded from OptimizeOptions::inner_pow. Required whenever the network
  /// has composite indices: empty does not fall back to the idxsz provider
  /// (k=1) -- inner_aware_volume throws instead, since such a fallback grossly
  /// mis-sizes multi-composite tensors, e.g. a 4-PAO integral. No default
  /// (matching OptimizeOptions::inner_pow and PeakBatchedModel::inner_pow);
  /// pass an explicit no-op only for composite-free work.
  std::function<double(Index const&, std::size_t)> inner_pow = {};
  /// When true, only persistent (volatile-leaf-free) subnetworks are batched;
  /// threaded from BatchPolicy::persistent_only. Batched objectives only.
  bool batch_persistent_only = false;
};

/// A type-erased provider mapping an Index to its extent. Used by the public
/// optimize() API. Callers reaching for the templated opt::single_term_opt
/// overloads (constrained by \ref opt::has_index_extent) should pass the
/// callable directly instead of wrapping it here — that keeps the cost
/// function's call site inlineable, whereas a value of this alias goes
/// through std::function's type-erased dispatch on every Index lookup.
using index_to_extent_t = std::function<std::size_t(Index const&)>;

/// Options that control behavior of \ref sequant::optimize.
struct OptimizeOptions {
  SEQUANT_DESIGNATED_INIT_ONLY;
  /// Objective function to minimize.
  ObjectiveFunction objective_function = ObjectiveFunction::DenseFLOPs;

  /// Whether to reorder summands so terms with shared intermediates appear
  /// closer to each other.
  ReorderSum reorder = ReorderSum::Reorder;

  /// Common-subexpression-elimination options. All disabled by default;
  /// enabling can reduce op counts at the cost of additional optimization time.
  CSEOptions CSE = {};

  /// Caller-supplied Index to extent provider. If empty, defaults to
  /// \c IndexSpace::approximate_size().
  index_to_extent_t idx_to_extent = {};

  /// Optional k-aware inner-composite extent for CSV/PNO tensor-of-tensor
  /// indices, used by the peak objectives' footprint accounting. For a group of
  /// \c k composites sharing a proto-index set, \c inner_pow(composite, k)
  /// should return the k-th power mean of the per-pair domain
  /// (\c (Sum_pairs d^k / nocc^N)^(1/k)), so that the product over the group
  /// times the outer \c nocc^N equals the true block-sparse volume
  /// \c Sum_pairs d^k. Required whenever the network has composite indices:
  /// leaving it empty does not fall back to \ref idx_to_extent (which would
  /// grossly mis-size multi-composite tensors and invert factorization
  /// choices, e.g. a 4-PAO integral) -- the sizing code throws instead. No
  /// default: pass an explicit no-op only for composite-free work.
  std::function<double(Index const&, std::size_t)> inner_pow = {};

  /// Batchability policy: bundles the per-index and per-leaf predicates
  /// that govern batched evaluation. All predicate fields default to empty (no
  /// batchable indices, no volatile leaves). The key sub-fields are:
  ///   - `is_batchable_contracted_index` / `is_batchable_external_index`: the
  ///     two role building blocks marking an Index as living in a batchable
  ///     space in the contracted / external role. The derived union
  ///     `is_batchable_index()` (= the eval cache's accept_aux) is the runtime
  ///     accept.
  ///   - `batch_target_size`: per-index slice size (an upper bound); a sliced
  ///     batchable index ix contributes min(extent, batch_target_size(ix)), a
  ///     conservative over-estimate of the realized tile-floored batch. Only
  ///     consulted by DenseSpaceTimeBatched.
  ///   - `is_volatile_leaf`: marks a LEAF tensor as volatile (its value
  ///     changes between replays). Empty => no tensor is volatile => cost
  ///     weighting is disabled and volatile_weight is ignored. CC callers pass
  ///     label==volatile_label, the same classification the runtime eval cache
  ///     uses, so the optimizer's cost model and the cache agree.
  BatchPolicy batch_policy = {};

  /// Real-valued weight on the cost of each volatile contraction (re-evaluated
  /// on every replay of the network), while persistent (volatile-independent)
  /// contractions are counted once. Conceptually the expected number of
  /// replays. Default 1.0 (no change). Consulted whenever
  /// batch_policy.is_volatile_leaf is non-empty, by DenseFLOPs and by the peak
  /// objectives: it scales the per-contraction cost in all of them (see
  /// AdditiveModel and PeakModel/PeakBatchedModel's `cflops`). That cost is the
  /// primary mode for the time-first objectives (so volatile_weight materially
  /// changes which factorization DenseTimeSpace* picks) but only the tie-break
  /// among equal-peak schedules for the space-first ones.
  double volatile_weight = 1.0;

  /// Relative peak tolerance for the peak objectives' final selection: among
  /// the Pareto frontier of (peak, flops) trade-offs, pick the fewest-flops
  /// schedule whose peak is within (1 + peak_flops_tolerance) of the minimum.
  /// 0 = strict peak-min (flop tie-break only on exact peak ties). The default
  /// 0.10 trades up to a 10% peak increase for a (often much larger) flop
  /// reduction -- e.g. forming a persistent composite integral instead of
  /// recomputing a ladder term. Only consulted by \c DenseSpaceTime. \c
  /// DenseSpaceTimeBatched's final selection is instead ceiling-gated by
  /// \ref BatchPolicy::peak_threshold.
  double peak_flops_tolerance = 0.10;

  /// Per-intermediate memory-footprint penalty added to the single-term
  /// optimization cost. For every binary contraction, the storage footprint of
  /// its RESULT intermediate (the product of the extents of the result's
  /// indices, i.e. its element count) is multiplied by this weight and added to
  /// the contraction cost. Unlike the FLOPs cost, this penalty is NOT scaled by
  /// \ref volatile_weight (peak footprint is a one-time materialization cost,
  /// not a per-replay one). 0 (default) disables the penalty, recovering the
  /// pure FLOPs/Size behavior.
  ///
  /// Rationale: the FLOPs cost is blind to the storage size of the
  /// intermediates it materializes, so it will happily pick an order (and thus
  /// expose, as a shareable subtree, a common subexpression) that carries a
  /// free large-space index -- e.g. a half-transformed DF integral that still
  /// carries a free projected-AO (PAO) index -- because forming it once is
  /// FLOPs-cheap. Such an intermediate can be enormous. A nonzero
  /// footprint_weight biases single-term optimization toward orders that defer
  /// or avoid materializing such large intermediates (e.g. transforming both
  /// large legs before exposing a shared subtree), trading a controlled amount
  /// of extra FLOPs for a lower peak footprint. Only consulted when
  /// objective_function == ObjectiveFunction::DenseFLOPs; the units are
  /// FLOPs-per-element, so a
  /// useful magnitude is on the order of the contracted-index extent that the
  /// offending intermediate would otherwise leave free.
  double footprint_weight = 0.0;

  /// Roofline parameters for the peak objectives' secondary (tie-break) cost;
  /// see \ref RooflineParams. machine_balance == 0 (default) => pure-flop
  /// tie-break (no behavior change). Consulted only by DenseSpaceTime /
  /// DenseSpaceTimeBatched.
  RooflineParams roofline = {};

  /// Optional out-channel: when non-null, optimize() records for each
  /// top-level summand's optimized ExprPtr its per-contraction-node
  /// sliced-sets (RPN / post-order, left-first, matching single_term_opt's
  /// Product construction). Consumed by binarize via BinarizationOptions to
  /// annotate eval nodes. Default null => no behavior change. Populated under
  /// either batched objective (\c DenseSpaceTimeBatched and \c
  /// DenseTimeSpaceBatched both thread it through as \c out_axes); every
  /// other objective leaves any pre-existing map entry for that summand
  /// untouched (it is simply never written).
  std::shared_ptr<std::unordered_map<sequant::Expr const*,
                                     container::vector<NodeBatchAnnotation>>>
      term_batch_axes = {};

  /// Prune disconnected (outer-product) subsets from the single-term
  /// contraction DP: a subset whose induced subgraph is disconnected under the
  /// "share a contractible (non-target) index" relation is an outer product the
  /// optimal tree never forms, so skipping it prunes search space without
  /// changing the result. \c true (default) prunes; \c false searches the full
  /// (unpruned) space and is exposed mainly for validation. (The environment
  /// variable \c SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING force-disables pruning
  /// regardless of this flag.)
  bool prune_outer_products = true;
};

}  // namespace sequant

#endif  // SEQUANT_CORE_OPTIMIZE_OPTIONS_HPP
