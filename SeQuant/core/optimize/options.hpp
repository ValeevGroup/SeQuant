#ifndef SEQUANT_CORE_OPTIMIZE_OPTIONS_HPP
#define SEQUANT_CORE_OPTIMIZE_OPTIONS_HPP

#include <SeQuant/core/batch_policy.hpp>

#include <cstddef>
#include <functional>

namespace sequant {

class Index;
class Tensor;

/// Objective function to minimize in single-term and top-level optimize
/// routines. The `Dense*` models assume dense tensors:
/// - `DenseFLOPs` counts floating-point operations.
/// - `DenseSize` counts result-tensor storage elements (summed over
///   intermediates) -- a gross-traffic proxy, not a peak.
/// - `DensePeakSize` minimizes peak memory: the maximum over the evaluation
///   schedule of the combined size of all simultaneously-live tensors
///   (intermediates AND resident input leaves, the all-co-resident model),
///   and it chooses the evaluation order that minimizes that peak. Unlike the
///   order-independent `DenseFLOPs`/`DenseSize`, the contraction order is a
///   real lever here. NOTE: `DensePeakSize` does not yet support
///   common-subexpression elimination (`CSEOptions::subnet` must be false).
/// - `DensePeakSizeBatched` extends `DensePeakSize` with a per-index
///   batchability model: each index satisfying
///   `OptimizeOptions::batch_policy.is_batchable_index` is treated as
///   independently sliced to
///   `min(extent, batch_policy.batch_target_size(ix))` elements per index. The
///   DP minimises peak over the worst-case sliced configuration. Only consulted
///   by the batched oracle and DP; requires
///   `batch_policy.is_batchable_index` and `batch_policy.batch_target_size`
///   to be set.
///
/// Leaves room for `Sparse*` models later.
enum class ObjectiveFunction {
  DenseFLOPs,
  DenseSize,
  DensePeakSize,
  DensePeakSizeBatched
};

/// Whether to reorder summands so terms with shared intermediates appear
/// closer to each other.
enum class ReorderSum { Reorder, NoReorder };

/// Common-subexpression-elimination (CSE) options for single-term
/// optimization. `subnet` recognizes equivalent subnetworks while searching for
/// an evaluation order, trading extra search time for potentially lower op
/// counts. (Room to grow, e.g. a maximum subnet size to consider.)
struct CSEOptions {
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
  /// Machine balance beta = 8*F/B in FLOPs per element of traffic. 0 = off.
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
  /// Marks a LEAF tensor as volatile (amplitude-dependent), so the contraction
  /// forming any subset that contains it is replayed every iteration. Empty =>
  /// nothing volatile (replay weighting off).
  std::function<bool(Tensor const&)> is_volatile_leaf = {};
  /// Replay weight on volatile contractions (conceptually the replay count).
  double volatile_weight = 1.0;
  /// Per-intermediate storage-footprint penalty (DenseFLOPs/DenseSize only; see
  /// OptimizeOptions::footprint_weight). Not used by the peak objectives.
  double footprint_weight = 0.0;
  /// Relative peak tolerance for the peak objectives' final selection; see
  /// OptimizeOptions::peak_flops_tolerance.
  double peak_flops_tolerance = 0.10;
  /// Roofline parameters for the peak objectives' secondary cost; see
  /// \ref RooflineParams. machine_balance == 0 => pure-flop tie-break.
  RooflineParams roofline = {};
  /// In-flight batch-contribution footprint multiplier for
  /// DensePeakSizeBatched; see BatchPolicy::accumulation_factor. 0 (default) =
  /// no penalty.
  double accumulation_factor = 0.0;
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
  /// \c Sum_pairs d^k. If empty, composites fall back to \ref idx_to_extent
  /// (k=1), which grossly under-sizes multi-composite tensors.
  std::function<double(Index const&, std::size_t)> inner_pow = {};

  /// Batchability policy: bundles the three per-index and per-leaf predicates
  /// that govern batched evaluation. All three fields default to empty (no
  /// batchable indices, no volatile leaves). The three sub-fields are:
  ///   - `is_batchable_index`: marks an Index as living in a batchable space
  ///     (e.g. DF/RI aux; = the eval cache's accept_aux).
  ///   - `batch_target_size`: per-index slice size; a sliced batchable index
  ///     ix contributes min(extent, batch_target_size(ix)). Only consulted by
  ///     DensePeakSizeBatched.
  ///   - `is_volatile_leaf`: marks a LEAF tensor as volatile (its value
  ///     changes between replays). Empty => no tensor is volatile => cost
  ///     weighting is disabled and volatile_weight is ignored. CC callers pass
  ///     label==volatile_label, the same classification the runtime eval cache
  ///     uses, so the optimizer's cost model and the cache agree.
  BatchPolicy batch_policy = {};

  /// Real-valued weight on the cost of each volatile contraction (re-evaluated
  /// on every replay of the network), while persistent (volatile-independent)
  /// contractions are counted once. Conceptually the expected number of
  /// replays. Default 1.0 (no change). Only consulted when
  /// batch_policy.is_volatile_leaf is non-empty and objective_function ==
  /// ObjectiveFunction::DenseFLOPs.
  double volatile_weight = 1.0;

  /// Relative peak tolerance for the peak objectives' final selection: among
  /// the Pareto frontier of (peak, flops) trade-offs, pick the fewest-flops
  /// schedule whose peak is within (1 + peak_flops_tolerance) of the minimum.
  /// 0 = strict peak-min (flop tie-break only on exact peak ties). The default
  /// 0.10 trades up to a 10% peak increase for a (often much larger) flop
  /// reduction -- e.g. forming a persistent 4-PNO integral instead of
  /// recomputing a particle-ladder. Only consulted by DensePeakSize /
  /// DensePeakSizeBatched.
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
  /// tie-break (no behavior change). Consulted only by DensePeakSize /
  /// DensePeakSizeBatched.
  RooflineParams roofline = {};
};

}  // namespace sequant

#endif  // SEQUANT_CORE_OPTIMIZE_OPTIONS_HPP
