#ifndef SEQUANT_CORE_OPTIMIZE_OPTIONS_HPP
#define SEQUANT_CORE_OPTIMIZE_OPTIONS_HPP

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
///   `OptimizeOptions::is_batchable_index` is treated as independently sliced
///   to `min(extent, batch_target_size)` elements. The DP minimises peak over
///   the worst-case sliced configuration. Only consulted by the batched oracle
///   and DP; requires `is_batchable_index` and `batch_target_size` to be set.
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

  /// Marks a LEAF tensor as volatile: its value changes between replays of the
  /// network, so any contraction depending on it is re-evaluated on every
  /// replay. Empty (default) ⇒ no tensor is volatile ⇒ cost weighting is
  /// disabled and volatile_weight is ignored (behavior identical to before this
  /// feature). CC callers pass label==volatile_label, the same classification
  /// the runtime eval cache uses, so the optimizer's cost model and the cache
  /// agree.
  std::function<bool(Tensor const&)> is_volatile_leaf = {};

  /// Real-valued weight on the cost of each volatile contraction (re-evaluated
  /// on every replay of the network), while persistent (volatile-independent)
  /// contractions are counted once. Conceptually the expected number of
  /// replays. Default 1.0 (no change). Only consulted when is_volatile_leaf is
  /// non-empty and objective_function == ObjectiveFunction::DenseFLOPs.
  double volatile_weight = 1.0;

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

  /// Predicate marking an Index as living in a batchable space the runtime
  /// slices over (e.g. DF/RI aux; = the eval cache's accept_aux). Each distinct
  /// batchable index is sliced independently. Only consulted by
  /// ObjectiveFunction::DensePeakSizeBatched.
  std::function<bool(Index const&)> is_batchable_index = {};

  /// Per-index slice size: a sliced batchable index contributes
  /// min(extent, batch_target_size(ix)). Empty (default nullptr/empty function)
  /// disables the batched discount. Only consulted by DensePeakSizeBatched.
  std::function<std::size_t(Index const&)> batch_target_size = {};
};

}  // namespace sequant

#endif  // SEQUANT_CORE_OPTIMIZE_OPTIONS_HPP
