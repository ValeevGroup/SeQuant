#ifndef SEQUANT_CORE_OPTIMIZE_OPTIONS_HPP
#define SEQUANT_CORE_OPTIMIZE_OPTIONS_HPP

#include <cstddef>
#include <functional>

namespace sequant {

class Index;
class Tensor;

/// Cost metric to optimize for in single-term and top-level optimize routines.
enum class OptFor { Flops, Memsize };

/// Whether to reorder summands so terms with shared intermediates appear
/// closer to each other.
enum class ReorderSum { Reorder, NoReorder };

/// Whether single-term optimization should recognize equivalent subnetworks
/// while searching for an evaluation order, trading extra search time for
/// potentially lower op counts.
enum class SubnetCSE { Enable, Disable };

/// A type-erased provider mapping an Index to its extent. Used by the public
/// optimize() API. Callers reaching for the templated opt::single_term_opt
/// overloads (constrained by \ref opt::has_index_extent) should pass the
/// callable directly instead of wrapping it here — that keeps the cost
/// function's call site inlineable, whereas a value of this alias goes
/// through std::function's type-erased dispatch on every Index lookup.
using index_to_extent_t = std::function<std::size_t(Index const&)>;

/// Options that control behavior of \ref sequant::optimize.
struct OptimizeOptions {
  /// Cost metric to minimize.
  OptFor opt_for = OptFor::Flops;

  /// Whether to reorder summands so terms with shared intermediates appear
  /// closer to each other.
  ReorderSum reorder = ReorderSum::Reorder;

  /// Whether single-term optimization should perform subnetwork
  /// common-subexpression recognition. Disabled by default; enabling can
  /// reduce op counts at the cost of additional optimization time.
  SubnetCSE subnet_cse = SubnetCSE::Disable;

  /// Caller-supplied Index to extent provider. If empty, defaults to
  /// \c IndexSpace::approximate_size().
  index_to_extent_t idx_to_extent = {};

  /// Marks a LEAF tensor as volatile: its value changes between replays of the
  /// network, so any contraction depending on it is re-evaluated on every
  /// replay. Empty (default) ⇒ no tensor is volatile ⇒ cost weighting is
  /// disabled and n_replay is ignored (behavior identical to before this
  /// feature). CC callers pass label==volatile_label, the same classification
  /// the runtime eval cache uses, so the optimizer's cost model and the cache
  /// agree.
  std::function<bool(Tensor const&)> is_volatile_leaf = {};

  /// Expected number of times the network is replayed with the volatile inputs
  /// mutated; the cost of each volatile contraction is multiplied by this,
  /// while persistent (volatile-independent) contractions are counted once.
  /// Default 1 (no change). Only consulted when is_volatile_leaf is non-empty
  /// and opt_for == Flops.
  unsigned n_replay = 1;

  /// Per-intermediate memory-footprint penalty added to the single-term
  /// optimization cost. For every binary contraction, the storage footprint of
  /// its RESULT intermediate (the product of the extents of the result's
  /// indices, i.e. its element count) is multiplied by this weight and added to
  /// the contraction cost. Unlike the FLOPs cost, this penalty is NOT scaled by
  /// \ref n_replay (peak footprint is a one-time materialization cost, not a
  /// per-replay one). 0 (default) disables the penalty, recovering the pure
  /// FLOPs/Memsize behavior.
  ///
  /// Rationale: the FLOPs cost is blind to the storage size of the
  /// intermediates it materializes, so it will happily pick an order (and thus
  /// expose, as a shareable subtree, a common subexpression) that carries a
  /// free large-space index — e.g. a half-transformed DF integral that still
  /// carries a free projected-AO (PAO) index — because forming it once is
  /// FLOPs-cheap. Such an intermediate can be enormous. A nonzero
  /// footprint_weight biases single-term optimization toward orders that defer
  /// or avoid materializing such large intermediates (e.g. transforming both
  /// large legs before exposing a shared subtree), trading a controlled amount
  /// of extra FLOPs for a lower peak footprint. Only consulted when
  /// opt_for == Flops; the units are FLOPs-per-element, so a useful magnitude
  /// is on the order of the contracted-index extent that the offending
  /// intermediate would otherwise leave free.
  double footprint_weight = 0.0;
};

}  // namespace sequant

#endif  // SEQUANT_CORE_OPTIMIZE_OPTIONS_HPP
