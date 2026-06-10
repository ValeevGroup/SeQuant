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
};

}  // namespace sequant

#endif  // SEQUANT_CORE_OPTIMIZE_OPTIONS_HPP
