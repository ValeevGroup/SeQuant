#ifndef SEQUANT_CORE_BATCH_POLICY_HPP
#define SEQUANT_CORE_BATCH_POLICY_HPP

#include <cstddef>
#include <functional>

namespace sequant {

class Index;
class Tensor;

/// One batchability policy shared by the single-term optimizer and the runtime
/// batched evaluator (make_evaluator, Task A3). All predicates default empty.
struct BatchPolicy {
  std::function<bool(Index const&)> is_batchable_index = {};
  std::function<std::size_t(Index const&)> batch_target_size = {};
  std::function<bool(Tensor const&)> is_volatile_leaf = {};

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
};

}  // namespace sequant

#endif  // SEQUANT_CORE_BATCH_POLICY_HPP
