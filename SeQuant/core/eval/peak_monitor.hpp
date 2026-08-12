#ifndef SEQUANT_CORE_EVAL_PEAK_MONITOR_HPP
#define SEQUANT_CORE_EVAL_PEAK_MONITOR_HPP

#include <cstddef>
#include <functional>
#include <vector>

namespace sequant::eval {

/// \brief Location of a (candidate) high-water event: the co-resident byte
///        count observed and the hash of the op that was being evaluated when
///        it was observed (0 if no op node was in scope at the call site).
struct PeakEvent {
  std::size_t bytes = 0;
  std::size_t op_hash = 0;
};

/// \brief DIAGNOSTIC (analysis-only): one alive cache entry snapshotted at a
///        high-water advance -- the value's CSE hash and its modeled byte size.
struct PeakLiveEntry {
  std::size_t hash = 0;
  std::size_t bytes = 0;
};

/// \brief Hierarchy-wide co-resident memory high-water tracker.
///
/// A single \c PeakMonitor can be wired (via \c CacheManager::set_peak_monitor)
/// onto the root of a scope-chain of \c CacheManager instances so every
/// \c note_working_set() call anywhere in the chain -- root cache, per-batch
/// scratch, nested scopes -- folds into ONE running high-water mark, with a
/// hook fired each time that mark advances (so a caller can capture *where*
/// the current peak was observed, e.g. for a schedule visualizer).
struct PeakMonitor {
  /// Running high-water mark (bytes), monotonically non-decreasing.
  std::size_t hwmark_bytes = 0;

  /// Location of the current high-water (the event that last advanced
  /// \c hwmark_bytes).
  PeakEvent peak{};

  /// Optional callback fired each time \c hwmark_bytes advances (i.e. exactly
  /// when \c peak is updated). Empty (default) => no-op.
  std::function<void(PeakEvent const&)> on_peak{};

  /// DIAGNOSTIC (analysis-only): optional callback fired -- by
  /// \c CacheManager::note_working_set, BEFORE \c observe advances the mark --
  /// each time a new global high-water is about to be recorded, carrying the
  /// total co-resident bytes and the full alive-entry set (chain-wide) at that
  /// instant. Empty (default) => \c note_working_set never enumerates the live
  /// set, so the fold is byte-identical to the pre-hook path. Used to decompose
  /// a realized peak into its co-resident values.
  std::function<void(std::size_t total_bytes,
                     std::vector<PeakLiveEntry> const&)>
      on_peak_liveset{};

  /// Observe one candidate high-water sample. Advances \c hwmark_bytes (and
  /// fires \c on_peak) only if \p bytes exceeds the current mark; a sample at
  /// or below the mark is a no-op.
  void observe(std::size_t bytes, std::size_t op_hash) noexcept {
    if (bytes > hwmark_bytes) {
      hwmark_bytes = bytes;
      peak = {bytes, op_hash};
      if (on_peak) on_peak(peak);
    }
  }
};

}  // namespace sequant::eval

#endif  // SEQUANT_CORE_EVAL_PEAK_MONITOR_HPP
