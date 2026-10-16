#ifndef SEQUANT_CORE_EXPORT_MEMORY_MODEL_HPP
#define SEQUANT_CORE_EXPORT_MEMORY_MODEL_HPP

namespace sequant {

/// Selects how ScheduleWalkingGenerationVisitor (schedule_export.hpp)
/// realizes a cell's live span as create/load/unload/destroy calls; also
/// doubles, via natural_memory_model.hpp, as the per-cell answer to "which
/// of these two models does this cell's own live span satisfy without any
/// rewriting."
///
/// As the visitor's operating mode, this is a fixed property of the target
/// backend (Generator::memory_model(), generator.hpp) -- queried from the
/// Generator being driven, not chosen by whoever constructs the visitor, so
/// that a backend without stack discipline requirements (MemoryModel::
/// RandomAccess) pays nothing beyond returning that constant.
enum class MemoryModel {
  /// The schedule's own life/order is passed through unmodified: a cell is
  /// created once, at its production, stays resident across all its reads,
  /// and is destroyed on the read that drains its life. Live spans may be
  /// non-nested (a shared value's span can outlast values created after it).
  /// This is what a backend without stack-like memory requirements uses.
  ///
  /// Every cell trivially satisfies this as its natural model, since it
  /// imposes no nesting requirement at all.
  RandomAccess,

  /// Every cell whose natural model (natural_memory_model.hpp) falls short
  /// of Stack -- because it is genuinely shared, or because its one read is
  /// separated from its production by other, unrelated allocation/
  /// deallocation activity -- is force-persisted: stored immediately after
  /// its own computation, then reloaded independently at each later read
  /// (unloaded after a non-draining read, destroyed after the draining
  /// one). Every cell whose natural model is already Stack nests normally,
  /// exactly as under RandomAccess. The result is that every live span
  /// nests properly (LIFO), matching a backend whose target language only
  /// supports stack-like memory allocation (e.g. ITF).
  Stack,
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_MEMORY_MODEL_HPP
