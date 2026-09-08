#ifndef SEQUANT_EVAL_DAG_SCOPE_HPP
#define SEQUANT_EVAL_DAG_SCOPE_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/index.hpp>

#include <cstddef>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace sequant {

/// \brief A realized DAG-scope loop's stable id: its position in the canonical
///        enumeration of every loop realized by the schedule (see
///        \c SlicedModeAssignment::levels, ordered_schedule.hpp). Defined here
///        (rather than in ordered_schedule.hpp) so the low-level cache seam
///        below -- consumed by \c CacheManager -- can name it without pulling
///        in the schedule machinery.
using LoopId = std::size_t;

/// \brief The STABLE identity of a batch loop: which loop-GROUP (\c depth) and
///        which member-SLOT within that group (\c loop_slot).
///
/// \details The cell table, the per-occurrence mode<->loop atlas, and
/// \c ordered_n_batches_by_loop (see ordered_executor.hpp) all key on this
/// identity; the layout (\c altitude_ordinal /
/// \c latitude_ordinal on \c DagScopeLevel) never enters it. \c depth
/// distinguishes even two groups of the SAME space (an "external" and a
/// "contracted" group of one space); \c loop_slot distinguishes the members of
/// one group. See doc/dev/specs/2026-08-28-batched-dag-loop-identity-design.md.
struct LoopKey {
  std::size_t depth;  //!< which loop-group
  int loop_slot;      //!< which member-slot within the group (0-based)

  /// \brief A single opaque color encoding the FULL loop identity (\c depth AND
  /// \c loop_slot), for use where one \c std::size_t must distinguish loops:
  /// \c ordered_n_batches_by_loop (one batch count PER LOOP, not per
  /// loop-group -- two members of one group are DISTINCT loops with distinct
  /// batch counts; see ordered_executor.hpp). Keying on \c depth alone would
  /// conflate same-group sibling loops. (\c loop_slot < 4096 in every
  /// realized schedule.)
  [[nodiscard]] std::size_t color() const {
    return (depth << 12) | static_cast<std::size_t>(loop_slot);
  }

  friend bool operator==(LoopKey const& a, LoopKey const& b) {
    return a.depth == b.depth && a.loop_slot == b.loop_slot;
  }
  friend bool operator!=(LoopKey const& a, LoopKey const& b) {
    return !(a == b);
  }
};

/// \brief A batch loop's realized placement: its identity (\c depth, \c
///        loop_slot) plus its LAYOUT (\c altitude_ordinal, \c
///        latitude_ordinal).
///
/// \details Identity coordinates (\c depth, \c loop_slot) name the loop; layout
/// coordinates say where a schedule placed it. \c altitude_ordinal is the
/// nesting rank the schedule assigned the slot within its group
/// (free/interchangeable);
/// \c latitude_ordinal is the PASS index within a forced-split nest (formerly
/// \c ordinal): a nest holding members of more than one pass emits one
/// sibling block per pass, and \c latitude_ordinal disambiguates them (see
/// \c forced_split_levels, ordered_schedule.hpp). \c space is a color for
/// fusion group-matching, never identity. See
/// doc/dev/specs/2026-08-28-batched-dag-loop-identity-design.md.
struct DagScopeLevel {
  std::size_t depth;   //!< which loop-group (identity)
  std::wstring space;  //!< color only, NOT identity
  int loop_slot = 0;   //!< which member-slot within the group (identity)
  int altitude_ordinal =
      0;  //!< layout: nesting rank of the slot within its group
  int latitude_ordinal = 0;  //!< layout: pass index (was: ordinal)

  [[nodiscard]] LoopKey key() const { return LoopKey{depth, loop_slot}; }

  friend bool operator==(DagScopeLevel const& lhs, DagScopeLevel const& rhs) {
    // Full-tuple comparison; identity migrates to key() alone in a later
    // task.
    return lhs.depth == rhs.depth && lhs.space == rhs.space &&
           lhs.loop_slot == rhs.loop_slot &&
           lhs.altitude_ordinal == rhs.altitude_ordinal &&
           lhs.latitude_ordinal == rhs.latitude_ordinal;
  }

  friend bool operator!=(DagScopeLevel const& lhs, DagScopeLevel const& rhs) {
    return !(lhs == rhs);
  }
};

/// \brief A node's mode (canon_indices position) -> DagScopeLevel map: for
///        each mode of a node's result, the DAG-scope loop it runs under (if
///        any).
struct ModeToLevel {
  /// indexed by canon_indices position; nullopt where that mode does not run
  /// under any DAG-scope loop.
  container::svector<std::optional<DagScopeLevel>> by_mode;
};

}  // namespace sequant

#endif  // SEQUANT_EVAL_DAG_SCOPE_HPP
