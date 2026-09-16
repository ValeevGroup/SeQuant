#ifndef SEQUANT_EVAL_DAG_SCOPE_HPP
#define SEQUANT_EVAL_DAG_SCOPE_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/index.hpp>

#include <cstddef>
#include <functional>
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

/// \brief The stable identity of a batch loop: which loop-group (\c depth) and
///        which member-slot within that group (\c loop_slot).
///
/// \details The cell table, the per-occurrence mode<->loop atlas, and
/// \c ordered_n_batches_by_loop (see ordered_executor.hpp) all key on this
/// identity; the layout (\c altitude_ordinal /
/// \c latitude_ordinal on \c DagScopeLevel) never enters it. \c depth
/// distinguishes even two groups of the same space (an "external" and a
/// "contracted" group of one space); \c loop_slot distinguishes the members of
/// one group. See doc/dev/specs/2026-09-12-batched-array-dag-eval-as-built.md,
/// section 5.1.
struct LoopKey {
  std::size_t depth;  //!< which loop-group
  int loop_slot;      //!< which member-slot within the group (0-based)

  friend bool operator==(LoopKey const& a, LoopKey const& b) {
    return a.depth == b.depth && a.loop_slot == b.loop_slot;
  }
  friend bool operator!=(LoopKey const& a, LoopKey const& b) {
    return !(a == b);
  }
};

/// \brief A batch loop's realized placement: its identity (\c depth, \c
///        loop_slot) plus its layout (\c altitude_ordinal, \c
///        latitude_ordinal).
///
/// \details Identity coordinates (\c depth, \c loop_slot) name the loop; layout
/// coordinates say where a schedule placed it. \c altitude_ordinal is the
/// nesting rank the schedule assigned the slot within its group
/// (free/interchangeable);
/// \c latitude_ordinal is the pass index within a forced-split nest: a nest
/// holding members of more than one pass emits one
/// sibling block per pass, and \c latitude_ordinal disambiguates them (see
/// \c forced_split_levels, ordered_schedule.hpp). \c space is a color for
/// fusion group-matching, never identity. See
/// doc/dev/specs/2026-09-12-batched-array-dag-eval-as-built.md, section 5.1.
struct DagScopeLevel {
  std::size_t depth;   //!< which loop-group (identity)
  std::wstring space;  //!< color only, not identity
  int loop_slot = 0;   //!< which member-slot within the group (identity)
  int altitude_ordinal =
      0;  //!< layout: nesting rank of the slot within its group
  int latitude_ordinal = 0;  //!< layout: pass index

  [[nodiscard]] LoopKey key() const { return LoopKey{depth, loop_slot}; }

  friend bool operator==(DagScopeLevel const& lhs, DagScopeLevel const& rhs) {
    // Full-tuple comparison; \c key() alone is the loop's identity.
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

/// \brief Hash of a \c LoopKey, so the full loop identity (\c depth and
///        \c loop_slot) can key an \c unordered_map directly.
///
/// \details Anything that counts or looks up something per loop (e.g.
/// \c ordered_n_batches_by_loop, one batch count per loop -- two members of
/// one loop-group are distinct loops with distinct batch counts, so keying on
/// \c depth alone would conflate them) keys on the pair itself through this
/// hash plus \c LoopKey::operator==. There is deliberately no packed
/// single-\c size_t "color": packing \c loop_slot into a fixed bit field
/// silently aliases two distinct loops once the slot numbering (an unbounded
/// per-space counter in peak_profile.hpp) exceeds the field, and a
/// hash+equality pair has no such bound.
template <>
struct std::hash<sequant::LoopKey> {
  std::size_t operator()(sequant::LoopKey const& k) const noexcept {
    std::size_t h = std::hash<std::size_t>{}(k.depth);
    // boost-style combine; the two fields are small and would otherwise
    // collide trivially under a plain xor.
    h ^= std::hash<int>{}(k.loop_slot) + 0x9e3779b97f4a7c15ULL + (h << 6) +
         (h >> 2);
    return h;
  }
};

#endif  // SEQUANT_EVAL_DAG_SCOPE_HPP
