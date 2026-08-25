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

/// \brief A DAG-scope loop's nest position, the value a \c ModeToLevel map
///        keys by mode to answer "which loop (if any) does mode m of this
///        node's result run under".
///
/// \details Mirrors \c ScopeBlock{axis.space, ordinal}: \c depth is the
/// nest position (outer loops have smaller \c depth), \c space is the
/// base_key of the space the loop iterates over (kept for diagnostics only,
/// not part of loop identity beyond what \c ordinal disambiguates), and
/// \c ordinal disambiguates sibling loops over the same space at the same
/// depth.
struct DagScopeLevel {
  std::size_t depth;
  std::wstring space;
  int ordinal;

  friend bool operator==(DagScopeLevel const& lhs, DagScopeLevel const& rhs) {
    return lhs.depth == rhs.depth && lhs.space == rhs.space &&
           lhs.ordinal == rhs.ordinal;
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

  /// \return the mode (position in \c by_mode) that runs under \p level, or
  ///         nullopt if no mode does.
  [[nodiscard]] std::optional<std::size_t> mode_of(
      DagScopeLevel const& level) const {
    for (std::size_t i = 0; i < by_mode.size(); ++i)
      if (by_mode[i] && *by_mode[i] == level) return i;
    return std::nullopt;
  }
};

/// \brief The runtime consumption view of the loop-colored canonical layout
///        (sliced-value canonical-layout / loop-coloring design, sec.4): the
///        per-VALUE datum \c slice_to_use reads to resolve WHICH physical mode
///        of a fetched value a batch loop slices.
///
/// \details This is the hash-keyed projection of \c SlicedModeAssignment (which
/// is value_id-keyed and lives in ordered_schedule.hpp): for each value's
/// \c EvalExpr::hash_value(), the list of (this value's OWN sliced-mode
/// \c Index, the \c LoopId that slices it) pairs the schedule assigned. The
/// executor holds it as a non-owning cache seam (\c CacheManager, mirroring
/// the \c ModeToLevel seam it supersedes as a slice-mode source) and consults
/// it per fetch: given the enclosing loop's \c DagScopeLevel, \c loop_of_level
/// yields the loop's id and \c mode_of yields the fetched value's own sliced
/// \c Index for it, which the executor maps to a physical slot on the fetched
/// node (space-mapped when a CSE-folded occurrence binds the space under a
/// relabeled physical index). An unwired seam (null) or a value with no entry
/// leaves the fetch UNSLICED -- exactly the participation the built-within gate
/// establishes.
struct LoopColoredSliceSeam {
  /// Every DISTINCT realized loop, in canonical order (copied verbatim from
  /// \c SlicedModeAssignment::levels): a \c LoopId indexes this list.
  container::vector<DagScopeLevel> levels;

  /// value hash -> (this value's own sliced-mode PHYSICAL POSITION, slicing
  /// LoopId) pairs. The position is the mode's index in the value's own carried
  /// (canon_indices) order -- computed at schedule time IN THE VALUE'S OWN
  /// index-frame, so it is used at runtime directly, never by re-matching a
  /// label against the fetched node (which is a different index-frame). This is
  /// the per-value (first-occurrence / non-divergent) map; a value whose
  /// occurrences physically differ (divergent CSE) is resolved per-occurrence
  /// via \c by_hash_consumer below.
  std::unordered_map<std::size_t,
                     container::svector<std::pair<std::size_t, LoopId>>>
      by_hash;

  /// CONSUMER-DISAMBIGUATED sliced-mode facts (sliced-value canonical-layout /
  /// loop-coloring design, PILLAR 2): value hash -> that value's
  /// (own sliced-mode Index, slicing LoopId, CONSUMER node hash) triples. The
  /// consumer hash is the hash of the eval-node whose fetch of this value binds
  /// \p loop to \p Index -- i.e. the member-root/use-site contraction the
  /// ordered executor is currently evaluating (\c
  /// CacheManager::current_consumer bracketing each \c evaluate_impl). This is
  /// the datum that tells apart the w8-symmetric case: one value's TWO free occ
  /// modes both stamped under the ONE occ loop, each bound by a DIFFERENT
  /// consumer (design sec.2, "pos0 here, pos1 there"). It is consulted by \c
  /// mode_of ONLY when \c by_hash records more than one mode under the queried
  /// loop; a value with a single mode per loop (aux, and every non-symmetric
  /// value) never reaches this map, keeping those fetches BYTE-IDENTICAL to the
  /// consumer-blind resolution.
  /// value hash -> (this occurrence's own sliced-mode PHYSICAL POSITION,
  /// slicing LoopId, CONSUMER node hash) triples. The position is computed at
  /// schedule time in THAT occurrence's own index-frame (index of the mode in
  /// occ.carried), so the runtime uses it directly for the fetch attributed to
  /// that consumer -- no cross-frame label match. This is the frame-correct
  /// resolution for a divergent (relabeled) CSE occurrence and for the
  /// symmetric case (one value's two free occ modes bound by two different
  /// consumers, at two different positions).
  std::unordered_map<
      std::size_t,
      container::svector<std::tuple<std::size_t, LoopId, std::size_t>>>
      by_hash_consumer;

  /// \return the \c LoopId whose realized \c DagScopeLevel equals \p level, or
  ///         nullopt if no realized loop matches (\p level is the root scope or
  ///         a level this seam does not enumerate).
  [[nodiscard]] std::optional<LoopId> loop_of_level(
      DagScopeLevel const& level) const {
    for (std::size_t i = 0; i < levels.size(); ++i)
      if (levels[i] == level) return i;
    return std::nullopt;
  }

  /// \return \p hash's own sliced-mode \c Index sliced by \p loop,
  /// disambiguated
  ///         by the fetching \p consumer, or nullopt if \p hash is unknown to
  ///         this seam or is not sliced by \p loop.
  ///
  /// \details When exactly ONE of \p hash's modes is sliced by \p loop (the
  /// common single-mode-per-loop case -- aux and every non-symmetric value),
  /// that mode is returned and \p consumer is IGNORED, so those fetches stay
  /// byte-identical to the old consumer-blind first-match. When MORE THAN ONE
  /// mode is sliced by \p loop (the w8-symmetric ambiguity: a value's two free
  /// occ modes both stamped under the one occ loop, each bound by a different
  /// use-site), \p consumer selects THIS use-site's own mode from
  /// \c by_hash_consumer. If \p consumer is null or has no recorded fact (a
  /// fetch outside any tracked consumer), the resolution falls back to the same
  /// deterministic first-match as before -- never a hard failure.
  [[nodiscard]] std::optional<std::size_t> mode_of(
      std::size_t hash, LoopId loop,
      std::optional<std::size_t> consumer = std::nullopt) const {
    // A consumer-attributed PER-OCCURRENCE fact is the frame-correct answer and
    // takes precedence: it is the sliced mode's physical POSITION in THIS
    // fetch's own occurrence frame, so it is correct for a divergent
    // (relabeled) occurrence and disambiguates the symmetric case (two
    // consumers binding one loop to two different positions).
    if (consumer) {
      auto const cit = by_hash_consumer.find(hash);
      if (cit != by_hash_consumer.end())
        for (auto const& [pos, lid, ch] : cit->second)
          if (lid == loop && ch == *consumer) return pos;
    }
    // Otherwise the per-value position (valid when the value's occurrences do
    // not physically diverge -- one position per loop). If MORE than one
    // position is recorded under this loop the value is ambiguous WITHOUT a
    // matching consumer fact: do NOT guess a position -- leave the fetch
    // unsliced (nullopt). A genuinely sliced fetch reaching this is a scheduler
    // gap that must be fixed by recording the occurrence fact, never papered
    // over with a first-match.
    auto const it = by_hash.find(hash);
    if (it == by_hash.end()) return std::nullopt;
    std::optional<std::size_t> first;
    std::size_t nmatch = 0;
    for (auto const& [pos, lid] : it->second)
      if (lid == loop) {
        if (nmatch == 0) first = pos;
        ++nmatch;
      }
    if (nmatch == 1) return first;
    return std::nullopt;  // 0 -> unsliced; >1 -> ambiguous, never guessed
  }

  /// \return whether \p hash has ANY recorded sliced-mode fact under \p loop --
  /// in either the per-occurrence (\c by_hash_consumer, for any consumer) or
  /// the per-value (\c by_hash) map. A value that participates in a loop (some
  /// occurrence of it is sliced there) but for which \c mode_of returns nullopt
  /// at a particular fetch is a SCHEDULER GAP, not an invariant value: the
  /// completeness guard in \c slice_to_use uses this to tell the two apart and
  /// fail loud rather than silently leave an operand unsliced. A value with NO
  /// fact under \p loop is genuinely invariant to it (correctly unsliced).
  [[nodiscard]] bool participates(std::size_t hash, LoopId loop) const {
    if (auto const cit = by_hash_consumer.find(hash);
        cit != by_hash_consumer.end())
      for (auto const& [pos, lid, ch] : cit->second)
        if (lid == loop) return true;
    if (auto const it = by_hash.find(hash); it != by_hash.end())
      for (auto const& [pos, lid] : it->second)
        if (lid == loop) return true;
    return false;
  }
};

}  // namespace sequant

#endif  // SEQUANT_EVAL_DAG_SCOPE_HPP
