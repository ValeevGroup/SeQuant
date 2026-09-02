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
/// \details The seam, value-id / occurrence-id coloring, and the per-occurrence
/// mode<->loop atlas all key on this identity; the layout (\c altitude_ordinal
/// /
/// \c latitude_ordinal on \c DagScopeLevel) never enters it. \c depth
/// distinguishes even two groups of the SAME space (an "external" and a
/// "contracted" group of one space); \c loop_slot distinguishes the members of
/// one group. See doc/dev/specs/2026-08-28-batched-dag-loop-identity-design.md.
struct LoopKey {
  std::size_t depth;  //!< which loop-group
  int loop_slot;      //!< which member-slot within the group (0-based)

  /// \brief A single opaque color encoding the FULL loop identity (\c depth AND
  /// \c loop_slot), for use where one \c std::size_t must distinguish loops:
  /// the value-id coloring (one cache PER LOOP, not per loop-group -- two
  /// members of one group are DISTINCT loops with distinct caches) and the
  /// home-scope filter that matches a home-sliced mode's loop against the
  /// enclosing scope. Keying on \c depth alone would conflate same-group
  /// sibling loops. (\c loop_slot < 4096 in every realized schedule.)
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
/// \c latitude_ordinal is the legality producer/consumer (PROCON) pass index
/// (formerly \c ordinal). \c space is a color for fusion group-matching, never
/// identity. See doc/dev/specs/2026-08-28-batched-dag-loop-identity-design.md.
struct DagScopeLevel {
  std::size_t depth;   //!< which loop-group (identity)
  std::wstring space;  //!< color only, NOT identity
  int loop_slot = 0;   //!< which member-slot within the group (identity)
  int altitude_ordinal =
      0;  //!< layout: nesting rank of the slot within its group
  int latitude_ordinal = 0;  //!< layout: PROCON pass index (was: ordinal)

  [[nodiscard]] LoopKey key() const { return LoopKey{depth, loop_slot}; }

  friend bool operator==(DagScopeLevel const& lhs, DagScopeLevel const& rhs) {
    // Task 1 (vocabulary landing): full-tuple comparison keeps behavior
    // byte-identical (loop_slot == altitude_ordinal == 0 until Task 3; passes
    // stay distinct via latitude_ordinal until Task 6). Identity migrates to
    // key() in later tasks.
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

  /// value hash -> loop colors (LoopKey::color) of the loops the value is
  /// PRODUCED-SLICED on (its residency home coloring, home_mode_depth). A
  /// fetch inside such a loop must NOT slice again (the stored form already
  /// is the batch); a fetch inside any OTHER enclosing loop the value carries
  /// a mode of is a use-induced slice decided by mode_of(). Empty when the
  /// caller did not publish residency (legacy hops-prefix slicing applies).
  std::unordered_map<std::size_t, container::svector<std::size_t>> home_colors;
  [[nodiscard]] bool produced_sliced_on(std::size_t hash,
                                        std::size_t color) const {
    auto const it = home_colors.find(hash);
    if (it == home_colors.end()) return false;
    return std::find(it->second.begin(), it->second.end(), color) !=
           it->second.end();
  }

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

  /// consumer node hash -> the loops that consumer is PRODUCED sliced on (it is
  /// building a batch of each). The use-induced-slicing oracle: a fetched
  /// operand carrying loop L's mode is sliced iff the CURRENT consumer is here
  /// under L -- i.e. building an L-batch. Drives both the slice decision for a
  /// whole-produced operand and the completeness guard (a consumer NOT sliced
  /// on L reads its operands whole there, never a gap). Non-transitive: keyed
  /// on the consumer's OWN production, never a use-induced slice it acquired
  /// for a different consumer.
  std::unordered_map<std::size_t, container::svector<LoopId>>
      consumer_sliced_loops;

  /// EXPLICIT per-occurrence INVARIANT facts: value hash -> (LoopId, CONSUMER
  /// hash) pairs recording that this consumer's fetch of the value is
  /// correctly UNSLICED on that loop (the value's occurrence in this consumer's
  /// frame does not carry the loop's mode). Consulted by the completeness
  /// guard: a fetch with no slice fact is a GAP only if it also has no
  /// invariant fact -- a CSE-shared value that carries the mode in one frame
  /// (sliced there) and not in another (invariant there) must not trip the
  /// guard on the latter.
  std::unordered_map<std::size_t,
                     container::svector<std::pair<LoopId, std::size_t>>>
      by_hash_consumer_invariant;

  /// \return whether \p hash's fetch by \p consumer is recorded as correctly
  /// UNSLICED on \p loop (an explicit invariant decision, not a gap).
  [[nodiscard]] bool invariant_for(std::size_t hash, LoopId loop,
                                   std::size_t consumer) const {
    auto const it = by_hash_consumer_invariant.find(hash);
    if (it == by_hash_consumer_invariant.end()) return false;
    for (auto const& [lid, ch] : it->second)
      if (lid == loop && ch == consumer) return true;
    return false;
  }

  /// \return whether \p consumer is produced sliced on \p loop (building an
  /// \p loop-batch), so its operands carrying that loop's mode must be sliced.
  [[nodiscard]] bool consumer_slices(std::size_t consumer, LoopId loop) const {
    auto const it = consumer_sliced_loops.find(consumer);
    if (it == consumer_sliced_loops.end()) return false;
    return std::find(it->second.begin(), it->second.end(), loop) !=
           it->second.end();
  }

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
