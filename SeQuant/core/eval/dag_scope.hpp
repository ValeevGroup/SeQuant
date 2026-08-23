#ifndef SEQUANT_EVAL_DAG_SCOPE_HPP
#define SEQUANT_EVAL_DAG_SCOPE_HPP

#include <SeQuant/core/container.hpp>

#include <cstddef>
#include <optional>
#include <string>

namespace sequant {

/// \brief A DAG-scope loop's nest position, keyed by \c mode_to_level to
///        answer "which loop (if any) does mode m of this node's result run
///        under".
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

}  // namespace sequant

#endif  // SEQUANT_EVAL_DAG_SCOPE_HPP
