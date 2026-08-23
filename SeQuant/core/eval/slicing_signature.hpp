#ifndef SEQUANT_EVAL_SLICING_SIGNATURE_HPP
#define SEQUANT_EVAL_SLICING_SIGNATURE_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/dag_scope.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <cstddef>
#include <optional>
#include <utility>

namespace sequant {

/// \return the position of index \p ix in \p node's canonical result indices
///         (i.e. the corresponding tensor mode), or nullopt if absent.
[[nodiscard]] inline std::optional<std::size_t> index_position(
    meta::eval_node auto const& node, Index const& ix) {
  auto const& idxs = node->canon_indices();
  for (std::size_t p = 0; p < idxs.size(); ++p)
    if (idxs[p] == ix) return p;
  return std::nullopt;
}

/// \brief The slicing signature of \p node over \p modes: the canonical
///        result-slot position of each mode (via \c index_position), or nullopt
///        where \p node does not carry that mode on its result.
///
/// \details The signature is a function of the CANONICAL node alone, so it is
/// identical across canonically-equal occurrences EXCEPT where the modes bind
/// different physical indices -- which is exactly the divergence the batched
/// runtime must not share (see \c make_batched_scratch in eval.hpp, and the
/// design spec doc/dev/specs/2026-08-07-remat-cse-aware-split-design.md). It is
/// the unified form of that path's per-node `sig` (the batch mode's position)
/// plus `ext_sig` (the external modes' positions): pass `{batch_mode}` followed
/// by the external axes as \p modes.
template <meta::eval_node Node>
[[nodiscard]] container::svector<std::optional<std::size_t>> slicing_signature(
    Node const& node, container::svector<Index> const& modes) {
  container::svector<std::optional<std::size_t>> sig;
  sig.reserve(modes.size());
  for (auto const& m : modes) sig.push_back(index_position(node, m));
  return sig;
}

/// \return true iff every occurrence in \p occurrences has the SAME slicing
///         signature over \p modes -- i.e. the value may be materialized once
///         and shared across all of them. False means at least one mode binds a
///         different physical slot across occurrences (a relabeled mode): the
///         value cannot be shared sliced and must be SPLIT (materialized per
///         occurrence).
template <typename Range>
[[nodiscard]] bool signatures_consistent(
    Range const& occurrences, container::svector<Index> const& modes) {
  std::optional<container::svector<std::optional<std::size_t>>> first;
  for (auto const& n : occurrences) {
    auto sig = slicing_signature(n, modes);
    if (!first)
      first = std::move(sig);
    else if (*first != sig)
      return false;
  }
  return true;
}

/// \brief Zips a \c slicing_signature (one entry per loop axis, in nest
///        order) with the matching \p levels (same order) into a
///        \c ModeToLevel sized to the node's \p rank -- the single point
///        turning "positions of loop axes on a node" into "mode->level".
///
/// \param rank the node's result rank (i.e. \c node->canon_indices().size()),
///        used to size \c ModeToLevel::by_mode.
/// \param sig a slicing signature (as returned by \c slicing_signature): one
///        optional mode-position per loop axis.
/// \param levels the DAG-scope loop level of each loop axis in \p sig, same
///        size and order as \p sig.
/// \pre \c sig.size() == levels.size()
[[nodiscard]] inline ModeToLevel mode_to_level_from_signature(
    std::size_t rank, container::svector<std::optional<std::size_t>> const& sig,
    container::svector<DagScopeLevel> const& levels) {
  SEQUANT_ASSERT(sig.size() == levels.size());
  ModeToLevel m2l;
  m2l.by_mode.resize(rank);
  for (std::size_t i = 0; i < sig.size(); ++i)
    if (sig[i]) m2l.by_mode[*sig[i]] = levels[i];
  return m2l;
}

}  // namespace sequant

#endif  // SEQUANT_EVAL_SLICING_SIGNATURE_HPP
