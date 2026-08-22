#ifndef SEQUANT_EVAL_MEMBER_AXIS_HPP
#define SEQUANT_EVAL_MEMBER_AXIS_HPP

#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/meta.hpp>

#include <algorithm>
#include <cstddef>
#include <optional>
#include <utility>

/// \file member_axis.hpp
/// \brief The DAG-space -> node-physical-mode map used by BOTH batched
/// executors.
///
/// The schedule and its loop nest are indexed by DAG-scope coordinates
/// (index SPACES, e.g. the aux space Κ); a node's own physical modes are its
/// tree-local index labels (Κ_1, Κ_2, ...), meaningful only within its tree.
/// A batch loop over a space must therefore be MAPPED, at each node boundary,
/// to that node's own physical index of the space before it is used to slice
/// or accumulate the node. These helpers are that map (node + space base_key ->
/// node's physical Index / result position). Formerly private to
/// scope_executor.hpp (\c walk_scope); hoisted here so \c ordered_executor.hpp
/// can apply the identical map without a circular include (scope_executor.hpp
/// includes ordered_executor.hpp).

namespace sequant::eval::detail {

/// \return the position on \p n's own result of the FIRST index of TYPE \p base
///         (\c IndexSpace::base_key()), or nullopt. Unlike \c index_position,
///         which matches an EXACT \c Index, this matches by space so it finds
///         whichever physical label \p n carries for that space -- the
///         result-side half of the DAG-space -> physical-mode map. \note If a
///         node carried more than one index of \p base this returns the first;
///         no current batch axis binds two physical modes of one space on a
///         single node.
template <meta::eval_node node_t>
[[nodiscard]] std::optional<std::size_t> result_position_type(
    node_t const& n, std::wstring const& base) {
  auto const& idxs = n->canon_indices();
  for (std::size_t p = 0; p < idxs.size(); ++p)
    if (idxs[p].space().base_key() == base) return p;
  return std::nullopt;
}

/// \return the position on \p n's own result of the first index of TYPE \p base
///         (\c IndexSpace::base_key()) that is ALSO one of \p n's STAMPED
///         sliced modes (\c sliced_modes(), the cross-occurrence batch meet),
///         or nullopt. Unlike \c result_position_type -- which matches ANY
///         same-space index and so would slice an index the optimizer never
///         batched (e.g. an internal / contracted occ that merely shares the
///         external occ's space) -- this honors the stamp: a space-mapped slice
///         must never touch a mode outside the batched set. Keeps the runtime
///         index-generic rather than space-generic (see the space-map fallback
///         in \c slice_to_use, eval.hpp).
template <meta::eval_node node_t>
[[nodiscard]] std::optional<std::size_t> sliced_result_position_type(
    node_t const& n, std::wstring const& base) {
  auto const& idxs = n->canon_indices();
  auto const& sliced = n->sliced_modes();
  for (std::size_t p = 0; p < idxs.size(); ++p)
    if (idxs[p].space().base_key() == base &&
        std::find(sliced.begin(), sliced.end(), idxs[p]) != sliced.end())
      return p;
  return std::nullopt;
}

/// \return (leaf, position) of the first leaf below \p n (or \p n itself) that
///         carries an index of TYPE \p base, or nullopt. TYPE-keyed counterpart
///         of \c sequant::find_leaf_carrying.
template <meta::eval_node node_t>
[[nodiscard]] std::optional<std::pair<node_t, std::size_t>>
find_leaf_carrying_type(node_t const& n, std::wstring const& base) {
  if (n.leaf()) {
    if (auto const p = result_position_type(n, base)) return std::pair{n, *p};
    return std::nullopt;
  }
  if (auto l = find_leaf_carrying_type(n.left(), base)) return l;
  return find_leaf_carrying_type(n.right(), base);
}

/// \return the PHYSICAL index \p root batches as an External (spectator) mode
/// of
///         TYPE \p base at its own root (from \c batched_here()), or nullopt.
///         An External loop is realized per-member over the member's own
///         physical index, so scatter uses this (not the schedule's canonical
///         representative) -- the schedule's mode only names the TYPE.
template <meta::eval_node node_t>
[[nodiscard]] std::optional<Index> member_external_axis(
    node_t const& root, std::wstring const& base) {
  if (root.leaf()) return std::nullopt;
  for (auto const& [ix, knd] : root->batched_here())
    if (knd == BatchModeType::External && ix.space().base_key() == base)
      return ix;
  return std::nullopt;
}

/// \return the PHYSICAL index \p root batches as a Contracted mode of TYPE \p
///         base (from \c batched_here() -- the authoritative source of WHICH
///         physical index is batched), or, failing that, the physical label a
///         leaf below \p root carries for that type.
///
/// \details Index labels are meaningful only within a single tree, so the
/// schedule's canonical scope mode (\c ScopeNode::mode) must be MAPPED to each
/// member's OWN physical axis before it is pushed into that member's batch
/// context / used to slice it -- exactly as \c member_external_axis does for
/// the external path. Reusing the schedule's canonical mode for every member
/// would (silently, since each tree today binds a single aux label) fail to
/// slice a member that binds the contracted mode under a different physical
/// label, building it full and mis-accumulating it. A contracted mode is summed
/// away below the root (never on its result), so this reads it off \c
/// batched_here() or a carrying leaf, not the root's own \c canon_indices.
template <meta::eval_node node_t>
[[nodiscard]] std::optional<Index> member_contracted_axis(
    node_t const& root, std::wstring const& base) {
  if (!root.leaf())
    for (auto const& [ix, knd] : root->batched_here())
      if (knd == BatchModeType::Contracted && ix.space().base_key() == base)
        return ix;
  if (auto const lf = find_leaf_carrying_type(root, base))
    return lf->first->canon_indices()[lf->second];
  return std::nullopt;
}

}  // namespace sequant::eval::detail

#endif  // SEQUANT_EVAL_MEMBER_AXIS_HPP
