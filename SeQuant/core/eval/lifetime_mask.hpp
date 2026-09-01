#ifndef SEQUANT_EVAL_LIFETIME_MASK_HPP
#define SEQUANT_EVAL_LIFETIME_MASK_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node_compare.hpp>
#include <SeQuant/core/index.hpp>

#include <algorithm>
#include <unordered_map>
#include <utility>

namespace sequant {

namespace detail {

/// Intersect \p acc in place with \p other by \c Index identity (operator==):
/// keep only elements of \p acc that also appear in \p other, preserving
/// \p acc's order. Used to accumulate the cross-occurrence meet.
inline void lifetime_mask_intersect_in_place(
    container::svector<Index>& acc, container::svector<Index> const& other) {
  container::svector<Index> keep;
  for (auto const& ix : acc)
    if (std::find(other.begin(), other.end(), ix) != other.end())
      keep.push_back(ix);
  acc = std::move(keep);
}

/// A node's OWN result-slot modes: its canonical indices, each taken as itself.
///
/// Batch loops are always over PLAIN occ/aux modes (the DP's batchable modes;
/// a PAO/PNO composite \c a<i,j> is the CSV inner dimension, never a loop
/// axis), so a batch mode lives on this node iff it appears here as a plain
/// slot. A composite slot \c a<i,j> therefore contributes only mode \c a, NOT
/// its proto pair \c i,j: in array land \c a<i,j> is just mode \c a over an \c
/// <i,j>-tied range, and slicing an \c i / \c j loop does nothing to mode \c a.
/// Callers (the residency meet in \c stamp_residency_impl, the occurrence key
/// in
/// \c occurrence_key.hpp) intersect the in-scope batch modes with this set to
/// keep only those on \p n's result; a node carrying none of a loop's mode is
/// left unsliced by it (loop-invariant).
template <typename Node>
container::svector<Index> slot_modes_of(Node const& n) {
  auto const& ci = n->canon_indices();
  return {ci.begin(), ci.end()};
}

/// Shared cross-occurrence meet walk underlying \c stamp_lifetime_masks.
/// \p modes_of extracts the batch loops a node OPENS (each physical loop once,
/// at its open site -- \c batch_loops_opened_here, NOT the per-carrying-node
/// \c node_slice_mask), which the top-down accumulation propagates down as a
/// set (no dedup); \p setter stamps the resulting per-canonical-node meet onto
/// the node (e.g. \c set_sliced_modes). Everything else -- the meet map, the
/// accumulation, the per-slot filter, the two-pass order -- is identical
/// between entry points; only the selector and setter differ.
template <meta::eval_node_range R, typename ModesOf, typename Setter>
void stamp_residency_impl(R const& forest, ModesOf const& modes_of,
                          Setter const& setter) noexcept {
  using Node = std::ranges::range_value_t<R>;

  // Running per-canonical intersection, keyed by a live pointer INTO the
  // forest (hashed/compared by pointee structure via the tree
  // hasher/comparator) so grouping does not deep-copy subtrees.
  std::unordered_map<Node const*, container::svector<Index>,
                     TreeNodeHasher<Node>, TreeNodeEqualityComparator<Node>>
      meet;
  // Every internal-node occurrence, in visitation order, for the stamping
  // pass.
  container::svector<Node const*> occ;

  // Pass 1: top-down walk accumulating the enclosing loops OPENED at or above
  // n. The full \p acc is passed DOWN to children (descendants must know which
  // loops enclose them), but a node's contribution to the meet is \p acc
  // filtered to the modes that live on that node's OWN result slots.
  auto walk = [&](auto&& self, Node const& n,
                  container::svector<Index> acc) -> void {
    if (n.leaf()) return;  // leaves are not stamped (they carry no meet)
    // acc = the enclosing loops opened at or above n, each appearing ONCE: a
    // physical loop is opened at a single site (\p modes_of reads opens, not
    // the per-carrying-node slice mask), so accumulating opens down every
    // root-to- node path visits each loop exactly once -- no dedup needed, acc
    // is a set by construction. (This is the whole point of sourcing opens: an
    // External loop reaches its carriers here by inheritance, and a Contracted
    // loop reaches its below-the-reduction carriers likewise, without either
    // being double-counted.)
    for (auto const& ix : modes_of(n)) acc.push_back(ix);
    // n's meet contribution: acc filtered to n's OWN result slots, in acc
    // order.
    auto const slots = slot_modes_of(n);
    container::svector<Index> node_modes;
    for (auto const& m : acc)
      if (std::find(slots.begin(), slots.end(), m) != slots.end())
        node_modes.push_back(m);
    occ.push_back(&n);
    if (auto it = meet.find(&n); it == meet.end())
      meet.emplace(&n, node_modes);  // first occurrence seeds the
                                     // intersection
    else
      lifetime_mask_intersect_in_place(it->second, node_modes);
    self(self, n.left(), acc);
    self(self, n.right(), acc);
  };
  for (auto const& tree : forest) walk(walk, tree, {});

  // Pass 2: stamp every occurrence with its canonical meet. The forest is
  // logically mutable (only the parameter binding is const); the setter
  // reaches the node payload to stamp it.
  for (Node const* n : occ)
    if (auto it = meet.find(n); it != meet.end()) setter(n, it->second);
}

}  // namespace detail

/// Stamp each canonical eval node's cross-occurrence sliced-mode mask
/// (\c EvalExpr::sliced_modes) -- the runtime residency \c place_at_this_level
/// consumes to home each value. A mode slices a canonical node iff it slices
/// EVERY occurrence of that node in \p forest (a *meet* / set-intersection over
/// occurrences). A node's occurrence-local sliced set is the union of the
/// enclosing loops OPENED at or above it -- \c batch_loops_opened_here (any
/// \c BatchModeType: External opened at its open site, Contracted opened at its
/// reduction site), each physical loop appearing once -- filtered to the modes
/// that live on the node's OWN result slots (\c slot_modes_of, i.e.
/// \c canon_indices() taken as-is -- see there for why a composite slot
/// \c a<i,j> is just mode \c a, never its \c i,j proto pair). A node invariant
/// to an outer batched loop -- it does not carry that loop's mode on any slot
/// -- is thus left all-full even under a batched ancestor, so it stays eligible
/// for loop-invariant reuse. (Sourcing OPENS, not the per-carrying-node
/// \c node_slice_mask, is what makes the accumulation a set: an External loop
/// reaches its carriers by inheritance rather than a redundant own-node stamp,
/// and a Contracted loop reaches its below-the-reduction carriers the same way
/// -- the case that genuinely needs the down-propagation.)
///
/// This all-batched-modes meet subsumes the former per-occurrence
/// \c contracted_modes bolt-on: a node variant to an outer contracted (aux)
/// loop carries that aux free on a result slot, so the aux mode survives the
/// meet and lands in \c sliced_modes directly.
///
/// Occurrences are grouped by canonical identity (\c hash_value plus structural
/// \c TreeNodeEqualityComparator equivalence). The meet is a set-intersection
/// by
/// \c Index identity: canonicalization gives consistent labels across
/// occurrences, so a genuinely sliced-everywhere mode survives, while a
/// block-agnostic node (e.g. \c s*C) whose occurrences bind disjoint concrete
/// modes intersects to empty (all-full).
///
/// Idempotent; a no-op on the OFF path: with no \c node_slice_mask() stamps
/// every occurrence set is empty, so every meet is empty and every mask is
/// all-full
/// (\c EvalExpr::sliced_modes_ is default-empty), leaving runtime behavior
/// unchanged.
template <meta::eval_node_range R>
void stamp_lifetime_masks(R const& forest) noexcept {
  using Node = std::ranges::range_value_t<R>;
  using Data = typename Node::value_type;

  // The batch loops OPENED at a node (EVERY kind: External at its open site,
  // Contracted at its reduction site), each physical loop named exactly ONCE --
  // NOT node_slice_mask(), which the DP stamps on every CARRYING node. Reading
  // opens is what makes the accumulation below a genuine set: each enclosing
  // loop reaches a node once, propagated down from its single open site, so an
  // External loop reaches its carriers by INHERITANCE (not a redundant own-node
  // stamp) and a Contracted loop reaches its below-the-reduction carriers the
  // SAME way -- the one case that genuinely needs the propagation.
  // (peak_profile reads opens for the identical reason.) The resulting
  // sliced_modes IS the runtime residency place_at_this_level consumes and \c
  // home_scope returns.
  auto opened_loops_of = [](Node const& n) {
    container::svector<Index> v;
    for (auto const& [ix, kind] : n->batch_loops_opened_here()) v.push_back(ix);
    return v;
  };

  detail::stamp_residency_impl(
      forest, opened_loops_of, [](Node const* n, container::svector<Index> m) {
        const_cast<Data&>(**n).set_sliced_modes(std::move(m));
      });
}

/// \brief The perfect-CSE seed home residency of \p n: the unified
/// all-batched-modes cross-occurrence meet on its own result slots -- i.e. the
/// runtime \c EvalExpr::sliced_modes (stamped by \c stamp_lifetime_masks). The
/// former separate \c seed_residency field was retired once the runtime mask
/// dropped its External-only filter and the two became identical. Returns the
/// residency mode-set; runtime depth resolution against a batch context is the
/// existing rl-walk, reused by the peak profile / remat.
template <meta::eval_node Node>
container::svector<Index> const& home_scope(Node const& n) noexcept {
  return n->sliced_modes();
}

}  // namespace sequant

#endif  // SEQUANT_EVAL_LIFETIME_MASK_HPP
