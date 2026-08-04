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

/// Append \p ix to \p v, proto-expanded: a composite index (e.g. a PNO tied
/// to an occ pair) contributes its \c proto_indices() instead of itself (a
/// composite index is domain-tied to its occ pair, so slicing the pair
/// slices it too); a plain index contributes itself. Shared by the batch-mode
/// side (\c stamp_residency_impl's \c modes_of selectors) and the slot side
/// (\c slot_modes_of) of the same proto expansion.
inline void proto_expand_into(container::svector<Index>& v, Index const& ix) {
  if (ix.has_proto_indices())
    for (auto const& p : ix.proto_indices()) v.push_back(p);
  else
    v.push_back(ix);
}

/// Modes exposed by node \p n's OWN canonical slots, proto-expanded (the SLOT
/// side of the same proto expansion the batch-mode selectors apply): a plain
/// slot exposes itself; a composite slot \c a<i,j> exposes its proto indices
/// \c i,j. A mode belongs in a node's residency mask only if it is one of
/// these -- i.e. it lives on \p n's own result. A node invariant to an outer
/// batched loop (it does not carry that loop's mode on any of its slots) is
/// thus NOT stamped by it.
template <typename Node>
container::svector<Index> slot_modes_of(Node const& n) {
  container::svector<Index> v;
  for (auto const& s : n->canon_indices()) proto_expand_into(v, s);
  return v;
}

/// Shared cross-occurrence meet walk underlying both \c stamp_lifetime_masks
/// and \c stamp_seed_residency. \p modes_of extracts a node's
/// occurrence-local batch modes (already proto-expanded); \p setter stamps
/// the resulting per-canonical-node meet onto the node (e.g. \c
/// set_sliced_modes / \c set_seed_residency). Everything else -- the meet
/// map, the top-down accumulation, the per-slot filter, the two-pass order --
/// is identical between the two entry points; only the selector and the
/// setter differ.
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

  // Pass 1: top-down walk accumulating ancestor+self batch modes. The full
  // \p acc is passed DOWN to children (descendants must know which modes are
  // batched above them), but a node's contribution to the meet is \p acc
  // filtered to the modes that slice one of that node's OWN slots.
  auto walk = [&](auto&& self, Node const& n,
                  container::svector<Index> acc) -> void {
    if (n.leaf()) return;  // leaves are not stamped (they carry no meet)
    for (auto const& ix : modes_of(n)) acc.push_back(ix);
    // Filter acc to modes that slice one of n's own slots (per-slot
    // semantics), preserving acc's order; store THAT as n's meet
    // contribution.
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
/// (\c EvalExpr::sliced_modes). A mode slices a canonical node iff it slices
/// EVERY occurrence of that node in \p forest (a *meet* / set-intersection over
/// occurrences). A node's occurrence-local sliced set is the union of the
/// External \c batched_here() stamps of the node and all its ancestors,
/// expanded proto-aware: a batched composite index contributes its
/// \c proto_indices() (a PNO/composite index is domain-tied to its occ pair, so
/// slicing the pair slices it too). That union is then filtered to the modes
/// that slice one of the node's OWN canonical slots (\c canon_indices(),
/// proto-expanded on the slot side): a mode belongs in the mask only if it
/// lives on this node's result. A node invariant to an outer external loop --
/// it does not carry that loop's mode on any slot -- is thus left all-full even
/// under a batched ancestor, so it stays eligible for loop-invariant reuse.
///
/// Occurrences are grouped by canonical identity (\c hash_value plus structural
/// \c TreeNodeEqualityComparator equivalence). The meet is a set-intersection
/// by
/// \c Index identity: canonicalization gives a shared proto pair consistent
/// labels across occurrences, so a genuinely sliced-everywhere mode survives,
/// while a block-agnostic node (e.g. \c s*C) whose occurrences bind disjoint
/// concrete modes intersects to empty (all-full).
///
/// Idempotent; a no-op on the OFF path: with no External \c batched_here()
/// stamps every occurrence set is empty, so every meet is empty and every mask
/// is all-full (\c EvalExpr::sliced_modes_ is default-empty), leaving runtime
/// behavior unchanged. Nothing consumes the mask yet; this pass is purely
/// additive.
template <meta::eval_node_range R>
void stamp_lifetime_masks(R const& forest) noexcept {
  using Node = std::ranges::range_value_t<R>;
  using Data = typename Node::value_type;

  // Occurrence-local External batch modes AT a node, proto-expanded.
  auto ext_modes_of = [](Node const& n) {
    container::svector<Index> v;
    for (auto const& [ix, kind] : n->batched_here())
      if (kind == BatchModeType::External) detail::proto_expand_into(v, ix);
    return v;
  };

  detail::stamp_residency_impl(
      forest, ext_modes_of, [](Node const* n, container::svector<Index> m) {
        const_cast<Data&>(**n).set_sliced_modes(std::move(m));
      });
}

///
/// \brief Stamp each canonical eval node's \c EvalExpr::seed_residency: the
/// cross-occurrence meet (see \c stamp_lifetime_masks for the full meet /
/// per-slot-filter semantics) of ALL batched modes -- any \c BatchModeType,
/// not just \c External -- on the node's own result slots. This is the
/// perfect-CSE \c home_scope seed residency used by the (Phase 3) static
/// predictor; it is a superset selector of \c stamp_lifetime_masks (dropping
/// the \c External-only filter on the batch-mode side), sharing the identical
/// walk, meet, and per-slot filter via \c detail::stamp_residency_impl.
///
/// Purely additive: writes only \c EvalExpr::seed_residency_, never \c
/// sliced_modes_, so it does not affect the runtime External-only masking
/// path (\c stamp_lifetime_masks / \c place_at_this_level) in any way.
///
template <meta::eval_node_range R>
void stamp_seed_residency(R const& forest) noexcept {
  using Node = std::ranges::range_value_t<R>;
  using Data = typename Node::value_type;

  // Occurrence-local batch modes AT a node, proto-expanded, EVERY kind (not
  // just External).
  auto all_batched_modes_of = [](Node const& n) {
    container::svector<Index> v;
    for (auto const& [ix, kind] : n->batched_here())
      detail::proto_expand_into(v, ix);
    return v;
  };

  detail::stamp_residency_impl(
      forest, all_batched_modes_of,
      [](Node const* n, container::svector<Index> m) {
        const_cast<Data&>(**n).set_seed_residency(std::move(m));
      });
}

}  // namespace sequant

#endif  // SEQUANT_EVAL_LIFETIME_MASK_HPP
