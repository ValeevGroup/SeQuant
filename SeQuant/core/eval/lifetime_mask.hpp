#ifndef SEQUANT_EVAL_LIFETIME_MASK_HPP
#define SEQUANT_EVAL_LIFETIME_MASK_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node_compare.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>

#include <algorithm>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sequant {

/// Lifetime-mask stamping: \c stamp_lifetime_masks and the detail helpers it
/// alone depends on (\c slot_modes_of is also read by \c stamp_occurrence_homes
/// below, qualified as \c eval::detail::slot_modes_of).
namespace eval {

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

/// A node's own result-slot modes: its canonical indices, each taken as itself.
///
/// Batch loops are always over plain occ/aux modes (the DP's batchable modes;
/// a PAO/PNO composite \c a<i,j> is the CSV inner dimension, never a loop
/// axis), so a batch mode lives on this node iff it appears here as a plain
/// slot. A composite slot \c a<i,j> therefore contributes only mode \c a, not
/// its proto pair \c i,j: in array land \c a<i,j> is just mode \c a over an \c
/// <i,j>-tied range, and slicing an \c i / \c j loop does nothing to mode \c a.
/// Callers (the residency meet in \c stamp_residency_impl, and \c
/// value_key_impl's per-value keying below) intersect the in-scope batch modes
/// with this set to keep only those on \p n's result; a node carrying none of
/// a loop's mode is left unsliced by it (loop-invariant).
template <typename Node>
container::svector<Index> slot_modes_of(Node const& n) {
  auto const& ci = n->canon_indices();
  return {ci.begin(), ci.end()};
}

/// Shared cross-occurrence meet walk underlying \c stamp_lifetime_masks.
/// \p modes_of extracts the batch loops a node opens (each physical loop once,
/// at its open site -- \c batch_loops_opened_here, not the per-carrying-node
/// \c node_slice_mask), which the top-down accumulation propagates down as a
/// set (no dedup); \p setter stamps the resulting per-canonical-node meet onto
/// the node (e.g. \c set_sliced_modes). Everything else -- the meet map, the
/// accumulation, the per-slot filter, the two-pass order -- is identical
/// between entry points; only the selector and setter differ.
template <meta::eval_node_range R, typename ModesOf, typename Setter>
void stamp_residency_impl(R const& forest, ModesOf const& modes_of,
                          Setter const& setter) {
  using Node = std::ranges::range_value_t<R>;

  // Running per-canonical intersection, keyed by a live pointer into the
  // forest (hashed/compared by pointee structure via the tree
  // hasher/comparator) so grouping does not deep-copy subtrees.
  std::unordered_map<Node const*, container::svector<Index>,
                     TreeNodeHasher<Node>, TreeNodeEqualityComparator<Node>>
      meet;
  // Every internal-node occurrence, in visitation order, for the stamping
  // pass, paired with a pointer to its meet entry. Pass 2 does not re-find an
  // occurrence in `meet`: a lookup there is a structural comparison that walks
  // the whole subtree, so re-finding would make stamping an in-place Sum tree
  // O(nodes x subtree), i.e. quadratic in forest size. Pass 1 already knows the
  // entry, so it records it. `std::unordered_map` keeps references to its
  // elements valid across rehashing, so these stay good while pass 1 inserts.
  container::svector<std::pair<Node const*, container::svector<Index>*>> occ;

  // Pass 1: top-down walk accumulating the enclosing loops opened at or above
  // n. The full \p acc is passed down to children (descendants must know which
  // loops enclose them), but a node's contribution to the meet is \p acc
  // filtered to the modes that live on that node's own result slots.
  // The descent is iterative (an explicit frame stack), not recursive: an
  // equation's residual/energy is a single in-place Sum tree with one node per
  // summand, so its left spine is as deep as the number of terms -- thousands
  // for a large equation -- and a recursive descent would overflow the call
  // stack. Pushing the right child before the left one makes the pop order
  // pre-order (node, left subtree, right subtree), which is the order `occ`
  // and the meet are built in.
  struct Frame {
    Node const* n;
    container::svector<Index> acc;
  };
  std::vector<Frame> stack;
  auto const walk_node = [&](Node const& n, container::svector<Index> acc) {
    if (n.leaf()) return;  // leaves are not stamped (they carry no meet)
    // acc = the enclosing loops opened at or above n, each appearing once: a
    // physical loop is opened at a single site (\p modes_of reads opens, not
    // the per-carrying-node slice mask), so accumulating opens down every
    // root-to- node path visits each loop exactly once -- no dedup needed, acc
    // is a set by construction. (This is the whole point of sourcing opens: an
    // External loop reaches its carriers here by inheritance, and a Contracted
    // loop reaches its below-the-reduction carriers likewise, without either
    // being double-counted.)
    for (auto const& ix : modes_of(n)) acc.push_back(ix);
    // n's meet contribution: acc filtered to n's own result slots, in acc
    // order.
    auto const slots = slot_modes_of(n);
    container::svector<Index> node_modes;
    for (auto const& m : acc)
      if (std::find(slots.begin(), slots.end(), m) != slots.end())
        node_modes.push_back(m);
    if (auto it = meet.find(&n); it == meet.end()) {
      // first occurrence seeds the intersection
      auto const ins = meet.emplace(&n, std::move(node_modes)).first;
      occ.push_back({&n, &ins->second});
    } else {
      lifetime_mask_intersect_in_place(it->second, node_modes);
      occ.push_back({&n, &it->second});
    }
    stack.push_back({&n.right(), acc});
    stack.push_back({&n.left(), std::move(acc)});
  };
  for (auto const& tree : forest) {
    stack.push_back({&tree, {}});
    while (!stack.empty()) {
      Frame f = std::move(stack.back());
      stack.pop_back();
      walk_node(*f.n, std::move(f.acc));
    }
  }

  // Pass 2: stamp every occurrence with its canonical meet. The forest is
  // logically mutable (only the parameter binding is const); the setter
  // reaches the node payload to stamp it. No lookup: pass 1 recorded which
  // entry each occurrence belongs to, and the entry it points at now holds the
  // final intersection (later occurrences narrowed it in place).
  for (auto const& [n, modes] : occ) setter(n, *modes);
}

}  // namespace detail

/// Stamp each canonical eval node's cross-occurrence sliced-mode mask
/// (\c EvalExpr::sliced_modes) -- the runtime residency \c place_at_this_level
/// consumes to home each value. A mode slices a canonical node iff it slices
/// every occurrence of that node in \p forest (a *meet* / set-intersection over
/// occurrences). A node's occurrence-local sliced set is the union of the
/// enclosing loops opened at or above it -- \c batch_loops_opened_here (any
/// \c BatchModeType: External opened at its open site, Contracted opened at its
/// reduction site), each physical loop appearing once -- filtered to the modes
/// that live on the node's own result slots (\c slot_modes_of, i.e.
/// \c canon_indices() taken as-is -- see there for why a composite slot
/// \c a<i,j> is just mode \c a, never its \c i,j proto pair). A node invariant
/// to an outer batched loop -- it does not carry that loop's mode on any slot
/// -- is thus left all-full even under a batched ancestor, so it stays eligible
/// for loop-invariant reuse. (Sourcing opens, not the per-carrying-node
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
/// Idempotent; a no-op on the off path: with no \c node_slice_mask() stamps
/// every occurrence set is empty, so every meet is empty and every mask is
/// all-full
/// (\c EvalExpr::sliced_modes_ is default-empty), leaving runtime behavior
/// unchanged.
template <meta::eval_node_range R>
void stamp_lifetime_masks(R const& forest) {
  using Node = std::ranges::range_value_t<R>;
  using Data = typename Node::value_type;

  // The batch loops opened at a node (every kind: External at its open site,
  // Contracted at its reduction site), each physical loop named exactly once --
  // not node_slice_mask(), which the DP stamps on every carrying node. Reading
  // opens is what makes the accumulation below a genuine set: each enclosing
  // loop reaches a node once, propagated down from its single open site, so an
  // External loop reaches its carriers by inheritance (not a redundant own-node
  // stamp) and a Contracted loop reaches its below-the-reduction carriers the
  // same way -- the one case that genuinely needs the propagation.
  // (peak_profile reads opens for the identical reason.) The resulting
  // sliced_modes is the runtime residency place_at_this_level consumes and \c
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

}  // namespace eval

/// \brief The home residency of \p n: the loops opened at or above this node
/// filtered to its own result slots, stamped per occurrence by \c
/// stamp_occurrence_homes. The body -- \c n->occurrence_home() -- is the
/// definition.
///
/// It is not a cross-occurrence meet. \c EvalExpr::sliced_modes (stamped by \c
/// eval::stamp_lifetime_masks) is the meet, and it is a different quantity,
/// read by the forest-descent route only; one node sliced along different
/// modes in different terms keeps each occurrence's own slicing here, and
/// value identity (\c value_key_of) tells those occurrences apart. See the
/// as-built design, \c
/// doc/dev/specs/2026-09-12-batched-array-dag-eval-as-built.md section 5.4.
template <meta::eval_node Node>
container::svector<Index> const& home_scope(Node const& n) noexcept {
  return n->occurrence_home();
}

/// \brief Stamps every occurrence's own home (\c EvalExpr::occurrence_home):
/// the loops opened at or above the node, filtered to its own result slots,
/// with no cross-occurrence meet. One node sliced along different modes in
/// different terms, or read whole in one term and sliced in another, keeps
/// each occurrence's slicing; value identity (\c value_key_of) then tells the
/// occurrences apart instead of the meet folding them to a whole home
/// (explicit-cells design section 11). The table-driven engine's home; the
/// forest-descent path keeps \c eval::stamp_lifetime_masks' meet.
template <meta::eval_node_range R>
void stamp_occurrence_homes(R const& forest) {
  using Node = std::ranges::range_value_t<R>;
  using Data = typename Node::value_type;
  // Iterative for the same reason as \c stamp_residency_impl above: the
  // residual's Sum spine is as deep as the number of terms, and a recursive
  // descent would overflow the stack. Right child pushed before left, so the
  // pop order is the recursion's pre-order.
  struct Frame {
    Node const* n;
    container::svector<Index> acc;
  };
  std::vector<Frame> stack;
  auto const walk_node = [&](Node const& n, container::svector<Index> acc) {
    if (n.leaf()) return;
    for (auto const& [ix, kind] : n->batch_loops_opened_here())
      acc.push_back(ix);
    auto const slots = eval::detail::slot_modes_of(n);
    container::svector<Index> home;
    for (auto const& m : acc)
      if (std::find(slots.begin(), slots.end(), m) != slots.end())
        home.push_back(m);
    const_cast<Data&>(*n).set_occurrence_home(std::move(home));
    stack.push_back({&n.right(), acc});
    stack.push_back({&n.left(), std::move(acc)});
  };
  for (auto const& tree : forest) {
    stack.push_back({&tree, {}});
    while (!stack.empty()) {
      Frame f = std::move(stack.back());
      stack.pop_back();
      walk_node(*f.n, std::move(f.acc));
    }
  }
}

/// \brief The value key of a node: its node id (\c hash_value, the canonical
/// colored graph of its tensor network, label-free) combined with the sorted
/// canonical positions it is home-sliced on (explicit-cells design section
/// 11). Equal to the node id when nothing is home-sliced, so every unbatched
/// value keeps the identity it has today. Two occurrences are one value iff
/// they are one node and are home-sliced on the same positions: one node
/// sliced along two different modes of its array in two terms is two values,
/// each loop-local in its own nest, never resident whole.
inline std::size_t value_key(std::size_t node_hash,
                             container::svector<std::size_t> positions) {
  if (positions.empty()) return node_hash;
  std::sort(positions.begin(), positions.end());
  std::size_t h = node_hash;
  hash::combine(h, positions.size());
  for (std::size_t p : positions) hash::combine(h, p);
  return h;
}

namespace detail {
/// (key, any node of the subtree is home-sliced). The key is over the
/// production subtree: a node's own id and home-sliced positions combined
/// with its operands' keys -- a value that reduces a differently-sliced
/// operand in another nest is a different production, hence a different
/// value (the recursive, whole-subtree form of the node's own canonical
/// identity). A subtree with nothing home-sliced keys to the node id.
template <meta::eval_node Node>
std::pair<std::size_t, bool> value_key_impl(Node const& n) {
  // The left spine is unwound iteratively (the `spine` vector and the
  // bottom-up loop below), exactly as TreeNodeEqualityComparator does: an
  // equation's residual/energy is one in-place Sum tree with a left spine as
  // deep as the number of terms, and recursing down it would overflow the
  // call stack. Right children (single terms) and Product operands are bounded
  // in depth and stay recursive. Aside from that unwinding this is a faithful
  // transcription of the recursive form.
  container::svector<Node const*> spine;
  for (Node const* c = &n;; c = &(*c).left()) {
    spine.push_back(c);
    if ((*c).leaf()) break;
  }
  std::pair<std::size_t, bool> below{};  // the current node's left child result
  for (auto sit = spine.rbegin(); sit != spine.rend(); ++sit) {
    Node const& nd = **sit;
    container::svector<std::size_t> pos;
    auto const& carried = nd->canon_indices();
    for (Index const& m : home_scope(nd))
      for (std::size_t p = 0; p < carried.size(); ++p)
        if (carried[p] == m) {
          pos.push_back(p);
          break;
        }
    bool sliced = !pos.empty();
    std::size_t h = value_key(nd->hash_value(), std::move(pos));
    if (!nd.leaf()) {
      auto [lk, ls] = below;  // the left child, already folded by this loop
      auto [rk, rs] = value_key_impl(nd.right());
      if (ls || rs) {
        sliced = true;
        // A Product (contraction) is commutative and the binarizer may emit the
        // same contraction as (X,Y) in one term and (Y,X) in another, so the
        // two operand keys are combined in a canonical order -- ascending --
        // and the two spellings key to one value. (The node id they are
        // combined into already folds the operand order; combining them as
        // emitted did not, which split swapped occurrences into two builds.) A
        // Sum's operands are not interchangeable -- its left child is the
        // in-place accumulator -- so they stay in the emitted order.
        //
        // The order is taken on the keys, deliberately not via
        // canonical_children (eval_node_compare.hpp), whose rule is different
        // (scalar operand last, then ascending node id). Do not "unify" the
        // two: a node id ties whenever the two operands are one value, and
        // their keys can still differ there -- one value home-sliced two
        // different ways is two keys -- so canonical_children would fall back
        // to the emitted order and leave exactly this combination
        // order-dependent. Ordering the keys themselves is order-independent
        // unconditionally. Both rules are value-determined, so the two
        // canonicalizations do not need to agree.
        if (nd->is_product() && lk > rk) std::swap(lk, rk);
        hash::combine(h, lk);
        hash::combine(h, rk);
      }
    }
    below = {sliced ? h : nd->hash_value(), sliced};
  }
  return below;
}
}  // namespace detail

/// \brief \c value_key of a forest node over its production subtree: its
/// node id combined with the canonical positions (indices in \c
/// canon_indices) of its \c home_scope modes and with its operands' keys;
/// equal to the node id when nothing in the subtree is home-sliced (see \c
/// detail::value_key_impl).
template <meta::eval_node Node>
std::size_t value_key_of(Node const& n) {
  if (std::size_t const k = n->value_key(); k != 0) return k;  // stamped
  return detail::value_key_impl(n).first;
}

}  // namespace sequant

#endif  // SEQUANT_EVAL_LIFETIME_MASK_HPP
