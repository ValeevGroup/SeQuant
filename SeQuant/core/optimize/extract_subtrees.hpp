#ifndef SEQUANT_OPTIMIZE_EXTRACT_SUBTREES_HPP
#define SEQUANT_OPTIMIZE_EXTRACT_SUBTREES_HPP

#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node_compare.hpp>

#include <concepts>
#include <cstddef>
#include <deque>
#include <ranges>
#include <unordered_set>
#include <utility>

// ideas for refactor
//
// - add a function overload (a different function name is fine as well) of
//   'extract_subtrees' that makes it easier to work with the pattern: leaf
//   tensors with labels in the list 'lbls' poison the intermediate (parent),
//   hence, extract out the sub-trees that are not poisoned and also not leaves.
//   This is the same pattern the example 'pred' shows in the doxx of
//   'extract_subtrees' (except more than one tensor labels to be generic)

namespace sequant::opt {

///
/// \brief Stores subtrees captured by extract_subtrees and provides a
///        bijection between an integer Id and the captured subtree.
///
/// Each insertion either reuses the Id of a tree-equal subtree already
/// present (under TreeNodeHasher/TreeNodeEqualityComparator) or assigns
/// a fresh Id equal to size() at the time of insertion. The id is the
/// single source of truth: it is stamped onto the captured subtree's
/// root (EvalExpr::extracted_id()) on insertion, and onto the synthetic
/// leaf that replaces the subtree in the forest (which is built by
/// copying the stamped root's payload). Callers can map from a synthetic
/// leaf in the rewritten forest back to the original subtree via at(id).
///
template <meta::eval_node Node>
class ExtractedSubtrees {
 public:
  using Id = std::size_t;
  using value_type = Node;

  /// Insert @p subtree. If a tree-equal subtree is already present,
  /// returns its existing Id (recovered from the stored root's stamp).
  /// Otherwise stamps a fresh Id equal to the previous size() onto the
  /// subtree's root, stores it, and returns the Id.
  Id intern(Node subtree) {
    if (auto it = dedup_.find(&subtree); it != dedup_.end())
      return (**it)->extracted_id().value();
    const Id id = by_id_.size();
    ExtractedIdSetter{}.set(static_cast<EvalExpr&>(*subtree), id);
    by_id_.push_back(std::move(subtree));
    dedup_.insert(&by_id_.back());
    return id;
  }

  [[nodiscard]] bool empty() const noexcept { return by_id_.empty(); }
  [[nodiscard]] std::size_t size() const noexcept { return by_id_.size(); }

  [[nodiscard]] Node const& at(Id id) const { return by_id_.at(id); }

  [[nodiscard]] auto begin() const noexcept { return by_id_.begin(); }
  [[nodiscard]] auto end() const noexcept { return by_id_.end(); }

 private:
  // deque so that pointers stored as dedup_ keys remain stable across inserts
  std::deque<Node> by_id_;
  std::unordered_set<Node const*, TreeNodeHasher<Node>,
                     TreeNodeEqualityComparator<Node>>
      dedup_;
};

namespace detail {

template <meta::eval_node Node>
void capture(Node& n, ExtractedSubtrees<Node>& result) {
  // intern moves n into the inventory on fresh insert (where it also
  // stamps the stored root with the id) or matches it against an existing
  // tree-equal entry (in which case n is consumed by the find/move but
  // the inventory is unchanged). The synthetic leaf payload is then a
  // copy of the stored root's payload, so it inherits the id stamp.
  const std::size_t id = result.intern(std::move(n));
  auto leaf_payload = *result.at(id);
  EvalOpSetter{}.reset(static_cast<EvalExpr&>(leaf_payload));
  n = Node{std::move(leaf_payload)};
}

template <meta::eval_node Node, typename Pred>
bool extract_subtrees_visit(Node& n, Pred& pred,
                            ExtractedSubtrees<Node>& result) {
  bool pl = false;
  bool pr = false;
  if (!n.leaf()) {
    pl = extract_subtrees_visit(n.left(), pred, result);
    pr = extract_subtrees_visit(n.right(), pred, result);
  }
  const bool pn = pred(static_cast<Node const&>(n), pl, pr);
  if (!pn) {
    if (pl) capture<Node>(n.left(), result);
    if (pr) capture<Node>(n.right(), result);
  }
  return pn;
}

}  // namespace detail

///
/// Extract maximal `pred`-matching subtrees from a forest of EvalNodes,
/// replacing each match in place with a leaf carrying a copy of the
/// original payload (expr, hash, canon_indices, phase, connectivity)
/// with `op_type_` reset to null so the synthetic leaf reads as a
/// primary expression, and `extracted_id_` set to the Id of the captured
/// subtree in the returned ExtractedSubtrees.
///
/// `pred(n, pl, pr) -> bool` is tested in post-order; `pl`/`pr` are the
/// predicate values of `n`'s children (both `false` if `n` is a leaf).
/// A node `n` is captured iff `pred(n) ∧ (n is a root ∨ ¬pred(parent(n)))`.
/// Nested pred-positive descendants of a captured subtree are inlined,
/// not re-extracted.
///
/// Tree-equal captured subtrees (under TreeNodeEqualityComparator) share
/// a single entry in the result and the corresponding synthetic leaves
/// all carry the same Id.
///
/// Example — extract intermediates that don't involve any tensor with
/// label "t". `pred(leaf) = false` keeps bare leaves out of the result
/// (and prevents whole-tree root collapse); the per-leaf test is folded
/// into the parent's check via `child_ok`:
///
///   auto pred = [](EvalNode<EvalExpr> const& n, bool pl, bool pr) -> bool {
///     if (n.leaf()) return false;
///     auto child_ok = [](EvalNode<EvalExpr> const& c, bool pc) {
///       return c.leaf() ? c->is_tensor() && c->as_tensor().label() != L"t"
///                       : pc;
///     };
///     return child_ok(n.left(), pl) && child_ok(n.right(), pr);
///   };
///   auto intermediates = sequant::opt::extract_subtrees(forest, pred);
///   // intermediates.at(leaf->extracted_id().value()) recovers the original
///   // subtree for any synthetic leaf left behind in `forest`.
///
template <meta::eval_node_range Range, typename Pred>
  requires std::indirectly_writable<std::ranges::iterator_t<Range>,
                                    std::ranges::range_value_t<Range>> &&
           std::predicate<Pred&, std::ranges::range_value_t<Range> const&, bool,
                          bool>
auto extract_subtrees(Range& nodes, Pred pred) {
  using Node = std::ranges::range_value_t<Range>;
  ExtractedSubtrees<Node> result;
  for (auto& node : nodes) {
    if (detail::extract_subtrees_visit(node, pred, result))
      detail::capture<Node>(node, result);
  }
  return result;
}

}  // namespace sequant::opt

#endif  // SEQUANT_OPTIMIZE_EXTRACT_SUBTREES_HPP
