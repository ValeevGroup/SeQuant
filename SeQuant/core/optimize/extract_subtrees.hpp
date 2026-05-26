#ifndef SEQUANT_OPTIMIZE_EXTRACT_SUBTREES_HPP
#define SEQUANT_OPTIMIZE_EXTRACT_SUBTREES_HPP

#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node_compare.hpp>

#include <concepts>
#include <ranges>
#include <unordered_set>

namespace sequant::opt {

namespace detail {

template <meta::eval_node Node, typename ResultSet>
void capture(Node& n, ResultSet& result) {
  result.insert(n);  // FullBinaryNode copy ctor is deep
  n = Node{*n};      // collapse to a leaf carrying a copy of n's payload
}

template <meta::eval_node Node, typename Pred, typename ResultSet>
bool extract_subtrees_visit(Node& n, Pred& pred, ResultSet& result) {
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
/// replacing each match in place with a leaf carrying a verbatim copy of
/// the original payload (expr, hash, canon_indices, phase, connectivity,
/// op_type).
///
/// `pred(n, pl, pr) -> bool` is tested in post-order; `pl`/`pr` are the
/// predicate values of `n`'s children (both `false` if `n` is a leaf).
/// A node `n` is captured iff `pred(n) ∧ (n is a root ∨ ¬pred(parent(n)))`.
/// Nested pred-positive descendants of a captured subtree are inlined,
/// not re-extracted. The synthetic leaf retains the original payload's
/// `op_type_`; consumers should treat the FullBinaryNode shape as
/// authoritative rather than branching on the payload's `op_type_`.
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
///
template <meta::eval_node_range Range, typename Pred>
  requires std::indirectly_writable<std::ranges::iterator_t<Range>,
                                    std::ranges::range_value_t<Range>> &&
           std::predicate<Pred&, std::ranges::range_value_t<Range> const&, bool,
                          bool>
auto extract_subtrees(Range& nodes, Pred pred) {
  using Node = std::ranges::range_value_t<Range>;
  using ResultSet = std::unordered_set<Node, TreeNodeHasher<Node>,
                                       TreeNodeEqualityComparator<Node>>;
  ResultSet result;
  for (auto& node : nodes) {
    if (detail::extract_subtrees_visit(node, pred, result))
      detail::capture<Node>(node, result);
  }
  return result;
}

}  // namespace sequant::opt

#endif  // SEQUANT_OPTIMIZE_EXTRACT_SUBTREES_HPP
