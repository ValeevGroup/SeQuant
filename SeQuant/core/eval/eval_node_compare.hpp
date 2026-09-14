#ifndef SEQUANT_EVAL_EVAL_NODE_COMPARE_HPP
#define SEQUANT_EVAL_EVAL_NODE_COMPARE_HPP

#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/tensor.hpp>
#include <SeQuant/external/bliss/graph.hh>

#include <cstddef>
#include <unordered_map>
#include <utility>

namespace sequant {

/// Functor to compute the hash of a given (evaluation) tree node.
///
/// This is the NODE-id hasher (\c hash::value): batching-blind, used by
/// \c compute_dag_boulevard / CSE and the top-level cache. The runtime cache
/// keys by \c CachedValue (\c value_id.hpp), whose hasher (\c
/// CachedValueHasher) and equality (\c CachedValueEqual) reduce to this
/// hasher and \c TreeNodeEqualityComparator on the wrapped node.
template <typename TreeNode, bool force_hash_collisions = false>
struct TreeNodeHasher {
  /// Trait used by the C++ STL allowing heterogenous lookups
  using is_transparent = void;

  std::size_t operator()(const TreeNode *node) const { return (*this)(*node); }

  std::size_t operator()(const TreeNode &node) const {
    if constexpr (force_hash_collisions) {
      return 0;
    }
    return hash::value(*node);
  }
};

/// \brief 3-way comparison establishing the CANONICAL order of a commutative
///        (Product) node's two operands -- a function of the operands' VALUES
///        only, never of the order the binarizer happened to emit them in.
///
/// \see canonical_children, the view that reads a node's children through it.
template <typename TreeNode>
int canonical_operand_cmp(TreeNode const &a, TreeNode const &b);

/// \brief A Product node's two children in CANONICAL order (the node itself
///        is NOT modified; its left/right stay exactly as emitted).
///
/// \see canonical_operand_cmp
template <typename TreeNode>
std::pair<TreeNode const &, TreeNode const &> canonical_children(
    TreeNode const &n);

/// Functor to compare two trees for equivalence
/// Explicit equivalence checking mitigates (accidental) hash collisions
template <typename TreeNode>
struct TreeNodeEqualityComparator {
  /// Trait used by the C++ STL allowing heterogenous lookups
  using is_transparent = void;

  TreeNodeEqualityComparator() = default;
  TreeNodeEqualityComparator(std::vector<Index> indices)
      : block_comparator_(std::move(indices)) {}

  bool operator()(const TreeNode *lhs, const TreeNode *rhs) const {
    return (*this)(*lhs, *rhs);
  }

  bool operator()(const TreeNode &lhs, const TreeNode *rhs) const {
    return (*this)(lhs, *rhs);
  }

  bool operator()(const TreeNode *lhs, const TreeNode &rhs) const {
    return (*this)(*lhs, rhs);
  }

  bool operator()(const TreeNode &lhs_in, const TreeNode &rhs_in) const {
    // The left-child descent is performed iteratively (this loop, not
    // recursion): an equation's residual/energy is kept as a single in-place
    // Sum-tree with one node per summand, so its left spine is as deep as the
    // number of terms -- thousands for e.g. a UCC BCH energy expansion -- and a
    // recursive descent down that spine would overflow the call stack. Right
    // children (individual terms) and Product operands are bounded in depth and
    // stay recursive; only the deep spine is unwound into this loop. Aside from
    // that unwinding this is a faithful transcription of the recursive
    // comparison -- same per-node checks, same ordered/unordered child logic.
    const TreeNode *lhsp = &lhs_in;
    const TreeNode *rhsp = &rhs_in;
    while (true) {
      const TreeNode &lhs = *lhsp;
      const TreeNode &rhs = *rhsp;

      if (lhs.leaf() != rhs.leaf()) {
        return false;
      }

      if (lhs.size() != rhs.size()) {
        return false;
      }

      if (hash::value(*lhs) != hash::value(*rhs)) {
        return false;
      }

      if (lhs->type_id() != rhs->type_id()) {
        return false;
      }

      if (lhs->is_constant() || lhs->is_variable() || lhs->is_power()) {
        if (*lhs->expr() != *rhs->expr()) {
          return false;
        }
      } else if (lhs->is_tensor() && !lhs->has_connectivity_graph()) {
        // Tensor nodes that carry a canonical connectivity graph
        // (tensor-network intermediates and proto-indexed leaves) are compared
        // exactly by that graph in the connectivity check below: it is the same
        // canonical colored graph the eval-node hash is derived from, its 3-way
        // cmp is a complete network identity, and it already folds bra<->ket
        // orientation (the orientation was canonicalized into the graph). Block
        // comparison here would instead compare on the stored bra/ket slot
        // order — a partial, orientation-sensitive signature — and wrongly
        // separate e.g. a CSV/PNO coefficient C{a;μ̃} from its equivalent
        // C{μ̃;a}, or the two g·C intermediates from transforming a real DF
        // factor g(μ̃,μ̃,Κ) on its bra vs its ket leg. So only graph-less tensor
        // nodes — protoindex-free leaves (block-canonicalized in place at
        // construction) and scalar*tensor results — are compared here by block.
        const Tensor &lhs_tensor = lhs->as_tensor();
        const Tensor &rhs_tensor = rhs->as_tensor();

        if (!block_comparator_(lhs_tensor, rhs_tensor)) {
          return false;
        }
      }

      if (lhs->has_connectivity_graph() != rhs->has_connectivity_graph()) {
        return false;
      }

      // Check connectivity in products / contractions
      if (lhs->has_connectivity_graph()) {
        SEQUANT_ASSERT(lhs->has_connectivity_graph());
        SEQUANT_ASSERT(rhs->has_connectivity_graph());

        if (bliss::ConstGraphCmp::cmp(lhs->connectivity_graph(),
                                      rhs->connectivity_graph()) != 0) {
          return false;
        }
      }

      if (lhs.leaf()) {
        return true;
      }

      // A Product (contraction) is commutative: the binarizer may emit the
      // same contraction as (X,Y) in one term and (Y,X) in another, and the
      // node id (hash) and the canonical connectivity graph both fold that
      // operand order -- so the child comparison must fold it too, else two
      // swapped occurrences of the same value are wrongly split into two cache
      // entries (built twice).
      //
      // It folds it by comparing the children through a CANONICAL VIEW
      // (canonical_children) rather than by trying both pairings: the view
      // picks which child comes first from the children's VALUES alone, so the
      // two spellings present the same ordered pair here and ONE recursion per
      // child suffices. Matching the children as an unordered pair would fold
      // them too, but it recurses twice per child and so is exponential in
      // product depth. Evaluation order is untouched: the node keeps the
      // left/right the DP emitted (cost, and peak in particular, is
      // order-dependent), only IDENTITY reads them through the view.
      //
      // The recursive child comparison is still REQUIRED and is NOT redundant
      // with the graph check above: the connectivity graph encodes only the
      // immediate two factors' connectivity, not each factor's recursive
      // build, so two products with the same immediate graph but different
      // sub-values must still be told apart. Product operand subtrees are
      // bounded in depth (a contraction of a fixed set of factors), so they
      // stay recursive.
      if (lhs->op_type() && *lhs->op_type() == EvalOp::Product) {
        auto const [lfirst, lsecond] = canonical_children(lhs);
        auto const [rfirst, rsecond] = canonical_children(rhs);
        return (*this)(lfirst, rfirst) && (*this)(lsecond, rsecond);
      }

      // Non-Product internal node (Sum, scalar*tensor product, adjoint): its
      // left/right assignment is canonical (e.g. the in-place Sum tree is
      // left-folded; an Adjoint's right child is a sentinel), so operand order
      // carries meaning and the children are compared in order. The right child
      // (a single summand / the scalar factor / the adjoint sentinel) is
      // bounded in depth and compared recursively; the left child is the deep
      // spine, so rather than recurse into it we loop back to the top with it
      // as the new (lhs, rhs) -- unwinding the spine iteratively.
      if (!(*this)(lhs.right(), rhs.right())) {
        return false;
      }
      lhsp = &lhs.left();
      rhsp = &rhs.left();
    }
  }

 private:
  IndexSpecificTensorBlockEqualComparator block_comparator_;
};

namespace detail {

/// 3-way compare of two canonical index lists: length first, then
/// \c Index::full_label lexicographically.
template <typename Indices>
int cmp_canon_indices(Indices const &a, Indices const &b) {
  if (a.size() != b.size()) return a.size() < b.size() ? -1 : 1;
  for (std::size_t i = 0; i < a.size(); ++i) {
    auto const &la = a[i].full_label();
    auto const &lb = b[i].full_label();
    if (la != lb) return la < lb ? -1 : 1;
  }
  return 0;
}

/// \brief Deterministic, value-determined 3-way tie-break used by
///        \c canonical_operand_cmp when two operands share a node id AND are
///        NOT equal -- i.e. a genuine 64-bit hash collision between distinct
///        values.
///
/// \details It inspects, in order, every node-local discriminator the equality
/// comparator itself uses -- subtree size, leaf-ness, expression type, result
/// type, the canonical connectivity graph, the canonical index list (which
/// subsumes the index-block signature used for graph-less tensor nodes), the
/// phase, the node label, and, for the scalar-valued nodes the comparator
/// compares symbolically, the expression itself -- and only then descends into
/// the children.
///
/// The descent goes through \c canonical_children at Product nodes, exactly as
/// the equality comparator does, so the result is a function of the subtrees'
/// VALUES and not of the order the binarizer emitted any node's operands in.
/// For the same reason it never compares an INTERNAL node's \c expr(): a
/// Product node's expression is assembled from its operands in emission order.
/// Like the equality comparator, the descent is iterative along the left
/// (potentially thousands-deep Sum) spine and recursive into bounded children.
///
/// A wrong answer here cannot make two distinct values compare equal -- the
/// equality comparator still compares both subtrees in full. The worst a
/// non-value-determined tie-break could do is order two colliding operands
/// differently in two terms, i.e. MISS a fold, which is the behaviour that
/// preceded any folding at all.
template <typename TreeNode>
int deep_tie_break(TreeNode const &lhs_in, TreeNode const &rhs_in) {
  const TreeNode *lhsp = &lhs_in;
  const TreeNode *rhsp = &rhs_in;
  while (true) {
    const TreeNode &lhs = *lhsp;
    const TreeNode &rhs = *rhsp;

    if (&lhs == &rhs) return 0;
    if (lhs.size() != rhs.size()) return lhs.size() < rhs.size() ? -1 : 1;
    if (lhs.leaf() != rhs.leaf()) return lhs.leaf() ? -1 : 1;
    if (lhs->type_id() != rhs->type_id())
      return lhs->type_id() < rhs->type_id() ? -1 : 1;
    if (lhs->result_type() != rhs->result_type())
      return lhs->result_type() < rhs->result_type() ? -1 : 1;
    if (lhs->has_connectivity_graph() != rhs->has_connectivity_graph())
      return lhs->has_connectivity_graph() ? -1 : 1;
    if (lhs->has_connectivity_graph()) {
      int const c = bliss::ConstGraphCmp::cmp(lhs->connectivity_graph(),
                                              rhs->connectivity_graph());
      if (c != 0) return c < 0 ? -1 : 1;
    }
    if (int const c =
            cmp_canon_indices(lhs->canon_indices(), rhs->canon_indices());
        c != 0)
      return c;
    if (lhs->canon_phase() != rhs->canon_phase())
      return lhs->canon_phase() < rhs->canon_phase() ? -1 : 1;
    if (lhs->label() != rhs->label())
      return lhs->label() < rhs->label() ? -1 : 1;
    // Symbolic form: only where it is emission-order-independent -- a leaf's
    // own tensor/constant/variable, and the scalar-valued nodes the equality
    // comparator compares by `*lhs->expr() != *rhs->expr()`.
    if (lhs.leaf() || lhs->is_constant() || lhs->is_variable() ||
        lhs->is_power()) {
      bool const le = static_cast<bool>(lhs->expr());
      bool const re = static_cast<bool>(rhs->expr());
      if (le != re) return le ? -1 : 1;
      if (le) {
        auto const la = lhs->expr()->to_latex();
        auto const lb = rhs->expr()->to_latex();
        if (la != lb) return la < lb ? -1 : 1;
      }
    }

    if (lhs.leaf()) return 0;

    if (lhs->op_type() && *lhs->op_type() == EvalOp::Product) {
      auto const [lfirst, lsecond] = canonical_children(lhs);
      auto const [rfirst, rsecond] = canonical_children(rhs);
      if (int const c = deep_tie_break(lfirst, rfirst); c != 0) return c;
      return deep_tie_break(lsecond, rsecond);
    }

    if (int const c = deep_tie_break(lhs.right(), rhs.right()); c != 0)
      return c;
    lhsp = &lhs.left();
    rhsp = &rhs.left();
  }
}

}  // namespace detail

/// \brief 3-way comparison establishing the CANONICAL order of a commutative
///        (Product) node's two operands.
///
/// \details A contraction is commutative, so the binarizer is free to emit its
/// two operands in either order -- and does, since the single-term DP's
/// contraction sequence decides which of the pair lands on the stack first.
/// The node id and the canonical connectivity graph both fold that choice, so
/// the two spellings ARE one value; this comparison lets the equality
/// comparator (and anything else deriving IDENTITY from the children) see them
/// as one, by reading the children through \c canonical_children instead of
/// through the emitted left/right. Nothing is reordered: the tree keeps the
/// evaluation order the DP chose, because cost -- peak in particular -- is
/// order-dependent.
///
/// The key is a function of the operands' VALUES only and is total:
///   1. a scalar operand sorts AFTER a non-scalar one, matching how
///      \c binarize already builds a `scalar * tensor` node (tensor left,
///      \c Constant right) and keeping `c * T` and `T * c` one value;
///   2. then ascending node id (\c hash::value, the operand-order-independent
///      \c EvalExpr hash) -- O(1), and the only step that runs in practice;
///   3. on a node-id tie, equal VALUES compare 0, established by the equality
///      comparator itself, so the order can never separate two operands the
///      comparator would fold (and either order is canonical for them);
///   4. and only a genuine 64-bit collision between DISTINCT values reaches
///      \c detail::deep_tie_break.
///
/// \note Steps 3 and 4 walk the operand subtrees, so a node whose two operands
///       collide on the node id costs O(subtree) here rather than O(1). That
///       is bounded by the subtree and cannot recur into itself (both steps
///       descend strictly), unlike the either-pairing match it replaces, whose
///       double recursion per child was exponential in product depth.
///
/// \return <0 if \p a belongs first, >0 if \p b does, 0 if the two are
///         interchangeable.
template <typename TreeNode>
int canonical_operand_cmp(TreeNode const &a, TreeNode const &b) {
  if (&a == &b) return 0;

  bool const a_scalar = a->is_scalar();
  bool const b_scalar = b->is_scalar();
  if (a_scalar != b_scalar) return a_scalar ? 1 : -1;

  std::size_t const ha = hash::value(*a);
  std::size_t const hb = hash::value(*b);
  if (ha != hb) return ha < hb ? -1 : 1;

  if (TreeNodeEqualityComparator<TreeNode>{}(a, b)) return 0;

  return detail::deep_tie_break(a, b);
}

/// \brief A node's two children in CANONICAL order -- the order-independent
///        VIEW that every derivation of node IDENTITY reads them through.
///
/// \details Returns the children swapped iff \c canonical_operand_cmp says the
/// right one belongs first. The node is NOT modified and nothing is cached on
/// it: the decision is one hash comparison in every case that occurs in
/// practice (see \c canonical_operand_cmp), and leaving the tree untouched is
/// the whole point -- evaluation, slicing positions, operand reads, dry-run op
/// costs and the peak sweep all keep the left/right the DP emitted.
///
/// Meaningful for a commutative (Product) node; callers use the emitted order
/// for everything else (a Sum's left child is the in-place accumulator, an
/// Adjoint's right child is a sentinel), so this is not applied there.
template <typename TreeNode>
std::pair<TreeNode const &, TreeNode const &> canonical_children(
    TreeNode const &n) {
  SEQUANT_ASSERT(!n.leaf());
  if (canonical_operand_cmp(n.left(), n.right()) > 0)
    return {n.right(), n.left()};
  return {n.left(), n.right()};
}

/// A map between (sub)tree hashes and how often they have been found
/// This is identical to SubexpressionUsageCounts except that we store node
/// pointers here (lower memory footprint but higher risk of dangling pointers)
template <typename TreeNode, bool force_hash_collisions = false>
using SubexpressionHashCollector =
    std::unordered_map<const TreeNode *, std::size_t,
                       TreeNodeHasher<TreeNode, force_hash_collisions>,
                       TreeNodeEqualityComparator<TreeNode>>;

/// A map between (sub)trees and how often they have been found
template <typename TreeNode, bool force_hash_collisions = false>
using SubexpressionUsageCounts =
    std::unordered_map<TreeNode, std::size_t,
                       TreeNodeHasher<TreeNode, force_hash_collisions>,
                       TreeNodeEqualityComparator<TreeNode>>;

/// A map between (sub)trees and the name chosen to represent the associated
/// intermediate
template <typename TreeNode, bool force_hash_collisions = false>
using SubexpressionNames =
    std::unordered_map<TreeNode, std::wstring,
                       TreeNodeHasher<TreeNode, force_hash_collisions>,
                       TreeNodeEqualityComparator<TreeNode>>;

}  // namespace sequant

#endif  // SEQUANT_EVAL_EVAL_NODE_COMPARE_HPP
