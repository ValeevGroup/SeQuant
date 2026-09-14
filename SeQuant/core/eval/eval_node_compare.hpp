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
/// The key is a function of the operands' VALUES only, and is O(1) -- it reads
/// nothing but the two nodes' own already-computed data, never their subtrees:
///   1. a scalar operand sorts AFTER a non-scalar one, matching how
///      \c binarize already builds a `scalar * tensor` node (tensor left,
///      \c Constant right) and keeping `c * T` and `T * c` one value;
///   2. then ascending node id (\c hash::value, the operand-order-independent
///      \c EvalExpr hash);
///   3. a tie resolves to 0, i.e. to the EMITTED order.
///
/// Step 3 is exactly right for the case that occurs: two operands that are ONE
/// value have one node id, either order is canonical for them, and an ordered
/// comparison of two swapped spellings succeeds whichever way each is read.
/// The only other way to reach it is a genuine 64-bit hash collision between
/// DISTINCT values, and there the cost is a MISSED fold (the two spellings
/// stay two cache entries, as they were before any folding existed) -- never a
/// wrong one, because \c TreeNodeEqualityComparator still compares both
/// subtrees in full and rejects a mismatch. Resolving such a collision instead
/// would mean comparing the operand subtrees here, and a full comparison per
/// visit turns the comparator quadratic on nested hash-equal operands (a
/// `t2 * t2` self-contraction is enough): T(n) = 4 T(n/2). Keeping this O(1)
/// is what makes the comparator provably linear.
///
/// \return <0 if \p a belongs first, >0 if \p b does, 0 if the two are
///         interchangeable (read them in the order they were emitted).
template <typename TreeNode>
int canonical_operand_cmp(TreeNode const &a, TreeNode const &b) {
  if (&a == &b) return 0;

  bool const a_scalar = a->is_scalar();
  bool const b_scalar = b->is_scalar();
  if (a_scalar != b_scalar) return a_scalar ? 1 : -1;

  std::size_t const ha = hash::value(*a);
  std::size_t const hb = hash::value(*b);
  if (ha != hb) return ha < hb ? -1 : 1;

  return 0;
}

/// \brief A node's two children in CANONICAL order -- the order-independent
///        VIEW that every derivation of node IDENTITY reads them through.
///
/// \details Returns the children swapped iff \c canonical_operand_cmp says the
/// right one belongs first. The node is NOT modified and nothing is cached on
/// it: the decision is one comparison of two already-computed hashes, and
/// leaving the tree untouched is the whole point -- evaluation, slicing
/// positions, operand reads, dry-run op costs and the peak sweep all keep the
/// left/right the DP emitted.
///
/// Meaningful for a commutative (Product) node; callers use the emitted order
/// for everything else (a Sum's left child is the in-place accumulator, an
/// Adjoint's right child is a sentinel), so this is not applied there.
///
/// \note The returned pair holds REFERENCES into \p n, so \p n must outlive
///       it; binding a temporary is rejected by the deleted overload below.
template <typename TreeNode>
std::pair<TreeNode const &, TreeNode const &> canonical_children(
    TreeNode const &n) {
  SEQUANT_ASSERT(!n.leaf());
  if (canonical_operand_cmp(n.left(), n.right()) > 0)
    return {n.right(), n.left()};
  return {n.left(), n.right()};
}

/// The result aliases \p n's children, so a temporary would dangle. (A
/// `TreeNode const&&` parameter is not a forwarding reference, so this
/// overload takes rvalues only and leaves lvalue calls to the one above.)
template <typename TreeNode>
std::pair<TreeNode const &, TreeNode const &> canonical_children(
    TreeNode const &&n) = delete;

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
    // comparison -- same per-node checks, same child logic.
    const TreeNode *lhsp = &lhs_in;
    const TreeNode *rhsp = &rhs_in;
    while (true) {
      const TreeNode &lhs = *lhsp;
      const TreeNode &rhs = *rhsp;

      if (lhs.leaf() != rhs.leaf()) {
        return false;
      }

      // Cheapest discriminator first: hash::value is the EvalExpr's cached
      // node id (a load), size() the subtree's cached node count (also a load
      // now, see FullBinaryNode::size_ -- it used to WALK the subtree, which
      // made every probe of a structurally-keyed map over a deep tree
      // quadratic in tree size). Both are necessary conditions, so the order
      // is semantically free; the node id is the more selective of the two.
      if (hash::value(*lhs) != hash::value(*rhs)) {
        return false;
      }

      if (lhs.size() != rhs.size()) {
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
      // picks which child comes first from the children's VALUES alone -- in
      // O(1), reading only the two nodes' own data -- so the two spellings
      // present the same ordered pair here and ONE recursion per child
      // suffices. The whole comparison is therefore LINEAR in the tree.
      // Matching the children as an unordered pair would fold them too, but it
      // recurses twice per child and so is exponential in product depth.
      // Evaluation order is untouched: the node keeps the left/right the DP
      // emitted (cost, and peak in particular, is order-dependent), only
      // IDENTITY reads them through the view.
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
