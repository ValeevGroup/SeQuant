#ifndef SEQUANT_EVAL_EVAL_NODE_COMPARE_HPP
#define SEQUANT_EVAL_EVAL_NODE_COMPARE_HPP

#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/tensor.hpp>
#include <SeQuant/external/bliss/graph.hh>

#include <compare>
#include <cstddef>
#include <cstdlib>
#include <unordered_map>
#include <utility>

namespace sequant {

namespace detail {
/// @return true (the default) unless SEQUANT_EVAL_LAX_LAYOUT is set: the
/// result mode layout is part of eval-node identity (see the use site in
/// TreeNodeEqualityComparator)
inline bool strict_layout_identity() {
  static const bool lax = std::getenv("SEQUANT_EVAL_LAX_LAYOUT") != nullptr;
  return !lax;
}
}  // namespace detail

/// Functor to compute the hash of a given (evaluation) tree node.
///
/// This is the node-id hasher (\c hash::value): batching-blind, used by
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

/// \brief 3-way comparison establishing the canonical order of a commutative
///        (Product) node's two operands.
///
/// \details A contraction is commutative, so the binarizer emits its two
/// operands in whichever order the single-term DP's contraction sequence
/// yields, i.e. in the order in which the pair lands on the DP's stack. The
/// node id and the canonical connectivity graph both fold that choice, so the
/// two spellings are one value; this comparison lets the equality comparator
/// (and anything else deriving identity from the children) see them as one, by
/// reading the children through \c canonical_children instead of through the
/// emitted left/right. Nothing is reordered: the tree keeps the evaluation
/// order the DP chose, because cost -- peak in particular -- is order
/// dependent.
///
/// The order is a function of the operands' values only, and is decided in
/// O(1), reading nothing but the two nodes' own already-computed data, never
/// their subtrees:
///   1. a scalar operand sorts after a non-scalar one, matching how
///      \c binarize builds a `scalar * tensor` node (tensor left, \c Constant
///      right) and keeping `c * T` and `T * c` one value;
///   2. then ascending node id (\c hash::value, the operand-order-independent
///      \c EvalExpr hash);
///   3. operands agreeing on both are equivalent, and are read in the order in
///      which they were emitted.
///
/// Case 3 is exactly right for the case that occurs: two operands that are one
/// value have one node id, either order is canonical for them, and an ordered
/// comparison of two swapped spellings succeeds whichever way each is read.
/// The only other way to reach it is a genuine 64-bit hash collision between
/// distinct values, and there the cost is a missed fold -- the two spellings
/// stay two cache entries -- never a wrong fold, because
/// \c TreeNodeEqualityComparator still compares both subtrees in full and
/// rejects a mismatch. Resolving such a collision instead would mean comparing
/// the operand subtrees here, and a full comparison per visit turns the
/// comparator quadratic on nested hash-equal operands (a `t2 * t2`
/// self-contraction is enough): T(n) = 4 T(n/2). Keeping this O(1) is what
/// makes the comparator provably linear.
///
/// \param a one operand of a commutative node
/// \param b the other operand of that node
/// \return \c std::strong_ordering::less if \p a belongs first,
///         \c std::strong_ordering::greater if \p b does, and
///         \c std::strong_ordering::equal if the two are interchangeable
///         (read them in the order in which they were emitted).
template <typename TreeNode>
std::strong_ordering canonical_operand_cmp(TreeNode const &a,
                                           TreeNode const &b) {
  if (&a == &b) return std::strong_ordering::equal;

  bool const a_scalar = a->is_scalar();
  bool const b_scalar = b->is_scalar();
  if (a_scalar != b_scalar)
    return a_scalar ? std::strong_ordering::greater
                    : std::strong_ordering::less;

  std::size_t const ha = hash::value(*a);
  std::size_t const hb = hash::value(*b);
  if (ha != hb) return ha <=> hb;

  return std::strong_ordering::equal;
}

/// \brief A node's two children in canonical order -- the order-independent
///        view that every derivation of node identity reads them through.
///
/// \details Returns the children swapped iff \c canonical_operand_cmp says the
/// right one belongs first. The node is left alone and nothing is cached on
/// it: the decision is one comparison of two already-computed hashes, and
/// leaving the tree untouched is the whole point -- evaluation, slicing
/// positions, operand reads, dry-run op costs and the peak sweep all keep the
/// left/right the DP emitted.
///
/// Meaningful for a commutative (Product) node; callers use the emitted order
/// for everything else (a Sum's left child is the in-place accumulator, an
/// Adjoint's right child is a sentinel), so this is not applied there.
///
/// \param n a non-leaf node
/// \pre \p n is not a leaf
/// \return \p n's children, the canonically first one first
/// \note The returned pair holds references into \p n, so \p n must outlive
///       it; binding a temporary is rejected by the deleted overload below.
template <typename TreeNode>
std::pair<TreeNode const &, TreeNode const &> canonical_children(
    TreeNode const &n) {
  SEQUANT_ASSERT(!n.leaf());
  if (std::is_gt(canonical_operand_cmp(n.left(), n.right())))
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
    // stay recursive; only the deep spine is unwound into this loop.
    const TreeNode *lhsp = &lhs_in;
    const TreeNode *rhsp = &rhs_in;
    while (true) {
      const TreeNode &lhs = *lhsp;
      const TreeNode &rhs = *rhsp;

      if (lhs.leaf() != rhs.leaf()) {
        return false;
      }

      // Cheapest discriminator first: hash::value is the EvalExpr's cached
      // node id (a load) and size() the subtree's cached node count (also a
      // load, see FullBinaryNode::size_). Both are necessary conditions, so
      // the order is semantically free; the node id is the more selective of
      // the two.
      if (hash::value(*lhs) != hash::value(*rhs)) {
        return false;
      }

      if (lhs.size() != rhs.size()) {
        return false;
      }

      if (lhs->type_id() != rhs->type_id()) {
        return false;
      }

      // NB the canonicalization transform (phase / conjugation / bra-ket
      // swap) is deliberately NOT part of the identity: a slot holds the
      // canonical value and every consumer applies its own transform on
      // retrieval (apply_canon_phase), so +T / -T / T* share one slot -- and
      // the hash-keyed value maps of the ordered (DAG) executor must agree
      // with this comparator on what is one value.

      // The two nodes must lay their result modes out the same way (default;
      // SEQUANT_EVAL_LAX_LAYOUT=1 opts out).
      //
      // The hash and the connectivity comparison identify nodes across index
      // renamings and across bra<->ket orientation -- that is what makes a
      // subexpression shareable -- and CanonTransform carries the residual
      // phase / conjugation / bra-ket swap, but nothing carries a permutation
      // of the result modes. Two same-space external indices of an isomorphic
      // network can be ordered either way, since bliss breaks an automorphic
      // orbit by input vertex order, and a cached buffer served under the
      // other ordering is a transposed value (measured: a residual block's
      // product node C+.(g.C) laid out (i_1,i_2;..) served to its twin laid
      // out (i_2,i_1;..) put a PNS-MP1 energy 7 % off).
      if (detail::strict_layout_identity() &&
          lhs->layout_fingerprint() != rhs->layout_fingerprint()) {
        return false;
      }

      if (lhs->is_constant() || lhs->is_variable() || lhs->is_power()) {
        if (*lhs->expr() != *rhs->expr()) {
          return false;
        }
      } else if (lhs->is_tensor() && !lhs->has_connectivity_graph()) {
        // Tensor nodes that carry a canonical connectivity graph
        // (tensor-network intermediates and proto-indexed leaves) are compared
        // exactly by that graph in the connectivity check below: it is the
        // same canonical colored graph the eval-node hash is derived from, its
        // 3-way cmp is a complete network identity, and it folds bra<->ket
        // orientation (the orientation is canonicalized into the graph). Block
        // comparison here would instead compare on the stored bra/ket slot
        // order -- a partial, orientation-sensitive signature -- and wrongly
        // separate e.g. a CSV/PNO coefficient C{a;m} from its equivalent
        // C{m;a} (m an aux/PAO index), or the two g.C intermediates obtained
        // by transforming a real density-fitting factor g(m,m,K) on its bra vs
        // its ket leg. Hence only graph-less tensor nodes -- protoindex-free
        // leaves (block-canonicalized in place at construction) and
        // scalar*tensor results -- are compared here by block.
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
      // operand order -- so the child comparison folds it too, else two
      // swapped occurrences of the same value are split into two cache entries
      // (and built twice).
      //
      // It folds it by comparing the children through a canonical view
      // (canonical_children) rather than by trying both pairings: the view
      // picks which child comes first from the children's values alone -- in
      // O(1), reading only the two nodes' own data -- so the two spellings
      // present the same ordered pair here and one recursion per child
      // suffices. The whole comparison is therefore linear in the tree.
      // Evaluation order is untouched: the node keeps the left/right the DP
      // emitted (cost, and peak in particular, is order-dependent); only
      // identity reads them through the view.
      //
      // The recursive child comparison is not redundant with the graph check
      // above: the connectivity graph encodes only the immediate two factors'
      // connectivity, not each factor's recursive build, so two products with
      // the same immediate graph but different sub-values are told apart here.
      // Product operand subtrees are bounded in depth (a contraction of a
      // fixed set of factors), so they stay recursive.
      if (lhs->op_type() && *lhs->op_type() == EvalOp::Product) {
        const auto &[lfirst, lsecond] = canonical_children(lhs);
        const auto &[rfirst, rsecond] = canonical_children(rhs);
        return (*this)(lfirst, rfirst) && (*this)(lsecond, rsecond);
      }

      // Non-Product internal node (Sum, scalar*tensor product, adjoint): its
      // left/right assignment is canonical (e.g. the in-place Sum tree is
      // left-folded; an Adjoint's right child is a sentinel), so operand order
      // carries meaning and the children are compared in order. The right
      // child (a single summand / the scalar factor / the adjoint sentinel) is
      // bounded in depth and compared recursively; the left child is the deep
      // spine, so rather than recurse into it the loop starts over with it as
      // the new (lhs, rhs), unwinding the spine iteratively.
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
