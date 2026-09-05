#ifndef SEQUANT_EVAL_EVAL_NODE_COMPARE_HPP
#define SEQUANT_EVAL_EVAL_NODE_COMPARE_HPP

#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/tensor.hpp>
#include <SeQuant/external/bliss/graph.hh>

#include <cstddef>
#include <unordered_map>

#include <atomic>

#include <cstdlib>

#include <iostream>

namespace sequant {

namespace detail {
/// @return true if SEQUANT_EVAL_STRICT_LAYOUT is set: make the result mode
/// layout part of eval-node identity (see the use site in
/// TreeNodeEqualityComparator)
inline bool strict_layout_identity() {
  static const bool on = std::getenv("SEQUANT_EVAL_STRICT_LAYOUT") != nullptr;
  return on;
}
}  // namespace detail

/// Functor to compute the hash of a given (evaluation) tree node
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

  bool operator()(const TreeNode &lhs, const TreeNode &rhs) const {
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

    // canon_phase (+1/-1) is part of the node's value identity: two nodes with
    // the same canonical graph/leaf but opposite antisymmetric-reorder parity
    // evaluate to negatives of each other (+T vs -T) and must not share a CSE
    // cache slot. It is already folded into hash_value() (so cross-phase pairs
    // normally hash apart and fail the check above), which is a no-op for real
    // closed-shell paths where every phase is +1; this guards the residual case
    // of a hash collision, mirroring how the graph is both hashed and compared.
    if (lhs->canon_phase() != rhs->canon_phase()) {
      return false;
    }

    // SEQUANT_EVAL_STRICT_LAYOUT=1 (diagnostic, OFF by default): additionally
    // require the two nodes to lay their result modes out the same way.
    //
    // The hash and the connectivity comparison identify nodes across index
    // renamings and across bra<->ket orientation -- that is what makes a
    // subexpression shareable -- and CanonTransform carries the residual
    // phase / conjugation / bra-ket swap, but nothing carries a PERMUTATION of
    // the result modes. Two same-space external indices of an isomorphic
    // network can be ordered either way, since bliss breaks an automorphic
    // orbit by input vertex order.
    //
    // Measured (h2o tpns=0 PNS-CCD, 2026-09-03): enforcing this separates 52
    // node pairs on the certified path and 59 on the CSV-then-DF one, changes
    // NEITHER energy, and costs ~35 % per iteration (HSeOH PNS-CCD 3.4 ->
    // 4.7 s). So the orderings it separates are value-compatible -- an
    // automorphic exchange of two externals is a genuine symmetry of the
    // network -- and the guard stays off. It remains available as a bisection
    // tool for layout-suspect wrong numbers.
    if (detail::strict_layout_identity() &&
        lhs->layout_fingerprint() != rhs->layout_fingerprint()) {
      return false;
    }

    if (lhs->is_constant() || lhs->is_variable() || lhs->is_power()) {
      if (*lhs->expr() != *rhs->expr()) {
        return false;
      }
    } else if (lhs->is_tensor() && !lhs->has_connectivity_graph()) {
      // Tensor nodes that carry a canonical connectivity graph (tensor-network
      // intermediates and proto-indexed leaves) are compared exactly by that
      // graph in the connectivity check below: it is the same canonical colored
      // graph the eval-node hash is derived from, its 3-way cmp is a complete
      // network identity, and it already folds bra<->ket orientation (the
      // orientation was canonicalized into the graph). Block comparison here
      // would instead compare on the stored bra/ket slot order — a partial,
      // orientation-sensitive signature — and wrongly separate e.g. a CSV/PNO
      // coefficient C{a;μ̃} from its equivalent C{μ̃;a}, or the two g·C
      // intermediates from transforming a real DF factor g(μ̃,μ̃,Κ) on its bra
      // vs its ket leg. So only graph-less tensor nodes — protoindex-free
      // leaves (block-canonicalized in place at construction) and scalar*tensor
      // results — are compared here by block.
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

    if (!lhs.leaf()) {
      // Note: We're assuming that the assignment of a subtree into left and
      // right is made consistently (canonical) and hence we don't compare
      // left with right
      if (!(*this)(lhs.left(), rhs.left())) {
        return false;
      }
      if (!(*this)(lhs.right(), rhs.right())) {
        return false;
      }
    }

    return true;
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
