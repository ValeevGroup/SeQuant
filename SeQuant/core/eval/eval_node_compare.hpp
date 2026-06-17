#ifndef SEQUANT_EVAL_EVAL_NODE_COMPARE_HPP
#define SEQUANT_EVAL_EVAL_NODE_COMPARE_HPP

#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/tensor.hpp>
#include <SeQuant/external/bliss/graph.hh>

#include <cstddef>
#include <unordered_map>

namespace sequant {

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

      TensorBlockEqualComparator cmp;
      if (!cmp(lhs_tensor, rhs_tensor)) {
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
