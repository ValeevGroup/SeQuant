#ifndef SEQUANT_EVAL_EVAL_NODE_COMPARE_HPP
#define SEQUANT_EVAL_EVAL_NODE_COMPARE_HPP

#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/tensor.hpp>
#include <SeQuant/external/bliss/graph.hh>

#include <cstddef>
#include <functional>
#include <optional>
#include <unordered_map>

namespace sequant {

/// Functor to compute the hash of a given (evaluation) tree node.
///
/// Pillar 1 (slice-colored value identity): an optional \c id_override lets a
/// sub-top cache key by the home-slice-colored VALUE-id instead of the plain
/// node-id. When set, it returns the value-id hash (\c
/// eval::value_id_hash) -- which is byte-identical to \c hash::value for an
/// unsliced value, so the override is safe to install everywhere; the main
/// (top-level) cache leaves it null for zero hot-path cost. Injected from
/// \c cache_manager.hpp so this low-level header does not depend on the
/// tensor-network coloring machinery.
template <typename TreeNode, bool force_hash_collisions = false>
struct TreeNodeHasher {
  /// Trait used by the C++ STL allowing heterogenous lookups
  using is_transparent = void;

  std::function<std::size_t(const TreeNode &)> id_override{};

  std::size_t operator()(const TreeNode *node) const { return (*this)(*node); }

  std::size_t operator()(const TreeNode &node) const {
    if constexpr (force_hash_collisions) {
      return 0;
    }
    if (id_override) return id_override(node);
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

  /// Pillar 1: optional slice-colored equality. When set, it is consulted FIRST
  /// on the top-level (lhs, rhs) compare: a returned value is the identity
  /// answer (two values sliced on different slots compare UNEQUAL even though
  /// their node structure is identical; symmetric slicings fold), and \c
  /// nullopt means "not slice-relevant -- use the structural comparison below"
  /// (both values unsliced), keeping the null/unsliced path byte-identical.
  /// Injected from \c cache_manager.hpp.
  std::function<std::optional<bool>(const TreeNode &, const TreeNode &)>
      colored_eq_override{};

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
    // Pillar 1: a slice-colored identity decision (if any) wins on the
    // top-level compare -- it distinguishes values sliced on different slots
    // that are structurally identical. nullopt => not slice-relevant, fall
    // through to the byte-identical structural comparison below.
    if (colored_eq_override)
      if (auto const c = colored_eq_override(lhs_in, rhs_in); c.has_value())
        return *c;
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

      // A Product (contraction) is commutative: the binarizer may emit the same
      // contraction as (X,Y) in one term and (Y,X) in another, and the hash and
      // canonical connectivity graph both fold that operand order -- so the
      // child comparison must fold it too, else two swapped occurrences of the
      // same value are wrongly split into two cache entries (built twice).
      // Match the children as an UNORDERED pair. The recursive child comparison
      // is still REQUIRED and is NOT redundant with the graph check above: the
      // connectivity graph encodes only the immediate two factors'
      // connectivity, not each factor's recursive build, so two products with
      // the same immediate graph but different sub-values must still be told
      // apart -- and are, because the unordered match fails when no child
      // pairing is equal. Product operand subtrees are bounded in depth (a
      // contraction of a fixed set of factors), so they stay recursive.
      if (lhs->op_type() && *lhs->op_type() == EvalOp::Product) {
        bool const in_order = (*this)(lhs.left(), rhs.left()) &&
                              (*this)(lhs.right(), rhs.right());
        bool const swapped = (*this)(lhs.left(), rhs.right()) &&
                             (*this)(lhs.right(), rhs.left());
        if (!in_order && !swapped) {
          return false;
        }
        return true;
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
