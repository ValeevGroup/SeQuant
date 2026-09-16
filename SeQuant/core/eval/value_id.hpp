#ifndef SEQUANT_EVAL_VALUE_ID_HPP
#define SEQUANT_EVAL_VALUE_ID_HPP

#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node_compare.hpp>
#include <SeQuant/core/hash.hpp>

#include <cstddef>
#include <utility>

namespace sequant::eval {

/// \brief The runtime-cache key: a forest node.
///
/// \details On the forest-descent path a cached value's identity is its
/// node's canonical hash (\c hash::value) -- exactly what \c
/// CachedValueHasher / \c CachedValueEqual key \c CacheManager's map by. The
/// batched (table-driven) executor never consults this map for values at
/// all: it keys by cell id instead (see \c cell_table.hpp), so no per-scope
/// slice coloring is threaded through here. Node-facing operations forward
/// to \c node via \c operator-> so \c CacheManager's internals retarget
/// mechanically.
template <meta::eval_node Node>
struct CachedValue {
  Node node;  //!< the (possibly shared) forest node

  /// Implicit from a bare node: every cache call site that passes a node
  /// keys the map by that node's identity.
  CachedValue(Node n) : node(std::move(n)) {}

  [[nodiscard]] auto operator->() const { return node.operator->(); }
  [[nodiscard]] Node const& operator*() const { return node; }
};

/// \brief Hasher for \c CachedValue: the node's canonical hash, i.e.
///        byte-identical to \c TreeNodeHasher.
template <meta::eval_node Node, bool force_hash_collisions = false>
struct CachedValueHasher {
  using is_transparent = void;
  [[nodiscard]] std::size_t operator()(CachedValue<Node> const& cv) const {
    return (*this)(cv.node);
  }
  /// Heterogeneous overload: probe the map with a bare node. Without it every
  /// lookup would convert the node to a CachedValue first, and that conversion
  /// copies the node -- which deep-copies its whole subtree (binary_node.hpp),
  /// so probing a deep left-leaning Sum-tree would cost O(terms^2) node
  /// allocations. The hash is the node's canonical hash either way.
  [[nodiscard]] std::size_t operator()(Node const& n) const {
    if constexpr (force_hash_collisions) return 0;
    return hash::value(*n);
  }
};

/// \brief Equality for \c CachedValue: the structural
///        \c TreeNodeEqualityComparator on the nodes.
template <meta::eval_node Node>
struct CachedValueEqual {
  using is_transparent = void;
  [[nodiscard]] bool operator()(CachedValue<Node> const& a,
                                CachedValue<Node> const& b) const {
    return (*this)(a.node, b.node);
  }
  /// Heterogeneous overloads: compare a stored key against a bare node
  /// probe, without materializing a CachedValue (which would deep-copy the
  /// probed node's subtree -- see CachedValueHasher).
  [[nodiscard]] bool operator()(CachedValue<Node> const& a,
                                Node const& b) const {
    return (*this)(a.node, b);
  }
  [[nodiscard]] bool operator()(Node const& a,
                                CachedValue<Node> const& b) const {
    return (*this)(a, b.node);
  }
  [[nodiscard]] bool operator()(Node const& a, Node const& b) const {
    return TreeNodeEqualityComparator<Node>{}(a, b);
  }
};

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_VALUE_ID_HPP
