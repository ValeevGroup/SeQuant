#ifndef SEQUANT_EVAL_VALUE_ID_HPP
#define SEQUANT_EVAL_VALUE_ID_HPP

#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node_compare.hpp>
#include <SeQuant/core/hash.hpp>

#include <cstddef>
#include <utility>

namespace sequant::eval {

/// \brief The runtime-cache key: a non-owning handle to a forest node.
///
/// \details On the forest-descent path a cached value's identity is its
/// node's canonical hash (\c hash::value) -- exactly what \c
/// CachedValueHasher / \c CachedValueEqual key \c CacheManager's map by. The
/// batched (table-driven) executor never consults this map for values at
/// all: it keys by cell id instead (see \c cell_table.hpp), so no per-scope
/// slice coloring is threaded through here. Node-facing operations forward
/// to the node via \c operator-> so \c CacheManager's internals retarget
/// mechanically.
///
/// The key does NOT own its node. Copying an eval node deep-copies its whole
/// subtree (binary_node.hpp), so an owning key made registering every internal
/// node of an equation's left-leaning Sum-tree (one spine node per summand,
/// each subtree larger than the last) cost O(terms^2) node copies. A key is
/// a pointer instead, so registration is O(1) per node at any \c min_repeats.
/// Whoever builds the \c CacheManager guarantees that every keyed node
/// outlives it and stays at its address (see \c cache_manager()). Node
/// identity (hash, structure, connectivity graph) is fixed at construction;
/// the batch annotations a node may receive later do not enter it.
template <meta::eval_node Node>
struct CachedValue {
  Node const* node;  //!< the forest node; not owned

  /// Keys the map by @p n's identity; @p n must outlive the key.
  explicit CachedValue(Node const& n) noexcept : node(&n) {}
  /// A temporary would dangle.
  CachedValue(Node const&&) = delete;

  /// The node's payload (its \c EvalExpr), as on a bare node.
  [[nodiscard]] auto operator->() const { return node->operator->(); }
  [[nodiscard]] Node const& operator*() const noexcept { return *node; }
};

/// \brief Hasher for \c CachedValue: the node's canonical hash, i.e.
///        byte-identical to \c TreeNodeHasher.
template <meta::eval_node Node, bool force_hash_collisions = false>
struct CachedValueHasher {
  using is_transparent = void;
  [[nodiscard]] std::size_t operator()(CachedValue<Node> const& cv) const {
    return (*this)(*cv);
  }
  /// Heterogeneous overload: probe the map with a bare node, which need not be
  /// one of the keyed nodes (e.g. a structurally equal node of another tree).
  /// The hash is the node's canonical hash either way.
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
    return (*this)(*a, *b);
  }
  /// Heterogeneous overloads: compare a stored key against a bare node
  /// probe.
  [[nodiscard]] bool operator()(CachedValue<Node> const& a,
                                Node const& b) const {
    return (*this)(*a, b);
  }
  [[nodiscard]] bool operator()(Node const& a,
                                CachedValue<Node> const& b) const {
    return (*this)(a, *b);
  }
  [[nodiscard]] bool operator()(Node const& a, Node const& b) const {
    return TreeNodeEqualityComparator<Node>{}(a, b);
  }
};

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_VALUE_ID_HPP
