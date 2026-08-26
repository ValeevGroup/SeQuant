#ifndef SEQUANT_EVAL_VALUE_ID_HPP
#define SEQUANT_EVAL_VALUE_ID_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/eval_node_compare.hpp>
#include <SeQuant/core/eval/occurrence_key.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor_network/typedefs.hpp>

#include <cstddef>
#include <utility>

namespace sequant::eval {

/// \brief Pillar 1 (slice-colored value identity): the per-scope coloring the
///        value-id keying reads.
///
/// \details `ctx_modes` are the sliced modes in scope (the value's home-sliced
/// modes for a value-id; a use-site's sliced modes for an occurrence-id), and
/// `colors` maps each to its DAG-scope DEPTH. Both are in the SAME canonical
/// index-frame as the nodes this coloring keys (value frame for a value-id) --
/// built directly from `OrderedSchedule::home_mode_depth`, whose keys are the
/// value's own `carried` indices, so no cross-frame match is ever performed
/// (the loop-open lesson). An empty coloring is the unsliced degenerate case:
/// keying reduces byte-for-byte to the plain node-id (`hash::value`).
struct ValueIdColoring {
  container::svector<Index> ctx_modes;        //!< sliced modes (frame keys)
  tensor_network::NamedIndexColorMap colors;  //!< mode -> DAG-depth color

  [[nodiscard]] bool empty() const noexcept { return ctx_modes.empty(); }
};

/// \brief Build a \c ValueIdColoring from an \c
/// OrderedSchedule::home_mode_depth
///        entry (a value's home-sliced mode -> depth pairs, value frame).
[[nodiscard]] inline ValueIdColoring value_id_coloring(
    container::svector<std::pair<Index, int>> const& mode_depth) {
  ValueIdColoring c;
  for (auto const& [m, d] : mode_depth) {
    c.ctx_modes.push_back(m);
    c.colors.emplace(m, static_cast<std::size_t>(d));
  }
  return c;
}

/// \brief The value-id HASH of \p node under \p coloring.
///
/// \details When \p node carries an in-scope sliced mode (so the value is
/// genuinely sliced in this scope), the id is the LOOP-COLORED
/// \c occurrence_key hash -- distinguishing e.g. `I(i,_)` from `I(_,i)` for a
/// non-symmetric `I`, folding symmetric ones. Otherwise (no coloring, empty
/// coloring, or a node with no in-scope sliced slot -- an unsliced value in a
/// sub-top scope, and every value in the main cache) the id is the plain
/// node-id \c hash::value, **byte-identical to \c TreeNodeHasher today**. The
/// equality partner is the graph comparison (\c RouterKeyEqual) for the colored
/// path and the structural \c TreeNodeEqualityComparator for the fallback --
/// see \c eval_node_compare.hpp.
/// \brief Whether \p node is genuinely sliced under \p coloring in this scope:
/// a
///        non-empty coloring that actually colors an in-scope slot of the node.
///
/// \details The negative case (no coloring / empty coloring / a node with no
/// in-scope sliced slot) is exactly where the value-id degenerates to the plain
/// node-id -- so this predicate is the single fork shared by \c value_id_hash
/// (which id to compute) and \c CachedValueEqual (which comparison to run).
template <meta::eval_node Node>
[[nodiscard]] bool value_is_sliced(Node const& node,
                                   ValueIdColoring const* coloring) {
  return coloring && !coloring->empty() &&
         !in_scope_batched_on_node(node, coloring->ctx_modes).empty();
}

template <meta::eval_node Node>
[[nodiscard]] std::size_t value_id_hash(Node const& node,
                                        ValueIdColoring const* coloring) {
  if (value_is_sliced(node, coloring))
    return occurrence_key(node, coloring->ctx_modes, &coloring->colors)
        .hash_value();
  return hash::value(*node);
}

/// \brief The value-keyed runtime-cache / remat-router key: a forest node
/// paired
///        with the home-slice coloring recorded for its value.
///
/// \details Its identity is the home-slice-colored value-id (\c
/// CachedValueHasher / \c CachedValueEqual), so two values of ONE node that
/// slice different slots are distinct keys, while an unsliced value (empty
/// coloring) is byte-identical to keying by the node itself. This is what lets
/// the cache tell apart split values that share a single forest node (\c
/// build_value_node_map is hash-keyed): the discriminating coloring is recorded
/// here beside the node, not re-derived from the shared node. Node-facing
/// operations forward to \c node via
/// \c operator-> so \c CacheManager's internals retarget mechanically.
template <meta::eval_node Node>
struct CachedValue {
  Node node;                 //!< the (possibly shared) forest node
  ValueIdColoring coloring;  //!< recorded home-slice coloring; empty => full

  [[nodiscard]] auto operator->() const { return node.operator->(); }
  [[nodiscard]] Node const& operator*() const { return node; }
};

/// \brief Hasher for \c CachedValue: the home-slice-colored value-id. Empty
///        coloring => \c hash::value(*node), byte-identical to \c
///        TreeNodeHasher.
template <meta::eval_node Node>
struct CachedValueHasher {
  using is_transparent = void;
  [[nodiscard]] std::size_t operator()(CachedValue<Node> const& cv) const {
    return value_id_hash(cv.node, &cv.coloring);
  }
};

/// \brief Equality for \c CachedValue. When both operands are genuinely sliced,
///        compares the colored \c occurrence_key graphs (symmetric slicings
///        fold, asymmetric ones separate); otherwise falls back to the
///        structural
///        \c TreeNodeEqualityComparator on the nodes -- byte-identical to
///        today. A slice-relevance mismatch is unequal (guards a hash collision
///        between a sliced and an unsliced value).
template <meta::eval_node Node>
struct CachedValueEqual {
  using is_transparent = void;
  [[nodiscard]] bool operator()(CachedValue<Node> const& a,
                                CachedValue<Node> const& b) const {
    bool const as = value_is_sliced(a.node, &a.coloring);
    bool const bs = value_is_sliced(b.node, &b.coloring);
    if (as != bs) return false;
    if (as)
      return RouterKeyEqual{}(
          occurrence_key(a.node, a.coloring.ctx_modes, &a.coloring.colors),
          occurrence_key(b.node, b.coloring.ctx_modes, &b.coloring.colors));
    return TreeNodeEqualityComparator<Node>{}(a.node, b.node);
  }
};

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_VALUE_ID_HPP
