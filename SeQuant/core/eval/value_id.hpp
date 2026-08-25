#ifndef SEQUANT_EVAL_VALUE_ID_HPP
#define SEQUANT_EVAL_VALUE_ID_HPP

#include <SeQuant/core/container.hpp>
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
template <meta::eval_node Node>
[[nodiscard]] std::size_t value_id_hash(Node const& node,
                                        ValueIdColoring const* coloring) {
  if (coloring && !coloring->empty() &&
      !in_scope_batched_on_node(node, coloring->ctx_modes).empty())
    return occurrence_key(node, coloring->ctx_modes, &coloring->colors)
        .hash_value();
  return hash::value(*node);
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_VALUE_ID_HPP
