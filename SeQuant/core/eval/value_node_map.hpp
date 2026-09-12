#ifndef SEQUANT_EVAL_VALUE_NODE_MAP_HPP
#define SEQUANT_EVAL_VALUE_NODE_MAP_HPP

#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/lifetime_mask.hpp>

#include <cstddef>
#include <unordered_map>

namespace sequant::eval {

///
/// \brief The value_id -> forest-node bridge (design integration point 1).
///
/// \details A schedule addresses a value by its \c ValueCell::value_id (e.g.
/// \c BuildStep::value_id); the value's eval node is recovered through the
/// \c ValueCell::hash that cell carries -- exactly the
/// \c EvalExpr::hash_value() identity \c CacheManager dedups by. This maps
/// every distinct forest node by that hash, so a \c value_id resolves as
/// `map[rich.cells[value_id].hash]`. Under perfect CSE many occurrences share
/// one hash; a single representative node (the first visited) is kept, which is
/// all a homed value needs (every occurrence is the same value). Pure lookup
/// construction -- no execution.
///
/// \note Lives in its own header so that `ordered_executor.hpp` -- which
/// needs this bridge to resolve a `BuildStep::value_id` to a forest node --
/// can use it without depending on an executor header.
///
template <meta::eval_node_range R>
[[nodiscard]] std::unordered_map<std::size_t, std::ranges::range_value_t<R>>
build_value_node_map(R const& forest) {
  using node_t = std::ranges::range_value_t<R>;
  std::unordered_map<std::size_t, node_t> out;
  // Two passes. Value keys first, over the WHOLE forest: a whole value's key
  // IS its node hash, and that entry must be the whole occurrence -- an
  // earlier sliced occurrence of the same node (key != hash) must not claim
  // the hash slot first. Then node hashes, only where no value claimed them
  // (a monitor's op hash, a leaf: any node of that hash).
  auto visit = [&out](auto&& self, node_t const& n, bool keys) -> void {
    if (keys)
      out.emplace(value_key_of(n), n);
    else
      out.emplace(n->hash_value(), n);
    if (!n.leaf()) {
      self(self, n.left(), keys);
      self(self, n.right(), keys);
    }
  };
  for (auto const& t : forest) visit(visit, t, true);
  for (auto const& t : forest) visit(visit, t, false);
  return out;
}

/// \brief Like \c build_value_node_map but keyed by VALUE id (\c
/// value_key_of: node id + home-sliced positions), one representative node
/// per value -- the map the ordered executor resolves a \c value_id through,
/// so a value is always built from one of its own occurrences' nodes (a node
/// of the same hash home-sliced on other positions is another value).
template <meta::eval_node_range R>
[[nodiscard]] std::unordered_map<std::size_t, std::ranges::range_value_t<R>>
build_value_key_node_map(R const& forest) {
  using node_t = std::ranges::range_value_t<R>;
  std::unordered_map<std::size_t, node_t> out;
  auto visit = [&out](auto&& self, node_t const& n) -> void {
    out.emplace(value_key_of(n), n);
    if (!n.leaf()) {
      self(self, n.left());
      self(self, n.right());
    }
  };
  for (auto const& t : forest) visit(visit, t);
  return out;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_VALUE_NODE_MAP_HPP
