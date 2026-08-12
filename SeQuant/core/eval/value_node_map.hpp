#ifndef SEQUANT_EVAL_VALUE_NODE_MAP_HPP
#define SEQUANT_EVAL_VALUE_NODE_MAP_HPP

#include <SeQuant/core/eval/eval_expr.hpp>

#include <cstddef>
#include <unordered_map>

namespace sequant::eval {

///
/// \brief The value_id -> forest-node bridge (design integration point 1).
///
/// \details \c ScopeNode::homed_values are \c ValueCell::value_id's; a value's
/// eval node is recovered through the \c ValueCell::hash it carries -- exactly
/// the \c EvalExpr::hash_value() identity \c CacheManager dedups by. This maps
/// every distinct forest node by that hash, so a \c value_id resolves as
/// `map[rich.cells[value_id].hash]`. Under perfect CSE many occurrences share
/// one hash; a single representative node (the first visited) is kept, which is
/// all a homed value needs (every occurrence is the same value). Pure lookup
/// construction -- no execution.
///
/// \note Factored into its own header (rather than living inline in
/// `scope_executor.hpp`, where it was originally introduced) so that
/// `ordered_executor.hpp` -- which needs this same bridge to resolve a
/// `BuildStep::value_id` to a forest node -- can use it without a circular
/// include: `scope_executor.hpp`'s coexistence driver dispatches into
/// `ordered_executor.hpp`'s `evaluate_ordered_schedule`, so the reverse
/// direction (`ordered_executor.hpp` including `scope_executor.hpp` for just
/// this helper) would cycle. Both headers include this one instead; the
/// symbol stays `sequant::eval::build_value_node_map`, unchanged for every
/// existing caller.
///
template <meta::eval_node_range R>
[[nodiscard]] std::unordered_map<std::size_t, std::ranges::range_value_t<R>>
build_value_node_map(R const& forest) {
  using node_t = std::ranges::range_value_t<R>;
  std::unordered_map<std::size_t, node_t> out;
  auto visit = [&out](auto&& self, node_t const& n) -> void {
    out.emplace(n->hash_value(), n);
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
