
#ifndef SEQUANT_ANALYSIS_ALGORITHM_HPP
#define SEQUANT_ANALYSIS_ALGORITHM_HPP

#include "imed.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_node.hpp>

#include <filesystem>

namespace sequant {
///
/// Maps a node position to a set of connected node positions.
///
using adjacency_map = container::map<size_t, container::set<size_t>>;

///
/// \param node An EvalNode.
/// \return Returns the storage requirement as AsyCost
///
template <typename NodeT, typename = std::enable_if_t<IsEvalNode<NodeT>>>
AsyCost single_node_storage(NodeT const& node) {
  using namespace sequant;
  auto total = AsyCost::zero();
  if (node->is_tensor()) {
    auto const& t = node->as_tensor();
    size_t no = ranges::count_if(t.const_braket(), [](Index const& i) {
      return i.space() == IndexSpace::active_occupied;
    });
    size_t nv = t.bra_rank() + t.ket_rank() - no;
    total += AsyCost{no, nv};
  }
  return total;
}

///
/// \param nodes An iterable of EvalNodes.
/// \return Given an iterable of EvalNodes, returns a set of Imeds sorted by
///         their keys.
///
template <typename NodesT,
          typename std::enable_if_t<IsIterableOfEvalNodes<NodesT>, bool> = true>
container::set<Imed> intermediates(NodesT const& nodes) {
  container::set<Imed> result;
  size_t p = 0;
  auto visitor = [&p, &result](auto const& n) -> bool {
    auto&& put = result.emplace(hash::value(*n));
    if (put.second) {
      // inserting this imed for the first time
      put.first->flops = asy_cost(n);
      put.first->memory = single_node_storage(n);
      put.first->pos.push_back(p);
      put.first->expr.append(to_expr(n));
      return true;  // keep visiting children nodes
    } else {
      put.first->pos.push_back(p);
      put.first->expr.append(to_expr(n));
      return false;  // no need to visit children nodes
    }
  };

  for (auto const& n : nodes) {
    n.visit_internal(visitor);
    ++p;
  }

  return result;
}

///
/// \return Given a sequant Sum, returns a set of Imeds sorted by their keys.
///
container::set<Imed> intermediates(Sum const&);

///
/// \brief Given a set of intermediates returns a map from the position of a
///        term to the set of positions of connected terms. Two terms are
///        connected if they share a common intermediate.
/// \param imeds The intermediates to analyze.
/// \return A map from the position of a term to the set of positions of the
///         connected terms.
template <typename ImedsT, typename = std::enable_if_t<
                               std::is_convertible_v<Imed, IteredT<ImedsT>>>>
adjacency_map adjacency_set(ImedsT const& imeds) {
  using ranges::views::tail;

  adjacency_map result;
  for (auto const& im : imeds) {
    result.emplace(im.pos[0],
                   tail(im.pos) | ranges::to<container::set<size_t>>);
  }
  return result;
}

// /////////////////////////////////////////////////////////////////////////////
// Deprecated functions below
// /////////////////////////////////////////////////////////////////////////////

///
/// \param sum A Sum object representing a many-body equation.
/// \return Returns a Sum object with terms that take part in the
///         common subexpression elimination.
///
template <typename Pred = std::function<bool(Imed const&)>,
          typename = std::enable_if_t<std::is_invocable_r_v<bool, Pred, Imed>>>
[[deprecated]] Sum connected_terms(
    Sum const& sum, Pred const& pred = [](auto const&) { return true; }) {
  using ranges::views::filter;
  using ranges::views::transform;

  auto const imeds = intermediates(sum);
  container::set<size_t> added;
  Sum result;
  for (auto const& i :
       imeds | filter(Imed::filter_count<std::greater<>>(1)) | filter(pred)) {
    for (auto p : i.pos)
      if (added.emplace(p).second) result.append(sum.summand(p)->clone());
  }
  return result;
}

///
/// \brief Get the adjacency matrix of the common sub-expressions between terms.
///
/// \param sum A Sum object representing terms that have common
///            sub-expressions.
/// \return Adjacency matrix of the common sub-expressions graph.
///
///
[[deprecated]] container::vector<container::vector<size_t>> cse_graph(
    Sum const& sum);

[[deprecated]] void write_cse_graph(
    std::filesystem::path const& file,
    container::vector<container::vector<size_t>> const& graph);

}  // namespace sequant

#endif  // SEQUANT_ANALYSIS_ALGORITHM_HPP
