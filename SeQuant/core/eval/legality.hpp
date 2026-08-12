#ifndef SEQUANT_EVAL_LEGALITY_HPP
#define SEQUANT_EVAL_LEGALITY_HPP

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/index.hpp>

#include <algorithm>
#include <cstddef>
#include <unordered_map>

namespace sequant::eval {

///
/// \brief Legality analysis output types (SP1 of the ordered-scope batched-
/// eval design): per-value classification of which batch-loop axes a value's
/// computation depends on, and (later tasks) whether it is legal to home the
/// value at each such axis.
///

///
/// \brief The role a batch-loop axis plays at one value (Task 2 fills \c
/// per_axis; this enum is declared now so \c AxisClass / \c CellLegality have
/// a complete type in Task 1).
///
enum class LoopRole {
  LoopLocal,      //!< the axis is fully consumed within the value's own build
  Reduction,      //!< the axis is summed (contracted) at or below the value
  LoopCarried,    //!< the axis survives into the value's own result indices
  LoopInvariant,  //!< the axis encloses the value but the value does not
                  //!< depend on it (a pure replication context)
};

///
/// \brief The classification of one (value, build-site axis) pair.
///
struct AxisClass {
  Index axis;
  LoopRole role;
};

///
/// \brief The legality record of one value (one \c ValueCell::hash): its
/// build-site (Task 1), per-axis role (Task 2), and home-placement legality
/// (Task 3/4).
///
struct CellLegality {
  std::size_t hash = 0;  //!< == ValueCell::hash

  //!< The batch-loop axes this value depends on AT ITS OWN NODE -- carried in
  //!< its own result OR contracted at its own node (NOT recursively over its
  //!< operand subtrees; each operand is itself a value with its own,
  //!< separately-analyzed \c CellLegality). Filled by \c build_site_of.
  container::svector<Index> build_site;

  //!< One entry per \c build_site axis, classifying its role. Filled in
  //!< Task 2.
  container::svector<AxisClass> per_axis;

  //!< The shallowest legal home mode-set for this value. Filled in Task 3.
  container::svector<Index> home_floor;

  //!< Loops that must re-enter (the value cannot be homed above them without
  //!< a re-materializing split). Filled in Task 3/4.
  container::svector<Index> forced_split_axes;
};

///
/// \brief A forest-wide legality analysis: one \c CellLegality per distinct
/// value. Assembled in Task 2 (\c analyze_legality); this header only
/// declares the container shape.
///
struct LegalitySchedule {
  container::svector<CellLegality> cells;
};

///
/// \brief The AT-NODE build-site of the value rooted at \p node: the
/// batch-loop axes (per \p policy.is_batchable_index()) that appear AT this
/// node -- carried in \p node's own result indices, OR contracted AT \p node.
///
/// \details Deliberately NOT a subtree union: in the ordered-scope model
/// every node is itself a VALUE, and a value is one contraction of already
/// -cached OPERAND values (its children), so its build-site is only the axes
/// it carries or contracts at its own node -- the axes below it belong to
/// the operand VALUES' own (separately-analyzed) build-sites, not to this
/// one. Concretely, the union of
///   - \p node's own \c canon_indices() (its result/free indices), and
///   - \p node's own \c contracted_indices(node) (see eval.hpp; empty for a
///     leaf or a non-product node),
/// filtered to the batchable subset. The result is de-duplicated by \c Index
/// identity (space + ordinal + proto-indices); order follows first
/// discovery (carried indices, then contracted indices).
///
[[nodiscard]] inline container::svector<Index> build_site_of(
    meta::eval_node auto const& node, BatchPolicy const& policy) {
  container::svector<Index> result;
  auto const batchable = policy.is_batchable_index();
  auto const add_if_new = [&](Index const& ix) {
    if (!batchable(ix)) return;
    if (std::find(result.begin(), result.end(), ix) == result.end())
      result.push_back(ix);
  };

  for (Index const& ix : node->canon_indices()) add_if_new(ix);
  for (Index const& ix : contracted_indices(node)) add_if_new(ix);
  return result;
}

///
/// \brief Classify one batch-loop axis \p axis against one value, given the
/// value's own result indices \p carried, the axes it contracts AT its own
/// node \p contracted_below (see \c contracted_indices in eval.hpp), and
/// every use-site \p occurrences of the value in the forest.
///
/// \details The four-way decision tree (SP1 design, Task 2):
///   - Q1: does \p carried hold an index of \p axis's \c IndexSpace (compared
///     by \c base_key(), i.e. TYPE not identity)?
///     - NO  -> Q2a: does \p contracted_below hold an index of that type
///       (the value reduces the axis at its own node)? -> \c Reduction;
///       otherwise the axis merely encloses the value without touching it
///       -> \c LoopInvariant.
///     - YES -> Q2b: for every occurrence that has an ENCLOSING loop of that
///       axis type in its \c OccurrenceRec::ectx, compare (via the
///       ordinal-and-proto-aware \c Index::operator==) that loop's own
///       \c Index against the occurrence's own carried index of the same
///       type. Bound to the SAME \c Index at every such occurrence (lockstep
///       with the enclosing loop) -> \c LoopLocal; bound to a DIFFERENT
///       \c Index at any occurrence (a free / cross-iteration read) ->
///       \c LoopCarried. If no occurrence has an enclosing loop of that type
///       at all, the axis is a plain free result index with nothing to lock
///       it to a loop iteration -> \c LoopCarried.
///
[[nodiscard]] inline LoopRole classify_axis(
    container::svector<Index> const& carried,
    container::svector<Index> const& contracted_below, Index const& axis,
    container::svector<OccurrenceRec> const& occurrences) {
  auto const same_type = [&](Index const& ix) {
    return ix.space().base_key() == axis.space().base_key();
  };

  bool const carries_axis =
      std::any_of(carried.begin(), carried.end(), same_type);

  if (!carries_axis) {
    bool const reduces_axis = std::any_of(contracted_below.begin(),
                                          contracted_below.end(), same_type);
    return reduces_axis ? LoopRole::Reduction : LoopRole::LoopInvariant;
  }

  bool found_enclosing = false;
  for (OccurrenceRec const& occ : occurrences) {
    auto const ectx_it =
        std::find_if(occ.ectx.begin(), occ.ectx.end(),
                     [&](auto const& e) { return same_type(e.first); });
    if (ectx_it == occ.ectx.end())
      continue;  // no enclosing loop of this
                 // type at this occurrence
    found_enclosing = true;

    auto const carried_it =
        std::find_if(occ.carried.begin(), occ.carried.end(), same_type);
    if (carried_it == occ.carried.end() || !(*carried_it == ectx_it->first))
      return LoopRole::LoopCarried;  // free / cross-iteration binding
  }
  return found_enclosing ? LoopRole::LoopLocal : LoopRole::LoopCarried;
}

///
/// \brief Fold a \c RichSchedule into the first real \c LegalitySchedule: one
/// \c CellLegality per \c ValueCell, with its at-node \c build_site (\c
/// build_site_of) and \c per_axis classification (\c classify_axis) filled
/// in. \c home_floor / \c forced_split_axes are left empty (Task 3/4).
///
/// \details \c ValueCell does not itself retain the forest node its \c hash
/// resolves to (only \c carried, a copy of the node's \c canon_indices()), so
/// \p forest is re-walked once up front into a hash -> representative-node
/// map (same pattern as \c placement_remat.hpp's router build and \c
/// scope_executor.hpp's \c build_value_node_map) to recover each cell's own
/// \c contracted_indices(node). Under perfect CSE every occurrence of a value
/// shares one node shape, so the first-visited representative suffices.
///
template <meta::eval_node_range R>
[[nodiscard]] inline LegalitySchedule analyze_legality(
    RichSchedule const& rich, R const& forest, BatchPolicy const& policy) {
  using Node = std::ranges::range_value_t<R>;

  std::unordered_map<std::size_t, Node> node_of;
  auto visit = [&](auto&& self, Node const& n) -> void {
    node_of.emplace(n->hash_value(), n);
    if (!n.leaf()) {
      self(self, n.left());
      self(self, n.right());
    }
  };
  for (auto const& tree : forest) visit(visit, tree);

  auto const batchable = policy.is_batchable_index();

  LegalitySchedule out;
  out.cells.reserve(rich.cells.size());
  for (ValueCell const& vc : rich.cells) {
    CellLegality cl;
    cl.hash = vc.hash;

    container::svector<Index> contracted_below;
    if (auto const it = node_of.find(vc.hash); it != node_of.end())
      for (Index const& ix : contracted_indices(it->second))
        contracted_below.push_back(ix);

    container::svector<Index> site;
    auto const add_if_new = [&](Index const& ix) {
      if (!batchable(ix)) return;
      if (std::find(site.begin(), site.end(), ix) == site.end())
        site.push_back(ix);
    };
    for (Index const& ix : vc.carried) add_if_new(ix);
    for (Index const& ix : contracted_below) add_if_new(ix);

    for (Index const& axis : site) {
      AxisClass ac;
      ac.axis = axis;
      ac.role =
          classify_axis(vc.carried, contracted_below, axis, vc.occurrences);
      cl.per_axis.push_back(std::move(ac));
    }
    cl.build_site = std::move(site);
    out.cells.push_back(std::move(cl));
  }
  return out;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_LEGALITY_HPP
