#ifndef SEQUANT_EVAL_LEGALITY_HPP
#define SEQUANT_EVAL_LEGALITY_HPP

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/index.hpp>

#include <algorithm>
#include <cstddef>

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

  //!< The batch-loop axes this value's computation depends on -- carried in
  //!< its own result OR contracted anywhere at or below it in the
  //!< contraction DAG. Filled by \c build_site_of (Task 1).
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
/// \brief The build-site of the value rooted at \p node: the batch-loop axes
/// (per \p policy.is_batchable_index()) that appear anywhere in \p node's
/// subtree -- carried in a descendant's (or \p node's own) result indices, OR
/// contracted at \p node or any node below it.
///
/// \details Derived bottom-up over the binary contraction DAG:
///   - a LEAF contributes its own \c canon_indices(), filtered to the
///     batchable subset;
///   - an INTERNAL node contributes the union of its children's build-sites
///     plus its own \c contracted_indices(node) (see eval.hpp), also
///     filtered to the batchable subset.
/// The result is de-duplicated by \c Index identity (space + ordinal +
/// proto-indices); order follows first discovery (left subtree, then right,
/// then this node's own contracted axes).
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

  if (node.leaf()) {
    for (Index const& ix : node->canon_indices()) add_if_new(ix);
    return result;
  }

  for (Index const& ix : build_site_of(node.left(), policy)) add_if_new(ix);
  for (Index const& ix : build_site_of(node.right(), policy)) add_if_new(ix);
  for (Index const& ix : contracted_indices(node)) add_if_new(ix);
  return result;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_LEGALITY_HPP
