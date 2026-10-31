#ifndef SEQUANT_CORE_EXPORT_NAME_RECYCLING_HPP
#define SEQUANT_CORE_EXPORT_NAME_RECYCLING_HPP

#include <SeQuant/core/export/context.hpp>
#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/export/export_node.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/irreps/any.hpp>
#include <SeQuant/core/irreps/eval_tree.hpp>
#include <SeQuant/core/optimize/name_recycling.hpp>
#include <SeQuant/core/utility/tensor.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <map>
#include <ranges>
#include <span>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sequant {

/// Determines ranges of consecutive expressions that might be executed
/// repeatedly as part of a loop over batch indices. The result is
/// conservative: for every depth k, every maximal range of consecutive
/// expressions whose (first k) batch indices agree is reported, which covers
/// every loop a generator may create when fusing the batching of consecutive
/// expressions.
///
/// @param expressions Range of expression trees in the order of export
/// @param ctx The context holding the batch indices of the expressions
/// @returns Ranges [begin, end) of positions in @p expressions
template <std::ranges::random_access_range Range>
std::vector<std::pair<std::size_t, std::size_t>> batch_loop_regions(
    const Range &expressions, const ExportContext &ctx) {
  std::vector<std::pair<std::size_t, std::size_t>> regions;

  const std::size_t size = std::ranges::size(expressions);
  std::vector<std::span<const Index>> batch_indices;
  batch_indices.reserve(size);
  std::size_t max_depth = 0;
  for (const auto &tree : expressions) {
    batch_indices.push_back(ctx.batch_indices(tree->id()));
    max_depth = std::max(max_depth, batch_indices.back().size());
  }

  auto shares_prefix = [&](std::size_t lhs, std::size_t rhs,
                           std::size_t depth) {
    return batch_indices[lhs].size() >= depth &&
           batch_indices[rhs].size() >= depth &&
           std::ranges::equal(batch_indices[lhs].first(depth),
                              batch_indices[rhs].first(depth));
  };

  for (std::size_t depth = 1; depth <= max_depth; ++depth) {
    std::size_t begin = 0;
    while (begin < size) {
      if (batch_indices[begin].size() < depth) {
        ++begin;
        continue;
      }

      std::size_t end = begin + 1;
      while (end < size && shares_prefix(begin, end, depth)) {
        ++end;
      }

      regions.emplace_back(begin, end);
      begin = end;
    }
  }

  return regions;
}

namespace detail {

/// IRREP context that resolves results of previously processed trees via the
/// given map (keyed by rewritten tensors) and everything else like
/// RegisteredIrrepContext
class ScannedIrrepContext {
 public:
  using set_type = AnyIrrepSet;
  using Registry = std::map<Tensor, AnyIrrepSet, TensorBlockLessThanComparator>;

  ScannedIrrepContext(const ExportContext &ctx, const Registry &scanned)
      : m_registered(ctx), m_ctx(ctx), m_scanned(scanned) {}

  AnyIrrepSet irreps_of(const Tensor &tensor) const {
    Tensor rewritten(tensor);
    m_ctx.rewrite(rewritten);
    if (auto it = m_scanned.find(rewritten); it != m_scanned.end()) {
      return it->second;
    }

    return m_registered.irreps_of(tensor);
  }
  AnyIrrepSet irreps_of(const Variable &variable) const {
    return m_registered.irreps_of(variable);
  }
  AnyIrrepSet totally_symmetric() const {
    return m_registered.totally_symmetric();
  }
  AnyIrrepSet empty_set() const { return m_registered.empty_set(); }
  AnyIrrepSet any_set() const { return m_registered.any_set(); }
  AnyIrrepSet product(const AnyIrrepSet &lhs, const AnyIrrepSet &rhs) const {
    return m_registered.product(lhs, rhs);
  }
  AnyIrrepSet conjugate(const AnyIrrepSet &set) const {
    return m_registered.conjugate(set);
  }

 private:
  RegisteredIrrepContext m_registered;
  const ExportContext &m_ctx;
  const Registry &m_scanned;
};

/// Determines the IRREPs of the results of the given trees (in order of
/// execution)
/// @param scanned The IRREPs of results of previously scanned trees (keyed by
/// rewritten tensors), which gets updated with the results of @p expressions
/// @returns The IRREPs of the results of the given trees keyed by the trees'
/// IDs
template <std::ranges::random_access_range Range, typename Predicate>
std::unordered_map<std::size_t, AnyIrrepSet> result_irreps(
    const Range &expressions, const ExportContext &ctx,
    const Predicate &overwrites, ScannedIrrepContext::Registry &scanned) {
  std::unordered_map<std::size_t, AnyIrrepSet> irreps;
  ScannedIrrepContext irrep_ctx(ctx, scanned);

  for (const auto &tree : expressions) {
    AnyIrrepSet result = irreps_of(tree, irrep_ctx);

    if (tree->is_tensor()) {
      Tensor rewritten = tree->as_tensor();
      ctx.rewrite(rewritten);

      auto it = scanned.find(rewritten);
      if (it == scanned.end()) {
        scanned.emplace(std::move(rewritten), result);
      } else if (overwrites(tree)) {
        it->second = result;
      } else {
        it->second.unite(result);
      }
    }

    irreps.emplace(tree->id(), std::move(result));
  }

  return irreps;
}

}  // namespace detail

/// Recycles names of intermediates in the given groups of expressions to be
/// exported (see opt::recycle_intermediate_names). Names are recycled within
/// each group, i.e. intermediates must be defined and used within the same
/// group. Two intermediates may only share a name if they are batched over the
/// same indices and, if @p ctx has an IRREP context, span the same IRREPs
/// (where IRREPs are propagated across all groups). Loops over batch indices
/// are taken into account when determining the lifetime of intermediates (see
/// batch_loop_regions). Tensors are identified the way the export does, i.e.
/// after being rewritten by @p ctx. Trees that need rebuilding are rebuilt via
/// to_export_tree, retaining their IDs such that information associated with
/// them (e.g. batch indices) stays valid.
///
/// @param groups Range of groups of expressions in the order of export, where
/// every group is a random-access range of expression trees in the order of
/// export
/// @param ctx The context used for the export
/// @param overwrites Predicate determining whether a given expression tree
/// overwrites rather than adds to the previous value of its result. Only the
/// results of such trees are considered for name recycling.
/// @param compatibility Shape compatibility required for recycling a name
/// @param label_gen Generator for new names of intermediates
/// @returns Statistics about the renaming (accumulated over all groups)
template <std::ranges::range Groups, typename Predicate>
  requires std::ranges::random_access_range<std::ranges::range_value_t<Groups>>
opt::NameRecyclingResult recycle_intermediate_names(
    Groups &groups, const ExportContext &ctx, const Predicate &overwrites,
    opt::ShapeCompatibility compatibility,
    opt::IntermediateNamePool::LabelGenerator label_gen =
        opt::IntermediateNamePool::default_label_generator()) {
  using TreeNode =
      std::ranges::range_value_t<std::ranges::range_value_t<Groups>>;

  auto overwriting = [&overwrites](std::size_t, const TreeNode &tree) {
    return overwrites(tree);
  };

  detail::ScannedIrrepContext::Registry scanned_irreps;
  opt::NameRecyclingResult stats;

  for (auto &expressions : groups) {
    std::unordered_map<std::size_t, AnyIrrepSet> irreps;
    if (ctx.has_irrep_context()) {
      irreps =
          detail::result_irreps(expressions, ctx, overwrites, scanned_irreps);
    }

    const opt::NameRecyclingResult group_stats =
        opt::recycle_intermediate_names(
            expressions,
            [](const ResultExpr &res) {
              return to_export_tree<typename TreeNode::value_type>(res);
            },
            opt::NameRecyclingOptions<TreeNode>{
                .is_recyclable = overwriting,
                .starts_new_value = overwriting,
                .shape_compatibility = compatibility,
                .extra_compatible =
                    [&](const TreeNode &previous, const TreeNode &current) {
                      if (!std::ranges::equal(
                              ctx.batch_indices(previous->id()),
                              ctx.batch_indices(current->id()))) {
                        return false;
                      }

                      return irreps.empty() || irreps.at(previous->id()) ==
                                                   irreps.at(current->id());
                    },
                .repeated_regions = batch_loop_regions(expressions, ctx),
                .label_gen = label_gen,
                // Tensors are identified by their rewritten form during export
                .normalize = [&ctx](Tensor &tensor) { ctx.rewrite(tensor); },
                .on_rebuilt =
                    [](const TreeNode &original, TreeNode &rebuilt) {
                      rebuilt->set_id(original->id());
                    }});

    stats.num_values += group_stats.num_values;
    stats.num_names += group_stats.num_names;
  }

  return stats;
}

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_NAME_RECYCLING_HPP
