#ifndef SEQUANT_COMMON_SUBEXPRESSION_ELIMINATION_HPP
#define SEQUANT_COMMON_SUBEXPRESSION_ELIMINATION_HPP

#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/eval/eval_node_compare.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/utility/aggregate.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <algorithm>
#include <concepts>
#include <format>
#include <functional>
#include <optional>
#include <ranges>
#include <span>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace sequant::opt {

namespace cse {

/// Functor that can be used to identify common subexpressions while iterating
/// expression trees
template <typename TreeNode, bool force_hash_collisions = false>
class SubexpressionIdentifier {
 public:
  SubexpressionIdentifier() = default;
  SubexpressionIdentifier(std::vector<Index> indices)
      : intermediate_hashs(
            SubexpressionHashCollector<TreeNode, force_hash_collisions>{}
                .bucket_count(),  // use default bucket count
            TreeNodeHasher<TreeNode, force_hash_collisions>{},
            TreeNodeEqualityComparator<TreeNode>(std::move(indices))) {}

  bool operator()(const TreeNode &tree) {
    using std::ranges::end;

    if (auto it = intermediate_hashs.find(&tree);
        it != end(intermediate_hashs)) {
      // The expression identified by tree has been used before -> stop
      // visiting of subtree
      // TODO: Early return might make us miss better subexpression choices that
      // may appear if we decide to not factor out the current one
      ++it->second;
      return false;
    }

    intermediate_hashs.emplace(&tree, 1);

    return true;
  }

  SubexpressionUsageCounts<TreeNode, force_hash_collisions>
  take_subexpression_map() {
    SubexpressionUsageCounts<TreeNode, force_hash_collisions> usages(
        SubexpressionUsageCounts<TreeNode, force_hash_collisions>{}
            .bucket_count(),  // use default bucket count
        intermediate_hashs.hash_function(), intermediate_hashs.key_eq());

    for (const auto &[node_ptr, usage_count] : intermediate_hashs) {
      if (usage_count < 2) {
        // Everything that is used less than 2 times is not a common
        // subexpression
        continue;
      }

      usages.emplace(*node_ptr, usage_count);
    }

    intermediate_hashs.clear();

    return usages;
  }

 private:
  // Note: we are using the hash collector rather than the
  // SubexpressionUsageCounts because we will collect many, many entries which
  // will only appear once and the hash collector only stores node pointers,
  // the latter stores full node objects which are more than an order of
  // magnitude larger, which can easily lead to memory bottlenecks for large
  // expressions.
  SubexpressionHashCollector<TreeNode, force_hash_collisions>
      intermediate_hashs;
};

/// Functor that will perform the elimination of previously found subexpressions
template <std::ranges::range VectorLike, typename TreeNode,
          typename Transformer, typename LabelGenerator,
          bool force_hash_collisions = false>
class SubexpressionReplacer {
 public:
  SubexpressionReplacer(
      VectorLike &expr_trees,
      const SubexpressionUsageCounts<TreeNode, force_hash_collisions> &map,
      const Transformer &transformer, const LabelGenerator &label_gen)
      : expr_trees(expr_trees),
        subexpressions(map),
        expr_to_tree(transformer),
        label_gen(label_gen) {}

  void perform_replacements(std::size_t begin, std::size_t end) {
    // Ensure we won't have to reallocate while iterating over the expressions
    expr_trees.reserve(expr_trees.size() + subexpressions.size());
    cse_definition_indices.reserve(subexpressions.size());

    for (std::size_t i = begin; i < end; ++i) {
      current_expr_idx = i;
      expr_trees.at(i + expr_offset).visit_internal(*this);
    }
  }

  bool operator()(TreeNode &tree) {
    using std::ranges::begin;
    using std::ranges::end;

    auto it = subexpressions.find(tree);

    if (it == end(subexpressions)) {
      return true;
    }

    auto label_it = cse_names.find(tree);

    std::wstring label = label_it == end(cse_names)
                             ? toUtf16(label_gen(tree, name_counter++))
                             : label_it->second;

    ExprPtr expr = [&]() {
      if (tree->is_tensor()) {
        Index::index_vector indices = tree->canon_indices();

        // Sort indices into distinct groups based on index basis (space and
        // basis instance) in order to create a deterministic set of indices.
        // Since there can't be any permutational symmetries between indices
        // of different bases, we can arrange them in arbitrary order as long
        // as we retain the relative order of indices that belong to the same
        // basis.
        std::ranges::stable_sort(indices, std::less<>{},
                                 [](const Index &idx) { return idx.basis(); });

        return ex<Tensor>(label, bra(), ket(), aux(std::move(indices)),
                          Symmetry::Nonsymm, BraKetSymmetry::Nonsymm,
                          ColumnSymmetry::Nonsymm);
      }

      return ex<Variable>(label);
    }();

    ExprPtr original = tree->expr();

    std::optional<TreeNode> intermediate_definition;
    if (label_it == end(cse_names)) {
      cse_names.emplace(tree, label);
      // Store tree as the way to compute the current CSE (under the
      // designated name and with canonical indexing)
      // Since tree node objects are immutable, we have to disassemble
      // and reassemble the tree in order to get the desired outcome.

      ResultExpr intermediate = [&]() {
        if (expr->is<Tensor>()) {
          return ResultExpr(expr->as<Tensor>(), to_expr(tree));
        }

        return ResultExpr(expr->as<Variable>(), to_expr(tree));
      }();

      intermediate_definition = expr_to_tree(intermediate);
    }

    if (tree->canon_phase() != 1) {
      expr *= ex<Constant>(tree->canon_phase());
    }

    if (tree.root()) {
      SEQUANT_ASSERT(original->is<Tensor>() || original->is<Variable>());

      ResultExpr result = [&]() {
        if (original->is<Tensor>()) {
          return ResultExpr(original->as<Tensor>(), std::move(expr));
        }

        return ResultExpr(original->as<Variable>(), std::move(expr));
      }();

      tree = expr_to_tree(result);
    } else {
      tree = expr_to_tree(expr);
    }

    if (intermediate_definition.has_value()) {
      SEQUANT_ASSERT(expr_trees.size() > current_expr_idx + expr_offset);
      auto expr_it = expr_trees.begin();
      std::advance(expr_it, current_expr_idx + expr_offset);

      cse_definition_indices.emplace_back(current_expr_idx + expr_offset);

      expr_trees.insert(expr_it, std::move(intermediate_definition.value()));
      expr_offset += 1;
    }

    return false;
  }

  const std::vector<std::size_t> &cse_indices() const {
    return cse_definition_indices;
  }

 private:
  VectorLike &expr_trees;
  std::size_t expr_offset = 0;
  std::size_t current_expr_idx = 0;
  const SubexpressionUsageCounts<TreeNode, force_hash_collisions>
      &subexpressions;
  const Transformer &expr_to_tree;
  SubexpressionNames<TreeNode, force_hash_collisions> cse_names;
  std::size_t name_counter = 1;
  const LabelGenerator &label_gen;
  std::vector<std::size_t> cse_definition_indices;
};

}  // namespace cse

template <typename TreeNode>
struct CSEOptions {
  SEQUANT_DESIGNATED_INIT_ONLY;
  // TODO: provide predicate that computes the FLOPs required to compute the CSE
  // and compares that with the expected speed for reading the CSE from disk.
  // Use typical FLOPs/s and I/O rates from e.g. https://ssd.userbenchmark.com/
  // https://openbenchmarking.org/test/pts/mt-dgemm
  // to determine whether or not pre-computing the CSE is worth it.
  // Note: of course the number of reusages of the CSE plays into this.
  std::function<bool(const TreeNode &, std::size_t)> filter_predicate =
      [](const TreeNode &, std::size_t) { return true; };
  std::function<std::string(const TreeNode &, std::size_t)> label_gen =
      [](const TreeNode &, std::size_t counter) {
        return std::format("CSE{}", counter);
      };
  std::vector<Index> batch_indices = {};
};

/// @return a CSEOptions::filter_predicate that accepts a subexpression only if
/// it is used at least @p min_usage times and contains at least
/// @p min_tensors tensors
/// @note subexpressions with fewer than two tensors are usually not worth
/// creating, storing and reusing; e.g. in the symmetrization
/// `1/2 R[abij] + 1/2 R[baji]`, `1/2 R` is such a subexpression
template <typename TreeNode>
std::function<bool(const TreeNode &, std::size_t)> min_usage_filter(
    std::size_t min_usage, std::size_t min_tensors = 2) {
  return
      [min_usage, min_tensors](const TreeNode &tree, std::size_t usage_count) {
        if (usage_count < min_usage) return false;

        std::size_t num_tensors = 0;
        tree.visit_leaf([&num_tensors](const TreeNode &node) {
          num_tensors += node->is_tensor();
        });
        return num_tensors >= min_tensors;
      };
}

/// Takes the range of expression trees and performs common subexpression
/// elimination on them. Evaluation trees for intermediates are inserted into
/// the given container as needed.
///
/// @param expr_trees Container (vector-like) of evaluation trees
/// @param expr_to_tree Functor that can convert ExprPtr and ResultExpr objects
/// into suitable evaluation trees
/// @param filter_predicate A predicate to filter out subexpressions that shall
/// not be eliminated
/// @returns The positions of the inserted intermediate definitions
///
/// @note Every intermediate gets its own name. Once the final order of
/// evaluation is known, opt::recycle_intermediate_names (see
/// name_recycling.hpp) can be used to reuse names of intermediates that are
/// no longer needed.
template <std::ranges::range VectorLike,
          std::regular_invocable<ExprPtr> Transformer,
          bool force_hash_collisions = false>
  requires std::regular_invocable<Transformer, ResultExpr> &&
           meta::eval_node<std::ranges::range_value_t<VectorLike>>
std::vector<std::size_t> eliminate_common_subexpressions(
    VectorLike &expr_trees, std::size_t expr_begin, std::size_t expr_end,
    const Transformer &expr_to_tree,
    const CSEOptions<std::ranges::range_value_t<VectorLike>> opts = {}) {
  using std::ranges::begin;
  using std::ranges::end;

  using TreeNode = std::ranges::range_value_t<VectorLike>;

  cse::SubexpressionIdentifier<TreeNode, force_hash_collisions> identifier(
      opts.batch_indices);

  for (const TreeNode &current : std::ranges::subrange(
           expr_trees.begin() + expr_begin, expr_trees.begin() + expr_end)) {
    current.visit_internal(identifier);
  }

  SubexpressionUsageCounts<TreeNode, force_hash_collisions>
      subexpression_usages = identifier.take_subexpression_map();

  // Apply filter
  std::erase_if(subexpression_usages, [&opts](const auto &pair) {
    return !opts.filter_predicate(pair.first, pair.second);
  });

  cse::SubexpressionReplacer<VectorLike, TreeNode, Transformer,
                             decltype(opts.label_gen), force_hash_collisions>
      replacer(expr_trees, subexpression_usages, expr_to_tree, opts.label_gen);

  replacer.perform_replacements(expr_begin, expr_end);

  return replacer.cse_indices();
}

template <std::ranges::range VectorLike,
          std::regular_invocable<ExprPtr> Transformer,
          bool force_hash_collisions = false>
  requires std::regular_invocable<Transformer, ResultExpr> &&
           meta::eval_node<std::ranges::range_value_t<VectorLike>>
std::vector<std::size_t> eliminate_common_subexpressions(
    VectorLike &expr_trees, const Transformer &expr_to_tree,
    const CSEOptions<std::ranges::range_value_t<VectorLike>> opts = {}) {
  return eliminate_common_subexpressions(expr_trees, 0, expr_trees.size(),
                                         expr_to_tree, opts);
}

/// Result of eliminate_common_subexpressions_batched
struct BatchedCSEResult {
  /// Positions of the inserted intermediate definitions
  std::vector<std::size_t> cse_positions;
  /// Ranges [begin, end) of the expression trees (including the inserted
  /// definitions) that were processed together
  std::vector<std::pair<std::size_t, std::size_t>> partitions;
  /// The batch indices that were used for each partition
  std::vector<std::vector<Index>> batch_indices;
};

/// Like eliminate_common_subexpressions, but takes into account that
/// expressions can be computed in batches over different sets of indices.
/// Consecutive expressions that share the same batch indices are processed
/// together, such that subexpressions are only shared among expressions that
/// are batched identically.
///
/// @param expr_trees Container (vector-like) of evaluation trees
/// @param expr_batch_indices The batch indices of each of the given trees. An
/// empty optional means that the tree has no well-defined set of batch indices,
/// which causes it to be processed on its own, without batching. Must have as
/// many elements as there are trees on input.
/// @param expr_to_tree See eliminate_common_subexpressions
/// @param opts See eliminate_common_subexpressions. The batch indices are taken
/// from @p expr_batch_indices instead of from the options.
template <std::ranges::range VectorLike,
          std::regular_invocable<ExprPtr> Transformer,
          bool force_hash_collisions = false>
  requires std::regular_invocable<Transformer, ResultExpr> &&
           meta::eval_node<std::ranges::range_value_t<VectorLike>>
BatchedCSEResult eliminate_common_subexpressions_batched(
    VectorLike &expr_trees,
    std::span<const std::optional<std::vector<Index>>> expr_batch_indices,
    const Transformer &expr_to_tree,
    CSEOptions<std::ranges::range_value_t<VectorLike>> opts = {}) {
  SEQUANT_ASSERT(expr_batch_indices.size() == expr_trees.size());

  BatchedCSEResult result;

  std::pair<std::size_t, std::size_t> partition = {0, 0};
  while (partition.second < expr_trees.size()) {
    partition.first = partition.second;

    // Positions of inserted definitions have to be skipped to get to the
    // index of the original expression
    SEQUANT_ASSERT(partition.first >= result.cse_positions.size());
    const auto &batch_indices =
        expr_batch_indices[partition.first - result.cse_positions.size()];

    if (batch_indices.has_value()) {
      for (; partition.second < expr_trees.size(); ++partition.second) {
        SEQUANT_ASSERT(partition.second >= result.cse_positions.size());
        if (expr_batch_indices[partition.second -
                               result.cse_positions.size()] != batch_indices) {
          break;
        }
      }
    }

    if (partition.first == partition.second) {
      ++partition.second;
      opts.batch_indices.clear();
    } else {
      opts.batch_indices = *batch_indices;
    }

    std::vector<std::size_t> new_cses = eliminate_common_subexpressions(
        expr_trees, partition.first, partition.second, expr_to_tree, opts);

    result.cse_positions.insert(result.cse_positions.end(), new_cses.begin(),
                                new_cses.end());

    partition.second += new_cses.size();

    result.batch_indices.push_back(std::move(opts.batch_indices));
    result.partitions.push_back(partition);
  }

  return result;
}

/// Orders trees by their batch indices, which are treated as having a
/// "priority": we start at the highest priority (first) batch index and compare
/// on the first batch index that only one of the trees has. Trees with the
/// lower priority batch indices (or even none) precede the ones with the higher
/// priority (or more) batch indices. Trees that are not tensors are treated as
/// having no batch indices.
///
/// @returns Whether @p lhs strictly precedes @p rhs
template <typename TreeNode>
bool batch_priority_less(const TreeNode &lhs, const TreeNode &rhs,
                         std::span<const Index> batch_indices) {
  auto contains = [](const TreeNode &tree, const Index &idx) {
    if (!tree->is_tensor()) {
      return false;
    }
    auto indices = tree->as_tensor().const_indices();
    return std::ranges::find(indices, idx) != indices.end();
  };

  for (const Index &idx : batch_indices) {
    const bool lhs_contains = contains(lhs, idx);
    const bool rhs_contains = contains(rhs, idx);

    if (lhs_contains != rhs_contains) {
      return rhs_contains;
    }
  }

  return false;
}

namespace cse {

/// Adds the number of tensor leaves of @p node to @p num_tensors and returns
/// whether none of its leaves satisfies @p is_volatile_leaf. Note that the
/// traversal stops at the first volatile leaf.
template <typename TreeNode>
bool count_tensors_if_persistent(
    const TreeNode &node,
    const std::function<bool(const TreeNode &)> &is_volatile_leaf,
    std::size_t &num_tensors) {
  if (node.leaf()) {
    num_tensors += node->is_tensor();
    return !is_volatile_leaf(node);
  }

  return count_tensors_if_persistent(node.left(), is_volatile_leaf,
                                     num_tensors) &&
         count_tensors_if_persistent(node.right(), is_volatile_leaf,
                                     num_tensors);
}

}  // namespace cse

/// Replaces all maximal subtrees that are "persistent" by intermediates, no
/// matter how often they are used. A subtree is persistent if none of its
/// leaves is volatile and if it contains at least @p min_tensors tensors.
/// Definitions of the intermediates are inserted into the given container
/// right before their first use.
///
/// @param expr_trees Container (vector-like) of evaluation trees
/// @param is_volatile_leaf Predicate that identifies leaves that are not
/// persistent
/// @param expr_to_tree See eliminate_common_subexpressions
/// @param label_gen Generates the label of an intermediate from its tree and
/// a counter
/// @param min_tensors The minimum number of tensors a persistent subtree must
/// contain in order to be extracted (see also min_usage_filter)
/// @returns The positions of the inserted intermediate definitions
template <std::ranges::range VectorLike,
          std::regular_invocable<ExprPtr> Transformer>
  requires std::regular_invocable<Transformer, ResultExpr> &&
           meta::eval_node<std::ranges::range_value_t<VectorLike>>
std::vector<std::size_t> extract_persistent_subexpressions(
    VectorLike &expr_trees,
    const std::function<bool(const std::ranges::range_value_t<VectorLike> &)>
        &is_volatile_leaf,
    const Transformer &expr_to_tree,
    const std::function<std::string(
        const std::ranges::range_value_t<VectorLike> &, std::size_t)>
        &label_gen,
    std::size_t min_tensors = 2) {
  using TreeNode = std::ranges::range_value_t<VectorLike>;

  auto is_persistent = [&](const TreeNode &node) {
    std::size_t num_tensors = 0;
    return cse::count_tensors_if_persistent(node, is_volatile_leaf,
                                            num_tensors) &&
           num_tensors >= min_tensors;
  };

  // Maximal persistent subtrees (their subtrees are not visited)
  SubexpressionUsageCounts<TreeNode> persistent;
  for (const TreeNode &tree : expr_trees) {
    tree.visit_internal([&](const TreeNode &node) {
      if (!is_persistent(node)) {
        return true;
      }

      ++persistent[node];
      return false;
    });
  }

  if (persistent.empty()) {
    return {};
  }

  cse::SubexpressionReplacer<VectorLike, TreeNode, Transformer,
                             std::remove_cvref_t<decltype(label_gen)>>
      replacer(expr_trees, persistent, expr_to_tree, label_gen);
  replacer.perform_replacements(0, expr_trees.size());

  return replacer.cse_indices();
}

}  // namespace sequant::opt

#endif  // SEQUANT_COMMON_SUBEXPRESSION_ELIMINATION_HPP
