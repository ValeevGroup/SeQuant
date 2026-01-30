#ifndef SEQUANT_COMMON_SUBEXPRESSION_ELIMINATION_HPP
#define SEQUANT_COMMON_SUBEXPRESSION_ELIMINATION_HPP

#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/tensor.hpp>

#include <concepts>
#include <format>
#include <optional>
#include <ranges>
#include <string>
#include <unordered_map>

namespace sequant {

namespace cse {

/// Functor to compute the hash of a given (evaluation) tree node
template <typename TreeNode, bool force_hash_collisions = false>
struct TreeNodeHasher {
  /// Trait used by the C++ STL allowing heterogenous lookups
  using is_transparent = void;

  std::size_t operator()(const TreeNode *node) const { return (*this)(*node); }

  std::size_t operator()(const TreeNode &node) const {
    if constexpr (force_hash_collisions) {
      return 0;
    }

    return hash::value(*node);
  }
};

// Functor to compare two trees for equivalence
// Explicit equivalence checking mitigates (accidental) hash collisions
template <typename TreeNode>
struct TreeNodeEqualityComparator {
  /// Trait used by the C++ STL allowing heterogenous lookups
  using is_transparent = void;

  bool operator()(const TreeNode *lhs, const TreeNode *rhs) const {
    return (*this)(*lhs, *rhs);
  }

  bool operator()(const TreeNode &lhs, const TreeNode *rhs) const {
    return (*this)(lhs, *rhs);
  }

  bool operator()(const TreeNode *lhs, const TreeNode &rhs) const {
    return (*this)(*lhs, rhs);
  }

  bool operator()(const TreeNode &lhs, const TreeNode &rhs) const {
    if (lhs.leaf() != rhs.leaf()) {
      return false;
    }

    if (lhs.size() != rhs.size()) {
      return false;
    }

    if (hash::value(*lhs) != hash::value(*rhs)) {
      return false;
    }

    if (lhs->is_scalar() != rhs->is_scalar() ||
        lhs->is_tensor() != rhs->is_tensor() ||
        lhs->is_constant() != rhs->is_constant() ||
        lhs->is_variable() != rhs->is_variable() ||
        lhs->is_product() != rhs->is_product() ||
        lhs->is_sum() != rhs->is_sum()) {
      return false;
    }

    if (lhs->is_constant() || lhs->is_variable()) {
      if (*lhs->expr() != *rhs->expr()) {
        return false;
      }
    } else if (lhs->is_tensor()) {
      const Tensor &lhs_tensor = lhs->as_tensor();
      const Tensor &rhs_tensor = rhs->as_tensor();

      TensorBlockEqualComparator cmp;
      if (!cmp(lhs_tensor, rhs_tensor)) {
        return false;
      }
    }

    if (lhs->has_connectivity_graph() != rhs->has_connectivity_graph()) {
      return false;
    }

    // Check connectivity in products / contractions
    if (lhs->has_connectivity_graph()) {
      SEQUANT_ASSERT(lhs->has_connectivity_graph());
      SEQUANT_ASSERT(rhs->has_connectivity_graph());

      if (bliss::ConstGraphCmp::cmp(lhs->connectivity_graph(),
                                    rhs->connectivity_graph()) != 0) {
        return false;
      }
    }

    if (!lhs.leaf()) {
      // Note: We're assuming that the assignment of a subtree into left and
      // right is made consistently (canonical) and hence we don't compare
      // left with right
      if (!(*this)(lhs.left(), rhs.left())) {
        return false;
      }
      if (!(*this)(lhs.right(), rhs.right())) {
        return false;
      }
    }

    return true;
  }
};

/// A map between (sub)tree hashes and how often they have been found
/// This is identical to SubexpressionUsageCounts except that we store node
/// pointers here (lower memory footprint but higher risk of dangling pointers)
template <typename TreeNode, bool force_hash_collisions = false>
using SubexpressionHashCollector =
    std::unordered_map<const TreeNode *, std::size_t,
                       TreeNodeHasher<TreeNode, force_hash_collisions>,
                       TreeNodeEqualityComparator<TreeNode>>;

/// A map between (sub)trees and how often they have been found
template <typename TreeNode, bool force_hash_collisions = false>
using SubexpressionUsageCounts =
    std::unordered_map<TreeNode, std::size_t,
                       TreeNodeHasher<TreeNode, force_hash_collisions>,
                       TreeNodeEqualityComparator<TreeNode>>;

/// A map between (sub)trees and the name chosen to represent the associated
/// intermediate
template <typename TreeNode, bool force_hash_collisions = false>
using SubexpressionNames =
    std::unordered_map<TreeNode, std::wstring,
                       TreeNodeHasher<TreeNode, force_hash_collisions>,
                       TreeNodeEqualityComparator<TreeNode>>;

/// Functor that can be used to identify common subexpressions while iterating
/// expression trees
template <typename TreeNode, bool force_hash_collisions = false>
class SubexpressionIdentifier {
 public:
  SubexpressionIdentifier() = default;

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
    SubexpressionUsageCounts<TreeNode, force_hash_collisions> usages;
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

  void perform_replacements() {
    // Ensure we won't have to reallocate while iterating over the expressions
    expr_trees.reserve(expr_trees.size() + subexpressions.size());

    const std::size_t orig_tree_count = expr_trees.size();

    for (std::size_t i = 0; i < orig_tree_count; ++i) {
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

    std::wstring label = label_it == end(cse_names) ? label_gen(name_counter++)
                                                    : label_it->second;

    ExprPtr expr = [&]() {
      if (tree->is_tensor()) {
        return ex<Tensor>(label, bra(), ket(), aux(tree->canon_indices()),
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

      expr_trees.insert(expr_it, std::move(intermediate_definition.value()));
      expr_offset += 1;
    }

    return false;
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
};

// TODO: provide predicate that computes the FLOPs required to compute the CSE
// and compares that with the expected speed for reading the CSE from disk.
// Use typical FLOPs/s and I/O rates from e.g. https://ssd.userbenchmark.com/
// https://openbenchmarking.org/test/pts/mt-dgemm
// to determine whether or not pre-computing the CSE is worth it.
// Note: of course the number of reusages of the CSE plays into this.
struct AcceptAllPredicate {
  template <typename... Ts>
  bool operator()(const Ts &...) const {
    return true;
  }
};

}  // namespace cse

/// Takes the range of expression trees and performs common subexpression
/// elimination on them. Evaluation trees for intermediates are inserted into
/// the given container as needed.
///
/// @param expr_trees Container (vector-like) of evaluation trees
/// @param expr_to_tree Functor that can convert ExprPtr and ResultExpr objects
/// into suitable evaluation trees
/// @param filter_predicate A predicate to filter out subexpressions that shall
/// not be eliminated
template <std::ranges::range VectorLike,
          std::regular_invocable<ExprPtr> Transformer,
          std::predicate<std::ranges::range_value_t<VectorLike>, std::size_t>
              Filter = cse::AcceptAllPredicate,
          bool force_hash_collisions = false>
  requires std::regular_invocable<Transformer, ResultExpr> &&
           meta::eval_node<std::ranges::range_value_t<VectorLike>>
void eliminate_common_subexpressions(VectorLike &expr_trees,
                                     const Transformer &expr_to_tree,
                                     const Filter &filter_predicate = {}) {
  using std::ranges::begin;
  using std::ranges::end;

  using TreeNode = std::ranges::range_value_t<VectorLike>;

  cse::SubexpressionIdentifier<TreeNode, force_hash_collisions> identifier;

  for (const TreeNode &current : expr_trees) {
    current.visit_internal(identifier);
  }

  cse::SubexpressionUsageCounts<TreeNode, force_hash_collisions>
      subexpression_usages = identifier.take_subexpression_map();

  // Apply filter
  std::erase_if(subexpression_usages, [&filter_predicate](const auto &pair) {
    return !filter_predicate(pair.first, pair.second);
  });

  // TODO: Make label generator be passable from the outside
  auto label_gen = [](std::size_t c) { return std::format(L"CSE{}", c); };

  cse::SubexpressionReplacer<VectorLike, TreeNode, Transformer,
                             decltype(label_gen), force_hash_collisions>
      replacer(expr_trees, subexpression_usages, expr_to_tree, label_gen);

  replacer.perform_replacements();
}

}  // namespace sequant

#endif  // SEQUANT_COMMON_SUBEXPRESSION_ELIMINATION_HPP
