#ifndef SEQUANT_CORE_EXPORT_EXPORT_HPP
#define SEQUANT_CORE_EXPORT_EXPORT_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_expr.hpp>
#include <SeQuant/core/export/compute_selection.hpp>
#include <SeQuant/core/export/context.hpp>
#include <SeQuant/core/export/export_expr.hpp>
#include <SeQuant/core/export/export_node.hpp>
#include <SeQuant/core/export/expression_group.hpp>
#include <SeQuant/core/export/generator.hpp>
#include <SeQuant/core/export/utils.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/result_expr.hpp>

#include <algorithm>
#include <array>
#include <iterator>
#include <optional>
#include <ranges>
#include <set>
#include <span>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <unordered_map>

namespace sequant {

namespace details {
template <typename NodeData, typename Context>
class GenerationVisitor {
 public:
  using NodeID = std::decay_t<decltype(std::declval<ExportExpr>().id())>;

  GenerationVisitor(Generator<Context> &generator, Context &ctx,
                    const std::unordered_map<NodeID, ExprPtr> &scalarFactors)
      : m_generator(generator), m_ctx(ctx), m_scalarFactors(scalarFactors) {}

  void operator()(const ExportNode<NodeData> &node, TreeTraversal context) {
    // Note the context for leaf nodes is always TreeTraversal::Any
    switch (context) {
      case TreeTraversal::PreOrder:
        load_or_create(node);
        break;
      case TreeTraversal::PostOrder:
        process_computation(node);
        break;
      case TreeTraversal::Any:
        load_or_create(node);
        break;
      case TreeTraversal::InOrder:
      case TreeTraversal::PreAndInOrder:
      case TreeTraversal::PreAndPostOrder:
      case TreeTraversal::PostAndInOrder:
      case TreeTraversal::None:
        // Should be unreachable
        assert(false);
        break;
    }

    if (!m_rootID.has_value()) {
      m_rootID = node->id();
    } else if (m_rootID.value() == node->id()) {
      if (node->is_tensor()) {
        m_generator.persist(node->as_tensor(), m_ctx);
      } else if (node->is_variable()) {
        m_generator.persist(node->as_variable(), m_ctx);
      } else {
        throw std::runtime_error(
            "Root nodes are expected to be of type tensor or variable");
      }

      m_rootID.reset();
    }
  }

  template <typename ExprType>
  void load_or_create(const ExprType &expr, const bool isLeaf) {
    static_assert(
        std::is_same_v<ExprType, Tensor> || std::is_same_v<ExprType, Variable>,
        "This function currently only works for tensors and variables");

    const LoadStrategy loadStrategy = m_ctx.loadStrategy(expr);
    const ZeroStrategy zeroStrategy = m_ctx.zeroStrategy(expr);

    const auto [usedExprBefore,
                alreadyLoaded] = [&]() -> std::tuple<bool, bool> {
      if constexpr (std::is_same_v<ExprType, Tensor>) {
        auto iter = m_tensorUses.find(expr);
        bool usedExprBefore = iter != m_tensorUses.end();
        bool alreadyLoaded = usedExprBefore && iter->second > 0;

        m_tensorUses[expr]++;

        return {usedExprBefore, alreadyLoaded};
      } else {
        auto iter = m_variableUses.find(expr);
        bool usedExprBefore = iter != m_variableUses.end();
        bool alreadyLoaded = usedExprBefore && iter->second > 0;

        m_variableUses[expr]++;

        return {usedExprBefore, alreadyLoaded};
      }
    }();

    const bool create =
        !isLeaf && !usedExprBefore && loadStrategy == LoadStrategy::Create;
    const bool load = !create && !alreadyLoaded;
    const bool zeroCreate =
        static_cast<bool>(zeroStrategy & ZeroStrategy::ZeroOnCreate);
    const bool zeroLoad =
        static_cast<bool>(zeroStrategy & ZeroStrategy::ZeroOnLoad);
    const bool zeroReuse =
        static_cast<bool>(zeroStrategy & ZeroStrategy::ZeroOnReuse);

    if (create) {
      m_generator.create(expr, zeroCreate, m_ctx);
    } else if (load) {
      m_generator.load(expr, zeroLoad, m_ctx);
    } else if (zeroReuse) {
      assert(alreadyLoaded);
      m_generator.set_to_zero(expr, m_ctx);
    }
  }

  void load_or_create(const ExportNode<NodeData> &node) {
    if (node->is_constant()) {
      // Constants don't need to be loaded/created
      return;
    }
    if (!node->is_tensor() && !node->is_variable()) {
      throw std::runtime_error(
          "Unexpected expression type in "
          "GenerationVisitor::load_or_create");
    }

    if (node->is_tensor()) {
      load_or_create<Tensor>(node->as_tensor(), node.leaf());
    } else {
      load_or_create<Variable>(node->as_variable(), node.leaf());
    }
  }

  void drop(const Expr &expr) {
    if (expr.is<Tensor>()) {
      const Tensor &tensor = expr.as<Tensor>();

      assert(m_tensorUses[tensor] > 0);
      m_tensorUses[tensor]--;

      if (m_tensorUses[tensor] == 0) {
        m_generator.unload(tensor, m_ctx);
      }
    } else if (expr.is<Variable>()) {
      const Variable &variable = expr.as<Variable>();

      assert(m_variableUses[variable] > 0);
      m_variableUses[variable]--;

      if (m_variableUses[variable] == 0) {
        m_generator.unload(variable, m_ctx);
      }
    }
  }

  void process_computation(const ExportNode<NodeData> &node) {
    assert(node->op_type() != EvalOp::Id);
    assert(!node.leaf());

    // Assemble the expression that should be evaluated
    container::svector<ExprPtr> expressions;
    switch (node->op_type()) {
      case EvalOp::Prod:
        expressions.push_back(
            ex<Product>(ExprPtrList{node.left()->expr(), node.right()->expr()},
                        Product::Flatten::No));
        break;
      case EvalOp::Sum: {
        switch (node->compute_selection()) {
          case ComputeSelection::None:
            // This node only exists for tree connectivity purposes. No explicit
            // computation that should be exported.
            return;
          case ComputeSelection::Left:
            expressions.push_back(node.left()->expr());
            break;
          case ComputeSelection::Right:
            expressions.push_back(node.right()->expr());
            break;
          case ComputeSelection::Both:
            expressions.push_back(node.left()->expr());
            expressions.push_back(node.right()->expr());
            break;
        }
        break;
      }
      case EvalOp::Id:
        throw std::runtime_error(
            "Computation must not consist of an ID operation");
    }

    assert(!expressions.empty());

    // Account for any scalar prefactor(s) and if there are
    // variables among them, make sure to load them first
    container::svector<Variable> variables;
    if (auto iter = m_scalarFactors.find(node->id());
        iter != m_scalarFactors.end()) {
      if (iter->second->template is<Product>()) {
        for (const ExprPtr &current : iter->second->template as<Product>()) {
          assert(current->is<Constant>() || current->is<Variable>());

          if (current->is<Variable>()) {
            variables.push_back(current->as<Variable>());
            if (m_variableUses[variables.back()] == 0) {
              m_generator.load(variables.back(), false, m_ctx);
            }

            m_variableUses[variables.back()]++;
          }
        }
      } else if (iter->second->template is<Variable>()) {
        variables.push_back(iter->second->template as<Variable>());

        if (m_variableUses[variables.back()] == 0) {
          m_generator.load(variables.back(), false, m_ctx);
        }

        m_variableUses[variables.back()]++;
      }

      for (ExprPtr &current : expressions) {
        current = iter->second * current;
      }
    }

    for (ExprPtr current : expressions) {
      switch (node->result_type()) {
        case ResultType::Tensor:
          m_generator.compute(*current, node->as_tensor(), m_ctx);
          break;
        case ResultType::Scalar:
          m_generator.compute(*current, node->as_variable(), m_ctx);
          break;
      }
    }

    // Drop loaded variables
    for (auto iter = variables.rbegin(); iter != variables.rend(); ++iter) {
      drop(*iter);
    }

    // Drop used leaf elements
    drop(*node.right()->expr());
    drop(*node.left()->expr());
  }

 private:
  Generator<Context> &m_generator;
  const std::unordered_map<NodeID, ExprPtr> &m_scalarFactors;
  Context &m_ctx;
  std::map<Tensor, std::size_t, TensorBlockCompare> m_tensorUses;
  std::map<Variable, std::size_t> m_variableUses;

  std::optional<NodeID> m_rootID;
};

struct PreprocessResult {
  std::unordered_map<std::size_t, ExprPtr> scalarFactors;
  std::set<Index> indices;
  std::set<Tensor, TensorBlockCompare> tensors;
  std::set<Variable> variables;

  std::map<Tensor, std::size_t, TensorBlockCompare> tensorReferences;
  std::map<Variable, std::size_t> variableReferences;

  std::size_t renameCounter = 2;
};

template <typename T>
bool prune_scalar(ExportNode<T> &node, PreprocessResult &result) {
  if (!node.leaf() || !node->is_scalar()) {
    return false;
  }

  assert(!node.root());

  const auto iter = result.scalarFactors.find(node.parent()->id());
  ExprPtr parentFactor =
      iter == result.scalarFactors.end() ? nullptr : iter->second;

  if (iter != result.scalarFactors.end()) {
    // The parent will be overwritten, so we don't have to keep
    // this entry around
    result.scalarFactors.erase(iter);
  }

  ExprPtr factor = std::move(node->expr());

  assert(factor);
  assert(factor->is<Constant>() || factor->is<Variable>());

  if (factor->is<Variable>()) {
    result.variables.insert(factor->as<Variable>());
  }

  if (parentFactor) {
    factor = parentFactor * factor;
  }

  ExportNode<T> fill_in = [&]() {
    if (node->id() == node.parent().left()->id()) {
      return std::move(node.parent().right());
    } else {
      assert(node->id() == node.parent().right()->id());
      return std::move(node.parent().left());
    }
  }();

  if (!fill_in.leaf()) {
    fill_in->set_expr(std::move(node.parent()->expr()));
  } else if (!node.parent().root()) {
    if (node.parent()->id() == node.parent().parent().left()->id()) {
      node.parent().parent()->select_left();
    } else {
      assert(node.parent()->id() == node.parent().parent().right()->id());
      node.parent().parent()->select_right();
    }
  }

  result.scalarFactors[fill_in->id()] = std::move(factor);

  // We can only do the moving after finishing access to node as
  // node might become unreferenced and thus deleted once its parent
  // is overwritten.
  node.parent() = std::move(fill_in);

  return true;
}

void rename(Tensor &tensor, std::size_t counter);
void rename(Variable &variable, std::size_t counter);

template <typename ExprType, typename Node>
void preprocess(ExprType expr, ExportContext &ctx, Node &node,
                PreprocessResult &result) {
  static_assert(
      std::is_same_v<ExprType, Tensor> || std::is_same_v<ExprType, Variable>,
      "This function currently only works for tensors and variables");

  bool storeExpr = false;

  storeExpr |= ctx.rewrite(expr);

  if (node.leaf()) {
    // Anything but loading doesn't make sense for leaf nodes
    ctx.setLoadStrategy(expr, LoadStrategy::Load);
    ctx.setZeroStrategy(expr, ZeroStrategy::NeverZero);
  } else {
    const auto [usedBefore, currentlyLoaded] = [&]() -> std::tuple<bool, bool> {
      if constexpr (std::is_same_v<ExprType, Tensor>) {
        auto iter = result.tensorReferences.find(expr);
        const bool usedBefore = iter != result.tensorReferences.end();
        const bool currentlyLoaded = usedBefore && iter->second > 0;

        return {usedBefore, currentlyLoaded};
      } else {
        auto iter = result.variableReferences.find(expr);
        const bool usedBefore = iter != result.variableReferences.end();
        const bool currentlyLoaded = usedBefore && iter->second > 0;

        return {usedBefore, currentlyLoaded};
      }
    }();

    // If this result is going to be summed, we don't do any special result
    // reuse handling as we instead rely on the += semantics that allows us to
    // keep adding sub-results to the final sum without the need for explicit
    // intermediates for each step of the summation.
    // Also, we don't want to mess with the final result tensor.
    const bool handleReuse =
        !node.root() && node.parent()->op_type() != EvalOp::Sum;

    // Special handling for result objects that have been used before
    // However, for any object that is the result of a sum,
    // we don't want this special behavior.
    if (handleReuse && usedBefore) {
      assert(!node.leaf());

      if (currentlyLoaded) {
        // This expr is currently in use -> can't use it as a result as that
        // would overwrite the currently used value -> need to rename expr
        rename(expr, result.renameCounter++);
        storeExpr = true;
      } else {
        // This expr is reused to store a new result -> set it to zero before
        // adding new result
        ctx.setZeroStrategy(expr, ZeroStrategy::AlwaysZero);
      }
    }
  }

  if (storeExpr) {
    node->set_expr(ex<ExprType>(expr));
  }

  if constexpr (std::is_same_v<ExprType, Tensor>) {
    result.tensors.insert(expr);
    const auto &indices = expr.const_indices();
    result.indices.insert(indices.begin(), indices.end());

    result.tensorReferences[expr]++;
  } else {
    result.variables.insert(expr);

    result.variableReferences[expr]++;
  }
}

///
/// Preprocesses the provided binary tree by
/// - removing explicit appearances of scalar leafs. We don't want them to
///   be represented in the tree. Instead, we keep track of them in a different
///   way in order to be able to give scalar factors alongside the actual
///   tensor contraction they are supposed to scale (this is necessary
///   as there are backends which only support scaling in this context)
/// - rebalance the tree such that for any given non-leaf node, its left
///   subtree is always larger (or equally large) than its right one.
///   This ensures that we have to have the least amount of tensors loaded
///   at the same time, when generating code for a backend which only supports
///   stack-like memory allocations (e.g. when A is allocated before B, B
///   must be deleted before A can be deleted).
/// - Rename intermediate tensors that have the same name and describe the same
///   tensor block, which are required as two separate entities at the same
///   time when evaluating the tree (thus a single tensor object is
///   insufficient).
/// While doing so, we also collect meta information about the expression(s)
/// encoded in this tree.
///
template <typename T>
void preprocess(ExportNode<T> &tree, PreprocessResult &result,
                ExportContext &ctx) {
  // Pruning can only be done if
  // - tree is not a leaf (if we prune that, the tree has vanished)
  // - pruning the tree does NOT turn it into a single-leaf tree
  // - if it does turn it into a leaf, then parent must be non-zero
  //   (so we can store the scalar factor on that)
  while (!tree.leaf() &&
         (!tree.root() || !tree.left().leaf() || !tree.right().leaf()) &&
         prune_scalar(tree.left(), result)) {
    // In case the pruning led to tree becoming a leaf, we have to move the
    // pruned scalar factor out to its parent in order to be properly accounted
    // for (as leafs only get loaded and never computed)
    if (auto iter = result.scalarFactors.find(tree->id());
        iter != result.scalarFactors.end() && tree.leaf()) {
      assert(!tree.root());
      ExprPtr factor = std::move(iter->second);
      result.scalarFactors.erase(iter);
      result.scalarFactors[tree.parent()->id()] = std::move(factor);
    }
  }
  while (!tree.leaf() &&
         (!tree.root() || !tree.left().leaf() || !tree.right().leaf()) &&
         prune_scalar(tree.right(), result)) {
    if (auto iter = result.scalarFactors.find(tree->id());
        iter != result.scalarFactors.end() && tree.leaf()) {
      assert(!tree.root());
      ExprPtr factor = std::move(iter->second);
      result.scalarFactors.erase(iter);
      result.scalarFactors[tree.parent()->id()] = std::move(factor);
    }
  }

  if (tree->is_tensor()) {
    preprocess<Tensor>(tree->as_tensor(), ctx, tree, result);
  } else if (tree->is_variable()) {
    preprocess<Variable>(tree->as_variable(), ctx, tree, result);
  }

  if (tree.leaf()) {
    return;
  }

  if (tree.left().size() < tree.right().size()) {
    std::swap(tree.left(), tree.right());
  }

  if (tree->op_type() == EvalOp::Sum) {
    // We don't want to explicitly encode addition of non-leaf
    // nodes in the tree as that could lead to unnecessary intermediates
    // being created. Instead, we flush the top-most result of the
    // addition downwards, making use of the += semantic that is assumed
    // for all computations.
    // Note the explicit flushing down of the result name is required in
    // case the top-level summation node has a different name than the
    // intermediate nodes.
    ComputeSelection selection = ComputeSelection::Both;
    if (!tree.left().leaf()) {
      tree.left()->set_expr(tree->expr());
      selection &= ~ComputeSelection::Left;
    }
    if (!tree.right().leaf()) {
      tree.right()->set_expr(tree->expr());
      selection &= ~ComputeSelection::Right;
    }

    tree->set_compute_selection(selection);
  }

  preprocess(tree.left(), result, ctx);
  preprocess(tree.right(), result, ctx);

  // Mark tensors/variables as no longer in use
  if (tree.left()->is_tensor()) {
    const Tensor &tensor = tree.left()->as_tensor();
    assert(result.tensorReferences[tensor] > 0);
    result.tensorReferences[tensor]--;
  } else if (tree.left()->is_variable()) {
    const Variable &variable = tree.left()->as_variable();
    assert(result.variableReferences[variable] > 0);
    result.variableReferences[variable]--;
  }
  if (tree.right()->is_tensor()) {
    const Tensor &tensor = tree.right()->as_tensor();
    assert(result.tensorReferences[tensor] > 0);
    result.tensorReferences[tensor]--;
  } else if (tree.right()->is_variable()) {
    const Variable &variable = tree.right()->as_variable();
    assert(result.variableReferences[variable] > 0);
    result.variableReferences[variable]--;
  }
}

template <typename T>
void preprocess_and_maybe_log(ExportNode<T> &tree, PreprocessResult &result,
                              ExportContext &ctx) {
  if (true || Logger::instance().export_equations) {
    std::cout << "Tree before preprocessing:\n"
              << tree.tikz(
                     [](const ExportNode<T> &node) {
                       return "$" + toUtf8(to_latex(node->expr())) + "$";
                     },
                     [](const ExportNode<T>) -> std::string { return ""; })
              << "\n";
  }

  preprocess(tree, result, ctx);

  if (true || Logger::instance().export_equations) {
    std::cout << "Tree after pre-processing:\n"
              << tree.tikz(
                     [](const ExportNode<T> &node) {
                       return "$" + toUtf8(to_latex(node->expr())) + "$";
                     },
                     [](const ExportNode<T>) -> std::string { return ""; })
              << "\n";
  }
}

template <typename T, typename Context>
void export_expression(ExportNode<T> &expression, Generator<Context> &generator,
                       Context &ctx, PreprocessResult &pp_result) {
  ctx.set_current_expression_id(expression->id());

  generator.begin_expression(ctx);

  details::GenerationVisitor<T, Context> visitor(generator, ctx,
                                                 pp_result.scalarFactors);
  expression.visit(
      [&visitor](const FullBinaryNode<T> &node, TreeTraversal context) {
        visitor(node, context);
      },
      TreeTraversal::PreAndPostOrder);

  generator.end_expression(ctx);

  ctx.clear_current_expression_id();
}

template <typename Range, typename Context>
  requires std::ranges::range<Range>
void declare_all(const Range &range, Generator<Context> &generator,
                 Context &ctx) {
  using RangeType = std::ranges::range_value_t<Range>;

  for (const RangeType &val : range) {
    generator.declare(val, ctx);
  }

  using std::ranges::size;
  if constexpr (std::is_same_v<RangeType, Tensor>) {
    generator.all_tensors_declared(size(range), ctx);
  } else if constexpr (std::is_same_v<RangeType, Variable>) {
    generator.all_variables_declared(size(range), ctx);
  } else {
    static_assert(std::is_same_v<RangeType, Index>);
    generator.all_indices_declared(size(range), ctx);
  }
}

template <typename T, typename Compare = std::less<T>, typename Range>
  requires std::ranges::range<Range> &&
           std::is_same_v<std::ranges::range_value_t<Range>, PreprocessResult>
std::set<T, Compare> combine_and_clear_pp_results(Range &&range) {
  std::set<T, Compare> combined;

  constexpr bool is_index = std::is_same_v<T, Index>;
  constexpr bool is_variable = std::is_same_v<T, Variable>;
  constexpr bool is_tensor = std::is_same_v<T, Tensor>;
  static_assert(is_index || is_variable || is_tensor);

  for (PreprocessResult &current : range) {
    if constexpr (is_index) {
      combined.insert(std::make_move_iterator(current.indices.begin()),
                      std::make_move_iterator(current.indices.end()));
      current.indices.clear();
    } else if constexpr (is_variable) {
      combined.insert(std::make_move_iterator(current.variables.begin()),
                      std::make_move_iterator(current.variables.end()));
      current.variables.clear();
    } else if constexpr (is_tensor) {
      combined.insert(std::make_move_iterator(current.tensors.begin()),
                      std::make_move_iterator(current.tensors.end()));
      current.tensors.clear();
    }
  }

  return combined;
}

template <DeclarationScope scope, typename Range, typename Context>
  requires std::ranges::range<Range> && (!std::is_const_v<Range>)
void handle_declarations(Range &&range, Generator<Context> &generator,
                         Context &ctx) {
  generator.begin_declarations(scope, ctx);

  if (generator.index_declaration_scope() == scope) {
    std::set<Index> indices =
        details::combine_and_clear_pp_results<Index>(range);
    details::declare_all(indices, generator, ctx);
  }

  if (generator.variable_declaration_scope() == scope) {
    std::set<Variable> variables =
        details::combine_and_clear_pp_results<Variable>(range);
    details::declare_all(variables, generator, ctx);
  }

  if (generator.tensor_declaration_scope() == scope) {
    std::set<Tensor, TensorBlockCompare> tensors =
        details::combine_and_clear_pp_results<Tensor, TensorBlockCompare>(
            range);
    details::declare_all(tensors, generator, ctx);
  }

  generator.end_declarations(scope, ctx);
}

}  // namespace details

template <typename T = ExportExpr, typename Context>
void export_group(ExpressionGroup<T> group, Generator<Context> &generator,
                  Context ctx = {}) {
  export_groups<T, Context>(std::array{std::move(group)}, generator,
                            std::move(ctx));
}

template <typename T = ExportExpr, typename Context>
void export_expression(ExportNode<T> expression, Generator<Context> &generator,
                       Context ctx = {}) {
  export_groups<T, Context>(
      std::array{ExpressionGroup<T>{std::move(expression)}}, generator,
      std::move(ctx));
}

template <typename T = ExportExpr, typename Context, typename Range>
  requires std::ranges::range<std::remove_cvref_t<Range>> &&
           std::is_same_v<std::ranges::range_value_t<Range>,
                          ExpressionGroup<T>> &&
           (!std::is_const_v<Range>)
void export_groups(Range groups, Generator<Context> &generator,
                   Context ctx = {}) {
  using std::ranges::size;

  if (size(groups) > 1) {
    if (!generator.supports_named_sections()) {
      throw std::runtime_error("The generator for '" +
                               generator.get_format_name() +
                               "' doesn't support named sections");
    }

    for (const ExpressionGroup<T> &current : groups) {
      if (!current.is_named()) {
        throw std::runtime_error(
            "Can't have unnamed groups when exporting multiple groups at once");
      }
    }
  }

  generator.begin_export(ctx);

  const DeclarationScope index_decl_scope = generator.index_declaration_scope();
  const DeclarationScope variable_decl_scope =
      generator.variable_declaration_scope();
  const DeclarationScope tensor_decl_scope =
      generator.tensor_declaration_scope();

  // First step: preprocessing of all expressions
  container::svector<details::PreprocessResult> pp_results;
  for (ExpressionGroup<T> &current_group : groups) {
    pp_results.reserve(pp_results.size() + size(groups));

    for (ExportNode<T> &current_tree : current_group) {
      pp_results.emplace_back();

      ctx.set_current_expression_id(current_tree->id());

      details::preprocess_and_maybe_log(current_tree, pp_results.back(), ctx);

      ctx.clear_current_expression_id();
    }
  }

  // Perform global declarations
  details::handle_declarations<DeclarationScope::Global>(pp_results, generator,
                                                         ctx);

  // Now initiate the actual code generation
  std::size_t pp_idx = 0;
  for (ExpressionGroup<T> &current_group : groups) {
    std::string section_name;
    if (generator.supports_named_sections()) {
      if (current_group.is_named()) {
        ctx.set_current_section_name(current_group.name());
        generator.begin_named_section(current_group.name(), ctx);
        section_name = current_group.name();
      } else if (generator.requires_named_sections()) {
        section_name = "unnamed";
        ctx.set_current_section_name(section_name);
        generator.begin_named_section(section_name, ctx);
      }
    }

    // Handle section-level declarations
    assert(pp_results.size() >= pp_idx + size(current_group));
    details::handle_declarations<DeclarationScope::Section>(
        std::span(&pp_results.at(pp_idx), size(current_group)), generator, ctx);

    for (ExportNode<T> &current_tree : current_group) {
      // Handle expression-level declarations
      details::handle_declarations<DeclarationScope::Expression>(
          std::span{&pp_results.at(pp_idx), 1}, generator, ctx);

      details::export_expression(current_tree, generator, ctx,
                                 pp_results.at(pp_idx));

      pp_idx++;
    }

    if (!section_name.empty()) {
      generator.end_named_section(section_name, ctx);
    }

    ctx.clear_current_section_name();
  }

  assert(pp_idx == pp_results.size());

  generator.end_export(ctx);
}

template <typename NodeData = ExportExpr>
ExportNode<NodeData> to_export_tree(const ExprPtr &expr) {
  return binarize<NodeData>(expr);
}

template <typename NodeData = ExportExpr>
ExportNode<NodeData> to_export_tree(const ResultExpr &expr) {
  return binarize<NodeData>(expr);
}

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_EXPORT_HPP
