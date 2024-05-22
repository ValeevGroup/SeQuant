#ifndef SEQUANT_CORE_EXPORT_EXPORT_HPP
#define SEQUANT_CORE_EXPORT_EXPORT_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/export/context.hpp>
#include <SeQuant/core/export/generator.hpp>
#include <SeQuant/core/export/utils.hpp>
#include <SeQuant/core/logger.hpp>

#include <functional>
#include <optional>
#include <stack>
#include <stdexcept>
#include <unordered_map>

namespace sequant {

namespace {
template <typename NodeData, typename Context>
class GenerationVisitor {
 public:
  using NodeID = std::decay_t<decltype(std::declval<EvalExpr>().id())>;

  GenerationVisitor(Generator<Context> &generator, Context &ctx,
                    const std::unordered_map<NodeID, ExprPtr> &scalarFactors)
      : m_generator(generator), m_ctx(ctx), m_scalarFactors(scalarFactors) {}

  void operator()(const EvalNode<NodeData> &node, TreeTraversal context) {
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

  void load_or_create(const EvalNode<NodeData> &node) {
    if (node->is_tensor()) {
      const Tensor &tensor = node->as_tensor();

      const LoadStrategy loadStrategy = m_ctx.loadStrategy(tensor);
      const ZeroStrategy zeroStrategy = m_ctx.zeroStrategy(tensor);

      const auto iter = m_tensorUses.find(tensor);
      const bool usedTensorBefore = iter != m_tensorUses.end();
      const bool isLeaf = node.leaf();
      const bool alreadyLoaded = usedTensorBefore && iter->second > 0;
      const bool create =
          !isLeaf && !usedTensorBefore && loadStrategy == LoadStrategy::Create;
      const bool load = !create && !alreadyLoaded;
      const bool zeroCreate = !!(zeroStrategy & ZeroStrategy::ZeroOnCreate);
      const bool zeroLoad = !!(zeroStrategy & ZeroStrategy::ZeroOnLoad);
      const bool zeroReuse = !!(zeroStrategy & ZeroStrategy::ZeroOnReuse);

      m_tensorUses[tensor]++;

      if (create) {
        m_generator.create(tensor, zeroCreate, m_ctx);
      } else if (load) {
        m_generator.load(tensor, zeroLoad, m_ctx);
      } else if (zeroReuse) {
        assert(alreadyLoaded);
        m_generator.set_to_zero(tensor, m_ctx);
      }
    } else if (node->is_variable()) {
      const Variable &variable = node->as_variable();

      if (m_variableUses[variable] == 0) {
        m_generator.load(variable, m_ctx);
      }

      m_variableUses[variable]++;
    } else if (node->is_constant()) {
      // Constants don't need to be loaded/created
    } else {
      throw std::runtime_error(
          "Unexpected expression type in "
          "GenerationVisitor::load_or_create");
    }
  }

  void process_computation(const EvalNode<NodeData> &node) {
    assert(node->op_type() != EvalOp::Id);
    assert(!node.leaf());

    // Assemble the expression that should be evaluated
    ExprPtr expression;
    switch (node->op_type()) {
      case EvalOp::Prod:
        expression =
            ex<Product>(ExprPtrList{node.left()->expr(), node.right()->expr()},
                        Product::Flatten::No);
        break;
      case EvalOp::Sum:
        expression =
            ex<Sum>(ExprPtrList{node.left()->expr(), node.right()->expr()});
        break;
      case EvalOp::Id:
        throw std::runtime_error(
            "Computation must not consist of an ID operation");
    }

    assert(expression);

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
            m_generator.load(variables.back(), m_ctx);
          }
        }
      } else if (iter->second->template is<Variable>()) {
        variables.push_back(iter->second->template as<Variable>());
        m_generator.load(variables.back(), m_ctx);
      }

      expression = iter->second * expression;
    }

    switch (node->result_type()) {
      case ResultType::Tensor:
        m_generator.compute(*expression, node->as_tensor(), m_ctx);
        break;
      case ResultType::Scalar:
        m_generator.compute(*expression, node->as_variable(), m_ctx);
        break;
    }

    // Drop loaded variables
    for (auto iter = variables.rbegin(); iter != variables.rend(); ++iter) {
      m_generator.unload(*iter, m_ctx);
    }

    // Drop used leaf tensors
    assert(node.left()->is_tensor());
    assert(node.right()->is_tensor());
    const Tensor &lhs = node.left()->as_tensor();
    const Tensor &rhs = node.right()->as_tensor();

    assert(m_tensorUses[rhs] > 0);
    m_tensorUses[rhs]--;

    if (m_tensorUses[rhs] == 0) {
      m_generator.unload(rhs, m_ctx);
    }

    assert(m_tensorUses[lhs] > 0);
    m_tensorUses[lhs]--;

    if (m_tensorUses[lhs] == 0) {
      m_generator.unload(lhs, m_ctx);
    }
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
  std::map<Variable, bool> usedVariables;
};

template <typename T>
bool prune_scalar(EvalNode<T> &node, EvalNode<T> &parent,
                  PreprocessResult &result) {
  if (!node.leaf() || !node->is_scalar()) {
    return false;
  }

  const auto iter = result.scalarFactors.find(parent->id());
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

  // We can only do the moving after finishing access to node as
  // node might become unreferenced and thus deleted once its parent
  // is overwritten.
  if (node->id() == parent.left()->id()) {
    parent = std::move(parent.right());
  } else {
    assert(node->id() == parent.right()->id());
    parent = std::move(parent.left());
  }

  result.scalarFactors[parent->id()] = std::move(factor);

  return true;
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
///   stack-like memory allocations (e.g. when A is allocated before B, A
///   must be deleted before B can be deleted).
/// - Rename intermediate tensors that have the same name and describe the same
///   tensor block, which are required as two separate entities at the same
///   time when evaluating the tree (thus a single tensor object is
///   insufficient).
/// While doing so, we also collect meta information about the expression(s)
/// encoded in this tree.
///
template <typename T>
void preprocess(EvalNode<T> &tree, PreprocessResult &result, ExportContext &ctx,
                EvalNode<T> *parent = nullptr) {
  while (!tree.leaf() && prune_scalar(tree.left(), tree, result)) {
    // In case the pruning led to tree becoming a leaf, we have to move the
    // pruned scalar factor out to its parent in order to be properly accounted
    // for (as leafs only get loaded and never computed)
    if (auto iter = result.scalarFactors.find(tree->id());
        iter != result.scalarFactors.end() && tree.leaf()) {
      assert(parent);
      ExprPtr factor = std::move(iter->second);
      result.scalarFactors.erase(iter);
      result.scalarFactors[(*parent)->id()] = std::move(factor);
    }
  }
  while (!tree.leaf() && prune_scalar(tree.right(), tree, result)) {
    if (auto iter = result.scalarFactors.find(tree->id());
        iter != result.scalarFactors.end() && tree.leaf()) {
      assert(parent);
      ExprPtr factor = std::move(iter->second);
      result.scalarFactors.erase(iter);
      result.scalarFactors[(*parent)->id()] = std::move(factor);
    }
  }

  static std::size_t renameCounter = 2;
  if (tree->is_tensor()) {
    Tensor tensor = tree->as_tensor();

    if (ctx.rewrite(tensor)) {
      tree->set_expr(tensor.shared_from_this());
    }

    if (tree.leaf()) {
      // Anything but loading doesn't make sense for leaf nodes
      ctx.setLoadStrategy(tensor, LoadStrategy::Load);
      ctx.setZeroStrategy(tensor, ZeroStrategy::NeverZero);
    } else {
      auto iter = result.tensorReferences.find(tensor);

      if (iter != result.tensorReferences.end()) {
        assert(!tree.leaf());
        if (iter->second > 0) {
          // This tensor is already in use -> can't use it as a result tensor as
          // that would overwrite the currently used instance -> need to rename
          // the result tensor
          tree->set_expr(ex<Tensor>(
              std::wstring(tensor.label()) + std::to_wstring(renameCounter++),
              tensor.bra(), tensor.ket(), tensor.symmetry(),
              tensor.braket_symmetry(), tensor.particle_symmetry()));
          tensor = tree->as_tensor();
        } else {
          // This tensor is reused -> need to zero out previous values when
          // reusing it to hold a new result
          ctx.setZeroStrategy(tensor, ZeroStrategy::ZeroOnLoad |
                                          ZeroStrategy::ZeroOnReuse |
                                          ZeroStrategy::ZeroOnCreate);
        }
      }
    }

    result.tensors.insert(tensor);
    const auto &indices = tensor.const_braket();
    result.indices.insert(indices.begin(), indices.end());

    result.tensorReferences[tensor]++;
  } else if (tree->is_variable()) {
    Variable variable = tree->as_variable();

    if (ctx.rewrite(variable)) {
      tree->set_expr(variable.shared_from_this());
    }

    if (result.usedVariables[variable]) {
      assert(!tree.leaf());

      tree->set_expr(ex<Variable>(std::wstring(variable.label()) +
                                  std::to_wstring(renameCounter++)));
      assert(tree->expr());

      variable = tree->as_variable();
    }

    result.variables.insert(variable);

    result.usedVariables[variable] = true;
  }

  if (tree.leaf()) {
    return;
  }

  if (tree.left().size() < tree.right().size()) {
    std::swap(tree.left(), tree.right());
  }

  preprocess(tree.left(), result, ctx, &tree);
  preprocess(tree.right(), result, ctx, &tree);

  // Mark tensors/variables as no longer in use
  if (tree.left()->is_tensor()) {
    const Tensor &tensor = tree.left()->as_tensor();
    assert(result.tensorReferences[tensor] > 0);
    result.tensorReferences[tensor]--;
  } else if (tree.left()->is_variable()) {
    result.usedVariables[tree.left()->as_variable()] = false;
  }
  if (tree.right()->is_tensor()) {
    const Tensor &tensor = tree.right()->as_tensor();
    assert(result.tensorReferences[tensor] > 0);
    result.tensorReferences[tensor]--;
  } else if (tree.right()->is_variable()) {
    result.usedVariables[tree.right()->as_variable()] = false;
  }
}

}  // namespace

template <typename T, typename Context>
void export_expression(EvalNode<T> expression, Generator<Context> &generator,
                       Context ctx = {}) {
  PreprocessResult result;
  preprocess(expression, result, ctx);

  if (Logger::instance().export_equations) {
    Logger &logger = Logger::instance();
    *logger.stream << "Pre-processed equation tree:\n"
                   << expression.template tikz<std::string>(
                          [](const EvalNode<T> &node) {
                            return "$" + toUtf8(to_latex(node->expr())) + "$";
                          },
                          [](const EvalNode<T>) -> std::string { return ""; })
                   << "\n";
  }

  for (const Index &idx : result.indices) {
    generator.declare(idx, ctx);
  }

  if (!result.indices.empty()) {
    generator.insert_blank_lines(1, ctx);
  }

  for (const Variable &var : result.variables) {
    generator.declare(var, ctx);
  }

  if (!result.variables.empty()) {
    generator.insert_blank_lines(1, ctx);
  }

  for (const Tensor &tensor : result.tensors) {
    generator.declare(tensor, ctx);
  }

  if (!result.tensors.empty()) {
    generator.insert_blank_lines(1, ctx);
  }

  GenerationVisitor<T, Context> visitor(generator, ctx, result.scalarFactors);
  expression.visit(
      [&visitor](const FullBinaryNode<T> &node, TreeTraversal context) {
        visitor(node, context);
      },
      TreeTraversal::PreAndPostOrder);
}

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_EXPORT_HPP
