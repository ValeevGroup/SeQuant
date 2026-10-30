#ifndef SEQUANT_CORE_IRREPS_EVAL_TREE_HPP
#define SEQUANT_CORE_IRREPS_EVAL_TREE_HPP

#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/irreps/propagation.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <concepts>

namespace sequant {

/// Makes (binary) evaluation trees (EvalNode, ExportNode, ...) accessible to
/// the generic IRREP propagation
template <typename T>
  requires std::derived_from<T, EvalExpr>
struct IrrepNodeTraits<FullBinaryNode<T>> {
  static IrrepNodeKind kind(const FullBinaryNode<T> &node) {
    if (node.leaf()) {
      return IrrepNodeKind::Leaf;
    }

    SEQUANT_ASSERT(node->op_type().has_value());
    switch (node->op_type().value()) {
      case EvalOp::Sum:
        return IrrepNodeKind::Sum;
      case EvalOp::Product:
        return IrrepNodeKind::Product;
      case EvalOp::Adjoint:
        return IrrepNodeKind::Adjoint;
      case EvalOp::Power:
        // an element-wise power spans its operand's irreps: a Sum over the
        // single child (the operand) yields exactly those. That holds only
        // for a positive exponent; otherwise every symmetry-forbidden zero
        // block becomes inf or 1 and the result spans all irreps, which is
        // harmless for the totally symmetric denominators this is meant for
        return IrrepNodeKind::Sum;
      case EvalOp::Broadcast:
        // a broadcast summand spans its operand's irreps, as the summand of
        // the symbolic Sum it stands for does
        return IrrepNodeKind::Sum;
    }

    SEQUANT_UNREACHABLE;
  }

  template <typename F>
  static void for_each_child(const FullBinaryNode<T> &node, F &&f) {
    f(node.left());
    if (!is_unary(*node->op_type())) {
      // The right child of a unary node is only a sentinel
      f(node.right());
    }
  }

  template <IrrepPropagationContext Ctx>
  static typename Ctx::set_type leaf_irreps(const FullBinaryNode<T> &node,
                                            const Ctx &ctx) {
    detail::NoOpIrrepCallback no_op;
    return detail::fold_irreps<Expr>(*node->expr(), ctx, no_op);
  }
};

/// @param tree The (evaluation) tree to process
/// @param ctx The IRREP context
/// @param on_node Callback that is invoked as on_node(node, irreps) for every
/// node in the tree in post-order
/// @returns The IRREPs spanned by the result of the given tree
template <typename T, IrrepPropagationContext Ctx,
          typename OnNode = detail::NoOpIrrepCallback>
  requires std::derived_from<T, EvalExpr>
typename Ctx::set_type irreps_of(const FullBinaryNode<T> &tree, const Ctx &ctx,
                                 OnNode &&on_node = {}) {
  return detail::fold_irreps(tree, ctx, on_node);
}

}  // namespace sequant

#endif  // SEQUANT_CORE_IRREPS_EVAL_TREE_HPP
