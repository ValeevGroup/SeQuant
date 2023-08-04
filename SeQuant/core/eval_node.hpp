//
// Created by Bimal Gaudel on 5/24/21.
//

#ifndef SEQUANT_EVAL_NODE_HPP
#define SEQUANT_EVAL_NODE_HPP

#include "asy_cost.hpp"
#include "binary_node.hpp"
#include "eval_expr.hpp"

#include <SeQuant/core/math.hpp>

namespace sequant {

namespace {

template <typename, typename = void>
constexpr bool IsEvalExpr{};

template <typename T>
constexpr bool
    IsEvalExpr<T, std::enable_if_t<std::is_convertible_v<T, EvalExpr>>>{true};

template <typename, typename = void>
constexpr bool IsEvalNode{};

template <typename T>
constexpr bool IsEvalNode<FullBinaryNode<T>, std::enable_if_t<IsEvalExpr<T>>>{
    true};

template <typename T>
constexpr bool
    IsEvalNode<const FullBinaryNode<T>, std::enable_if_t<IsEvalExpr<T>>>{true};

}  // namespace

template <typename T,
          typename = std::enable_if_t<std::is_convertible_v<T, EvalExpr>>>
using EvalNode = FullBinaryNode<T>;

///
/// \brief Creates an evaluation tree from @c ExprPtr.
///
/// \tparam ExprT Can be @c EvalExpr or derived from it.
///
/// \param expr The expression to be converted to an evaluation tree.
///
/// \return A full-binary tree whose nodes are @c EvalExpr (or derived from it)
///         types.
///
template <typename ExprT,
          std::enable_if_t<IsEvalExpr<std::decay_t<ExprT>>, bool> = true>
EvalNode<ExprT> eval_node(ExprPtr const& expr) {
  using ranges::accumulate;
  using ranges::views::tail;
  using ranges::views::transform;

  if (expr->is<Tensor>()) return EvalNode<ExprT>{ExprT{expr->as<Tensor>()}};
  if (expr->is<Constant>()) return EvalNode<ExprT>{ExprT{expr->as<Constant>()}};
  if (expr->is<Variable>())
    return EvalNode<ExprT>{ExprT{expr->as<Variable>()}};
  assert(expr->is<Sum>() || expr->is<Product>());

  auto subxprs = *expr | ranges::views::transform([](auto const& x) {
    return eval_node<ExprT>(x);
  }) | ranges::to_vector;

  if (expr->is<Product>()) {
    auto const& prod = expr->as<Product>();
    if (prod.scalar() != 1)
      subxprs.emplace_back(eval_node<ExprT>(ex<Constant>(prod.scalar())));
  }

  auto const op = expr->is<Sum>() ? EvalOp::Sum : EvalOp::Prod;

  auto bnode = ranges::accumulate(
      ranges::views::tail(subxprs), std::move(*subxprs.begin()),
      [op](auto& lnode, auto& rnode) {
        auto pxpr = ExprT{*lnode, *rnode, op};
        return EvalNode<ExprT>(std::move(pxpr), std::move(lnode),
                               std::move(rnode));
      });

  return bnode;
}

///
/// \brief Creates an evaluation tree from Expr.
///
/// \tparam EvalNodeT is @c EvalNode<ExprT>, where @c ExprT is @c EvalExpr or
/// derived.
///
/// \param expr The expression to be converted to an evaluation tree.
///
/// \return A full-binary tree whose type is @c EvalNodeT.
///
template <typename EvalNodeT,
          std::enable_if_t<IsEvalNode<std::decay_t<EvalNodeT>>, bool> = true>
EvalNodeT eval_node(ExprPtr const& expr) {
  return eval_node<typename EvalNodeT::value_type>(expr);
}

template <typename ExprT>
ExprPtr to_expr(EvalNode<ExprT> const& node) {
  auto const op = node->op_type();
  auto const& evxpr = *node;

  if (node.leaf()) return evxpr.expr();

  if (op == EvalOp::Prod) {
    auto prod = Product{};

    ExprPtr lexpr = to_expr(node.left());
    ExprPtr rexpr = to_expr(node.right());

    prod.append(1, lexpr, Product::Flatten::No);
    prod.append(1, rexpr, Product::Flatten::No);

    assert(!prod.empty());

    if (prod.size() == 1 && !prod.factor(0)->is<Tensor>()) {
      return ex<Product>(Product{prod.scalar(), prod.factor(0)->begin(),
                                 prod.factor(0)->end(), Product::Flatten::No});
    } else {
      return ex<Product>(std::move(prod));
    }

  } else {
    assert(op == EvalOp::Sum && "unsupported operation type");
    return ex<Sum>(Sum{to_expr(node.left()), to_expr(node.right())});
  }
}

template <typename ExprT>
ExprPtr linearize_eval_node(EvalNode<ExprT> const& node) {
  if (node.leaf()) return to_expr(node);

  ExprPtr lres = linearize_eval_node(node.left());
  ExprPtr rres = linearize_eval_node(node.right());

  assert(lres);
  assert(rres);

  if (node->op_type() == EvalOp::Sum) return ex<Sum>(ExprPtrList{lres, rres});
  assert(node->op_type() == EvalOp::Prod);
  return ex<Product>(
      Product{1, ExprPtrList{lres, rres}, Product::Flatten::Yes});
}

namespace detail {
///
///  \tparam F function type that takes EvalNode const& argument and returns
///  bool.
///
/// \param node Node to compute asymptotic cost on.
///
/// \param exploit_symmetry Whether to use symmetry properties of an
/// intermediate to get reduced cost.
///
/// \param pred pred is called on every node and only those nodes that return
/// true will be used to compute cost.
///
/// \return Asymptotic cost of evaluation
/// in terms of number of occupied and virtual orbitals.
///
template <
    typename ExprT, typename F = std::function<bool(EvalNode<ExprT> const&)>,
    std::enable_if_t<std::is_invocable_r_v<bool, F, EvalNode<ExprT> const&>,
                     bool> = true>
AsyCost asy_cost_impl(EvalNode<ExprT> const& node, bool exploit_symmetry,
                      F const& pred) {
  if (node.leaf() || !pred(node)) return AsyCost::zero();

  return AsyCost{exploit_symmetry ? asy_cost_single_node(node)
                                  : asy_cost_single_node_symm_off(node)}  //
         +                                                                //
         asy_cost_impl(node.left(), exploit_symmetry, pred)               //
         +                                                                //
         asy_cost_impl(node.right(), exploit_symmetry, pred);
}
}  // namespace detail

template <typename ExprT>
AsyCost asy_cost_single_node_symm_off(EvalNode<ExprT> const& node) {
  if (node.leaf()) return AsyCost::zero();

  auto bks = ranges::views::concat(
      node.left()->expr()->template as<Tensor>().const_braket(),
      node.right()->expr()->template as<Tensor>().const_braket(),
      node->expr()->template as<Tensor>().const_braket());
  auto const uniques =
      bks | ranges::to<container::set<Index, Index::LabelCompare>>;

  size_t const nocc = ranges::count_if(uniques, [](auto&& idx) {
    return idx.space() == IndexSpace::active_occupied;
  });

  size_t const nvirt = uniques.size() - nocc;

  return AsyCost{node->op_type() == EvalOp::Prod ? 2 : 1, nocc, nvirt};
}

template <typename ExprT>
AsyCost asy_cost_single_node(EvalNode<ExprT> const& node) {
  auto cost = asy_cost_single_node_symm_off(node);
  auto const& ptensor = node->expr()->template as<Tensor>();
  // parent node symmetry
  auto const psym = ptensor.symmetry();
  // parent node bra symmetry
  auto const pbrank = ptensor.bra_rank();
  // parent node ket symmetry
  auto const pkrank = ptensor.ket_rank();

  if (psym == Symmetry::nonsymm || psym == Symmetry::invalid) {
    // do nothing
  } else {
    // ------
    // psym is Symmetry::symm or Symmetry::antisymm
    //
    // the rules of cost reduction are taken from
    //   doi:10.1016/j.procs.2012.04.044
    // ------

    auto const op = node->op_type();
    if (op == EvalOp::Sum) {
      cost = psym == Symmetry::symm
                 ? cost / (factorial(pbrank) * factorial(pkrank))
                 : cost / factorial(pbrank);
    } else if (op == EvalOp::Prod) {
      auto const lsym = node.left()->expr()->template as<Tensor>().symmetry();
      auto const rsym = node.right()->expr()->template as<Tensor>().symmetry();
      cost = (lsym == rsym && lsym == Symmetry::nonsymm)
                 ? cost / factorial(pbrank)
                 : cost / (factorial(pbrank) * factorial(pkrank));
    } else {
      assert(
          false &&
          "Unsupported evaluation operation for asymptotic cost computation.");
    }
  }

  return cost;
}

///
/// \param pred pred is called on every node and only those nodes that return
/// true will be used to compute cost. Default function: returns true.
///
template <
    typename ExprT, typename F = std::function<bool(EvalNode<ExprT> const&)>,
    std::enable_if_t<std::is_invocable_r_v<bool, F, EvalNode<ExprT> const&>,
                     bool> = true>
AsyCost asy_cost_symm_off(
    EvalNode<ExprT> const& node,
    F const& pred = [](EvalNode<ExprT> const&) { return true; }) {
  return detail::asy_cost_impl(node, false, pred);
}

///
/// \param pred pred is called on every node and only those nodes that return
/// true will be used to compute cost. Default function: returns true.
///
template <
    typename ExprT, typename F = std::function<bool(EvalNode<ExprT> const&)>,
    std::enable_if_t<std::is_invocable_r_v<bool, F, EvalNode<ExprT> const&>,
                     bool> = true>
AsyCost asy_cost(
    EvalNode<ExprT> const& node,
    F const& pred = [](EvalNode<ExprT> const&) { return true; }) {
  return detail::asy_cost_impl(node, true, pred);
}

}  // namespace sequant

#endif  // SEQUANT_EVAL_NODE_HPP
