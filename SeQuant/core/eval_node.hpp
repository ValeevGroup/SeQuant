//
// Created by Bimal Gaudel on 5/24/21.
//

#ifndef SEQUANT_EVAL_NODE_HPP
#define SEQUANT_EVAL_NODE_HPP

#include "asy_cost.hpp"
#include "binary_node.hpp"
#include "eval_expr.hpp"

#include <boost/math/special_functions/factorials.hpp>

namespace sequant {

template <typename T = EvalExpr,
          typename = std::enable_if_t<std::is_convertible_v<T, EvalExpr>>>
using EvalNode = FullBinaryNode<T>;

template <typename ExprT>
EvalNode<ExprT> to_eval_node(ExprPtr const& expr) {
  using ranges::accumulate;
  using ranges::views::tail;
  using ranges::views::transform;

  if (expr->is<Tensor>()) return EvalNode<ExprT>{ExprT{expr->as<Tensor>()}};
  if (expr->is<Constant>()) return EvalNode<ExprT>{ExprT{expr->as<Constant>()}};
  assert(expr->is<Sum>() || expr->is<Product>());

  auto subxprs = *expr | ranges::views::transform([](auto const& x) {
    return to_eval_node<ExprT>(x);
  }) | ranges::to_vector;

  if (expr->is<Product>()) {
    auto const& prod = expr->as<Product>();
    if (prod.scalar() != 1)
      subxprs.emplace_back(to_eval_node<ExprT>(ex<Constant>(prod.scalar())));
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

template <typename ExprT>
ExprPtr to_expr(EvalNode<ExprT> const& node) {
  auto const op = node->op_type();
  auto const& evxpr = *node;

  if (node.leaf()) return evxpr.expr()->clone();

  if (op == EvalOp::Prod) {
    ExprPtr lexpr = to_expr(node.left());
    ExprPtr rexpr = to_expr(node.right());

    if (lexpr->is<Constant>() && rexpr->is<Product>()) {
      auto& p = rexpr->as<Product>();
      p.scale(p.scalar() * lexpr->as<Constant>().value());
      return rexpr;
    }

    if (lexpr->is<Product>() && rexpr->is<Constant>()) {
      auto& p = lexpr->as<Product>();
      p.scale(p.scalar() * rexpr->as<Constant>().value());
      return lexpr;
    }

    auto prod = Product{};

    if (lexpr->is<Tensor>())
      prod.append(1, lexpr);
    else
      prod.append(lexpr);

    if (rexpr->is<Tensor>())
      prod.append(1, rexpr);
    else
      prod.append(rexpr);

    return ex<Product>(std::move(prod));
  } else {
    assert(op == EvalOp::Sum && "unsupported operation type");
    return ex<Sum>(Sum{to_expr(node.left()), to_expr(node.right())});
  }
}

template <typename ExprT>
ExprPtr linearize_eval_node(EvalNode<ExprT> const& node) {
  auto extend_product = [](Product& prod, ExprPtr const& factors) {
    for (auto&& f : *factors)
      if (f->is<Tensor>())
        prod.append(1, f);
      else
        prod.append(f);
  };

  if (node.leaf()) return to_expr(node);

  ExprPtr lres = linearize_eval_node(node.left());
  ExprPtr rres = linearize_eval_node(node.right());

  assert(lres);
  assert(rres);

  if (node->op_type() == EvalOp::Sum) return ex<Sum>(ExprPtrList{lres, rres});

  if (lres->is<Constant>() && rres->is<Product>()) {
    auto& prod = rres->as<Product>();
    prod.scale(prod.scalar() * lres->as<Constant>().value());
    return rres;
  }

  if (lres->is<Product>() && rres->is<Constant>()) {
    auto& prod = lres->as<Product>();
    prod.scale(prod.scalar() * rres->as<Constant>().value());
    return lres;
  }

  if (lres->is<Product>() && rres->is<Product>()) {
    auto& lprod = lres->as<Product>();
    auto& rprod = rres->as<Product>();

    auto scal = lprod.scalar() * rprod.scalar();

    auto result = Product{};
    extend_product(result, lres);
    extend_product(result, rres);
    return result.clone();
  }

  if (lres->is<Tensor>() && rres->is<Product>()) {
    auto result = Product{};
    result.scale(rres->as<Product>().scalar());
    result.append(1, lres);
    extend_product(result, rres);
    return result.clone();
  }

  if (lres->is<Product>() && rres->is<Tensor>()) {
    auto& prod = lres->as<Product>();
    prod.append(1, rres);
    return lres;
  }

  auto prod = Product{};
  prod.append(lres);
  prod.append(rres);
  return prod.clone();
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
                      F&& pred) {
  if (node.leaf() || !std::invoke(std::forward<F>(pred), node))
    return AsyCost::zero();

  return AsyCost{exploit_symmetry ? asy_cost_single_node(node)
                                  : asy_cost_single_node_symm_off(node)} +  //
         asy_cost_impl(node.left(), exploit_symmetry,
                       std::forward<F>(pred)) +                             //
         asy_cost_impl(node.right(), exploit_symmetry, std::forward<F>(pred));
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
  auto factorial = [](auto x) {
    return static_cast<int>(boost::math::factorial<double>(x));
  };
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
    F&& pred = [](EvalNode<ExprT> const&) { return true; }) {
  return detail::asy_cost_impl(node, false, std::forward<F>(pred));
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
    F&& pred = [](EvalNode<ExprT> const&) { return true; }) {
  return detail::asy_cost_impl(node, true, std::forward<F>(pred));
}

}  // namespace sequant

#endif  // SEQUANT_EVAL_NODE_HPP
