#include "eval_node.hpp"
#include "expr.hpp"

namespace sequant {

EvalNode to_eval_node(ExprPtr const& expr) {
  if (expr->is<Tensor>()) return EvalNode{EvalExpr{expr->as<Tensor>()}};

  auto subxprs = *expr | ranges::views::transform([](auto const& x) {
    return to_eval_node(x);
  }) | ranges::to_vector;

  auto bnode = ranges::accumulate(
      ranges::views::tail(subxprs), std::move(*subxprs.begin()),
      [](auto& lnode, auto& rnode) {
        auto pxpr = EvalExpr{*lnode, *rnode};
        if (pxpr.op() == EvalExpr::EvalOp::Prod) {
          pxpr *= lnode->scalar();
          pxpr *= rnode->scalar();

          lnode->scale(1.0);
          rnode->scale(1.0);
        }

        return EvalNode(std::move(pxpr), std::move(lnode), std::move(rnode));
      });

  if (expr->is<Product>()) *bnode *= expr->as<Product>().scalar();

  return bnode;
}

ExprPtr to_expr(EvalNode const& node) {
  auto const op = node->op();
  auto const& evxpr = *node;

  if (node.leaf()) {
    return evxpr.scalar() == Constant{1}
               ? evxpr.tensor().clone()
               : ex<Constant>(evxpr.scalar()) * evxpr.tensor().clone();
  }

  if (op == EvalExpr::EvalOp::Prod) {
    auto prod = Product{};
    prod.scale(evxpr.scalar().value());

    auto lexpr = to_expr(node.left());
    auto rexpr = to_expr(node.right());

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
    assert(op == EvalExpr::EvalOp::Sum && "unsupported operation type");
    return ex<Sum>(Sum{to_expr(node.left()), to_expr(node.right())});
  }
}

ExprPtr linearize_eval_node(EvalNode const& node) {
  if (node.leaf()) return to_expr(node);

  auto lres = to_expr(node.left());
  auto rres = to_expr(node.right());
  if (node->op() == EvalExpr::EvalOp::Sum)
    return ex<Sum>(ExprPtrList{lres, rres});

  if (node.left().leaf() && node.right().leaf()) {
    return ex<Product>(node->scalar().value(), ExprPtrList{lres, rres});
  } else if (!node.left().leaf() && !node.right().leaf()) {
    auto prod = Product(node->scalar().value(), ExprPtrList{});
    prod.append(lres);
    prod.append(rres);
    return ex<Product>(std::move(prod));
  } else if (node.left().leaf() && !node.right().leaf()) {
    auto prod = Product(node->scalar().value(), ExprPtrList{lres});
    prod.append(rres);
    return ex<Product>(std::move(prod));
  } else {  // (!node.left().leaf() && node.right().leaf())
    auto& res = lres->as<Product>();
    res.scale(node->scalar().value());
    res.append(rres);
    return lres;
  }
}

AsyCost asy_cost_single_node(EvalNode const& node) {
  if (node.leaf()) return AsyCost::zero();
  if (node->op() == EvalExpr::EvalOp::Sum) {
    size_t nocc =
        ranges::count_if(node->tensor().const_braket(), [](Index const& idx) {
          return idx.space() == IndexSpace::active_occupied;
        });

    return AsyCost{nocc, node->tensor().const_braket().size() - nocc};
  }

  auto bks = ranges::views::concat(node.left()->tensor().const_braket(),
                                   node.right()->tensor().const_braket(),
                                   node->tensor().const_braket());
  auto uniques = bks | ranges::to<container::set<Index, Index::LabelCompare>>;

  size_t nocc = ranges::count_if(uniques, [](auto&& idx) {
    return idx.space() == IndexSpace::active_occupied;
  });
  size_t nvirt = uniques.size() - nocc;

  return AsyCost{nocc, nvirt};
}

}  // namespace sequant
