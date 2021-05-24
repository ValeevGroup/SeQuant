#include "binarize_expr.hpp"

namespace sequant::utils {

binary_node<eval_expr> binarize_expr(const ExprPtr& expr) {
  if (expr->is<Tensor>())
    return binary_node<eval_expr>{eval_expr{expr->as<Tensor>()}};

  auto subxprs = *expr | ranges::views::transform([](auto const& x) {
    return binarize_expr(x);
  }) | ranges::to_vector;

  auto bnode = ranges::accumulate(
      ranges::views::tail(subxprs), std::move(*subxprs.begin()),
      [](auto& lnode, auto& rnode) {
        auto pxpr = eval_expr{*lnode, *rnode};
        if (pxpr.op() == eval_expr::eval_op::Prod) {
          pxpr *= lnode->scalar();
          pxpr *= rnode->scalar();

          lnode->scale(1.0);
          rnode->scale(1.0);
        }

        return binary_node<eval_expr>(std::move(pxpr), std::move(lnode),
                                      std::move(rnode));
      });

  if (expr->is<Product>()) *bnode *= expr->as<Product>().scalar();

  return bnode;
}

ExprPtr binarize_expr(binary_node<eval_expr> const& node) {
  auto const op = node->op();
  auto const& evxpr = *node;

  if (node.leaf()) {
    return evxpr.scalar() == Constant{1}
               ? evxpr.tensor().clone()
               : ex<Constant>(evxpr.scalar()) * evxpr.tensor().clone();
  }

  if (op == eval_expr::eval_op::Prod) {
    auto prod = Product{};
    prod.scale(evxpr.scalar().value());

    auto lexpr = binarize_expr(node.left());
    auto rexpr = binarize_expr(node.right());

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
    assert(op == eval_expr::eval_op::Sum && "unsupported operation type");
    return ex<Sum>(
        Sum{binarize_expr(node.left()), binarize_expr(node.right())});
  }
}

ExprPtr linearize_expr(binary_node<eval_expr> const& node) {
  if (node.leaf()) return binarize_expr(node);

  auto lres = linearize_expr(node.left());
  auto rres = linearize_expr(node.right());
  if (node->op() == eval_expr::eval_op::Sum)
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

ExprPtr linearize_expr(ExprPtr const& expr) {
  return linearize_expr(binarize_expr(expr));
}

}  // namespace sequant::utils
