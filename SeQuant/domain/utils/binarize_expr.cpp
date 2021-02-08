#include "binarize_expr.hpp"

namespace sequant::utils {

ExprPtr debinarize_eval_expr(binary_node<eval_expr> const& node) {
  auto op = node->op();
  auto const& evxpr = *node;

  switch (op) {
    case eval_expr::eval_op::Id:
      return evxpr.scalar().value() == 1.0
                 ? evxpr.tensor().clone()
                 : ex<Constant>(evxpr.scalar()) * evxpr.tensor().clone();

    case eval_expr::eval_op::Prod:
      return ex<Constant>(evxpr.scalar()) *
             ex<Product>(Product{debinarize_eval_expr(node.left()),
                                 debinarize_eval_expr(node.right())});

    case eval_expr::eval_op::Sum:
      return ex<Sum>(Sum{debinarize_eval_expr(node.left()),
                         debinarize_eval_expr(node.right())});
  }
}

}  // namespace sequant::utils
