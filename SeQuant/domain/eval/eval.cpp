#include "eval.hpp"

#include <TiledArray/expressions/index_list.h>
#include <TiledArray/expressions/permopt.h>

namespace sequant::eval {

std::string const& EvalExprTA::annot() const { return annot_; }

EvalExprTA::EvalExprTA(Tensor const& tnsr)
    : EvalExpr(tnsr), annot_{braket_to_annot(tnsr.braket())} {}

EvalExprTA::EvalExprTA(const EvalExprTA& left, const EvalExprTA& right,
                       EvalOp op)
    : EvalExpr{left, right, op} {
  using Tidxs = TA::expressions::IndexList;
  using TA::expressions::BipartiteIndexList;
  using TA::expressions::GEMMPermutationOptimizer;

  if (result_type() == EvalResult::Tensor) {
    if (tot()) {
      annot_ = braket_to_annot(expr()->as<Tensor>().const_braket());
    } else {
      annot_ =
          GEMMPermutationOptimizer(Tidxs{left.annot()}, Tidxs{right.annot()})
              .target_result_indices()
              .string();
    }
  }
}

//
// EvalNodeTA to_eval_node_ta(EvalNode const& node) {
//  if (node.leaf()) {
//    auto result = EvalNodeTA{EvalExprTA{node->tensor()}};
//    result->scale(node->scalar());
//    return result;
//  }
//
//  auto left = to_eval_node_ta(node.left());
//  left->scale(node.left()->scalar());
//
//  auto right = to_eval_node_ta(node.right());
//  right->scale(node.right()->scalar());
//
//  auto curr = EvalExprTA{*left, *right, node->op()};
//  curr.scale(node->scalar());
//  return EvalNodeTA{std::move(curr), std::move(left), std::move(right)};
//}
//
// EvalNodeTA to_eval_node_ta(ExprPtr const& expr) {
//  return to_eval_node_ta(to_eval_node(expr));
//}

}  // namespace sequant::eval
