#include "eval.hpp"

#include <TiledArray/expressions/index_list.h>
#include <TiledArray/expressions/permopt.h>

namespace sequant::eval {

std::string const& EvalExprTA::annot() const { return annot_; }

EvalExprTA::EvalExprTA(Tensor const& tnsr)
    : EvalExpr(tnsr), annot_{braket_to_annot(tnsr.braket())} {}

EvalExprTA::EvalExprTA(Constant const& c) : EvalExpr(c), annot_{} {}

EvalExprTA::EvalExprTA(const EvalExprTA& left, const EvalExprTA& right,
                       EvalOp op)
    : EvalExpr{left, right, op} {
  using Tidxs = TA::expressions::IndexList;
  using TA::expressions::BipartiteIndexList;
  using TA::expressions::GEMMPermutationOptimizer;

  if (result_type() == ResultType::Tensor) {
    annot_ = braket_to_annot(as_tensor().const_braket());
    if (left.result_type() == right.result_type() &&
        op_type() == EvalOp::Prod && !tot()) {
      // tensor x tensor confirmed
      auto lh = hash::value(left);
      auto rh = hash::value(right);
      auto const& left_ = lh <= rh ? left : right;
      auto const& right_ = lh <= rh ? right : left;
      annot_ = GEMMPermutationOptimizer(Tidxs{left_.annot()},  //
                                        Tidxs{right_.annot()})
                   .target_result_indices()
                   .string();
    }
  }
}

EvalExprBTAS::annot_t const& EvalExprBTAS::annot() const noexcept {
  return annot_;
}

EvalExprBTAS::EvalExprBTAS(Tensor const& t) noexcept
    : EvalExpr{t}, annot_{index_hash(t.const_braket()) | ranges::to<annot_t>} {}

EvalExprBTAS::EvalExprBTAS(Constant const& c) noexcept : EvalExpr{c} {}

EvalExprBTAS::EvalExprBTAS(EvalExprBTAS const& left,   //
                           EvalExprBTAS const& right,  //
                           EvalOp op) noexcept
    : EvalExpr{left, right, op} {
  if (result_type() == ResultType::Tensor) {
    assert(!tot() && "Tensor of tensor not supported in BTAS");
    annot_ = index_hash(as_tensor().const_braket()) | ranges::to<annot_t>;
  }
}

}  // namespace sequant::eval
