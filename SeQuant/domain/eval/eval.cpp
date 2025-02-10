#include <SeQuant/domain/eval/eval.hpp>

#include <TiledArray/expressions/index_list.h>
#include <TiledArray/expressions/permopt.h>

namespace sequant {

std::string const& EvalExprTA::annot() const { return annot_; }

EvalExprTA::EvalExprTA(Tensor const& tnsr) : EvalExpr(tnsr) {
  annot_ = indices_annot();
}

EvalExprTA::EvalExprTA(Constant const& c) : EvalExpr(c), annot_{} {}

EvalExprTA::EvalExprTA(Variable const& v) : EvalExpr(v), annot_{} {}

EvalExprTA::EvalExprTA(const EvalExprTA& left, const EvalExprTA& right,
                       EvalOp op)
    : EvalExpr{left, right, op} {
  using TA::expressions::BipartiteIndexList;
  using TA::expressions::GEMMPermutationOptimizer;

  if (result_type() == ResultType::Tensor) {
    annot_ = indices_annot();
    // clang-format off
// TODO: fix the following so that it works for ToT x ToT -> T
//    using Tidxs = TA::expressions::IndexList;
//    if (left.result_type() == right.result_type() &&
//        op_type() == EvalOp::Prod && !tot()) {
//      // tensor x tensor confirmed
//      auto lh = hash::value(left);
//      auto rh = hash::value(right);
//      auto const& left_ = lh <= rh ? left : right;
//      auto const& right_ = lh <= rh ? right : left;
//      annot_ = GEMMPermutationOptimizer(Tidxs{left_.annot()},  //
//                                        Tidxs{right_.annot()})
//                   .target_result_indices()
//                   .string();
//    }
    // clang-format on
  }
}

EvalExprBTAS::annot_t const& EvalExprBTAS::annot() const noexcept {
  return annot_;
}

EvalExprBTAS::EvalExprBTAS(Tensor const& t) noexcept
    : EvalExpr{t},
      annot_{index_hash(t.const_indices()) | ranges::to<annot_t>} {}

EvalExprBTAS::EvalExprBTAS(Constant const& c) noexcept : EvalExpr{c} {}

EvalExprBTAS::EvalExprBTAS(Variable const& v) noexcept : EvalExpr{v} {}

EvalExprBTAS::EvalExprBTAS(EvalExprBTAS const& left,   //
                           EvalExprBTAS const& right,  //
                           EvalOp op) noexcept
    : EvalExpr{left, right, op} {
  if (result_type() == ResultType::Tensor) {
    assert(!tot() && "Tensor of tensor not supported in BTAS");
    annot_ = index_hash(as_tensor().const_indices()) | ranges::to<annot_t>;
  }
}

}  // namespace sequant
