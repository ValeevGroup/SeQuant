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

}  // namespace sequant::eval
