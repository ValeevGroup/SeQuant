#include "eval.hpp"

#include <TiledArray/expressions/index_list.h>
#include <TiledArray/expressions/permopt.h>

namespace sequant {

namespace {

template <
    typename Iterable,
    std::enable_if_t<!std::is_same_v<InnerOuterIndices, std::decay_t<Iterable>>,
                     bool> = true>
std::string indices_to_annot(Iterable const& indices) noexcept {
  using ranges::views::intersperse;
  using ranges::views::join;
  using ranges::views::transform;

  auto idx_label = [](Index const& idx) { return to_string(idx.label()); };

  return indices | transform(idx_label) | intersperse(",") | join |
         ranges::to<std::string>;
}

std::string indices_to_annot(InnerOuterIndices const& inout) noexcept {
  auto const& in = inout.inner;
  auto const& out = inout.outer;
  if (out.empty()) {
    return indices_to_annot(in);
  } else {
    return indices_to_annot(in) + ";" + indices_to_annot(out);
  }
}

}  // namespace

std::string const& EvalExprTA::annot() const { return annot_; }

EvalExprTA::EvalExprTA(Tensor const& tnsr)
    : EvalExpr(tnsr), annot_{indices_to_annot(inner_outer_indices())} {}

EvalExprTA::EvalExprTA(Constant const& c) : EvalExpr(c), annot_{} {}

EvalExprTA::EvalExprTA(const EvalExprTA& left, const EvalExprTA& right,
                       EvalOp op)
    : EvalExpr{left, right, op} {
  using Tidxs = TA::expressions::IndexList;
  using TA::expressions::BipartiteIndexList;
  using TA::expressions::GEMMPermutationOptimizer;

  if (result_type() == ResultType::Tensor) {
    annot_ = indices_to_annot(inner_outer_indices());
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

}  // namespace sequant
