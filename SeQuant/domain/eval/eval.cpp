#include "eval.hpp"

#include <TiledArray/expressions/index_list.h>
#include <TiledArray/expressions/permopt.h>

namespace sequant::eval {

namespace {

///
/// Given an iterable of Index objects, generate a string annotation
/// that can be used for TiledArray tensor expressions.
/// Tensor-of-tensors also supported.
template <typename Indices>
std::string braket_to_annot(Indices const& indices) {
  using ranges::find;
  using ranges::views::filter;
  using ranges::views::intersperse;
  using ranges::views::join;

  // make a comma-separated string out of an iterable of strings
  auto add_commas = [](auto const& strs) -> std::string {
    return strs | intersperse(",") | join | ranges::to<std::string>;
  };

  container::svector<std::string> idxs{}, pidxs{};
  for (auto&& idx : indices) {
    idxs.emplace_back(sequant::to_string(idx.label()));
    for (auto&& pidx : idx.proto_indices())
      pidxs.emplace_back(sequant::to_string(pidx.label()));
  }

  if (pidxs.empty()) {
    // not a tensor-of-tensor type expression
    return add_commas(idxs);
  } else {
    ranges::stable_sort(pidxs);
    ranges::actions::unique(pidxs);
    auto not_in_pidx = [&pidxs](auto&& l) {
      return find(pidxs, l) == pidxs.end();
    };
    return add_commas(pidxs) + ";" +
           add_commas(idxs | filter(not_in_pidx) | ranges::to<decltype(idxs)>);
  }
}
}  // namespace

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
