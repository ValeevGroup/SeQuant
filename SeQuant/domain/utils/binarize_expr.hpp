#ifndef SEQUANT_UTILS_BINARIZE_EXPR_HPP
#define SEQUANT_UTILS_BINARIZE_EXPR_HPP

#include "binary_node.hpp"
#include "eval_expr.hpp"
#include "eval_seq.hpp"

namespace sequant::utils {

const struct {
 public:
  eval_expr operator()(const eval_expr& x) const { return x; }

  eval_expr operator()(const eval_expr& x, const eval_expr& y) const {
    auto result = eval_expr{x, y};

    if (result.op() == eval_expr::eval_op::Prod) {
      result *= x.scalar();
      result *= y.scalar();

      const_cast<eval_expr&>(x).scale(1.0);
      const_cast<eval_expr&>(y).scale(1.0);
    }
    return result;
  }
} binarize_eval_expr;

struct binarize_flat_prod {
  explicit binarize_flat_prod(const Product& p) : prod{p} {
    for (const auto& f : p) assert(f->is<Tensor>());
  }

  binary_node<eval_expr> operator()(const eval_seq<size_t>& seq) {
    auto const xpr_seq = seq.transform([this](auto const& x) {
      return eval_expr{prod.at(x)->template as<Tensor>()};
    });

    auto result = xpr_seq.binarize(binarize_eval_expr);
    *result *= prod.scalar();

    return result;
  }

  binary_node<eval_expr> operator()() {
    using ranges::views::iota;
    using ranges::views::take;
    using ranges::views::transform;
    auto seq = eval_seq<size_t>{
        0,                                                           //
        iota((size_t)1)                                              //
            | take(prod.size() - 1)                                  //
            | transform([](auto x) { return eval_seq<size_t>{x}; })  //
            | ranges::to_vector};
    return operator()(seq);
  }

 private:
  const Product& prod;
};

/**
 * @tparam Cont type of a container.
 *
 * @param container Cont type container of eval_expr objects.
 */
template <typename Cont>
binary_node<eval_expr> binarize_evxpr_range(Cont&& container) {
  using ranges::begin;
  using ranges::to;
  using ranges::views::tail;
  using ranges::views::transform;

  const auto eseq = eval_seq<eval_expr>{
      *begin(container),                                                     //
      tail(container)                                                        //
          | transform([](const auto& x) { return eval_seq<eval_expr>{x}; })  //
          | to<std::vector<eval_seq<eval_expr>>>};

  return eseq.binarize(binarize_eval_expr);
}

/**
 * @tparam Cont type of a container.
 *
 * @param scaler the final result will be scaled by this.
 * @param container Cont type container of eval_expr objects.
 */
template <typename Cont, typename Scalar = std::complex<double>>
binary_node<eval_expr> binarize_evxpr_range(Scalar scal, Cont&& container) {
  auto result = binarize_evxpr_range(std::forward<Cont>(container));

  *result *= scal;
  return result;
}

ExprPtr debinarize_eval_expr(binary_node<eval_expr> const& node);

}  // namespace sequant::utils

#endif  // SEQUANT_UTILS_BINARIZE_EXPR_HPP
