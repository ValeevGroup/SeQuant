#ifndef SEQUANT_UTILS_BINARIZE_EXPR_HPP
#define SEQUANT_UTILS_BINARIZE_EXPR_HPP

#include "binary_expr.hpp"
#include "eval_expr.hpp"
#include "eval_sequence.hpp"

#include <iostream>

namespace sequant::utils {

const struct {
 public:
  eval_expr operator()(const eval_expr& x) const { return x; }

  eval_expr operator()(const eval_expr& x, const eval_expr& y) const {
    return eval_expr{x, y};
  }
} binarize_eval_expr;

struct binarize_flat_prod {
  binarize_flat_prod(const Product& p) : prod{p} {
    for (const auto& f : p) assert(f->is<Tensor>());
  }

  binary_expr<eval_expr>::node_ptr operator()(
      const eval_sequence<size_t>& seq) {
    auto xpr_seq =
        transform_eval_sequence<size_t, eval_expr>(seq, [this](auto x) {
          return eval_expr{prod.at(x)->template as<Tensor>()};
        });

    auto&& pred = binarize_eval_expr;

    auto without_scalar =
        binarize_eval_sequence<eval_expr, eval_expr>(xpr_seq, pred);

    if (prod.scalar() == 1.0) return std::move(without_scalar);

    auto scal = Constant{prod.scalar()};

    return make_binary_expr(pred(eval_expr{scal}, without_scalar->data()),
                            make_binary_expr(eval_expr{scal}),
                            std::move(without_scalar));
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
binary_expr<eval_expr>::node_ptr binarize_evxpr_range(Cont&& container) {
  const auto eseq = eval_sequence<eval_expr>{
      (*ranges::begin(container)),
      ranges::views::tail(container) |
          ranges::views::transform(
              [](const auto& x) { return eval_sequence<eval_expr>{x}; }) |
          ranges::to<std::vector<eval_sequence<eval_expr>>>};

  return binarize_eval_sequence<eval_expr, eval_expr>(eseq, binarize_eval_expr);
}

}  // namespace sequant::utils

#endif  // SEQUANT_UTILS_BINARIZE_EXPR_HPP
