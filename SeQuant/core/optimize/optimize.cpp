#include "SeQuant/core/optimize.hpp"

namespace sequant {

namespace opt {

ExprPtr tail_factor(ExprPtr const& expr) noexcept {
  if (expr->is<Tensor>())
    return expr->clone();

  else if (expr->is<Product>()) {
    auto scalar = expr->as<Product>().scalar();
    auto facs = ranges::views::tail(*expr);
    return ex<Product>(Product{scalar, ranges::begin(facs), ranges::end(facs)});
  } else {
    // sum
    auto summands = *expr | ranges::views::transform(
                                [](auto const& x) { return tail_factor(x); });
    return ex<Sum>(Sum{ranges::begin(summands), ranges::end(summands)});
  }
}

void pull_scalar(sequant::ExprPtr expr) noexcept {
  using sequant::Product;
  if (!expr->is<Product>()) return;
  auto& prod = expr->as<Product>();

  auto scal = prod.scalar();
  for (auto&& x : *expr)
    if (x->is<Product>()) {
      auto& p = x->as<Product>();
      scal *= p.scalar();
      p.scale(1.0 / p.scalar());
    }

  prod.scale(1.0 / prod.scalar());
  prod.scale(scal);
}

}  // namespace opt

EvalNode optimize(const ExprPtr& expr) {
  using ranges::views::transform;
  if (expr->is<Tensor>())
    return to_eval_node(expr);
  else if (expr->is<Product>()) {
    return *(opt::single_term_opt(expr->as<Product>()).optimal_seqs.begin());
  } else if (expr->is<Sum>()) {
    auto smands = *expr | transform([](auto const& s) {
      return to_expr(optimize(s));
    }) | ranges::to_vector;

    return to_eval_node(ex<Sum>(Sum{smands.begin(), smands.end()}));
  } else
    throw std::runtime_error{"optimization attempted on unsupported Expr type"};
}

}  // namespace sequant
