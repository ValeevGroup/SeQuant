//
// Created by Eduard Valeyev on 9/25/25.
//

#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/expressions/variable.hpp>

namespace sequant {

bool has_tensor(const ExprPtr& expr, std::wstring label) {
  if (expr->is<Constant>() || expr->is<Variable>()) return false;

  auto check_product = [&label](const Product& p) {
    return ranges::any_of(p.factors(), [&label](const auto& t) {
      return t->template is<AbstractTensor>() &&
             (t->template as<AbstractTensor>())._label() == label;
    });
  };

  if (expr->is<AbstractTensor>()) {
    return expr->as<AbstractTensor>()._label() == label;
  } else if (expr->is<Product>()) {
    return check_product(expr->as<Product>());
  } else if (expr->is<Sum>()) {
    return ranges::any_of(
        *expr, [&label](const auto& term) { return has_tensor(term, label); });
  } else
    return false;
}

ExprPtr remove_tensor(const ExprPtr& expr, std::wstring label) {
  if (expr->is<Sum>()) {
    Sum result{};
    for (auto& term : *expr) {
      result.append(remove_tensor(term, label));
    }
    return ex<Sum>(result);
  } else if (expr->is<Product>()) {
    return [](const Product& product, std::wstring label) {
      // filter out tensors with specified label
      auto new_product = std::make_shared<Product>();
      new_product->scale(product.scalar());
      for (auto&& term : product) {
        if (term->is<AbstractTensor>()) {
          if (term->as<AbstractTensor>()._label() != label)
            new_product->append(1, term.clone());
        } else
          new_product->append(1, term);
      }
      return new_product;
    }(expr->as<Product>(), label);
  } else if (expr->is<AbstractTensor>())
    return expr->as<AbstractTensor>()._label() == label ? ex<Constant>(1)
                                                        : expr;
  else if (expr->is<Constant>() || expr->is<Variable>())
    return expr;
  else
    throw std::runtime_error("Invalid Expr type in remove_tensor");
}

}  // namespace sequant
