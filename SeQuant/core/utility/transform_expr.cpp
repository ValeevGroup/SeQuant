#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/transform_expr.hpp>

namespace sequant {

ExprPtr transform_expr(const ExprPtr& expr,
                       const container::map<Index, Index>& index_replacements,
                       Constant::scalar_type scaling_factor) {
  if (expr->is<Constant>() || expr->is<Variable>()) {
    return ex<Constant>(scaling_factor) * expr;
  }

  auto transform_tensor = [&index_replacements](const Tensor& tensor) {
    auto result = std::make_shared<Tensor>(tensor);
    result->transform_indices(index_replacements);
    result->reset_tags();
    return result;
  };

  auto transform_product = [&transform_tensor,
                            &scaling_factor](const Product& product) {
    auto result = std::make_shared<Product>();
    result->scale(product.scalar());
    for (auto&& term : product) {
      if (term->is<Tensor>()) {
        auto tensor = term->as<Tensor>();
        result->append(1, transform_tensor(tensor));
      } else if (term->is<Variable>() || term->is<Constant>()) {
        result->append(1, term->clone());
      } else {
        throw std::runtime_error("Invalid Expr type in transform_product");
      }
    }
    result->scale(scaling_factor);
    return result;
  };

  if (expr->is<Tensor>()) {
    auto result =
        ex<Constant>(scaling_factor) * transform_tensor(expr->as<Tensor>());
    return result;
  } else if (expr->is<Product>()) {
    auto result = transform_product(expr->as<Product>());
    return result;
  } else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto& term : *expr) {
      result->append(transform_expr(term, index_replacements, scaling_factor));
    }
    return result;
  } else {
    throw std::runtime_error("Invalid Expr type in transform_expr");
  }
}

}  // namespace sequant
