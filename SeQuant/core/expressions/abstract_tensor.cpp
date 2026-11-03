//
// Created by Eduard Valeyev on 9/25/25.
//

#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/expressions/variable.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/algorithm/any_of.hpp>

#include <typeinfo>

namespace sequant {

void reset_tags(const ExprPtr& expr) {
  expr->visit(
      [](const ExprPtr& current) {
        if (current.is<AbstractTensor>()) {
          reset_tags(current->as<AbstractTensor>());
        }
      },
      /* atoms_only = */ true);
}

bool cross_type_less(const AbstractTensor& a, const AbstractTensor& b) {
  if (a._label() != b._label()) return a._label() < b._label();
  if (a._bra_rank() != b._bra_rank()) return a._bra_rank() < b._bra_rank();
  if (a._ket_rank() != b._ket_rank()) return a._ket_rank() < b._ket_rank();
  if (a._aux_rank() != b._aux_rank()) return a._aux_rank() < b._aux_rank();
  return expr_type_id_less(a, b);
}

bool expr_type_id_less(const AbstractTensor& a, const AbstractTensor& b) {
  if (typeid(a) == typeid(b)) return false;
  const auto* lhs_expr = dynamic_cast<const Expr*>(&a);
  const auto* rhs_expr = dynamic_cast<const Expr*>(&b);
  SEQUANT_ASSERT(lhs_expr && rhs_expr);
  return lhs_expr->type_id() < rhs_expr->type_id();
}

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
    throw Exception("Invalid Expr type in remove_tensor");
}

}  // namespace sequant
