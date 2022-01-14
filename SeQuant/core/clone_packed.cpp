//
// Created by Bimal Gaudel on 9/22/21.
//

#include "clone_packed.hpp"
#include "expr.hpp"
#include "tensor.hpp"

namespace sequant {

using ranges::views::transform;

ExprPtr clone_packed(Tensor const& t) {
  return t.clone();
}

ExprPtr clone_packed(Sum const& s) {
  auto const smands = s
                      | transform([](ExprPtr x){return clone_packed(x);})
                      | ranges::to_vector;
  return ex<Sum>(smands.begin(), smands.end());
}

ExprPtr clone_packed(Product const& p) {
  auto const facs = p
                    | transform([](ExprPtr x){return clone_packed(x);});
  auto const scal = p.scalar();

  auto result = ex<Product>(scal, ExprPtrList{});
  for (auto&& f: facs)
    result->as<Product>().append(f);
  return result;
}

ExprPtr clone_packed(ExprPtr expr) {

  if (!expr) return nullptr;
  else if (expr->is<Sum>()){
    return clone_packed(expr->as<Sum>());
  }
  else if (expr->is<Product>()){
    return clone_packed(expr->as<Product>());
  }
  else return expr->clone();
}

} //
