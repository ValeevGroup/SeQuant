//
// Created by Bimal Gaudel on 9/22/21.
//

#include "clone.hpp"
#include "expr.hpp"

sequant::ExprPtr sequant::clone(sequant::ExprPtr expr) {
  using ranges::views::transform;

  if (!expr) return nullptr;
  else if (expr->is<Sum>()){
    auto const smands = *expr
                  | transform([](ExprPtr x){return clone(x);})
                  | ranges::to_vector;
    return ex<Sum>(smands.begin(), smands.end());
  }
  else if (expr->is<Product>()){
    auto const facs = *expr
                | transform([](ExprPtr x){return clone(x);});
    auto const scal = expr->as<Product>().scalar();

    auto result = ex<Product>(scal, ExprPtrList{});
    for (auto&& f: facs)
      result->as<Product>().append(f);
    return result;
  }
  else return expr->clone();
}
