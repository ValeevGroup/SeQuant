//
// Created by Eduard Valeyev on 8/4/23.
//

#ifndef SEQUANT_CORE_VISITORS_HPP
#define SEQUANT_CORE_VISITORS_HPP

#include "SeQuant/core/abstract_tensor.hpp"
#include "SeQuant/core/expr.hpp"

namespace sequant {

/// Removes tags from an Expr
inline ExprPtr& remove_tags(ExprPtr& expr) {
  auto reset_tags = [](ExprPtr& expr) {
    auto expr_cast_to_tensor =
        std::dynamic_pointer_cast<AbstractTensor>(expr.as_shared_ptr());
    if (expr_cast_to_tensor) {
      expr_cast_to_tensor->_reset_tags();
      return;
    }

    auto expr_cast_to_taggable =
        std::dynamic_pointer_cast<Taggable>(expr.as_shared_ptr());
    if (expr_cast_to_taggable) {
      expr_cast_to_taggable->reset();
      return;
    }
  };
  expr->visit(reset_tags);
  return expr;
}

}  // namespace sequant

#endif  // SEQUANT_CORE_VISITORS_HPP
