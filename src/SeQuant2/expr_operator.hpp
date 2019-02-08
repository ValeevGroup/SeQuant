//
// Created by Eduard Valeyev on 4/2/18.
//

#ifndef SEQUANT2_EXPR_OPERATOR_HPP
#define SEQUANT2_EXPR_OPERATOR_HPP

namespace sequant2 {

inline ExprPtr
operator*(const ExprPtr &left, const ExprPtr &right) {
  auto left_is_product = left->is<Product>();
  auto right_is_product = right->is<Product>();
  if (!left_is_product && !right_is_product) {
    return ex<Product>(ExprPtrList{left, right});
  } else if (left_is_product) {
    auto left_product = std::static_pointer_cast<Product>(left);
    auto result = std::make_shared<Product>(*left_product);
    result->append(1, right);
    return result;
  } else {  // right_is_product
    auto right_product = std::static_pointer_cast<Product>(right);
    auto result = std::make_shared<Product>(*right_product);
    result->prepend(1, left);
    return result;
  }
  abort();  // unreachable
}

inline ExprPtr
operator+(const ExprPtr &left, const ExprPtr &right) {
  auto left_is_sum = left->is<Sum>();
  auto right_is_sum = right->is<Sum>();
  if (!left_is_sum && !right_is_sum) {
    return ex<Sum>(ExprPtrList{left, right});
  } else if (left_is_sum) {
    auto left_sum = std::static_pointer_cast<Sum>(left);
    auto result = std::make_shared<Sum>(*left_sum);
    result->append(right);
    return result;
  } else {  // right_is_sum
    auto right_sum = std::static_pointer_cast<Sum>(right);
    auto result = std::make_shared<Sum>(*right_sum);
    result->prepend(left);
    return result;
  }
  abort();  // unreachable
}

inline ExprPtr
operator-(const ExprPtr &left, const ExprPtr &right) {
  auto left_is_sum = left->is<Sum>();
  if (!left_is_sum) {
    return ex<Sum>(ExprPtrList{left, ex<Product>(-1.0, ExprPtrList{right})});
  } else if (left_is_sum) {
    auto left_sum = std::static_pointer_cast<Sum>(left);
    auto result = std::make_shared<Sum>(*left_sum);
    result->append(ex<Product>(-1.0, ExprPtrList{right}));
    return result;
  }
  abort();  // unreachable
}

}  // namespace sequant2

#endif //SEQUANT2_EXPR_OPERATOR_HPP
