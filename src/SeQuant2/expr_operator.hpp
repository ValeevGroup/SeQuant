//
// Created by Eduard Valeyev on 4/2/18.
//

#ifndef SEQUANT2_EXPR_OPERATOR_HPP
#define SEQUANT2_EXPR_OPERATOR_HPP

namespace sequant2 {

inline ExprPtr
operator*(const ExprPtr &left, const ExprPtr &right) {
  // naive version is to just make a Product
  // TODO why is ExprPtrList needed?
  auto result = std::make_shared<Product>(ExprPtrList{left, right});
  return result;
}

inline ExprPtr
operator+(const ExprPtr &left, const ExprPtr &right) {
  // naive version is to just make a Sum
  // TODO why is ExprPtrList needed?
  auto result = std::make_shared<Sum>(ExprPtrList{left, right});
  return result;
}

inline ExprPtr
operator-(const ExprPtr &left, const ExprPtr &right) {
  auto result = std::make_shared<Sum>(ExprPtrList{left, std::make_shared<Product>(-1.0, ExprPtrList{right})});
  return result;
}

}  // namespace sequant2

#endif //SEQUANT2_EXPR_OPERATOR_HPP
