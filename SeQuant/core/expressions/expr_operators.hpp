//
// Created by Eduard Valeyev on 4/2/18.
//

#ifndef SEQUANT_EXPRESSIONS_OPERATORS_HPP
#define SEQUANT_EXPRESSIONS_OPERATORS_HPP

#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <memory>

namespace sequant {

inline bool operator==(const ExprPtr &left, const ExprPtr &right) {
  return *left == *right;
}

inline ExprPtr operator*(const ExprPtr &left, const ExprPtr &right) {
  if (left.is<Constant>() && right.is<Constant>()) {
    auto c_ = left->clone();
    auto &c = c_.as<Constant>();
    c *= right.as<Constant>();
    return c_;
  }

  auto left_is_product = left->is<Product>();
  auto right_is_product = right->is<Product>();
  if (!left_is_product && !right_is_product) {
    return ex<Product>(ExprPtrList{left, right});
  } else if (left_is_product) {
    auto result = std::static_pointer_cast<Product>(left->clone());
    result->append(1, right);
    return result;
  } else {  // right_is_product
    auto result = std::static_pointer_cast<Product>(right->clone());
    result->prepend(1, left);
    return result;
  }

  SEQUANT_UNREACHABLE;
}

/// Unlike @code operator*(const ExprPtr&, const ExprPtr&) @endcode this
/// produces a non-commutative product (i.e. NCProduct)
inline ExprPtr operator^(const ExprPtr &left, const ExprPtr &right) {
  auto left_is_product = left->is<Product>();
  auto right_is_product = right->is<Product>();
  if (!left_is_product && !right_is_product) {
    return ex<NCProduct>(ExprPtrList{left, right});
  } else if (left_is_product) {
    auto result = std::make_shared<NCProduct>(left->clone().as<Product>());
    result->append(1, right);
    return result;
  } else {  // right_is_product
    auto result = std::make_shared<NCProduct>(right->clone().as<Product>());
    result->prepend(1, left);
    return result;
  }

  SEQUANT_UNREACHABLE;
}

template <typename T>
  requires(std::constructible_from<Constant, T>)
inline ExprPtr operator*(T left, const ExprPtr &right) {
  return ex<Constant>(std::move(left)) * right;
}

template <typename T>
  requires(std::constructible_from<Constant, T>)
inline ExprPtr operator*(const ExprPtr &left, T right) {
  return left * ex<Constant>(std::move(right));
}

inline ExprPtr operator+(const ExprPtr &left, const ExprPtr &right) {
  auto left_is_sum = left->is<Sum>();
  auto right_is_sum = right->is<Sum>();
  if (!left_is_sum && !right_is_sum) {
    return ex<Sum>(ExprPtrList{left, right});
  } else if (left_is_sum) {
    auto result = std::static_pointer_cast<Sum>(left->clone());
    result->append(right);
    return result;
  } else {  // right_is_sum
    auto result = std::static_pointer_cast<Sum>(right->clone());
    result->prepend(left);
    return result;
  }

  SEQUANT_UNREACHABLE;
}

inline ExprPtr operator-(const ExprPtr &left, const ExprPtr &right) {
  auto left_is_sum = left->is<Sum>();
  if (!left_is_sum) {
    return ex<Sum>(ExprPtrList{
        left,
        (right->is<Constant>() ? ex<Constant>(-right->as<Constant>().value())
                               : ex<Product>(-1, ExprPtrList{right}))});
  } else if (left_is_sum) {
    auto result = std::static_pointer_cast<Sum>(left->clone());
    if (right->is<Constant>())
      result->append(ex<Constant>(-right->as<Constant>().value()));
    else
      result->append(ex<Product>(-1, ExprPtrList{right}));
    return result;
  }

  SEQUANT_UNREACHABLE;
}

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_OPERATORS_HPP
