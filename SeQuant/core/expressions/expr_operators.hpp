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
#include <SeQuant/core/expressions/variable.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <concepts>

namespace sequant {

template <typename T>
  requires(std::constructible_from<Constant, T>)
ExprPtr operator+(const ExprPtr &lhs, T &&rhs) {
  return lhs + ex<Constant>(std::forward<T>(rhs));
}

template <typename T>
  requires(std::constructible_from<Constant, T>)
ExprPtr operator+(T &&lhs, const ExprPtr &rhs) {
  return ex<Constant>(std::forward<T>(lhs)) + rhs;
}

template <typename T>
  requires(std::constructible_from<Constant, T>)
ExprPtr operator-(const ExprPtr &lhs, T &&rhs) {
  return lhs - ex<Constant>(std::forward<T>(rhs));
}

template <typename T>
  requires(std::constructible_from<Constant, T>)
ExprPtr operator-(T &&lhs, const ExprPtr &rhs) {
  return ex<Constant>(std::forward<T>(lhs)) - rhs;
}

template <typename T>
  requires(std::constructible_from<Constant, T>)
ExprPtr operator*(const ExprPtr &lhs, T &&rhs) {
  return lhs * ex<Constant>(std::forward<T>(rhs));
}

template <typename T>
  requires(std::constructible_from<Constant, T>)
ExprPtr operator*(T &&lhs, const ExprPtr &rhs) {
  return ex<Constant>(std::forward<T>(lhs)) * rhs;
}

template <typename T>
  requires(std::is_arithmetic_v<T>)
ExprPtr operator/(const ExprPtr &lhs, T &&rhs) {
  return lhs * ex<Constant>(rational(1, std::forward<T>(rhs)));
}

inline ExprPtr operator/(const ExprPtr &lhs, const Constant &rhs) {
  return lhs * ex<Constant>(1.0 / rhs.value());
}

template <typename T>
  requires(std::constructible_from<Variable, T>)
ExprPtr operator+(T &&lhs, const ExprPtr &rhs) {
  return ex<Variable>(std::forward<T>(lhs)) + rhs;
}

template <typename T>
  requires(std::constructible_from<Variable, T>)
ExprPtr operator+(const ExprPtr &lhs, T &&rhs) {
  return lhs + ex<Variable>(std::forward<T>(rhs));
}

template <typename T>
  requires(std::constructible_from<Variable, T>)
ExprPtr operator-(T &&lhs, const ExprPtr &rhs) {
  return ex<Variable>(std::forward<T>(lhs)) - rhs;
}

template <typename T>
  requires(std::constructible_from<Variable, T>)
ExprPtr operator-(const ExprPtr &lhs, T &&rhs) {
  return lhs - ex<Variable>(std::forward<T>(rhs));
}

template <typename T>
  requires(std::constructible_from<Variable, T>)
ExprPtr operator*(T &&lhs, const ExprPtr &rhs) {
  return ex<Variable>(std::forward<T>(lhs)) * rhs;
}

template <typename T>
  requires(std::constructible_from<Variable, T>)
ExprPtr operator*(const ExprPtr &lhs, T &&rhs) {
  return lhs * ex<Variable>(std::forward<T>(rhs));
}

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_OPERATORS_HPP
