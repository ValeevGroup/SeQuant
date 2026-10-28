//
// Created by Eduard Valeyev on 4/2/18.
//

#ifndef SEQUANT_EXPRESSIONS_OPERATORS_HPP
#define SEQUANT_EXPRESSIONS_OPERATORS_HPP

#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_holder.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/expressions/traits.hpp>
#include <SeQuant/core/expressions/variable.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <concepts>
#include <type_traits>
#include <utility>

namespace sequant {

template <typename T>
  requires(std::constructible_from<Constant, T>)
ExprPtr operator+(const std::same_as<ExprPtr> auto &lhs, T &&rhs) {
  return lhs + ex<Constant>(std::forward<T>(rhs));
}

template <typename T>
  requires(std::constructible_from<Constant, T>)
ExprPtr operator+(T &&lhs, const std::same_as<ExprPtr> auto &rhs) {
  return ex<Constant>(std::forward<T>(lhs)) + rhs;
}

template <typename T>
  requires(std::constructible_from<Constant, T>)
ExprPtr operator-(const std::same_as<ExprPtr> auto &lhs, T &&rhs) {
  return lhs - ex<Constant>(std::forward<T>(rhs));
}

template <typename T>
  requires(std::constructible_from<Constant, T>)
ExprPtr operator-(T &&lhs, const std::same_as<ExprPtr> auto &rhs) {
  return ex<Constant>(std::forward<T>(lhs)) - rhs;
}

template <typename T>
  requires(std::constructible_from<Constant, T>)
ExprPtr operator*(const std::same_as<ExprPtr> auto &lhs, T &&rhs) {
  return lhs * ex<Constant>(std::forward<T>(rhs));
}

template <typename T>
  requires(std::constructible_from<Constant, T>)
ExprPtr operator*(T &&lhs, const std::same_as<ExprPtr> auto &rhs) {
  return ex<Constant>(std::forward<T>(lhs)) * rhs;
}

template <typename T>
  requires(std::is_arithmetic_v<std::remove_cvref_t<T>>)
ExprPtr operator/(const std::same_as<ExprPtr> auto &lhs, T &&rhs) {
  return lhs * ex<Constant>(rational(1, std::forward<T>(rhs)));
}

inline ExprPtr operator/(const std::same_as<ExprPtr> auto &lhs,
                         const Constant &rhs) {
  return lhs * ex<Constant>(1.0 / rhs.value());
}

template <typename T>
  requires(std::constructible_from<Variable, T>)
ExprPtr operator+(T &&lhs, const std::same_as<ExprPtr> auto &rhs) {
  return ex<Variable>(std::forward<T>(lhs)) + rhs;
}

template <typename T>
  requires(std::constructible_from<Variable, T>)
ExprPtr operator+(const std::same_as<ExprPtr> auto &lhs, T &&rhs) {
  return lhs + ex<Variable>(std::forward<T>(rhs));
}

template <typename T>
  requires(std::constructible_from<Variable, T>)
ExprPtr operator-(T &&lhs, const std::same_as<ExprPtr> auto &rhs) {
  return ex<Variable>(std::forward<T>(lhs)) - rhs;
}

template <typename T>
  requires(std::constructible_from<Variable, T>)
ExprPtr operator-(const std::same_as<ExprPtr> auto &lhs, T &&rhs) {
  return lhs - ex<Variable>(std::forward<T>(rhs));
}

template <typename T>
  requires(std::constructible_from<Variable, T>)
ExprPtr operator*(T &&lhs, const std::same_as<ExprPtr> auto &rhs) {
  return ex<Variable>(std::forward<T>(lhs)) * rhs;
}

template <typename T>
  requires(std::constructible_from<Variable, T>)
ExprPtr operator*(const std::same_as<ExprPtr> auto &lhs, T &&rhs) {
  return lhs * ex<Variable>(std::forward<T>(rhs));
}

// N.B. arguments that are expressions themselves are handled by the
// ExprHolder operators taking Expr
template <typename T>
  requires(std::constructible_from<Constant, T> &&
           !is_an_expr_v<std::remove_cvref_t<T>>)
ExprHolder operator+(const std::same_as<ExprHolder> auto &lhs, T &&rhs) {
  return lhs + Constant(std::forward<T>(rhs));
}

template <typename T>
  requires(std::constructible_from<Constant, T> &&
           !is_an_expr_v<std::remove_cvref_t<T>>)
ExprHolder operator+(T &&lhs, const std::same_as<ExprHolder> auto &rhs) {
  return Constant(std::forward<T>(lhs)) + rhs;
}

template <typename T>
  requires(std::constructible_from<Constant, T> &&
           !is_an_expr_v<std::remove_cvref_t<T>>)
ExprHolder operator-(const std::same_as<ExprHolder> auto &lhs, T &&rhs) {
  return lhs - Constant(std::forward<T>(rhs));
}

template <typename T>
  requires(std::constructible_from<Constant, T> &&
           !is_an_expr_v<std::remove_cvref_t<T>>)
ExprHolder operator-(T &&lhs, const std::same_as<ExprHolder> auto &rhs) {
  return Constant(std::forward<T>(lhs)) - rhs;
}

template <typename T>
  requires(std::constructible_from<Constant, T> &&
           !is_an_expr_v<std::remove_cvref_t<T>>)
ExprHolder operator*(const std::same_as<ExprHolder> auto &lhs, T &&rhs) {
  return lhs * Constant(std::forward<T>(rhs));
}

template <typename T>
  requires(std::constructible_from<Constant, T> &&
           !is_an_expr_v<std::remove_cvref_t<T>>)
ExprHolder operator*(T &&lhs, const std::same_as<ExprHolder> auto &rhs) {
  return Constant(std::forward<T>(lhs)) * rhs;
}

template <typename T>
  requires(std::is_arithmetic_v<std::remove_cvref_t<T>>)
ExprHolder operator/(const std::same_as<ExprHolder> auto &lhs, T &&rhs) {
  return lhs * Constant(rational(1, std::forward<T>(rhs)));
}

inline ExprHolder operator/(const std::same_as<ExprHolder> auto &lhs,
                            const Constant &rhs) {
  return lhs * Constant(1.0 / rhs.value());
}

template <typename T>
  requires(std::constructible_from<Variable, T> &&
           !is_an_expr_v<std::remove_cvref_t<T>>)
ExprHolder operator+(T &&lhs, const std::same_as<ExprHolder> auto &rhs) {
  return Variable(std::forward<T>(lhs)) + rhs;
}

template <typename T>
  requires(std::constructible_from<Variable, T> &&
           !is_an_expr_v<std::remove_cvref_t<T>>)
ExprHolder operator+(const std::same_as<ExprHolder> auto &lhs, T &&rhs) {
  return lhs + Variable(std::forward<T>(rhs));
}

template <typename T>
  requires(std::constructible_from<Variable, T> &&
           !is_an_expr_v<std::remove_cvref_t<T>>)
ExprHolder operator-(T &&lhs, const std::same_as<ExprHolder> auto &rhs) {
  return Variable(std::forward<T>(lhs)) - rhs;
}

template <typename T>
  requires(std::constructible_from<Variable, T> &&
           !is_an_expr_v<std::remove_cvref_t<T>>)
ExprHolder operator-(const std::same_as<ExprHolder> auto &lhs, T &&rhs) {
  return lhs - Variable(std::forward<T>(rhs));
}

template <typename T>
  requires(std::constructible_from<Variable, T> &&
           !is_an_expr_v<std::remove_cvref_t<T>>)
ExprHolder operator*(T &&lhs, const std::same_as<ExprHolder> auto &rhs) {
  return Variable(std::forward<T>(lhs)) * rhs;
}

template <typename T>
  requires(std::constructible_from<Variable, T> &&
           !is_an_expr_v<std::remove_cvref_t<T>>)
ExprHolder operator*(const std::same_as<ExprHolder> auto &lhs, T &&rhs) {
  return lhs * Variable(std::forward<T>(rhs));
}

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_OPERATORS_HPP
