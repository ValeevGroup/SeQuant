#ifndef SEQUANT_EXPRESSIONS_TRAITS_HPP
#define SEQUANT_EXPRESSIONS_TRAITS_HPP

#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/meta.hpp>

namespace sequant {

template <typename T>
constexpr bool is_an_expr_v = meta::is_base_of_v<Expr, T>;
template <typename T>
constexpr bool is_expr_v = meta::is_same_v<Expr, T>;

template <typename T>
constexpr bool is_a_constant_v = meta::is_base_of_v<Constant, T>;
template <typename T>
constexpr bool is_constant_v = meta::is_same_v<Constant, T>;

template <typename T>
constexpr bool is_a_sum_v = meta::is_base_of_v<Sum, T>;
template <typename T>
constexpr bool is_sum_v = meta::is_same_v<Sum, T>;

template <typename T>
constexpr bool is_a_product_v = meta::is_base_of_v<Product, T>;
template <typename T>
constexpr bool is_product_v = meta::is_same_v<Product, T>;

template <typename T>
constexpr bool is_a_variable_v = meta::is_base_of_v<Variable, T>;
template <typename T>
constexpr bool is_variable_v = meta::is_same_v<Variable, T>;

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_TRAITS_HPP
