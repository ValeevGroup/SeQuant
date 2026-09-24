//
// Created by Ajay Melekamburath on 4/27/26.
//

#ifndef SEQUANT_CORE_EXPORT_UTILS_HPP
#define SEQUANT_CORE_EXPORT_UTILS_HPP

#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/power.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/expressions/variable.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <functional>
#include <string>
#include <string_view>

namespace sequant::detail {

/// Formats a Power exponent for export framework
/// @param exponent the rational exponent
/// @param double_slash if true, use Julia's `//` rational syntax; otherwise
///        use `/` (Python style)
/// @return a string such as `2`, `(-3)`, `(1/2)`, `(-1//3)`
std::string format_power_exponent(const Power::exponent_type &exponent,
                                  bool double_slash);

/// Parenthesizes an already-stringified Power base when needed to keep
/// exponentiation precedence unambiguous in the target language.
/// @param base the Power base expression
/// @param base_str @p base already rendered to a string by the caller
/// @return @p base_str, wrapped in parens iff @p base is a Constant whose
///         value is a non-integer or negative real
std::string format_power_base(const ExprPtr &base, std::string base_str);

/// Formats a Power for export framework
/// @param power the Power to format
/// @param base_str the base of @p power already rendered to a string
/// @param pow_op the exponentiation operator of the target language
/// @param double_slash see format_power_exponent
/// @param wrap_conj renders the complex conjugate of its argument, applied to
///        a conjugated Variable base and to a conjugated @p power
/// @return the formatted power
std::string format_power(
    const Power &power, std::string base_str, std::string_view pow_op,
    bool double_slash,
    const std::function<std::string(std::string)> &wrap_conj);

/// Renders @p expr in infix notation: sums as `a + b`, products as their
/// scalar prefactor (if not 1) and factors joined by @p product_separator
/// @param expr the expression to render
/// @param product_separator the separator between the factors of a Product
/// @param represent renders a Tensor, Variable, Constant or Power; invoked
///        with the expression cast to its concrete type
/// @param who names the caller in the exception message
/// @throw Exception if @p expr contains any other type of expression
template <typename Represent>
std::string stringify_expr(const Expr &expr, std::string_view product_separator,
                           Represent &&represent, std::string_view who) {
  if (expr.is<Tensor>()) {
    return represent(expr.as<Tensor>());
  } else if (expr.is<Variable>()) {
    return represent(expr.as<Variable>());
  } else if (expr.is<Constant>()) {
    return represent(expr.as<Constant>());
  } else if (expr.is<Power>()) {
    return represent(expr.as<Power>());
  } else if (expr.is<Product>()) {
    const Product &product = expr.as<Product>();
    std::string repr;

    if (!product.scalar().is_identity()) {
      repr += represent(Constant(product.scalar()));
      repr += product_separator;
    }

    for (std::size_t i = 0; i < product.size(); ++i) {
      repr +=
          stringify_expr(*product.factor(i), product_separator, represent, who);

      if (i + 1 < product.size()) {
        repr += product_separator;
      }
    }

    return repr;
  } else if (expr.is<Sum>()) {
    const Sum &sum = expr.as<Sum>();
    std::string repr;

    for (std::size_t i = 0; i < sum.size(); ++i) {
      repr +=
          stringify_expr(*sum.summand(i), product_separator, represent, who);

      if (i + 1 < sum.size()) {
        repr += " + ";
      }
    }

    return repr;
  }

  throw Exception("Unsupported expression type in " + std::string(who));
}

}  // namespace sequant::detail

#endif  // SEQUANT_CORE_EXPORT_UTILS_HPP
