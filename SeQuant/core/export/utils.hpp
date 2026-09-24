//
// Created by Ajay Melekamburath on 4/27/26.
//

#ifndef SEQUANT_CORE_EXPORT_UTILS_HPP
#define SEQUANT_CORE_EXPORT_UTILS_HPP

#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/power.hpp>

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

}  // namespace sequant::detail

#endif  // SEQUANT_CORE_EXPORT_UTILS_HPP
