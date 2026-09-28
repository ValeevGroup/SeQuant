//
// Created by Ajay Melekamburath on 4/27/26.
//

#ifndef SEQUANT_CORE_EXPORT_UTILS_HPP
#define SEQUANT_CORE_EXPORT_UTILS_HPP

#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/power.hpp>

#include <string>
#include <string_view>

namespace sequant::detail {

/// How sanitize_identifier treats non-ASCII characters
enum class NonAsciiPolicy {
  /// Keep them as they are (for languages accepting Unicode identifiers)
  Keep,
  /// Replace each by `u` followed by its code point in hex (e.g. `u3b3`)
  Escape,
};

/// Turns a label into a valid identifier: ASCII characters other than
/// alphanumerics and `_` are replaced by `_`, and a leading digit is prefixed
/// with `_`.
/// @param label The label to sanitize
/// @param policy What to do with non-ASCII characters
/// @return The UTF-8 encoded identifier
std::string sanitize_identifier(std::wstring_view label, NonAsciiPolicy policy);

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

}  // namespace sequant::detail

#endif  // SEQUANT_CORE_EXPORT_UTILS_HPP
