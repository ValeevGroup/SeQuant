#ifndef SEQUANT_PARSE_EXPR_HPP
#define SEQUANT_PARSE_EXPR_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <regex>
#include <string>

///
///  Create SeQuant expression from string input.
///
/// @author: Bimal Gaudel
/// version: 21 July, 2021
///

namespace sequant {

ExprPtr parse_expr(std::wstring_view raw, Symmetry tensor_sym);

ExprPtr parse_expr_asymm(std::wstring_view raw);

}  // namespace sequant::utils

#endif  // SEQUANT_PARSE_EXPR_HPP
