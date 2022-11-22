#ifndef SEQUANT_PARSE_EXPR_HPP
#define SEQUANT_PARSE_EXPR_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <string_view>

///
///  Create SeQuant expression from string input.
///
/// @author: Bimal Gaudel
/// version: 21 July, 2021
///

namespace sequant {

///
/// \param raw A tensor algebra expression. A valid expression is
///            a product, a sum or a mix of those including parentheses.
///            eg. 'A{i1, a1; i2, a2}' A is the label (alphanumeric), i1, a1 are bra indices, i2, a2 are ket indices.
///                'A{i_1, a_1; i_2, a_2}' same as above with alternate notation for indices
///                'A_{i1, a1}^{i2, a2}' same as above with alternate notation for bra and ket
///                'A^{i2, a2}_{i1, a1}' same as above with alternate notation for bra and ket
///                'A{i1,i2; a1,a2} + B{i2,i1; a1,a2}' a sum of tensors
///                'A{i1,i2; a3,a4} * B{a3,a4; a1,a2}' a product of tensors
///                'A{i1,i2; a3,a4} * B{a3,a4; a1,a2} + C{i1,i2;a1,a2}' a sum and a product of tensors
///                'A{i1,i2; a3,a4} * (B{a3,a4; a1,a2} + C{a3,a4; a1,a2}) a parenthesized expression
/// \param tensor_sym The symmetry of all atomic tensors in the
///                   @c raw expression. Explicit tensor symmetry can
///                   be annotated in the expression itself. In which case, the
///                   annotated symmetry will be used.
///                   eg. 'g{i1, a1; i2, a2}:A' tensor with 'sequant::Symmetry::antisymm' annotation
///                       'g{i1, a1; i2, a2}:S' tensor with 'sequant::Symmetry::symm' annotation
///                       'g{i1, a1; i2, a2}:N' tensor with 'sequant::Symmetry::nonsymm' annotation
/// \return SeQuant expression.
ExprPtr parse_expr(std::wstring_view raw, Symmetry tensor_sym);

}  // namespace sequant::utils

#endif  // SEQUANT_PARSE_EXPR_HPP
