#ifndef SEQUANT_PARSE_EXPR_HPP
#define SEQUANT_PARSE_EXPR_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expr_fwd.hpp>

#include <string_view>
#include <stdexcept>
#include <string>

///
///  Create SeQuant expression from string input.
///
/// @author: Bimal Gaudel
/// version: 21 July, 2021
///

namespace sequant {

struct ParseError : std::runtime_error {
	std::size_t offset;
	std::size_t length;

	ParseError(std::size_t offset, std::size_t length, std::string message);
};

// clang-format off
///
/// \param raw A tensor algebra expression. A valid expression is
///            a product, a sum or a mix of those including parentheses.
///            eg. 'A{i1, a1; i2, a2}' A is the label (non-space), i1, a1 are bra indices, i2, a2 are ket indices.
///                'A{i_1, a_1; i_2, a_2}' same as above with alternate notation for indices
///                'A_{i1, a1}^{i2, a2}' same as above with alternate notation for bra and ket
///                'A^{i2, a2}_{i1, a1}' same as above with alternate notation for bra and ket
///                'A{i1,i2; a1,a2} + B{i2,i1; a1,a2}' a sum of tensors
///                'A{i1,i2; a3,a4} * B{a3,a4; a1,a2}' a product of tensors
///                'A{i1,i2; a3,a4} * B{a3,a4; a1,a2} + C{i1,i2;a1,a2}' a sum and a product of tensors
///                'A{i1,i2; a3,a4} * (B{a3,a4; a1,a2} + C{a3,a4; a1,a2}) a parenthesized expression
///                '0.5 * t{i1;a1} * f{i1; a1}' tensor product with a scalar
///                '1/2 * t{i1;a1} * f{i1; a1}' same as above (fractions supported)
///                '1./2. * t{i1;a1} * f{i1; a1}' same as above num. and denom. are automatically cast to double
///                '1.0/2.0 * t{i1;a1} * f{i1; a1}' same as above
///                't{i1,i2; a1<i1,i2>, a2<i1,i2>}' a tensor having indices with proto indices.
///                                                a1<i1,i2> is an index with i1 and i2 as proto-indices.
/// \param tensor_sym The symmetry of all atomic tensors in the
///                   @c raw expression. Explicit tensor symmetry can
///                   be annotated in the expression itself. In that case, the
///                   annotated symmetry will be used.
///                   eg. 'g{i1, a1; i2, a2}:A' tensor with 'sequant::Symmetry::antisymm' annotation
///                       'g{i1, a1; i2, a2}:S' tensor with 'sequant::Symmetry::symm' annotation
///                       'g{i1, a1; i2, a2}:N' tensor with 'sequant::Symmetry::nonsymm' annotation
/// \return SeQuant expression.
// clang-format on
ExprPtr parse_expr(std::wstring_view raw,
                   Symmetry tensor_sym = Symmetry::nonsymm);

///
/// Get a parsable string from an expression.
///
/// \param expr Expression to stringify that can be re-parsed to itself.
/// \param annot_sym Whether to add sequant::Symmetry annotation
///                  to each Tensor string.
/// \return wstring of the expression.
std::wstring deparse_expr(ExprPtr expr, bool annot_sym = true);

}  // namespace sequant

#endif  // SEQUANT_PARSE_EXPR_HPP
