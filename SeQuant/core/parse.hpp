#ifndef SEQUANT_PARSE_HPP
#define SEQUANT_PARSE_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/result_expr.hpp>
#include <SeQuant/core/tensor.hpp>

#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

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
///             Every tensor may optionally be annoted with index symmetry specifications. The general syntax is
///             <tensorSpec> [:<perm symm> [-<braket symm> [-<particle symm>]]]
///             (no whitespace is allowed at this place). Examples are
///             't{i1;i2}:A', 't{i1;i2}:A-S', 't{i1;i2}:N-C-S'
///             Possible values for <perm symm> are
///             - 'A' for antisymmetry (sequant::Symmetry::antisymm)
///             - 'S' for symmetric (sequant::Symmetry::symm)
///             - 'N' for non-symmetric (sequant::Symmetry::nonsymm)
///             Possible values for <braket symm> are
///             - 'C' for antisymmetry (sequant::BraKetSymmetry::conjugate)
///             - 'S' for symmetric (sequant::BraKetSymmetry::symm)
///             - 'N' for non-symmetric (sequant::BraKetSymmetry::nonsymm)
///             Possible values for <particle symm> are
///             - 'S' for symmetric (sequant::ParticleSymmetry::symm)
///             - 'N' for non-symmetric (sequant::ParticleSymmetry::nonsymm)
/// \param perm_symm Default index permutation symmetry to be used if tensors don't specify a permutation
///                  symmetry explicitly.
/// \param braket_symm Default BraKet symmetry to be used if tensors don't specify a BraKet symmetry explicitly.
/// \param particle_symm Default particle symmetry to be used if tensors don't specify a particle symmetry explicitly.
///                   @c raw expression. Explicit tensor symmetry can
///                   be annotated in the expression itself. In that case, the
///                   annotated symmetry will be used.
///                   eg. 'g{i1, a1; i2, a2}:A' tensor with 'sequant::Symmetry::antisymm' annotation
///                       'g{i1, a1; i2, a2}:S' tensor with 'sequant::Symmetry::symm' annotation
///                       'g{i1, a1; i2, a2}:N' tensor with 'sequant::Symmetry::nonsymm' annotation
/// \return SeQuant expression.
// clang-format on
ExprPtr parse_expr(std::wstring_view raw,
                   std::optional<Symmetry> perm_symm = {},
                   std::optional<BraKetSymmetry> braket_symm = {},
                   std::optional<ParticleSymmetry> particle_symm = {});

ResultExpr parse_result_expr(
    std::wstring_view raw, std::optional<Symmetry> perm_symm = {},
    std::optional<BraKetSymmetry> braket_symm = {},
    std::optional<ParticleSymmetry> particle_symm = {});

///
/// Get a parsable string from an expression.
///
/// An expression that has a flat structure (ie. product does not contain
/// sum or product subexpressions) is guaranteed to be reparsable to itself so
/// that equality comparison holds: x == parse_expr(deparse_expr(x)). For nested
/// expressions the equality comparison is not guaranteed, however, for such
/// an expression x, its corresponding evaluation node, eval_node(x), and
/// the parsed evaluation node, eval_node(parse_expr(deparse_expr(x))) are
/// equivalent.
///
/// \param expr Expression to stringify that can be re-parsed to itself.
/// \param annot_sym Whether to add sequant::Symmetry annotation
///                  to each Tensor string.
/// \return wstring of the expression.
std::wstring deparse(const ResultExpr &expr, bool annot_sym = true);
std::wstring deparse(const ExprPtr &expr, bool annot_sym = true);
std::wstring deparse(const Expr &expr, bool annot_sym = true);
std::wstring deparse(const Product &product, bool annot_sym);
std::wstring deparse(const Sum &sum, bool annot_sym);
std::wstring deparse(const Tensor &tensor, bool annot_sym = true);
std::wstring deparse(const AbstractTensor &tensor, bool annot_sym = true);
template <Statistics S>
std::wstring deparse(const NormalOperator<S> &nop);
std::wstring deparse(const Variable &variable);
std::wstring deparse(const Constant &constant);
std::wstring deparse(const Index &index);

}  // namespace sequant

#endif  // SEQUANT_PARSE_HPP
