//
// Created by Eduard Valeyev on 3/8/25.
//

#ifndef SEQUANT_DOMAIN_MBPT_RULES_CSV_HPP
#define SEQUANT_DOMAIN_MBPT_RULES_CSV_HPP

#include <SeQuant/core/basis.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/reserved.hpp>

#include <concepts>

namespace sequant::mbpt {

///
/// expands CSVs in an expression in terms of a basis (standard unoccupieds,
/// PAOs, AOs, etc.)
///
/// \param expr The expression to be CSV-transformed.
/// \param csv_basis the basis in terms of which the CSVs are expanded: a
///                  space (its own basis, orthonormal) or a basis instance
///                  registered under a name, whose registry entry says whether
///                  it is orthonormal (IndexBasis::metric()); the rank-1/1 CSV
///                  overlap of an orthonormal basis expands into `C·C` instead
///                  of `C·s·C`
/// \param coeff_tensor_label The label of the CSV-transformation tensors that
///                           will be introduced.
/// \param tensor_labels The labels of the tensors that will be
///                    transformed
/// \return The CSV-transformed expression if CSV-tensors with labels present
///         in @c tensor_labels appear in @c expr. Otherwise returns the input
///         expression itself.
/// \throw Exception if @p csv_basis has a basis instance that the default
///        context's registry does not name
ExprPtr csv_transform(ExprPtr const& expr, const IndexBasis& csv_basis,
                      std::wstring const& coeff_tensor_label = L"C",
                      container::svector<std::wstring> const& tensor_labels = {
                          L"f", L"g", sequant::reserved::overlap_label()});

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_RULES_CSV_HPP
