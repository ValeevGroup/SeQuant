//
// Created by Eduard Valeyev on 3/8/25.
//

#ifndef SEQUANT_DOMAIN_MBPT_RULES_CSV_HPP
#define SEQUANT_DOMAIN_MBPT_RULES_CSV_HPP

#include <SeQuant/core/expr.hpp>

namespace sequant::mbpt {

///
/// expands CSVs in an expression in terms of a basis (standard unoccupieds,
/// PAOs, AOs, etc.)
///
/// \param expr The expression to be CSV-transformed.
/// \param csv_basis the basis in terms of which the CSVs are expanded
/// \param coeff_tensor_label The label of the CSV-tranformation tensors that
///                           will be introduced.
/// \param tensor_labels The labels of the tensors that will be
///                    transformed
/// \return The CSV-transformed expression if CSV-tensors with labels present
///         in @c csv_tensors appear in @c expr. Otherwise returns the input
///         expression itself.
ExprPtr csv_transform(ExprPtr const& expr, const IndexSpace& csv_basis,
                      std::wstring const& coeff_tensor_label = L"C",
                      container::svector<std::wstring> const& tensor_labels = {
                          L"f", L"g", overlap_label()});

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_RULES_CSV_HPP
