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

namespace detail {
ExprPtr csv_transform(ExprPtr const& expr, const IndexBasis& csv_basis,
                      bool orthonormal, std::wstring const& coeff_tensor_label,
                      container::svector<std::wstring> const& tensor_labels);
}  // namespace detail

///
/// expands CSVs in an expression in terms of a basis (standard unoccupieds,
/// PAOs, AOs, etc.)
///
/// \param expr The expression to be CSV-transformed.
/// \param csv_basis the basis in terms of which the CSVs are expanded; a
///                  basis instance must be registered under a name. The
///                  registry's entry says whether the basis is orthonormal
///                  (IndexBasis::metric()); the rank-1/1 CSV overlap then
///                  expands into `C·C` instead of `C·s·C`
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

///
/// expands CSVs in an expression in terms of a basis space: its own basis,
/// orthonormal unless the space carries the LCAOQNS::ao or LCAOQNS::pao bit
///
/// \param expr The expression to be CSV-transformed.
/// \param csv_basis the basis in terms of which the CSVs are expanded
/// \param coeff_tensor_label The label of the CSV-transformation tensors that
///                           will be introduced.
/// \param tensor_labels The labels of the tensors that will be
///                    transformed
/// \return The CSV-transformed expression if CSV-tensors with labels present
///         in @c tensor_labels appear in @c expr. Otherwise returns the input
///         expression itself.
ExprPtr csv_transform(ExprPtr const& expr, const IndexSpace& csv_basis,
                      std::wstring const& coeff_tensor_label = L"C",
                      container::svector<std::wstring> const& tensor_labels = {
                          L"f", L"g", sequant::reserved::overlap_label()});

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_RULES_CSV_HPP
