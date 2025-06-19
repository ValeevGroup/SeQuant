//
// Created by Eduard Valeyev on 3/8/25.
//

#ifndef SEQUANT_DOMAIN_MBPT_RULES_DF_HPP
#define SEQUANT_DOMAIN_MBPT_RULES_DF_HPP

#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/space.hpp>

#include <string_view>

namespace sequant::mbpt {

///
/// Converts the 4-center (2-electron integral) tensors into a (antisymmetrized)
/// product of two rank-3 tensors.
///
/// \param expr The expression to be density-fit.
/// \param aux_space The index space representing the auxiliary indices
/// introduced through the decomposition.
/// \param tensor_label The label off the tensor that shall be decomposed
/// \param factor_label The label of the rank-3 tensors used in the
/// decomposition.
/// \return The density-fitted expression (potentially unchanged,
/// if the target tensor was not contained in the given expression)
///
[[nodiscard]] ExprPtr density_fit(ExprPtr const& expr, IndexSpace aux_space,
                                  std::wstring_view tensor_label,
                                  std::wstring_view factor_label);

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_RULES_DF_HPP
