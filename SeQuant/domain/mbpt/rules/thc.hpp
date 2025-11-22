//
// Create by Oliver Backhouse on 22/11/2025.
//

#ifndef SEQUANT_DOMAIN_MBPT_RULES_THC_HPP
#define SEQUANT_DOMAIN_MBPT_RULES_THC_HPP

#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/space.hpp>

#include <string_view>

namespace sequant::mbpt {

///
/// Converts the 4-center (2-electron integral) tensors into a tensor hyper-contracted
/// product of five rank-2 tensors.
///
/// \param expr The expression to be tensor-hyper-contracted.
/// \param aux_space_1 The index space representing the auxiliary indices
/// introduced through the decomposition.
/// \param tensor_label The label off the tensor that shall be decomposed
/// \param factor_label The label of the rank-2 mixed tensors used in the
/// decomposition.
/// \param aux_label The label of the rank-2 auxiliary tensor used in the decomposition.
[[nodiscard]] ExprPtr tensor_hypercontract(ExprPtr const& expr, IndexSpace aux_space,
                                          std::wstring_view tensor_label,
                                          std::wstring_view factor_label,
                                          std::wstring_view aux_label);

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_RULES_THC_HPP
