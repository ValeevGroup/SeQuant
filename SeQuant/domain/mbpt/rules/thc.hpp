//
// Created by Oliver Backhouse on 22/11/2025.
//

#ifndef SEQUANT_DOMAIN_MBPT_RULES_THC_HPP
#define SEQUANT_DOMAIN_MBPT_RULES_THC_HPP

#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/space.hpp>

#include <string_view>

namespace sequant::mbpt {

// clang-format off
///
/// Factorizes 2-particle tensors into five rank-2 tensors using
/// the ``tensor hypercontraction'' topology (see [DOI 10.1063/1.4732310](https://doi.org/10.1063/1.4732310);
/// also see [DOI 10.1021/acs.jctc.0c01310](https://doi.org/10.1021/acs.jctc.0c01310) for discussion
/// in the context of related pseudospectral and CP factorizations).
/// Namely, \f$ g_{b_1 b_2}^{k_1 k_2} \f$ is factorized into
/// \f$ B_{b_1}[r_1] B^{k_1}[r_1] B_{b_2}[r_2] B^{k_2}[r_2] C[r_1,r_2] \f$
///
/// \param expr The expression to be tensor-hyper-contracted.
/// \param aux_space The index space representing the auxiliary indices (\f$ r_1, r_2 \f$ in the example above)
/// introduced through the decomposition.
/// \param tensor_label The label of the tensor that shall be decomposed
/// \param outer_factor_label The label of the outer factor tensors (\f$ B \f$ in the example above).
/// \param core_tensor_label The label of the core tensor (\f$ C \f$ in the example above).
// clang-format on
[[nodiscard]] ExprPtr tensor_hypercontract(ExprPtr const& expr,
                                           IndexSpace aux_space,
                                           std::wstring_view tensor_label,
                                           std::wstring_view outer_factor_label,
                                           std::wstring_view core_tensor_label);

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_RULES_THC_HPP
