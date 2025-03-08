//
// Created by Eduard Valeyev on 3/8/25.
//

#ifndef MPQC_DOMAIN_MBPT_RULES_DF_HPP
#define MPQC_DOMAIN_MBPT_RULES_DF_HPP

#include "SeQuant/core/expr.hpp"

namespace sequant::mbpt {

///
/// Converts the 4-center 'g' tensors into a product of two rank-3 tensors.
///
/// \param expr The expression to be density-fit.
/// \param aux_label The label of the introduced auxilliary index. eg. 'x', 'p'.
/// \return The density-fit expression if 'g' of rank-4 present, otherwise the
///         input expression itself will be returned.
///
ExprPtr density_fit(ExprPtr const& expr, std::wstring const& aux_label);

}  // namespace sequant::mbpt

#endif  // MPQC_DOMAIN_MBPT_RULES_DF_HPP
