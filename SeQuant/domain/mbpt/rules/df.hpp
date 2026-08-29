//
// Created by Eduard Valeyev on 3/8/25.
//

#ifndef SEQUANT_DOMAIN_MBPT_RULES_DF_HPP
#define SEQUANT_DOMAIN_MBPT_RULES_DF_HPP

#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/space.hpp>

#include <functional>
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
/// \param should_split Optional predicate consulted for every tensor that
/// would otherwise be decomposed (label matches, bra/ket net rank 2, no aux
/// indices); returning false keeps that tensor as a 4-center leaf. Electron
/// \c k of the tensor is the index pair (bra[k], ket[k]), i.e. the pair that
/// the decomposition assigns to the k-th rank-3 factor. Empty (default) =
/// decompose every matching tensor.
/// \return The density-fitted expression (potentially unchanged,
/// if the target tensor was not contained in the given expression)
///
[[nodiscard]] ExprPtr density_fit(
    ExprPtr const& expr, IndexSpace aux_space, std::wstring_view tensor_label,
    std::wstring_view factor_label,
    std::function<bool(Tensor const&)> const& should_split = {});

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_RULES_DF_HPP
