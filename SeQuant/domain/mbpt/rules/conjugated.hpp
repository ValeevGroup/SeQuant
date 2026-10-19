#ifndef SEQUANT_DOMAIN_MBPT_RULES_CONJUGATED_HPP
#define SEQUANT_DOMAIN_MBPT_RULES_CONJUGATED_HPP

#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <utility>

namespace sequant::mbpt::detail {

/// @brief applies a factorization rule to a tensor that may carry a core
/// state.
///
/// The rule builds its factors from the tensor's slots, so it is given the
/// bare array: for `⁺` the array over the exchanged bundles, since `t⁺{q;p}`
/// denotes `conj t{p;q}`. Its factorization is then conjugated back, the
/// adjoint and the K-conjugation distributing over a product of c-number
/// tensors; what either does to a factor is decided by the factor's own
/// traits, so e.g. a Hermitian factor just exchanges its bundles.
/// @param tnsr the tensor to factorize
/// @param rule the rule, applied to a Tensor without states
/// @return the factorization of the value @p tnsr denotes
template <typename Rule>
[[nodiscard]] ExprPtr factorize_conjugated(const Tensor &tnsr, Rule &&rule) {
  if (!tnsr.adjointed() && !tnsr.kconjugated())
    return std::forward<Rule>(rule)(tnsr);
  Tensor bare = tnsr;
  const auto sign = bare.set_states(false, false);
  SEQUANT_ENFORCE(sign == 1, "clearing the states consumes no sign");
  if (tnsr.adjointed()) static_cast<AbstractTensor &>(bare)._swap_bra_ket();
  ExprPtr fit = std::forward<Rule>(rule)(bare);
  if (tnsr.adjointed()) fit = sequant::adjoint(fit);
  if (tnsr.kconjugated()) fit = sequant::kconjugate(fit);
  return fit;
}

}  // namespace sequant::mbpt::detail

#endif  // SEQUANT_DOMAIN_MBPT_RULES_CONJUGATED_HPP
