//
// Amplitude tensors and the check that they carry their basis grants.
//

#ifndef SEQUANT_DOMAIN_MBPT_BASIS_GRANTS_HPP
#define SEQUANT_DOMAIN_MBPT_BASIS_GRANTS_HPP

#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/domain/mbpt/op_registry.hpp>

namespace sequant::mbpt {

/// @return true if the label of @p t, without its adjoint marker, is a
/// registered Ex or Deex operator of @p reg (perturbation-order decorated
/// labels such as `t¹` must be registered as such)
bool is_amplitude_tensor(const AbstractTensor& t, const OpRegistry& reg);

/// Checks that every granted bra/ket slot of every amplitude tensor of
/// @p expr (see is_amplitude_tensor, default mbpt::Context's registry) carries
/// its grant, `basis_grant(label, slot.space())` with the adjoint marker
/// stripped: a granted leg is minted in a specific basis, which nothing
/// downstream may replace by another basis or by the space's own one. An
/// ungranted slot is unchecked: Wick lets it take a partner's instance, which
/// is exact.
///
/// The grant is looked up by the slot's exact IndexSpace, so @p expr must be
/// spin-free (or closed-shell spin-traced), as the granted spaces are.
/// @throw Exception naming the offending tensor and slot
void assert_amplitudes_carry_granted_basis(const ExprPtr& expr);

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_BASIS_GRANTS_HPP
