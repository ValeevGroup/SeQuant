#ifndef SEQUANT_DOMAIN_MBPT_MODELS_HF_HPP
#define SEQUANT_DOMAIN_MBPT_MODELS_HF_HPP

#include <SeQuant/core/expr_fwd.hpp>

#include <cstddef>

namespace sequant::mbpt::hf {

/// @return the Hartree-Fock energy as a function of an orbital rotation,
/// ⟨0|e^(−κ̂) H e^(κ̂)|0⟩ − ⟨0|H|0⟩ through @p order in κ, simplified, where
/// κ̂ = K̂ − K̂⁺ and K̂ = Σ κₒ{a;i} a†_a a_i takes every orbital of the
/// reference-occupied space (i) into every orbital of its complement (a),
/// which in an F12 context includes the CABS orbitals.
/// H is op::H(), normal ordered, so ⟨0|H|0⟩ is not part of the result. Its
/// derivatives at κ = 0 are the electronic gradient and Hessian
/// @throw Exception if the default mbpt context's operator registry lacks
/// κₒ (make_legacy_registry() and make_minimal_registry() have it), the
/// vacuum is not Vacuum::SingleProduct, the reference is not the vacuum, the
/// unoccupied space is not registered, or the mbpt context has CSV::Yes
/// (orbital rotations with pair-specific orbitals are not supported)
[[nodiscard]] ExprPtr rotated_energy(std::size_t order);

}  // namespace sequant::mbpt::hf

#endif  // SEQUANT_DOMAIN_MBPT_MODELS_HF_HPP
