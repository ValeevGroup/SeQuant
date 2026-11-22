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

/// @return the Hartree-Fock energy through @p order in κ, written in h and
/// g: ⟨0|H|0⟩ + rotated_energy(order) with every f{p;q} replaced by
/// h{p;q} + Σ_i g{p,i;q,i} (i new in the reference-occupied space),
/// simplified; h has the symmetries of f, g those of op::H()'s g
/// @throw Exception what rotated_energy() throws, and for a spin-free basis
[[nodiscard]] ExprPtr energy(std::size_t order);

/// @return @p expr in the integrals of symmetrically orthonormalized
/// orbitals, a series in @p strength through @p order : every h and g gets
/// its series() and, on each leg, a factor X = S^(−1/2), S(0) = 1, written
/// Σ_k C(−1/2,k) Δ^k where Δ is the series of the overlap without its first
/// term; X{leg;n} on a bra leg, X{n;leg} on a ket leg, n new in the complete
/// space. Each leaf is cut at total order @p order in @p strength , so a
/// derivative in @p strength of order above @p order is not that of the
/// full function. The overlap coefficients are labelled like the coefficients
/// of series() of the overlap (reserved::overlap_label(): s¹, s², …), have one
/// slot per side, Symmetry::Nonsymm and the bra-ket and column symmetries of
/// h, i.e. of op::H()'s f; to differentiate by one, spell it so, e.g.
/// s¹{p_1;p_2}:N-C-S. Unlike rotated_energy(), needs neither κₒ nor the
/// spaces of the rotation
/// @throw Exception if @p expr holds an mbpt::Operator, f (unless it is
/// labelled h), h⁺ or g⁺ (also under a TensorPower), or an h or g with a leg
/// that carries protoindices or a spin, and what series() throws for @p expr ,
/// @p strength and @p order with h, g or the overlap, but for a Variable
/// labelled like the overlap or an overlap under a TensorPower
[[nodiscard]] ExprPtr orthonormalized(const ExprPtr& expr,
                                      const ExprPtr& strength,
                                      std::size_t order);

}  // namespace sequant::mbpt::hf

#endif  // SEQUANT_DOMAIN_MBPT_MODELS_HF_HPP
