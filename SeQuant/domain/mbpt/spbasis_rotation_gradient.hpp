#ifndef SEQUANT_DOMAIN_MBPT_SPBASIS_ROTATION_GRADIENT_HPP
#define SEQUANT_DOMAIN_MBPT_SPBASIS_ROTATION_GRADIENT_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/index.hpp>

#include <string>

namespace sequant::mbpt {

/// @return the derivative of the scalar function @p L of orbital-indexed
/// tensors by a rotation of the orbitals: with every bra and ket slot of
/// every leaf labelled as one of @p labels transformed as
/// T{…u…} → T{…u…} + Σ_q κ{q;u} T{…q…}, the result is ∂L/∂κ{q;p} at κ = 0,
///
///   X{q;p} = Σ_T Σ_slots (∂L/∂T{…p…}) T{…q…},
///
/// a sum, over the slots of the leaves, of the term of @p L with @p q in the
/// slot and @p p in place of the slot's index elsewhere. @p p and @p q are
/// the free indices of the result. A slot contributes if its index runs over
/// the space of @p p or one that includes it; a slot whose space is disjoint
/// from that of @p p does not. @p q may run over any space: it decides the
/// blocks of the tensors the result is written in.
///
/// Every slot transforms alike, as for real orbitals; auxiliary slots are not
/// orbitals and do not transform. The antisymmetric part X{q;p} − X{p;q} is
/// the orbital gradient of @p L , and X is the right-hand side that the
/// orbital-response (Z-vector) equations and the energy-weighted density of a
/// nuclear gradient are built from.
/// A leaf spelled as the adjoint of a tensor labelled as one of @p labels
/// rotates too.
/// @throw Exception for a null @p L , one that holds operators, or one that
/// is not a scalar (an index in one slot of a term, or an Â or Ŝ); if @p p
/// or @p q has protoindices or they are one index; if @p L uses the label of
/// @p p or @p q ; and for a leaf to rotate that stands under a TensorPower,
/// has a slot whose index has protoindices or is a protoindex in its term,
/// or runs over a space that overlaps that of @p p without including it
[[nodiscard]] ExprPtr spbasis_rotation_gradient(
    const ExprPtr& L, const container::svector<std::wstring>& labels,
    const Index& p, const Index& q);

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_SPBASIS_ROTATION_GRADIENT_HPP
