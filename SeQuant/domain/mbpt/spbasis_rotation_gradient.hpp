#ifndef SEQUANT_DOMAIN_MBPT_SPBASIS_ROTATION_GRADIENT_HPP
#define SEQUANT_DOMAIN_MBPT_SPBASIS_ROTATION_GRADIENT_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/index.hpp>

#include <string>

namespace sequant::mbpt {

/// @return the derivative of the function @p L of tensors over the
/// single-particle basis by a transformation of that basis: with every bra
/// and ket slot of every leaf that rotates transformed as
/// T{…u…} → T{…u…} + Σ_q κ{q;u} T{…q…}, the result is ∂L/∂κ{q;p} at κ = 0,
///
///   X{q;p} = Σ_T Σ_slots (∂L/∂T{…p…}) T{…q…},
///
/// a sum, over the slots of the leaves, of the term of @p L with @p q in the
/// slot and @p p in place of the slot's index elsewhere. @p p and @p q are
/// free indices of the result. A slot contributes if its index runs over the
/// space of @p p or one that includes it; a slot whose space is disjoint from
/// that of @p p does not. @p q may run over any space: it decides the blocks
/// of the tensors the result is written in.
///
/// κ is any matrix, not only an antisymmetric one. The antisymmetric part
/// X{q;p} − X{p;q} is the gradient of @p L by a rotation of the basis (the
/// orbital gradient), and X is what the equations for the orbital multipliers
/// of a Lagrangian and the energy-weighted density of a nuclear gradient are
/// built from. Every slot transforms alike, as for real orbitals; auxiliary
/// slots are no slots of the basis and do not transform.
///
/// The leaves that rotate are those labelled as one of @p labels , those
/// spelled as the adjoint of such a tensor, and those with a definition:
///
/// - A leaf with one of @p definitions (the definition whose head has its
///   label and ranks) is a function of the basis also through the summed
///   indices of its expression, e.g. the Fock operator through the occupied
///   orbitals: f{p_1;p_2} = h{p_1;p_2} + g{p_1,i_1;p_2,i_1}. Its own slots
///   rotate as those of a leaf do, which keeps the result in terms of it, and
///   the leaves that rotate in its expression contribute through the slots
///   that hold a summed index. Every tensor of the expression that holds an
///   index of the head in a bra or ket slot must itself rotate: only then
///   does the leaf transform as a tensor. Definitions may use others, but
///   not themselves. A leaf without a definition is a tensor in a fixed
///   basis otherwise: f held as a leaf gives the derivative at a fixed Fock
///   operator.
///
/// @p L may have free indices, each in one slot of a term. The result is the
/// derivative of the tensor of those indices in the transformed basis: a slot
/// that holds a free index u contributes its term with @p q in the slot times
/// the Kronecker delta of u and @p p . The slots that hold an index of
/// @p fixed do not transform. With these the function applies to its own
/// result: X is the first derivative for every κ if the slots that hold
/// @p q are held fixed, since κ enters each slot linearly, so
///
///   spbasis_rotation_gradient(X, labels, r, s, definitions, {q})
///
/// is ∂²L/∂κ{q;p}∂κ{s;r} at κ = 0, and so on to any order, each with the
/// second indices of all earlier ones held fixed ({q, s} for the third).
/// The derivatives by the generator of a unitary transformation exp(κ)
/// differ from these, from the second order on, by terms of lower order.
/// @throw Exception for a null @p L , one that holds operators, or one with
/// an Â or Ŝ; if @p p or @p q has protoindices or they are one index; if
/// @p L or a definition uses the label of @p p or @p q ; for a leaf that
/// rotates that stands under a TensorPower, has a slot whose index has
/// protoindices, is a protoindex in its term or stands in more than two
/// slots of it, or runs over a space that overlaps that of @p p without
/// including it; for a leaf spelled as the adjoint of one with a definition;
/// for a tensor of the expression of a definition that holds an index of the
/// head and does not rotate; and for a definition that depends on itself
[[nodiscard]] ExprPtr spbasis_rotation_gradient(
    const ExprPtr& L, const container::svector<std::wstring>& labels,
    const Index& p, const Index& q,
    const container::svector<ResultExpr>& definitions,
    const container::svector<Index>& fixed = {});

/// spbasis_rotation_gradient() with no definitions and no fixed indices
[[nodiscard]] ExprPtr spbasis_rotation_gradient(
    const ExprPtr& L, const container::svector<std::wstring>& labels,
    const Index& p, const Index& q);

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_SPBASIS_ROTATION_GRADIENT_HPP
