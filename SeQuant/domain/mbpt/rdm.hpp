//
// Created by Conner Masteran on 7/1/21.
//

#ifndef SEQUANT_DOMAIN_MBPT_RDM_HPP
#define SEQUANT_DOMAIN_MBPT_RDM_HPP

#include <SeQuant/domain/mbpt/antisymmetrizer.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

namespace sequant {
namespace mbpt {
/// decompositions of reference densities and normal-ordered operators; they
/// are spin-orbital: every density they build is a spin-orbital γ
/// (antisymmetric if multi-body), never a spin-free Γ
namespace decompositions {

/// @return the expansion of the cumulant κ_k @p ex_ (any k ≥ 1) in densities,
/// κ_k = Σ_λ (-1)^(m-1) (m-1)! A[γ_λ₁ ⋯ γ_λₘ] over the integer partitions λ of
/// k into m parts, where A[⋯] is the sum of the distinct antisymmetrized terms
/// that mbpt::antisymmetrize generates, e.g.
/// κ₃ = γ₃ - A[γ₂γ₁] + 2 A[γ₁γ₁γ₁]; the result is expanded, not simplified
ExprPtr cumulant_to_density(ExprPtr ex_);

/// cumulant_to_density for a κ₂
ExprPtr cumulant2_to_density(ExprPtr ex_);

/// cumulant_to_density for a κ₃
ExprPtr cumulant3_to_density(ExprPtr ex_);

/// replaces every cumulant κ_k in @p expr by its expansion in densities, then
/// expands and simplifies
/// @note a κ is recognized by its label alone, whatever its symmetries
ExprPtr cumulants_to_densities(ExprPtr expr);

ExprPtr one_body_sub(ExprPtr ex_);

ExprPtr two_body_decomp(ExprPtr ex_, bool approx = false);

// express 3-body term as sums of 1 and 2-body term. as described in J. Chem.
// Phys. 132, 234107 (2010); https://doi.org/10.1063/1.3439395 eqn 17.
std::pair<ExprPtr, std::pair<std::vector<Index>, std::vector<Index>>>
three_body_decomp(ExprPtr ex_, bool approx = true);

std::pair<ExprPtr, std::pair<std::vector<Index>, std::vector<Index>>>
three_body_decomposition(ExprPtr ex_, int rank, bool fast = false);

// in general a three body substitution can be approximated with 1, 2, or 3 body
// terms(3 body has no approximation). this is achieved by replacing densities
// with with particle number > rank by the each successive cumulant
// approximation followed by neglect of the particle rank sized term.
// TODO this implementation is ambitious and currently we only support rank 2
// decompositions.
//
// fast implementation represent non-constant solution interms of like terms and
// permutation operators.
ExprPtr three_body_substitution(ExprPtr& input, int rank, bool fast = false);

}  // namespace decompositions
}  // namespace mbpt
}  // namespace sequant

#endif  // SEQUANT_DOMAIN_MBPT_RDM_HPP
