//
// Created by Conner Masteran on 7/1/21.
//

#ifndef SEQUANT_DOMAIN_MBPT_RDM_HPP
#define SEQUANT_DOMAIN_MBPT_RDM_HPP

#include <SeQuant/domain/mbpt/antisymmetrizer.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

namespace sequant {
namespace mbpt {
namespace decompositions {

ExprPtr cumu_to_density(ExprPtr ex_);

ExprPtr cumu2_to_density(ExprPtr ex_);

ExprPtr cumu3_to_density(ExprPtr ex_);

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
