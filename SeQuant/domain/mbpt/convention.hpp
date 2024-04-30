//
// Created by Eduard Valeyev on 2019-04-01.
//

#ifndef SEQUANT_CONVENTION_HPP
#define SEQUANT_CONVENTION_HPP
#include "SeQuant/core/index_space_registry.h"

namespace sequant {
namespace mbpt {

enum class Convention { QCiFS };

// clang-format off
// most "simple" models only need reference to 2 base spaces, occupied and unoccupied.
// This is minimal partitioning is sufficient for computing expectation values of Coupled-Cluster type operators.
IndexSpaceRegistry make_minimal_single_reference_subspaces();

//A registry containing a common partitioning for single reference F12 calculations.
// notably, this set contains an other_unoccupied space, α', commonly used to construct an approximately complete representation
IndexSpaceRegistry make_F12_single_reference_subspaces();

//This multireference partitioning contains an active space, x, which is assumed to have partial density
// although it is considered unoccupied with respect to a SingleProduct Vacuum. This leads to a variety of additional
// composite spaces with may or may not be occupied.
IndexSpaceRegistry make_standard_multireference_subspaces();

// This paritioning is a somewhat 'standard' choice of partitioning orbitals in a single reference.
// contains frozen_core, active_occupied, active_unoccupied, and inactive_unoccupied orbitals as base spaces.
IndexSpaceRegistry make_standard_single_reference_subspaces();

// This is a legacy partitioning similar to previous versions of SeQuant which had compile time hard coded partitioning.
// This is useful when verifying previously obtained results which have been canonicalized in this context.
IndexSpaceRegistry make_standard_single_reference_subspaces_v1();
}  // namespace mbpt
}  // namespace sequant

#endif  // SEQUANT_CONVENTION_HPP
