//
// Created by Eduard Valeyev on 2019-04-01.
//

#ifndef SEQUANT_CONVENTION_HPP
#define SEQUANT_CONVENTION_HPP

#include "SeQuant/core/index_space_registry.hpp"

namespace sequant {
namespace mbpt {

/// @brief Conventions for partitioning the single-particle Hilbert space
enum class Convention {
  Minimal,  //!< occupied/hole + unoccupied/unoccupied + their union
  SR,       //!< single determinant reference: occupied (frozen + active) +
            //!< unoccupied (active + frozen)
  MR,    //!< multi determinant reference: occupied (frozen + active) + active +
         //!< unoccupied (active + frozen)
  F12,   //!< SR + complement from complete basis, used for F12 methods
  QCiFS  //!< ``Quantum Chemistry in Fock Space'' = superset of above
};

void load(Convention conv = Convention::Minimal);

///@brief decorate IndexSpace labels with spin
std::wstring decorate_label(std::wstring label, bool up);

///@brief add fermionic spin spaces to registry
void add_fermi_spin(IndexSpaceRegistry &isr);

/// @name built-in definitions of IndexSpace
/// @{

/// Most standard models only need 2 base spaces, occupied and unoccupied.
/// This is minimal partitioning is sufficient for computing expectation values
/// of Coupled-Cluster type operators.
IndexSpaceRegistry make_min_sr_subspaces();

/// Common partitioning for single reference F12 calculations.
/// notably, this set contains an other_unoccupied space, α', commonly used to
/// construct an approximately complete representation
IndexSpaceRegistry make_F12_sr_subspaces();

/// Multireference partitioning contains an active space, x, which is assumed to
/// have partial density although it is considered unoccupied with respect to a
/// SingleProduct Vacuum. This leads to a variety of additional composite spaces
/// with may or may not be occupied.
IndexSpaceRegistry make_mr_subspaces();

/// 'Standard' choice of partitioning orbitals in a single reference.
/// Includes frozen_core, active_occupied, active_unoccupied, and
/// inactive_unoccupied orbitals as base spaces.
IndexSpaceRegistry make_sr_subspaces();

/// Legacy partitioning similar to previous versions of SeQuant which had
/// compile time hard coded partitioning. This is useful when verifying
/// previously obtained results which have been canonicalized in this context.
/// @param ignore_spin if true, do not add spin-specific spaces, and do not use
/// Spin::any as IndexSpace::QuantumNumbers for spin-free space
IndexSpaceRegistry make_legacy_subspaces(bool ignore_spin = false);

/// @}

}  // namespace mbpt
}  // namespace sequant

#endif  // SEQUANT_CONVENTION_HPP
