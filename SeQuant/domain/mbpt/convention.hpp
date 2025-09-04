//
// Created by Eduard Valeyev on 2019-04-01.
//

#ifndef SEQUANT_CONVENTION_HPP
#define SEQUANT_CONVENTION_HPP

#include <SeQuant/domain/mbpt/fwd.hpp>

#include <SeQuant/core/index_space_registry.hpp>

namespace sequant {
namespace mbpt {

/// @brief Conventions for partitioning the single-particle Hilbert space
enum class Convention {
  Minimal,  //!< occupied/hole + unoccupied/particle + their union
  SR,       //!< single determinant reference: occupied (frozen + active) +
            //!< unoccupied (active + frozen)
  MR,    //!< multi determinant reference: occupied (frozen + active) + active +
         //!< unoccupied (active + frozen)
  F12,   //!< SR + complement from complete basis, used for F12 methods
  QCiFS  //!< ``Quantum Chemistry in Fock Space'' = superset of above
};

/// @brief Conventions for representing spin quantum numbers
enum class SpinConvention {
  None,  //!< particles are assumed spin-free, spin bits are set to Spin::none
  Default,  //!< fermions are assumed spin-1/2 (Spin::any, Spin::up,
            //!< Spin::down), bosons are spin-free (Spin::none)
  Legacy,  //!< all particles are assumed spin-free, spin bits set to Spin::null
};

void load(Convention conv = Convention::Minimal,
          SpinConvention spconv = SpinConvention::Default);

/// @brief decorate IndexSpace labels with spin
std::wstring decorate_label(std::wstring label, bool up);

/// @brief add fermionic spin spaces to registry
void add_fermi_spin(IndexSpaceRegistry& isr);

/// @brief add AO spaces to registry

/// @param isr the IndexSpaceRegistry to which add the AO spaces
/// @param vbs if true, have separate virtual basis
void add_ao_spaces(std::shared_ptr<IndexSpaceRegistry>& isr, bool vbs = false,
                   bool abs = false);

/// @brief add DF spaces to registry
void add_df_spaces(std::shared_ptr<IndexSpaceRegistry>& isr);

/// @brief add PAO spaces to registry

/// expects \p isr to have a defined particle space
void add_pao_spaces(std::shared_ptr<IndexSpaceRegistry>& isr);

/// @brief add batching spaces to registry
void add_batching_spaces(std::shared_ptr<IndexSpaceRegistry>& isr);

/// @name built-in definitions of IndexSpace
/// @{

/// Most standard models only need 2 base spaces, occupied and unoccupied.
/// This is minimal partitioning sufficient for computing expectation values
/// in context of single-reference MBPT.
std::shared_ptr<IndexSpaceRegistry> make_min_sr_spaces(
    SpinConvention scv = SpinConvention::Default);

/// Common partitioning for single reference F12 calculations.
/// notably, this set contains an other_unoccupied space, α', commonly used to
/// construct an approximately complete representation
std::shared_ptr<IndexSpaceRegistry> make_F12_sr_spaces(
    SpinConvention spconv = SpinConvention::Default);

/// Multireference partitioning contains an active space, x, which is assumed to
/// have partial occupancy although it is considered unoccupied with respect to
/// a SingleProduct Vacuum. This leads to a variety of additional composite
/// spaces with may or may not be occupied.
std::shared_ptr<IndexSpaceRegistry> make_mr_spaces(
    SpinConvention spconv = SpinConvention::Default);

/// 'Standard' choice of partitioning orbitals in a single reference.
/// Includes frozen_core, active_occupied, active_unoccupied, and
/// inactive_unoccupied orbitals as base spaces.
std::shared_ptr<IndexSpaceRegistry> make_sr_spaces(
    SpinConvention spconv = SpinConvention::Default);

/// Legacy partitioning similar to previous versions of SeQuant which had
/// compile time hard coded partitioning. This is useful when verifying
/// previously obtained results which have been canonicalized in this context.
std::shared_ptr<IndexSpaceRegistry> make_legacy_spaces(
    SpinConvention spconv = SpinConvention::Default);

/// make fermi and bose space registries for multicomponent models
std::pair<std::shared_ptr<IndexSpaceRegistry>,
          std::shared_ptr<IndexSpaceRegistry>>
make_fermi_and_bose_spaces(SpinConvention spconv = SpinConvention::Default);

/// @}

}  // namespace mbpt
}  // namespace sequant

#endif  // SEQUANT_CONVENTION_HPP
