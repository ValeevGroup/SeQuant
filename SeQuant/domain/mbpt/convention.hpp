//
// Created by Eduard Valeyev on 2019-04-01.
//

#ifndef SEQUANT_CONVENTION_HPP
#define SEQUANT_CONVENTION_HPP

#include <SeQuant/domain/mbpt/fwd.hpp>

#include <SeQuant/core/basis.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/index_basis_registry.hpp>

#include <limits>
#include <string_view>

namespace sequant {
namespace mbpt {

/// @brief Conventions for partitioning the single-particle Hilbert space
enum class Convention {
  Minimal,  //!< occupied/hole + unoccupied/particle + their union
  SR,       //!< single determinant reference: occupied (frozen + active) +
            //!< unoccupied (active + frozen)
  MR,  //!< multi determinant reference: occupied (frozen + active) + active +
       //!< unoccupied (active + frozen)
  MinimalMR,  //!< MR, without frozen spaces
  F12,        //!< SR + complement from complete basis, used for F12 methods
  QCiFS       //!< ``Quantum Chemistry in Fock Space'' = superset of above
};

/// @brief Conventions for representing spin quantum numbers
enum class SpinConvention {
  None,  //!< particles are assumed spin-free, spin bits are set to Spin::none
  Default,  //!< fermions are assumed spin-1/2 (Spin::any, Spin::up,
            //!< Spin::down), bosons are spin-free (Spin::none)
  Legacy,  //!< all particles are assumed spin-free, spin bits set to Spin::null
};

/// @brief installs the index basis registry of a convention and the
///        single-product vacuum into the default Context
///
/// Every other setting of the current default Context, including the
/// canonicalizer configuration, is kept.
/// @throw Exception if the calling thread has an active scoped context (see
/// set_scoped_default_context()), whose settings would otherwise be
/// installed process-wide
void load(Convention conv = Convention::Minimal,
          SpinConvention spconv = SpinConvention::Default);

/// @brief decorate IndexSpace labels with spin
std::wstring decorate_label(std::wstring label, bool up);

/// @brief add fermionic spin spaces to registry
void add_fermi_spin(IndexBasisRegistry& isr);

/// the basis instance add_ao_basis registers the AO bases as: above every CSV
/// instance, below the PAO one, so an AO basis sorts after every other basis
/// of its space (IndexBasis orders by space, then instance)
inline constexpr IndexBasis::instance_type default_ao_basis_instance =
    std::numeric_limits<IndexBasis::instance_type>::max() - 1;

/// @brief registers the AO bases as named instances of the orbital spaces
/// they span, with a general metric (the AOs are not orthonormal)

/// The OBS AO basis `μ` spans the complete space (`p`, or `m` with \p vbs);
/// with \p vbs the VBS AO basis `Α` spans `e` and `Γ` = `μ` + `Α` spans their
/// union; with \p abs the ABS AO basis `σ` spans `α'`, `ρ` = `μ` + `σ`
/// spans the union of the OBS and `α'`, and with \p vbs also `Ρ` = `Γ` + `σ`
/// spans the union of all three. Every space a basis spans must be
/// registered. The extents default to the spaces' dimensions (populate them
/// with IndexBasisRegistry::extent(label, n) before the registry is
/// given to a Context, which holds it immutable). No spin-cased counterparts
/// (`μ↑`, `μ↓`, ...) are registered, unlike by add_pao_basis(): the AO bases
/// are the target of csv_transform() after spin tracing, and spin-casing an
/// index in one throws (see make_spinalpha()) unless they are registered by
/// hand.
/// @param isr the IndexBasisRegistry to which add the AO bases
/// @param spin_any the quantum numbers of the spin-agnostic orbital spaces
///        in the target convention (Spin::null for SpinConvention::Legacy,
///        else Spin::any)
/// @param vbs if true, have separate virtual basis
/// @param abs if true, have an auxiliary (F12) basis
/// @param instance the basis instance of the AO bases
/// @throw Exception if a label is already registered or a space an AO basis
///        spans is not; the registry is then left untouched
void add_ao_basis(
    std::shared_ptr<IndexBasisRegistry>& isr,
    IndexSpace::QuantumNumbers spin_any, bool vbs = false, bool abs = false,
    IndexBasis::instance_type instance = default_ao_basis_instance);

/// @deprecated the AO bases are named instances of the orbital spaces, see
/// add_ao_basis(), to which this forwards; unlike the AO spaces this used to
/// register, the bases need every space they span registered (with \p vbs
/// and \p abs together, also the union of the OBS and the ABS)
[[deprecated("the AO bases are named basis instances; use add_ao_basis")]] void
add_ao_spaces(std::shared_ptr<IndexBasisRegistry>& isr,
              IndexSpace::QuantumNumbers spin_any, bool vbs = false,
              bool abs = false);

/// @brief add DF spaces to registry
void add_df_spaces(std::shared_ptr<IndexBasisRegistry>& isr);

/// @brief add THC spaces to registry
void add_thc_spaces(std::shared_ptr<IndexBasisRegistry>& isr);

/// @deprecated the PAO basis is a named instance of the particle space, see
/// add_pao_basis(), to which this forwards; like it, this registers the OBS
/// AO basis `μ` if \p isr has none, which add_ao_spaces() cannot then
/// register, so call that one first
[[deprecated(
    "the PAO basis is a named basis instance; use add_pao_basis")]] void
add_pao_spaces(std::shared_ptr<IndexBasisRegistry>& isr,
               IndexSpace::QuantumNumbers spin_any);

/// the basis instance add_pao_basis registers by default: above every CSV
/// instance, so a PAO basis sorts after every other basis of the particle
/// space (IndexBasis orders by space, then instance)
inline constexpr IndexBasis::instance_type default_pao_basis_instance =
    std::numeric_limits<IndexBasis::instance_type>::max();

/// @brief registers the PAO basis as a named instance of the particle space

/// expects \p isr to have a defined particle space. The PAOs are the AOs
/// projected on the particle space, so the entry follows the OBS AO basis
/// `μ` (IndexBasisRegistry::follow()), which add_ao_basis() registers first
/// if \p isr has no `μ`, as the OBS AO basis alone; for the VBS or ABS AO
/// bases call add_ao_basis() first, since it cannot register `μ` twice. The
/// PAO bases' extent, metric (general) and field are those
/// of `μ`, set through `μ`'s label (populate the extent with
/// IndexBasisRegistry::extent(L"μ", n) before the registry is given to a
/// Context, which holds it immutable). The α- and β-spin PAO bases are
/// registered alongside, as the same instance of the spin-cased particle
/// spaces (if \p isr has them) under the spin-annotated label (`μ̃↑`,
/// `μ̃↓`), so that a PAO index can be spin-cased (see make_spinalpha());
/// they follow `μ` too, since the AOs do not depend on spin.
/// @param spin_any the quantum numbers of the spin-agnostic particle space
/// @param instance the basis instance of the PAO basis
/// @param label the label the PAO basis is registered under
/// @throw Exception if \p label is `μ` or it or a spin-annotated version of
///        it is already registered, if \p instance of a particle space is
///        already named, if `μ` is missing and cannot be registered (see
///        add_ao_basis()) or is registered but is a space or follows an
///        entry itself; the registry is then left untouched
void add_pao_basis(
    std::shared_ptr<IndexBasisRegistry>& isr,
    IndexSpace::QuantumNumbers spin_any,
    IndexBasis::instance_type instance = default_pao_basis_instance,
    std::wstring_view label = L"μ̃");

/// @brief add batching spaces to registry
void add_batching_spaces(std::shared_ptr<IndexBasisRegistry>& isr);

/// @name built-in definitions of IndexSpace
/// @{

/// Most standard models only need 2 base spaces, occupied and unoccupied.
/// This is minimal partitioning sufficient for computing expectation values
/// in context of single-reference MBPT.
std::shared_ptr<IndexBasisRegistry> make_min_sr_spaces(
    SpinConvention scv = SpinConvention::Default);

/// Common partitioning for single reference F12 calculations.
/// notably, this set contains an other_unoccupied space, α', commonly used to
/// construct an approximately complete representation
std::shared_ptr<IndexBasisRegistry> make_F12_sr_spaces(
    SpinConvention spconv = SpinConvention::Default);

/// Multireference partitioning contains an active space, x, which is assumed to
/// have partial occupancy although it is considered unoccupied with respect to
/// a SingleProduct Vacuum. This leads to a variety of additional composite
/// spaces with may or may not be occupied.
std::shared_ptr<IndexBasisRegistry> make_mr_spaces(
    SpinConvention spconv = SpinConvention::Default);

/// like make_mr_spaces, but without frozen orbitals
std::shared_ptr<IndexBasisRegistry> make_min_mr_spaces(
    SpinConvention spconv = SpinConvention::Default);

/// 'Standard' choice of partitioning orbitals in a single reference.
/// Includes frozen_core, active_occupied, active_unoccupied, and
/// inactive_unoccupied orbitals as base spaces.
std::shared_ptr<IndexBasisRegistry> make_sr_spaces(
    SpinConvention spconv = SpinConvention::Default);

/// Legacy partitioning similar to previous versions of SeQuant which had
/// compile time hard coded partitioning. This is useful when verifying
/// previously obtained results which have been canonicalized in this context.
std::shared_ptr<IndexBasisRegistry> make_legacy_spaces(
    SpinConvention spconv = SpinConvention::Default);

/// make fermi and bose space registries for multicomponent models
std::pair<std::shared_ptr<IndexBasisRegistry>,
          std::shared_ptr<IndexBasisRegistry>>
make_fermi_and_bose_spaces(SpinConvention spconv = SpinConvention::Default);

/// @}

/// @brief Checks whether ISR has batching space labelled by "z"
/// @throws Assertion failure if batching space "z" is not found in the registry
inline void check_for_batching_space() {
  SEQUANT_ASSERT(
      sequant::get_default_context().index_basis_registry()->contains(L"z"));
}

}  // namespace mbpt
}  // namespace sequant

#endif  // SEQUANT_CONVENTION_HPP
