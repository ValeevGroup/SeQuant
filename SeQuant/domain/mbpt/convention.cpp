//
// Created by Eduard Valeyev on 2019-04-01.
//

#include "convention.hpp"
#include "SeQuant/core/context.hpp"
#include "SeQuant/core/index.hpp"
#include "SeQuant/core/tensor.hpp"
#include "op.hpp"

#include <SeQuant/domain/mbpt/spin.hpp>

namespace sequant {
namespace mbpt {

void add_fermi_spin(IndexSpaceRegistry& isr) {
  IndexSpaceRegistry result;
  for (auto&& space : isr) {
    if (space.base_key() != L"") {
      IndexSpace spin_any(space.base_key(), space.type(), Spin::any,
                          space.approximate_size());
      IndexSpace spin_up(spinannotation_add(space.base_key(), Spin::alpha),
                         space.type(), Spin::alpha, space.approximate_size());
      IndexSpace spin_down(spinannotation_add(space.base_key(), Spin::beta),
                           space.type(), Spin::beta, space.approximate_size());
      result.add(spin_any);
      result.add(spin_up);
      result.add(spin_down);
    }
  }
  result.reference_occupied_space(isr.reference_occupied_space());
  result.vacuum_occupied_space(isr.vacuum_occupied_space());
  result.active_particle_space(isr.active_particle_space());
  result.active_hole_space(isr.active_hole_space());
  result.complete_space(isr.complete_space());
  isr = std::move(result);
}

IndexSpaceRegistry make_min_sr_subspaces() {
  IndexSpaceRegistry minimal_reference_registry;

  IndexSpace occupied(L"i", 0b01);
  minimal_reference_registry.add(occupied);

  IndexSpace unoccupied(L"a", 0b10);
  minimal_reference_registry.add(unoccupied);

  IndexSpace all(L"p", 0b11);
  minimal_reference_registry.add(all);

  minimal_reference_registry.reference_occupied_space(occupied);
  minimal_reference_registry.complete_space(all);
  minimal_reference_registry.vacuum_occupied_space(occupied);
  minimal_reference_registry.active_particle_space(occupied);
  minimal_reference_registry.active_hole_space(unoccupied);

  add_fermi_spin(minimal_reference_registry);

  return minimal_reference_registry;
}
IndexSpaceRegistry make_F12_sr_subspaces() {
  IndexSpaceRegistry standard_reference_registry;

  // add irreducible representations first
  IndexSpace frozen(L"o", 0b00001);
  standard_reference_registry.add(frozen);

  IndexSpace active_occ(L"i", 0b00010);
  standard_reference_registry.add(active_occ);

  IndexSpace active_uocc(L"a", 0b00100);
  standard_reference_registry.add(active_uocc);

  IndexSpace inactive_uocc(L"g", 0b01000);
  standard_reference_registry.add(inactive_uocc);

  IndexSpace other_uocc(L"α'", 0b10000);
  standard_reference_registry.add(other_uocc);

  IndexSpace occupied(L"m", frozen.type().unIon(active_occ.attr()).to_int32());
  standard_reference_registry.add(occupied);

  IndexSpace unoccupied(
      L"e", active_uocc.type().unIon(inactive_uocc.type()).to_int32());
  standard_reference_registry.add(unoccupied);

  IndexSpace active(L"x",
                    active_occ.type().unIon(active_uocc.type()).to_int32());
  standard_reference_registry.add(active);

  IndexSpace complete_unoccupied(
      L"α", unoccupied.type().unIon(other_uocc.type()).to_int32());
  standard_reference_registry.add(complete_unoccupied);

  IndexSpace all(L"p", unoccupied.type().unIon(occupied.type()).to_int32());
  standard_reference_registry.add(all);

  IndexSpace complete(
      L"κ", complete_unoccupied.type().unIon(occupied.type()).to_int32());
  standard_reference_registry.add(complete);

  // only necessary for some Operator definitions
  standard_reference_registry.reference_occupied_space(occupied);
  standard_reference_registry.complete_space(complete);
  standard_reference_registry.active_particle_space(active_occ);
  standard_reference_registry.active_hole_space(active_uocc);
  // necessary for SR wick algebra
  standard_reference_registry.vacuum_occupied_space(occupied);

  add_fermi_spin(standard_reference_registry);

  return standard_reference_registry;
}

// Multireference supspace uses a subset of its occupied orbitals to define a
// vacuum occupied subspace.
//  this leaves an active space which is partially occupied/unoccupied. This
//  definition is convenient when coupled with SR vacuum.
IndexSpaceRegistry make_mr_subspaces() {
  IndexSpaceRegistry multireference_registry;

  IndexSpace frozen(L"o", 0b00001);
  multireference_registry.add(frozen);

  IndexSpace unfrozen_vacuum_occupied(L"i", 0b00010);
  multireference_registry.add(unfrozen_vacuum_occupied);

  IndexSpace active(L"x", 0b00100);
  multireference_registry.add(active);

  IndexSpace active_unoccupied(L"a", 0b01000);
  multireference_registry.add(active_unoccupied);

  IndexSpace inactive_unoccupied(L"g", 0b10000);
  multireference_registry.add(inactive_unoccupied);

  IndexSpace vacuum_occupied(
      L"O", frozen.type().unIon(unfrozen_vacuum_occupied.type()));
  multireference_registry.add(vacuum_occupied);

  IndexSpace maybe_occupied(L"M", active.type()
                                      .unIon(unfrozen_vacuum_occupied.type())
                                      .unIon(frozen.type()));
  multireference_registry.add(maybe_occupied);

  IndexSpace maybe_active_occupied(
      L"I", active.type().unIon(unfrozen_vacuum_occupied.type()));
  multireference_registry.add(maybe_active_occupied);

  IndexSpace maybe_unoccupied(L"E", active.type()
                                        .unIon(active_unoccupied.type())
                                        .unIon(inactive_unoccupied.type()));
  multireference_registry.add(maybe_unoccupied);

  IndexSpace maybe_active_unoccupied(
      L"A", active.type().unIon(active_unoccupied.type()));
  multireference_registry.add(maybe_active_unoccupied);

  IndexSpace complete(L"p",
                      vacuum_occupied.type().unIon(maybe_unoccupied.type()));
  multireference_registry.add(complete);

  // only necessary for some Operator definitions
  multireference_registry.complete_space(complete);
  multireference_registry.reference_occupied_space(maybe_occupied);
  multireference_registry.active_hole_space(maybe_active_unoccupied);
  multireference_registry.active_particle_space(maybe_active_occupied);
  // needed for SR wick algebra
  multireference_registry.vacuum_occupied_space(vacuum_occupied);

  add_fermi_spin(multireference_registry);

  return multireference_registry;
}

IndexSpaceRegistry make_sr_subspaces() {
  IndexSpaceRegistry standard_reference_registry;
  IndexSpace frozen(L"o", 0b0001);
  standard_reference_registry.add(frozen);

  IndexSpace active_occ(L"i", 0b0010);
  standard_reference_registry.add(active_occ);

  IndexSpace active_uocc(L"a", 0b0100);
  standard_reference_registry.add(active_uocc);

  IndexSpace inactive_uocc(L"g", 0b1000);
  standard_reference_registry.add(inactive_uocc);

  IndexSpace occupied(L"m", frozen.type().unIon(active_occ.attr()).to_int32());
  standard_reference_registry.add(occupied);

  IndexSpace unoccupied(
      L"e", active_uocc.type().unIon(inactive_uocc.type()).to_int32());
  standard_reference_registry.add(unoccupied);

  IndexSpace complete(L"p",
                      unoccupied.type().unIon(occupied.type()).to_int32());
  standard_reference_registry.add(complete);

  // only necessary for some Operator definitions
  standard_reference_registry.reference_occupied_space(occupied);
  standard_reference_registry.complete_space(complete);
  standard_reference_registry.active_particle_space(active_occ);
  standard_reference_registry.active_hole_space(active_uocc);
  // necessary for SR wick algebra
  standard_reference_registry.vacuum_occupied_space(occupied);

  add_fermi_spin(standard_reference_registry);

  return standard_reference_registry;
}

IndexSpaceRegistry make_legacy_subspaces(bool ignore_spin) {
  IndexSpaceRegistry standard_reference_registry;
  IndexSpace frozen(L"o", 0b0000001);
  standard_reference_registry.add(frozen);

  IndexSpace active_occ(L"i", 0b0000100);
  standard_reference_registry.add(active_occ);

  IndexSpace inactive_occupied(L"n", 0b0000010);
  standard_reference_registry.add(inactive_occupied);

  IndexSpace active{L"u", 0b0001000};
  standard_reference_registry.add(active);

  IndexSpace active_uocc(L"a", 0b0010000);
  standard_reference_registry.add(active_uocc);

  IndexSpace inactive_uocc(L"g", 0b0100000);
  standard_reference_registry.add(inactive_uocc);

  IndexSpace occupied(L"m", 0b0000111);
  standard_reference_registry.add(occupied);

  IndexSpace unoccupied(L"e", 0b0110000);
  standard_reference_registry.add(unoccupied);

  IndexSpace other_unoccupied(L"α'", 0b1000000);
  standard_reference_registry.add(other_unoccupied);

  IndexSpace all_active(
      L"x", active_occ.type().unIon(active.type()).unIon(active_uocc.type()));
  standard_reference_registry.add(all_active);

  IndexSpace all(L"p", unoccupied.type().unIon(occupied.type()));
  standard_reference_registry.add(all);

  IndexSpace complete(L"κ", 0b1111111);
  standard_reference_registry.add(complete);

  IndexSpace maybe_occupied(L"M", active.type().unIon(occupied.type()));
  standard_reference_registry.add(maybe_occupied);

  IndexSpace maybe_unoccupied(L"E", active.type().unIon(unoccupied.type()));
  standard_reference_registry.add(maybe_unoccupied);

  // only necessary for some Operator definitions
  standard_reference_registry.reference_occupied_space(occupied);
  standard_reference_registry.complete_space(complete);
  standard_reference_registry.active_particle_space(active_occ);
  standard_reference_registry.active_hole_space(active_uocc);
  // necessary for SR wick algebra
  standard_reference_registry.vacuum_occupied_space(occupied);

  if (!ignore_spin) add_fermi_spin(standard_reference_registry);

  return standard_reference_registry;
}

}  // namespace mbpt
}  // namespace sequant
