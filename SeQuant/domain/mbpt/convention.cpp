//
// Created by Eduard Valeyev on 2019-04-01.
//

#include "convention.hpp"
#include "SeQuant/core/index.hpp"
#include "SeQuant/core/tensor.hpp"
#include "SeQuant/core/context.hpp"
#include "op.hpp"

namespace sequant {
namespace mbpt {


namespace {
enum class qns { none = 0, alpha = 1, beta = 2 };
auto qndecorate(qns qn, std::wstring_view label) {
  switch (static_cast<int>(qn)) {
    case 0:
      return std::wstring(label);
    case 1:
      return std::wstring(label) + L"↑";
    case 2:
      return std::wstring(label) + L"↓";
    default:
      assert(false && "invalid quantum number");
  }
  abort();  // unreachable
};
};  // namespace
IndexSpaceRegistry make_F12_single_reference_subspaces(){
  IndexSpaceRegistry standard_reference_registry;

  // add irreducible representations first
  IndexSpace frozen(L"o",{0b00001});
  standard_reference_registry.add(frozen);

  IndexSpace active_occ(L"i", {0b00010});
  standard_reference_registry.add(active_occ);

  IndexSpace active_uocc(L"a", {0b00100});
  standard_reference_registry.add(active_uocc);

  IndexSpace inactive_uocc(L"g", {0b01000});
  standard_reference_registry.add(inactive_uocc);

  IndexSpace other_uocc(L"α'", {0b10000});
  standard_reference_registry.add(other_uocc);

  IndexSpace occupied(L"m", frozen.type().unIon(active_occ.attr()).to_int32());
  standard_reference_registry.add(occupied);

  IndexSpace unoccupied(L"e",active_uocc.type().unIon(inactive_uocc.type()).to_int32());
  standard_reference_registry.add(unoccupied);

  IndexSpace complete_unoccupied(L"α",unoccupied.type().unIon(other_uocc.type()).to_int32());
  standard_reference_registry.add(complete_unoccupied);

  IndexSpace all(L"p",unoccupied.type().unIon(occupied.type()).to_int32());
  standard_reference_registry.add(all);

  IndexSpace complete(L"κ",complete_unoccupied.type().unIon(occupied.type()).to_int32());
  standard_reference_registry.add(complete);

  //only neccessary for some Operator definitions
  standard_reference_registry.assign_density_occupied(L"m");
  standard_reference_registry.assign_complete(L"κ");
  standard_reference_registry.assign_active_particle_space(L"i");
  standard_reference_registry.assign_active_hole_space(L"a");
  // neccessary for SR wick algebra
  standard_reference_registry.assign_vacuum_occupied(L"m");

  return standard_reference_registry;
}

//Multireference supspace uses a subset of its occupied orbitals to define a vacuum occupied subspace.
// this leaves an active space which is partially occupied/unoccupied. This definition is convenient when coupled with SR vacuum.
IndexSpaceRegistry make_standard_multireference_subspaces(){
  IndexSpaceRegistry multireference_registry;

  IndexSpace frozen(L"o",{0b00001});
  multireference_registry.add(frozen);

  IndexSpace unfrozen_vacuum_occupied(L"i",{0b00010});
  multireference_registry.add(unfrozen_vacuum_occupied);

  IndexSpace active(L"x",{0b00100});
  multireference_registry.add(active);

  IndexSpace active_unoccupied(L"a",{0b01000});
  multireference_registry.add(active_unoccupied);

  IndexSpace inactive_unoccupied(L"g",{0b10000});
  multireference_registry.add(inactive_unoccupied);

  IndexSpace vacuum_occupied(L"O", frozen.type().unIon(unfrozen_vacuum_occupied.type()));
  multireference_registry.add(vacuum_occupied);

  IndexSpace maybe_occupied(L"M",active.type().unIon(unfrozen_vacuum_occupied.type()).unIon(frozen.type()));
  multireference_registry.add(maybe_occupied);

  IndexSpace maybe_active_occupied(L"I", active.type().unIon(unfrozen_vacuum_occupied.type()));
  multireference_registry.add(maybe_active_occupied);

  IndexSpace maybe_unoccupied(L"E", active.type().unIon(inactive_unoccupied.type()).unIon(inactive_unoccupied.type()));
  multireference_registry.add(maybe_unoccupied);

  IndexSpace maybe_active_unoccupied(L"A", active.type().unIon(active_unoccupied.type()));
  multireference_registry.add(maybe_active_unoccupied);

  IndexSpace complete(L"p",vacuum_occupied.type().unIon(maybe_unoccupied.type()));
  multireference_registry.add(complete);

  //only neccessary for some Operator definitions
  multireference_registry.assign_complete(L"p");
  multireference_registry.assign_density_occupied(L"M");
  multireference_registry.assign_active_hole_space(L"A");
  multireference_registry.assign_active_particle_space(L"I");
  // needed for SR wick algebra
  multireference_registry.assign_vacuum_occupied(L"O");

  return multireference_registry;
}

IndexSpaceRegistry make_standard_single_reference_subspaces(){
  IndexSpaceRegistry standard_reference_registry;
  IndexSpace frozen(L"o",{0b0001});
  standard_reference_registry.add(frozen);

  IndexSpace active_occ(L"i", {0b0010});
  standard_reference_registry.add(active_occ);

  IndexSpace active_uocc(L"a", {0b0100});
  standard_reference_registry.add(active_uocc);

  IndexSpace inactive_uocc(L"g", {0b1000});
  standard_reference_registry.add(inactive_uocc);

  IndexSpace occupied(L"m", frozen.type().unIon(active_occ.attr()).to_int32());
  standard_reference_registry.add(occupied);

  IndexSpace unoccupied(L"e",active_uocc.type().unIon(inactive_uocc.type()).to_int32());
  standard_reference_registry.add(unoccupied);

  IndexSpace complete(L"p",unoccupied.type().unIon(occupied.type()).to_int32());
  standard_reference_registry.add(complete);

  //only neccessary for some Operator definitions
  standard_reference_registry.assign_density_occupied(L"m");
  standard_reference_registry.assign_complete(L"κ");
  standard_reference_registry.assign_active_particle_space(L"i");
  standard_reference_registry.assign_active_hole_space(L"a");
  // neccessary for SR wick algebra
  standard_reference_registry.assign_vacuum_occupied(L"m");

  return standard_reference_registry;
}



}  // namespace mbpt
}  // namespace sequant
