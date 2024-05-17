//
// Created by Eduard Valeyev on 2019-04-01.
//

#include "SeQuant/domain/mbpt/convention.hpp"
#include <SeQuant/domain/mbpt/spin.hpp>
#include "SeQuant/domain/mbpt/op.hpp"

#include "SeQuant/core/context.hpp"
#include "SeQuant/core/index.hpp"
#include "SeQuant/core/tensor.hpp"

namespace sequant {
namespace mbpt {

void load(Convention conv) {
  std::shared_ptr<IndexSpaceRegistry> isr;
  switch (conv) {
    case Convention::Minimal:
      isr = make_min_sr_spaces();
      break;
    case Convention::SR:
      isr = make_sr_spaces();
      break;
    case Convention::MR:
      isr = make_mr_spaces();
      break;
    case Convention::F12:
      isr = make_F12_sr_spaces();
      break;
    case Convention::QCiFS:
      isr = make_legacy_spaces();
      break;
  }
  set_default_context(sequant::Context(isr, Vacuum::SingleProduct));
}

void add_fermi_spin(std::shared_ptr<IndexSpaceRegistry>& isr) {
  auto result = std::make_shared<IndexSpaceRegistry>();
  for (auto&& space : *isr) {
    if (space.base_key() != L"") {
      IndexSpace spin_any(space.base_key(), space.type(), Spin::any,
                          space.approximate_size());
      IndexSpace spin_up(spinannotation_add(space.base_key(), Spin::alpha),
                         space.type(), Spin::alpha, space.approximate_size());
      IndexSpace spin_down(spinannotation_add(space.base_key(), Spin::beta),
                           space.type(), Spin::beta, space.approximate_size());
      result->add(spin_any);
      result->add(spin_up);
      result->add(spin_down);
    }
  }
  const bool nulltype_ok = true;
  result->reference_occupied_space(isr->reference_occupied_space(nulltype_ok));
  result->vacuum_occupied_space(isr->vacuum_occupied_space(nulltype_ok));
  result->particle_space(isr->particle_space(nulltype_ok));
  result->hole_space(isr->hole_space(nulltype_ok));
  result->complete_space(isr->complete_space(nulltype_ok));
  isr = std::move(result);
}

std::shared_ptr<IndexSpaceRegistry> make_min_sr_spaces() {
  auto isr = std::make_shared<IndexSpaceRegistry>();

  isr->add(L"i", 0b01, is_vacuum_occupied, is_reference_occupied, is_hole)
      .add(L"a", 0b10, is_particle)
      .add_union(L"p", {L"i", L"a"}, is_complete);
  add_fermi_spin(isr);

  return isr;
}

// Multireference supspace uses a subset of its occupied orbitals to define a
// vacuum occupied subspace.
//  this leaves an active space which is partially occupied/unoccupied. This
//  definition is convenient when coupled with SR vacuum.
std::shared_ptr<IndexSpaceRegistry> make_mr_spaces() {
  auto isr = std::make_shared<IndexSpaceRegistry>();

  isr->add(L"o", 0b00001)
      .add(L"i", 0b00010)
      .add(L"u", 0b00100)
      .add(L"a", 0b01000)
      .add(L"g", 0b10000)
      .add_union(L"O", {L"o", L"i"}, is_vacuum_occupied)
      .add_union(L"M", {L"o", L"i", L"u"}, is_reference_occupied)
      .add_union(L"I", {L"i", L"u"}, is_hole)
      .add_union(L"E", {L"u", L"a", L"g"})
      .add_union(L"A", {L"u", L"a"}, is_particle)
      .add_union(L"p", {L"M", L"E"}, is_complete);

  add_fermi_spin(isr);

  return isr;
}

std::shared_ptr<IndexSpaceRegistry> make_sr_spaces() {
  auto isr = std::make_shared<IndexSpaceRegistry>();

  isr->add(L"o", 0b0001)
      .add(L"i", 0b0010, is_hole)
      .add(L"a", 0b0100, is_particle)
      .add(L"g", 0b1000)
      .add_union(L"m", {L"o", L"i"}, is_vacuum_occupied, is_reference_occupied)
      .add_union(L"e", {L"a", L"g"})
      .add_union(L"x", {L"i", L"a"})
      .add_union(L"p", {L"m", L"e"}, is_complete);
  add_fermi_spin(isr);

  return isr;
}

std::shared_ptr<IndexSpaceRegistry> make_F12_sr_spaces() {
  auto isr = std::make_shared<IndexSpaceRegistry>();

  isr->add(L"o", 0b00001)
      .add(L"i", 0b00010, is_hole)
      .add(L"a", 0b00100, is_particle)
      .add(L"g", 0b01000)
      .add(L"α'", 0b10000)
      .add_union(L"m", {L"o", L"i"}, is_vacuum_occupied, is_reference_occupied)
      .add_union(L"e", {L"a", L"g"})
      .add_union(L"x", {L"i", L"a"})
      .add_union(L"p", {L"m", L"e"})
      .add_union(L"α", {L"e", L"α'"})
      .add_union(L"κ", {L"p", L"α'"}, is_complete);
  add_fermi_spin(isr);

  return isr;
}

std::shared_ptr<IndexSpaceRegistry> make_legacy_spaces(bool ignore_spin) {
  auto isr = std::make_shared<IndexSpaceRegistry>();

  isr->add(L"o", 0b0000001)
      .add(L"n", 0b0000010)
      .add(L"i", 0b0000100, is_hole)
      .add(L"u", 0b0001000)
      .add(L"a", 0b0010000, is_particle)
      .add(L"g", 0b0100000)
      .add(L"α'", 0b1000000)
      .add_union(L"m", {L"o", L"n", L"i"}, is_vacuum_occupied,
                 is_reference_occupied)
      .add_union(L"M", {L"m", L"u"})
      .add_union(L"e", {L"a", L"g"})
      .add_union(L"E", {L"u", L"e"})
      .add_union(L"x", {L"i", L"u", L"a"})
      .add_union(L"p", {L"m", L"x", L"e"})
      .add_union(L"κ", {L"p", L"α'"}, is_complete);

  if (!ignore_spin) add_fermi_spin(isr);

  return isr;
}

std::pair<std::shared_ptr<IndexSpaceRegistry>,
          std::shared_ptr<IndexSpaceRegistry>>
make_fermi_and_bose_spaces() {
  auto isr = std::make_shared<IndexSpaceRegistry>();

  isr->add(L"i", 0b001)               // fermi occupied
      .add(L"a", 0b010)               // fermi unoccupied
      .add_union(L"p", {L"i", L"a"})  // fermi all
      ;
  add_fermi_spin(isr);
  isr->add(L"β", 0b100);  // bose

  auto fermi_isr = std::make_shared<IndexSpaceRegistry>(isr->spaces());
  fermi_isr->vacuum_occupied_space(L"i");
  fermi_isr->reference_occupied_space(L"i");
  fermi_isr->hole_space(L"i");
  fermi_isr->particle_space(L"a");
  fermi_isr->complete_space(L"p");

  auto bose_isr = std::make_shared<IndexSpaceRegistry>(isr->spaces());
  bose_isr->vacuum_occupied_space(IndexSpace::null);
  bose_isr->reference_occupied_space(IndexSpace::null);
  bose_isr->hole_space(IndexSpace::null);
  bose_isr->particle_space(L"β");
  bose_isr->complete_space(L"β");

  return std::make_pair(std::move(fermi_isr), std::move(bose_isr));
}

}  // namespace mbpt
}  // namespace sequant
