//
// Created by Eduard Valeyev on 2019-04-01.
//

#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/rules/df.hpp>
#include <SeQuant/domain/mbpt/space_qns.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>

#include <cassert>
#include <cstdlib>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace sequant {
namespace mbpt {

void load(Convention conv, SpinConvention spconv) {
  std::shared_ptr<IndexSpaceRegistry> isr;
  switch (conv) {
    case Convention::Minimal:
      isr = make_min_sr_spaces(spconv);
      break;
    case Convention::SR:
      isr = make_sr_spaces(spconv);
      break;
    case Convention::MR:
      isr = make_mr_spaces(spconv);
      break;
    case Convention::F12:
      isr = make_F12_sr_spaces(spconv);
      break;
    case Convention::QCiFS:
      isr = make_legacy_spaces(spconv);
      break;
  }
  set_default_context({.index_space_registry_shared_ptr = isr,
                       .vacuum = Vacuum::SingleProduct});
}

void add_fermi_spin(IndexSpaceRegistry& isr) {
  IndexSpaceRegistry result = isr.clone();

  for (auto&& space : isr) {
    if (space.base_key() != L"") {
      IndexSpace spin_up(spinannotation_add(space.base_key(), Spin::alpha),
                         space.type(), Spin::alpha, space.approximate_size());
      IndexSpace spin_down(spinannotation_add(space.base_key(), Spin::beta),
                           space.type(), Spin::beta, space.approximate_size());
      result.add(spin_up);
      result.add(spin_down);
    }
  }
  const bool nulltype_ok = true;
  result.reference_occupied_space(isr.reference_occupied_space(nulltype_ok));
  result.vacuum_occupied_space(isr.vacuum_occupied_space(nulltype_ok));
  result.particle_space(isr.particle_space(nulltype_ok));
  result.hole_space(isr.hole_space(nulltype_ok));
  result.complete_space(isr.complete_space(nulltype_ok));

  isr = std::move(result);
}

void add_ao_spaces(std::shared_ptr<IndexSpaceRegistry>& isr, bool vbs,
                   bool abs) {
  // matches the MPQC layout, see spindex.h
  // this will not work for MR
  auto obs_lcao = isr->retrieve(vbs ? L"m" : L"p");
  isr->add(IndexSpace{L"μ", obs_lcao.type(), LCAOQNS::ao});  // OBS AO
  if (vbs) {
    auto vbs_lcao = isr->retrieve(L"e");
    isr->add(IndexSpace{L"Α", vbs_lcao.type(), LCAOQNS::ao})  // VBS AO
        .add_union(L"Γ", {L"μ", L"Α"});  // VBS+ = OBS + VBS
  }
  if (abs) {
    auto abs_lcao = isr->retrieve(L"α'");
    isr->add(IndexSpace{L"σ", abs_lcao.type(),
                        LCAOQNS::ao})  // ABS AO in F12 methods
        .add_union(L"ρ", {L"μ", L"σ"});
    if (vbs)                               // ABS+ = OBS + ABS
      isr->add_union(L"Ρ", {L"Γ", L"σ"});  // VABS+ = VBS+ + ABS
  }
}

void add_pao_spaces(std::shared_ptr<IndexSpaceRegistry>& isr) {
  auto uocc_space = isr->particle_space(/* nulltype_ok = */ false);
  isr->add(IndexSpace{L"μ̃", uocc_space, LCAOQNS::pao})  // OBS PAO
      ;
}

void add_df_spaces(std::shared_ptr<IndexSpaceRegistry>& isr) {
  // matches the MPQC layout, see spindex.h
  isr->add(IndexSpace{L"Κ", 0b00001, TensorFactorizationQNS::df})  // DFBS AO
      ;
}

std::shared_ptr<IndexSpaceRegistry> make_min_sr_spaces(SpinConvention spconv) {
  auto isr = std::make_shared<IndexSpaceRegistry>();

  const auto spin_any = IndexSpace::QuantumNumbers{
      spconv == SpinConvention::Legacy ? Spin::null : Spin::any};
  isr->add(L"i", 0b01, spin_any, is_vacuum_occupied, is_reference_occupied,
           is_hole)
      .add(L"a", 0b10, spin_any, is_particle)
      .add_union(L"p", {L"i", L"a"}, is_complete);
  if (spconv == SpinConvention::Default) add_fermi_spin(*isr);
  isr->physical_particle_attribute_mask(bitset_t(spin_any));

  return isr;
}

// Multireference supspace uses a subset of its occupied orbitals to define a
// vacuum occupied subspace.
//  this leaves an active space which is partially occupied/unoccupied. This
//  definition is convenient when coupled with SR vacuum.
std::shared_ptr<IndexSpaceRegistry> make_mr_spaces(SpinConvention spconv) {
  auto isr = std::make_shared<IndexSpaceRegistry>();

  const auto spin_any = IndexSpace::QuantumNumbers{
      spconv == SpinConvention::Legacy ? Spin::null : Spin::any};
  isr->add(L"o", 0b00001, spin_any)
      .add(L"i", 0b00010, spin_any)
      .add(L"u", 0b00100, spin_any)
      .add(L"a", 0b01000, spin_any)
      .add(L"g", 0b10000, spin_any)
      .add_union(L"O", {L"o", L"i"}, is_vacuum_occupied)
      .add_union(L"M", {L"o", L"i", L"u"}, is_reference_occupied)
      .add_union(L"I", {L"i", L"u"}, is_hole)
      .add_union(L"E", {L"u", L"a", L"g"})
      .add_union(L"A", {L"u", L"a"}, is_particle)
      .add_union(L"p", {L"M", L"E"}, is_complete);

  if (spconv == SpinConvention::Default) add_fermi_spin(*isr);
  isr->physical_particle_attribute_mask(bitset_t(spin_any));

  return isr;
}

std::shared_ptr<IndexSpaceRegistry> make_sr_spaces(SpinConvention spconv) {
  auto isr = std::make_shared<IndexSpaceRegistry>();

  const auto spin_any = IndexSpace::QuantumNumbers{
      spconv == SpinConvention::Legacy ? Spin::null : Spin::any};
  isr->add(L"o", 0b0001, spin_any)
      .add(L"i", 0b0010, spin_any, is_hole)
      .add(L"a", 0b0100, spin_any, is_particle)
      .add(L"g", 0b1000, spin_any)
      .add_union(L"m", {L"o", L"i"}, is_vacuum_occupied, is_reference_occupied)
      .add_union(L"e", {L"a", L"g"})
      .add_union(L"x", {L"i", L"a"})
      .add_union(L"p", {L"m", L"e"}, is_complete);
  if (spconv == SpinConvention::Default) add_fermi_spin(*isr);
  isr->physical_particle_attribute_mask(bitset_t(spin_any));

  return isr;
}

std::shared_ptr<IndexSpaceRegistry> make_F12_sr_spaces(SpinConvention spconv) {
  auto isr = std::make_shared<IndexSpaceRegistry>();

  const auto spin_any = IndexSpace::QuantumNumbers{
      spconv == SpinConvention::Legacy ? Spin::null : Spin::any};
  isr->add(L"o", 0b00001, spin_any)
      .add(L"i", 0b00010, spin_any, is_hole)
      .add(L"a", 0b00100, spin_any, is_particle)
      .add(L"g", 0b01000, spin_any)
      .add(L"α'", 0b10000, spin_any)
      .add_union(L"m", {L"o", L"i"}, is_vacuum_occupied, is_reference_occupied)
      .add_union(L"e", {L"a", L"g"})
      .add_union(L"x", {L"i", L"a"})
      .add_union(L"p", {L"m", L"e"})
      .add_unIon(L"h", {L"x", L"g"})
      .add_unIon(L"c", {L"g", L"α'"})
      .add_union(L"α", {L"e", L"α'"})
      .add_union(L"H", {L"i", L"α"})
      .add_union(L"κ", {L"p", L"α'"}, is_complete);
  if (spconv == SpinConvention::Default) add_fermi_spin(*isr);
  isr->physical_particle_attribute_mask(bitset_t(spin_any));

  return isr;
}

std::shared_ptr<IndexSpaceRegistry> make_legacy_spaces(SpinConvention spconv) {
  auto isr = std::make_shared<IndexSpaceRegistry>();

  const auto spin_any = IndexSpace::QuantumNumbers{
      spconv == SpinConvention::Legacy ? Spin::null : Spin::any};
  isr->add(L"o", 0b0000001, spin_any)
      .add(L"n", 0b0000010, spin_any)
      .add(L"i", 0b0000100, spin_any, is_hole)
      .add(L"u", 0b0001000, spin_any)
      .add(L"a", 0b0010000, spin_any, is_particle)
      .add(L"g", 0b0100000, spin_any)
      .add(L"α'", 0b1000000, spin_any)
      .add_union(L"m", {L"o", L"n", L"i"}, is_vacuum_occupied,
                 is_reference_occupied)
      .add_union(L"M", {L"m", L"u"})
      .add_union(L"e", {L"a", L"g"})
      .add_union(L"E", {L"u", L"e"})
      .add_union(L"x", {L"i", L"u", L"a"})
      .add_union(L"p", {L"m", L"x", L"e"})
      .add_union(L"κ", {L"p", L"α'"}, is_complete);

  if (spconv == SpinConvention::Default) add_fermi_spin(*isr);
  isr->physical_particle_attribute_mask(bitset_t(spin_any));

  return isr;
}

std::pair<std::shared_ptr<IndexSpaceRegistry>,
          std::shared_ptr<IndexSpaceRegistry>>
make_fermi_and_bose_spaces(SpinConvention spconv) {
  auto isr = std::make_shared<IndexSpaceRegistry>();

  const auto fspin_any = IndexSpace::QuantumNumbers{
      spconv == SpinConvention::Legacy ? Spin::null : Spin::any};
  isr->add(L"i", 0b001, fspin_any)    // fermi occupied
      .add(L"a", 0b010, fspin_any)    // fermi unoccupied
      .add_union(L"p", {L"i", L"a"})  // fermi all
      ;
  if (spconv == SpinConvention::Default) add_fermi_spin(*isr);
  const auto bspin_any = IndexSpace::QuantumNumbers{Spin::any};
  isr->add(L"β", 0b100, bspin_any);  // bose

  auto fermi_isr = std::make_shared<IndexSpaceRegistry>(isr->spaces());
  fermi_isr->vacuum_occupied_space(L"i");
  fermi_isr->reference_occupied_space(L"i");
  fermi_isr->hole_space(L"i");
  fermi_isr->particle_space(L"a");
  fermi_isr->complete_space(L"p");
  fermi_isr->physical_particle_attribute_mask(bitset_t(fspin_any));

  auto bose_isr = std::make_shared<IndexSpaceRegistry>(isr->spaces());
  bose_isr->vacuum_occupied_space(IndexSpace::null);
  bose_isr->reference_occupied_space(IndexSpace::null);
  bose_isr->hole_space(IndexSpace::null);
  bose_isr->particle_space(L"β");
  bose_isr->complete_space(L"β");
  bose_isr->physical_particle_attribute_mask(bitset_t(bspin_any));

  return std::make_pair(std::move(fermi_isr), std::move(bose_isr));
}

}  // namespace mbpt
}  // namespace sequant
