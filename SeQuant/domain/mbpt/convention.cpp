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
#include <SeQuant/core/utility/context.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <cassert>
#include <cstdlib>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sequant {
namespace mbpt {

void load(Convention conv, SpinConvention spconv) {
  if (sequant::detail::implicit_context_overlay<
          container::map<Statistics, sequant::Context>>())
    throw Exception(
        "mbpt::load: cannot install a process-wide context while a scoped "
        "context is active on this thread");
  std::shared_ptr<IndexBasisRegistry> isr;
  switch (conv) {
    case Convention::Minimal:
      isr = make_min_sr_spaces(spconv);
      break;
    case Convention::SR:
      isr = make_sr_spaces(spconv);
      break;
    case Convention::MinimalMR:
      isr = make_min_mr_spaces(spconv);
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
  sequant::Context ctx = get_default_context_snapshot();
  ctx.set(std::move(isr));
  ctx.set(Vacuum::SingleProduct);
  set_default_context(std::move(ctx));
}

void add_fermi_spin(IndexBasisRegistry& isr) {
  IndexBasisRegistry result = isr;

  for (auto&& space : isr) {
    if (space.base_key() != L"") {
      IndexSpace spin_up(spinannotation_add(space.base_key(), Spin::alpha),
                         space.type(), Spin::alpha, space.dimension());
      IndexSpace spin_down(spinannotation_add(space.base_key(), Spin::beta),
                           space.type(), Spin::beta, space.dimension());
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

void add_ao_basis(std::shared_ptr<IndexBasisRegistry>& isr,
                  IndexSpace::QuantumNumbers spin_any, bool vbs, bool abs,
                  IndexBasis::instance_type instance) {
  // matches the MPQC layout, see spindex.h
  // this will not work for MR
  // the AOs of a basis set span an orbital space; the union of two AO bases
  // spans the union of their spaces, which must be registered
  auto add = [&](std::wstring_view label, const IndexSpace& space) {
    isr->add(label, IndexBasis{space, instance}, IndexSpaceMetric::General);
  };
  auto union_space = [&](std::wstring_view label, const IndexSpace& s1,
                         const IndexSpace& s2) -> const IndexSpace& {
    const auto* space = isr->retrieve_ptr(s1.type() | s2.type(), spin_any);
    if (!space)
      throw Exception("add_ao_basis: the AO basis '" + toUtf8(label) +
                      "' spans the union of the spaces " +
                      toUtf8(s1.base_key()) + " and " + toUtf8(s2.base_key()) +
                      ", which is not registered");
    return *space;
  };
  const auto& obs = isr->retrieve(vbs ? L"m" : L"p");
  add(L"μ", obs);  // OBS AO
  if (vbs) {
    const auto& vbs_space = isr->retrieve(L"e");
    add(L"Α", vbs_space);                          // VBS AO
    add(L"Γ", union_space(L"Γ", obs, vbs_space));  // VBS+ = OBS + VBS
  }
  if (abs) {
    const auto& abs_space = isr->retrieve(L"α'");
    add(L"σ", abs_space);                          // ABS AO in F12 methods
    add(L"ρ", union_space(L"ρ", obs, abs_space));  // ABS+ = OBS + ABS
    if (vbs)                                       // VABS+ = VBS+ + ABS
      add(L"Ρ", union_space(L"Ρ", union_space(L"Γ", obs, isr->retrieve(L"e")),
                            abs_space));
  }
}

void add_ao_spaces(std::shared_ptr<IndexBasisRegistry>& isr,
                   IndexSpace::QuantumNumbers spin_any, bool vbs, bool abs) {
  add_ao_basis(isr, spin_any, vbs, abs);
}

void add_pao_spaces(std::shared_ptr<IndexBasisRegistry>& isr,
                    IndexSpace::QuantumNumbers spin_any) {
  add_pao_basis(isr, spin_any);
}

void add_pao_basis(std::shared_ptr<IndexBasisRegistry>& isr,
                   IndexSpace::QuantumNumbers spin_any,
                   IndexBasis::instance_type instance,
                   std::wstring_view label) {
  const auto& uocc =
      isr->retrieve(isr->particle_space(/* nulltype_ok = */ false), spin_any);
  isr->add(label, IndexBasis{uocc, instance}, IndexSpaceMetric::General);
}

void add_df_spaces(std::shared_ptr<IndexBasisRegistry>& isr) {
  // matches the MPQC layout, see spindex.h
  isr->add(IndexSpace{L"Κ", 0b00001, TensorFactorizationQNS::df});  // DFBS AO
}

void add_batching_spaces(std::shared_ptr<IndexBasisRegistry>& isr) {
  isr->add(IndexSpace{L"z", 0b100000, BatchingQNS::batch});  // Batching Space
}

void add_thc_spaces(std::shared_ptr<IndexBasisRegistry>& isr) {
  isr->add(IndexSpace{L"L", 0b000001, TensorFactorizationQNS::thc})  // THC AO
      ;
}

std::shared_ptr<IndexBasisRegistry> make_min_sr_spaces(SpinConvention spconv) {
  auto isr = std::make_shared<IndexBasisRegistry>();

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
std::shared_ptr<IndexBasisRegistry> make_mr_spaces(SpinConvention spconv) {
  auto isr = std::make_shared<IndexBasisRegistry>();

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

std::shared_ptr<IndexBasisRegistry> make_min_mr_spaces(SpinConvention spconv) {
  auto isr = std::make_shared<IndexBasisRegistry>();

  const auto spin_any = IndexSpace::QuantumNumbers{
      spconv == SpinConvention::Legacy ? Spin::null : Spin::any};
  isr->add(L"i", 0b0001, spin_any, is_vacuum_occupied)
      .add(L"u", 0b0010, spin_any)
      .add(L"a", 0b0100, spin_any)
      .add_union(L"I", {L"i", L"u"}, is_reference_occupied, is_hole)
      .add_union(L"A", {L"u", L"a"}, is_particle)
      .add_union(L"p", {L"I", L"a"}, is_complete);

  if (spconv == SpinConvention::Default) add_fermi_spin(*isr);
  isr->physical_particle_attribute_mask(bitset_t(spin_any));

  return isr;
}

std::shared_ptr<IndexBasisRegistry> make_sr_spaces(SpinConvention spconv) {
  auto isr = std::make_shared<IndexBasisRegistry>();

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

std::shared_ptr<IndexBasisRegistry> make_F12_sr_spaces(SpinConvention spconv) {
  auto isr = std::make_shared<IndexBasisRegistry>();

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

std::shared_ptr<IndexBasisRegistry> make_legacy_spaces(SpinConvention spconv) {
  auto isr = std::make_shared<IndexBasisRegistry>();

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

std::pair<std::shared_ptr<IndexBasisRegistry>,
          std::shared_ptr<IndexBasisRegistry>>
make_fermi_and_bose_spaces(SpinConvention spconv) {
  auto isr = std::make_shared<IndexBasisRegistry>();

  const auto fspin_any = IndexSpace::QuantumNumbers{
      spconv == SpinConvention::Legacy ? Spin::null : Spin::any};
  isr->add(L"i", 0b001, fspin_any)    // fermi occupied
      .add(L"a", 0b010, fspin_any)    // fermi unoccupied
      .add_union(L"p", {L"i", L"a"})  // fermi all
      ;
  if (spconv == SpinConvention::Default) add_fermi_spin(*isr);
  const auto bspin_any = IndexSpace::QuantumNumbers{Spin::any};
  isr->add(L"β", 0b100, bspin_any);  // bose

  auto fermi_isr = std::make_shared<IndexBasisRegistry>(isr->bases());
  fermi_isr->vacuum_occupied_space(L"i");
  fermi_isr->reference_occupied_space(L"i");
  fermi_isr->hole_space(L"i");
  fermi_isr->particle_space(L"a");
  fermi_isr->complete_space(L"p");
  fermi_isr->physical_particle_attribute_mask(bitset_t(fspin_any));

  auto bose_isr = std::make_shared<IndexBasisRegistry>(isr->bases());
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
