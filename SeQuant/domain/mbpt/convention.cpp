//
// Created by Eduard Valeyev on 2019-04-01.
//

#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/rules/df.hpp>
#include <SeQuant/domain/mbpt/space_qns.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/context.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>

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
  std::shared_ptr<IndexBasisRegistry> ibr;
  switch (conv) {
    case Convention::Minimal:
      ibr = make_min_sr_spaces(spconv);
      break;
    case Convention::SR:
      ibr = make_sr_spaces(spconv);
      break;
    case Convention::MinimalMR:
      ibr = make_min_mr_spaces(spconv);
      break;
    case Convention::MR:
      ibr = make_mr_spaces(spconv);
      break;
    case Convention::F12:
      ibr = make_F12_sr_spaces(spconv);
      break;
    case Convention::QCiFS:
      ibr = make_legacy_spaces(spconv);
      break;
  }
  sequant::Context ctx = get_default_context_snapshot();
  ctx.set(std::move(ibr));
  ctx.set(Vacuum::SingleProduct);
  set_default_context(std::move(ctx));
}

void add_fermi_spin(IndexBasisRegistry& ibr) {
  IndexBasisRegistry result = ibr;

  for (auto&& space : ibr) {
    if (space.base_key() != "") {
      IndexSpace spin_up(spinannotation_add(space.base_key(), Spin::alpha),
                         space.type(), Spin::alpha, space.dimension());
      IndexSpace spin_down(spinannotation_add(space.base_key(), Spin::beta),
                           space.type(), Spin::beta, space.dimension());
      result.add(spin_up);
      result.add(spin_down);
    }
  }
  const bool nulltype_ok = true;
  result.reference_occupied_space(ibr.reference_occupied_space(nulltype_ok));
  result.vacuum_occupied_space(ibr.vacuum_occupied_space(nulltype_ok));
  result.particle_space(ibr.particle_space(nulltype_ok));
  result.hole_space(ibr.hole_space(nulltype_ok));
  result.complete_space(ibr.complete_space(nulltype_ok));

  ibr = std::move(result);
}

namespace {

/// @throw Exception, naming @p caller, if a label of @p bases is registered
/// in @p ibr or @p instance of its space is named there
void require_unregistered(
    const IndexBasisRegistry& ibr, const char* caller,
    const container::svector<std::pair<std::string, IndexSpace>>& bases,
    IndexBasis::instance_type instance) {
  for (const auto& [label, space] : bases) {
    if (ibr.contains(label))
      throw Exception(std::string(caller) + ": the label '" +
                      std::string(label) + "' is already registered");
    if (auto taken = ibr.basis_label(IndexBasis{space, instance}))
      throw Exception(std::string(caller) + ": instance " +
                      std::to_string(instance) + " of space " +
                      std::string(space.base_key()) + " is already named '" +
                      std::string(*taken) + "'");
  }
}

}  // namespace

void add_ao_basis(std::shared_ptr<IndexBasisRegistry>& ibr,
                  IndexSpace::QuantumNumbers spin_any) {
  // the OBS AOs span the complete space less the CABS `α'`, if one is
  // registered as a space: the CABS is the F12 complement of the OBS, and
  // has no role. Every lookup and check precedes the registration, so a
  // failure leaves the registry untouched; the space is taken by value, since
  // add() may reallocate the registry's table
  if (!ibr->complete_space(/* nulltype_ok = */ true))
    throw Exception(
        "add_ao_basis: the registry designates no complete space, which the "
        "OBS AO basis 'μ' spans");
  IndexSpace::Type obs_type = ibr->complete_space(spin_any).type();
  if (const auto* cabs = ibr->retrieve_space_ptr("α'"))
    obs_type = obs_type.xOr(obs_type.intersection(cabs->type()));
  const auto* obs_ptr = ibr->retrieve_space_ptr(obs_type, spin_any);
  if (!obs_ptr)
    throw Exception(
        "add_ao_basis: the AO basis 'μ' spans the complete space less the "
        "CABS, which is not registered");
  const IndexSpace obs = *obs_ptr;
  require_unregistered(*ibr, "add_ao_basis", {{"μ", obs}},
                       default_ao_basis_instance);
  ibr->add("μ", obs, default_ao_basis_instance,
           {.metric = IndexBasisMetric::General});
}

void add_ao_basis(std::shared_ptr<IndexBasisRegistry>& ibr,
                  IndexSpace::QuantumNumbers spin_any, bool vbs, bool abs,
                  IndexBasis::instance_type instance) {
  // matches the MPQC layout, see spindex.h
  // this will not work for MR
  // the AOs of a basis set span an orbital space; the union of two AO bases
  // spans the union of their spaces, which must be registered. Every space is
  // looked up and every label checked before anything is registered, so a
  // failure leaves the registry untouched. The spaces are taken by value:
  // add() may reallocate the registry's table
  auto union_space = [&](std::string_view label, const IndexSpace& s1,
                         const IndexSpace& s2) -> IndexSpace {
    const auto* space =
        ibr->retrieve_space_ptr(s1.type() | s2.type(), spin_any);
    if (!space)
      throw Exception("add_ao_basis: the AO basis '" + std::string(label) +
                      "' spans the union of the spaces " +
                      std::string(s1.base_key()) + " and " +
                      std::string(s2.base_key()) + ", which is not registered");
    return *space;
  };
  container::svector<std::pair<std::string, IndexSpace>> bases;
  const IndexSpace obs = ibr->retrieve_space(vbs ? "m" : "p");
  bases.emplace_back("μ", obs);  // OBS AO
  IndexSpace vbs_plus;
  if (vbs) {
    const IndexSpace vbs_space = ibr->retrieve_space("e");
    vbs_plus = union_space("Γ", obs, vbs_space);
    bases.emplace_back("Α", vbs_space);  // VBS AO
    bases.emplace_back("Γ", vbs_plus);   // VBS+ = OBS + VBS
  }
  if (abs) {
    const IndexSpace abs_space = ibr->retrieve_space("α'");
    bases.emplace_back("σ", abs_space);  // ABS AO in F12 methods
    bases.emplace_back("ρ", union_space("ρ", obs, abs_space));  // OBS + ABS
    if (vbs)  // VABS+ = VBS+ + ABS
      bases.emplace_back("Ρ", union_space("Ρ", vbs_plus, abs_space));
  }
  require_unregistered(*ibr, "add_ao_basis", bases, instance);
  for (const auto& [label, space] : bases)
    ibr->add(label, space, instance, {.metric = IndexBasisMetric::General});
}

void add_ao_spaces(std::shared_ptr<IndexBasisRegistry>& ibr,
                   IndexSpace::QuantumNumbers spin_any, bool vbs, bool abs) {
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  add_ao_basis(ibr, spin_any, vbs, abs);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
}

void add_pao_spaces(std::shared_ptr<IndexBasisRegistry>& ibr,
                    IndexSpace::QuantumNumbers spin_any) {
  add_pao_basis(ibr, spin_any);
}

void add_pao_basis(std::shared_ptr<IndexBasisRegistry>& ibr,
                   IndexSpace::QuantumNumbers spin_any,
                   const PAOBasisSpec& spec) {
  // the PAOs are the AOs projected on the particle space (of either spin),
  // so every PAO basis follows the AO basis for its extent, metric and
  // field. Everything is checked before anything is registered, so a failure
  // leaves the registry untouched. The spaces are taken by value: add() may
  // reallocate the table
  const std::string label(spec.label), followed(spec.ao_basis);
  if (label == followed)
    throw Exception("add_pao_basis: '" + std::string(label) +
                    "' is the label of the AO basis the PAO bases follow");
  const auto uocc_type = ibr->particle_space(/* nulltype_ok = */ false);
  container::svector<std::pair<std::string, IndexSpace>> bases;
  bases.emplace_back(label, ibr->retrieve_space(uocc_type, spin_any));
  for (const auto spin : {Spin::alpha, Spin::beta})
    if (const auto* uocc_spin = ibr->retrieve_space_ptr(uocc_type, spin))
      bases.emplace_back(spinannotation_add(label, spin), *uocc_spin);
  require_unregistered(*ibr, "add_pao_basis", bases, spec.instance);
  bool registered_ao = false;
  if (!ibr->contains(followed)) {
    if (followed != "μ")
      throw Exception("add_pao_basis: the AO basis '" + std::string(followed) +
                      "' the PAO bases follow is not registered");
    add_ao_basis(ibr, spin_any);
    registered_ao = true;
  } else if (!ibr->retrieve(followed).has_basis_instance())
    throw Exception("add_pao_basis: '" + std::string(followed) +
                    "' is a space, not the AO basis the PAO bases follow");
  else if (auto source = ibr->follows(followed))
    throw Exception("add_pao_basis: the AO basis '" + std::string(followed) +
                    "' follows '" + *source +
                    "', so the PAO bases cannot follow it");
  // the AOs must span the PAOs they are projected on
  const IndexSpace ao_space = ibr->retrieve(followed).space();
  for (const auto& [l, space] : bases)
    if (!includes(ao_space, space)) {
      if (registered_ao) ibr->remove(followed);
      throw Exception("add_pao_basis: the AO basis '" + std::string(followed) +
                      "' spans " + std::string(ao_space.base_key()) +
                      ", which does not include the particle space " +
                      std::string(space.base_key()) + " of '" + std::string(l) +
                      "'");
    }
  for (const auto& [l, space] : bases)
    ibr->add(l, space, spec.instance, {.follows = followed});
}

void add_df_spaces(std::shared_ptr<IndexBasisRegistry>& ibr) {
  // matches the MPQC layout, see spindex.h
  ibr->add(IndexSpace{"Κ", 0b00001, TensorFactorizationQNS::df});  // DFBS AO
}

void add_batching_spaces(std::shared_ptr<IndexBasisRegistry>& ibr) {
  ibr->add(IndexSpace{"z", 0b100000, BatchingQNS::batch});  // Batching Space
}

void add_thc_spaces(std::shared_ptr<IndexBasisRegistry>& ibr) {
  ibr->add(IndexSpace{"L", 0b000001, TensorFactorizationQNS::thc})  // THC AO
      ;
}

std::shared_ptr<IndexBasisRegistry> make_min_sr_spaces(SpinConvention spconv) {
  auto ibr = std::make_shared<IndexBasisRegistry>();

  const auto spin_any = IndexSpace::QuantumNumbers{
      spconv == SpinConvention::Legacy ? Spin::null : Spin::any};
  ibr->add("i", 0b01,
           {.qns = spin_any,
            .vacuum_occupied = true,
            .reference_occupied = true,
            .hole = true})
      .add("a", 0b10, {.qns = spin_any, .particle = true})
      .add_union("p", {"i", "a"}, {.complete = true});
  if (spconv == SpinConvention::Default) add_fermi_spin(*ibr);
  ibr->physical_particle_attribute_mask(bitset_t(spin_any));

  return ibr;
}

// Multireference supspace uses a subset of its occupied orbitals to define a
// vacuum occupied subspace.
//  this leaves an active space which is partially occupied/unoccupied. This
//  definition is convenient when coupled with SR vacuum.
std::shared_ptr<IndexBasisRegistry> make_mr_spaces(SpinConvention spconv) {
  auto ibr = std::make_shared<IndexBasisRegistry>();

  const auto spin_any = IndexSpace::QuantumNumbers{
      spconv == SpinConvention::Legacy ? Spin::null : Spin::any};
  ibr->add("o", 0b00001, {.qns = spin_any})
      .add("i", 0b00010, {.qns = spin_any})
      .add("u", 0b00100, {.qns = spin_any})
      .add("a", 0b01000, {.qns = spin_any})
      .add("g", 0b10000, {.qns = spin_any})
      .add_union("O", {"o", "i"}, {.vacuum_occupied = true})
      .add_union("M", {"o", "i", "u"}, {.reference_occupied = true})
      .add_union("I", {"i", "u"}, {.hole = true})
      .add_union("E", {"u", "a", "g"})
      .add_union("A", {"u", "a"}, {.particle = true})
      .add_union("p", {"M", "E"}, {.complete = true});

  if (spconv == SpinConvention::Default) add_fermi_spin(*ibr);
  ibr->physical_particle_attribute_mask(bitset_t(spin_any));

  return ibr;
}

std::shared_ptr<IndexBasisRegistry> make_min_mr_spaces(SpinConvention spconv) {
  auto ibr = std::make_shared<IndexBasisRegistry>();

  const auto spin_any = IndexSpace::QuantumNumbers{
      spconv == SpinConvention::Legacy ? Spin::null : Spin::any};
  ibr->add("i", 0b0001, {.qns = spin_any, .vacuum_occupied = true})
      .add("u", 0b0010, {.qns = spin_any})
      .add("a", 0b0100, {.qns = spin_any})
      .add_union("I", {"i", "u"}, {.reference_occupied = true, .hole = true})
      .add_union("A", {"u", "a"}, {.particle = true})
      .add_union("p", {"I", "a"}, {.complete = true});

  if (spconv == SpinConvention::Default) add_fermi_spin(*ibr);
  ibr->physical_particle_attribute_mask(bitset_t(spin_any));

  return ibr;
}

std::shared_ptr<IndexBasisRegistry> make_sr_spaces(SpinConvention spconv) {
  auto ibr = std::make_shared<IndexBasisRegistry>();

  const auto spin_any = IndexSpace::QuantumNumbers{
      spconv == SpinConvention::Legacy ? Spin::null : Spin::any};
  ibr->add("o", 0b0001, {.qns = spin_any})
      .add("i", 0b0010, {.qns = spin_any, .hole = true})
      .add("a", 0b0100, {.qns = spin_any, .particle = true})
      .add("g", 0b1000, {.qns = spin_any})
      .add_union("m", {"o", "i"},
                 {.vacuum_occupied = true, .reference_occupied = true})
      .add_union("e", {"a", "g"})
      .add_union("x", {"i", "a"})
      .add_union("p", {"m", "e"}, {.complete = true});
  if (spconv == SpinConvention::Default) add_fermi_spin(*ibr);
  ibr->physical_particle_attribute_mask(bitset_t(spin_any));

  return ibr;
}

std::shared_ptr<IndexBasisRegistry> make_F12_sr_spaces(SpinConvention spconv) {
  auto ibr = std::make_shared<IndexBasisRegistry>();

  const auto spin_any = IndexSpace::QuantumNumbers{
      spconv == SpinConvention::Legacy ? Spin::null : Spin::any};
  ibr->add("o", 0b00001, {.qns = spin_any})
      .add("i", 0b00010, {.qns = spin_any, .hole = true})
      .add("a", 0b00100, {.qns = spin_any, .particle = true})
      .add("g", 0b01000, {.qns = spin_any})
      .add("α'", 0b10000, {.qns = spin_any})
      .add_union("m", {"o", "i"},
                 {.vacuum_occupied = true, .reference_occupied = true})
      .add_union("e", {"a", "g"})
      .add_union("x", {"i", "a"})
      .add_union("p", {"m", "e"})
      .add_unIon("h", {"x", "g"})
      .add_unIon("c", {"g", "α'"})
      .add_union("α", {"e", "α'"})
      .add_union("H", {"i", "α"})
      .add_union("κ", {"p", "α'"}, {.complete = true});
  if (spconv == SpinConvention::Default) add_fermi_spin(*ibr);
  ibr->physical_particle_attribute_mask(bitset_t(spin_any));

  return ibr;
}

std::shared_ptr<IndexBasisRegistry> make_legacy_spaces(SpinConvention spconv) {
  auto ibr = std::make_shared<IndexBasisRegistry>();

  const auto spin_any = IndexSpace::QuantumNumbers{
      spconv == SpinConvention::Legacy ? Spin::null : Spin::any};
  ibr->add("o", 0b0000001, {.qns = spin_any})
      .add("n", 0b0000010, {.qns = spin_any})
      .add("i", 0b0000100, {.qns = spin_any, .hole = true})
      .add("u", 0b0001000, {.qns = spin_any})
      .add("a", 0b0010000, {.qns = spin_any, .particle = true})
      .add("g", 0b0100000, {.qns = spin_any})
      .add("α'", 0b1000000, {.qns = spin_any})
      .add_union("m", {"o", "n", "i"},
                 {.vacuum_occupied = true, .reference_occupied = true})
      .add_union("M", {"m", "u"})
      .add_union("e", {"a", "g"})
      .add_union("E", {"u", "e"})
      .add_union("x", {"i", "u", "a"})
      .add_union("p", {"m", "x", "e"})
      .add_union("κ", {"p", "α'"}, {.complete = true});

  if (spconv == SpinConvention::Default) add_fermi_spin(*ibr);
  ibr->physical_particle_attribute_mask(bitset_t(spin_any));

  return ibr;
}

std::pair<std::shared_ptr<IndexBasisRegistry>,
          std::shared_ptr<IndexBasisRegistry>>
make_fermi_and_bose_spaces(SpinConvention spconv) {
  auto ibr = std::make_shared<IndexBasisRegistry>();

  const auto fspin_any = IndexSpace::QuantumNumbers{
      spconv == SpinConvention::Legacy ? Spin::null : Spin::any};
  ibr->add("i", 0b001, {.qns = fspin_any})  // fermi occupied
      .add("a", 0b010, {.qns = fspin_any})  // fermi unoccupied
      .add_union("p", {"i", "a"})           // fermi all
      ;
  if (spconv == SpinConvention::Default) add_fermi_spin(*ibr);
  const auto bspin_any = IndexSpace::QuantumNumbers{Spin::any};
  ibr->add("β", 0b100, {.qns = bspin_any});  // bose

  auto fermi_ibr = std::make_shared<IndexBasisRegistry>(ibr->bases());
  fermi_ibr->vacuum_occupied_space("i");
  fermi_ibr->reference_occupied_space("i");
  fermi_ibr->hole_space("i");
  fermi_ibr->particle_space("a");
  fermi_ibr->complete_space("p");
  fermi_ibr->physical_particle_attribute_mask(bitset_t(fspin_any));

  auto bose_ibr = std::make_shared<IndexBasisRegistry>(ibr->bases());
  bose_ibr->vacuum_occupied_space(IndexSpace::null);
  bose_ibr->reference_occupied_space(IndexSpace::null);
  bose_ibr->hole_space(IndexSpace::null);
  bose_ibr->particle_space("β");
  bose_ibr->complete_space("β");
  bose_ibr->physical_particle_attribute_mask(bitset_t(bspin_any));

  return std::make_pair(std::move(fermi_ibr), std::move(bose_ibr));
}

}  // namespace mbpt
}  // namespace sequant
