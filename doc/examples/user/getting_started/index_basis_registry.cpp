#include <SeQuant/core/basis.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/space_qns.hpp>

#include <algorithm>

void v1() {
  // start-snippet-1
  using namespace sequant;
  IndexBasisRegistry ibr;

  // base spaces
  ibr.add(L"i", 0b01).add(L"a", 0b10);
  // union of 2 base spaces
  // can create manually, as ibr.add(L"p", 0b11) , or explicitly ...
  ibr.add_union(L"p", {L"i", L"a"});  // union of i and a

  // can access unions and intersections of base and composite spaces
  SEQUANT_ASSERT(ibr.unIon(L"i", L"a") == ibr.retrieve_space(L"p"));
  SEQUANT_ASSERT(ibr.intersection(L"p", L"i") == ibr.retrieve_space(L"i"));

  // to use the vocabulary defined by ibr use it to make a Context object and
  // make it the default
  set_default_context({.index_basis_registry = std::move(ibr)});

  // now can use space labels to construct Index objects representing said
  // spaces
  Index i1(L"i_1");
  Index a1(L"a_1");
  Index p1(L"p_1");

  // set theoretic operations on spaces
  SEQUANT_ASSERT(i1.space().type().includes(a1.space().type()) == false);
  // end-snippet-1
}

void v2() {
  // start-snippet-2
  using namespace sequant;
  using namespace sequant::mbpt;

  // makes 2 base spaces, i and a, and their union
  set_default_context(
      {.index_basis_registry_shared_ptr = make_min_sr_spaces()});

  // set theoretic operations on spaces
  auto i1 = Index(L"i_1");
  auto a1 = Index(L"a_1");
  SEQUANT_ASSERT(i1.space().attr().intersection(a1.space().attr()).type() ==
                 IndexSpace::Type::null);
  SEQUANT_ASSERT(i1.space().attr().intersection(a1.space().attr()).qns() ==
                 mbpt::Spin::any);
  // end-snippet-2
}

void v3() {
  // start-snippet-3
  using namespace sequant;
  auto ibr = mbpt::make_min_sr_spaces();
  // a basis instance of a registered space can be registered under its own
  // label: here a localized basis of the space i, as its basis instance 1
  const IndexSpace i = ibr->retrieve_space(L"i");
  ibr->add(L"ĩ", i, 1);
  // a Context owns its registry: register everything first, then hand it over
  // (it is moved in)
  set_default_context({.index_basis_registry_shared_ptr = std::move(ibr)});
  const auto& registry = *get_default_context().index_basis_registry();
  Index loc1(registry.retrieve(L"ĩ"), 1);
  // same space ...
  SEQUANT_ASSERT(loc1.space() == i);
  // ... a specific basis of it, which carries its name: the name is part of
  // the basis, as a space's label is of the space, so the registry's entry
  // is the basis and a bare {space, instance} pair is another, unnamed one
  SEQUANT_ASSERT(loc1.basis() == registry.retrieve(L"ĩ"));
  SEQUANT_ASSERT(loc1.basis() != IndexBasis(i, 1));
  SEQUANT_ASSERT(same_instance(loc1.basis(), IndexBasis(i, 1)));
  SEQUANT_ASSERT(registry.resolve(IndexBasis(i, 1)) == loc1.basis());
  // a name is not a space
  SEQUANT_ASSERT(registry.contains(L"ĩ") &&
                 registry.retrieve_space_ptr(L"ĩ") == nullptr);
  // the space views see i once: the named entry is not a space
  SEQUANT_ASSERT(
      std::ranges::count_if(registry.spaces(), [](const IndexSpace& s) {
        return s.base_key() == L"i";
      }) == 1);
  // indices in a named basis are printed, serialized and deserialized by that
  // label, and constructed from it
  SEQUANT_ASSERT(Index(L"ĩ_1") == loc1);
  SEQUANT_ASSERT(loc1.full_label() == L"ĩ_1");  // not i_1<;1>
  SEQUANT_ASSERT(loc1.basis().base_key() == L"ĩ");
  SEQUANT_ASSERT(serialize(deserialize(L"t{ĩ_1;i_1}"), {.annot_symm = false}) ==
                 L"t{ĩ_1;i_1}");
  // end-snippet-3
}

int main() {
  v1();
  v2();
  v3();
  return 0;
}
