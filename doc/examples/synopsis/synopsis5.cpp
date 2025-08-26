//
// Created by Conner Masteran on 5/1/24.
//

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

void test0() {
  using namespace sequant;
  IndexSpaceRegistry isr;

  // base spaces
  isr.add(L"i", 0b01).add(L"a", 0b10);
  // union of 2 base spaces
  // can create manually, as isr.add(L"p", 0b11) , or explicitly ...
  isr.add_union(L"p", {L"i", L"a"});  // union of i and a

  // can access unions and intersections of base and composite spaces
  assert(isr.unIon(L"i", L"a") == isr.retrieve(L"p"));
  assert(isr.intersection(L"p", L"i") == isr.retrieve(L"i"));

  // to use the vocabulary defined by isr use it to make a Context object and
  // make it the default
  set_default_context({.index_space_registry = std::move(isr)});

  // now can use space labels to construct Index objects representing said
  // spaces
  Index i1(L"i_1");
  Index a1(L"a_1");
  Index p1(L"p_1");

  // set theoretic operations on spaces
  assert(i1.space().attr().intersection(a1.space().attr()) ==
         IndexSpace::Attr::null);
}

void test1() {
  using namespace sequant;
  IndexSpaceRegistry isr;

  isr.add(L"x", 0b001)
      .add(L"y", 0b010)
      .add(L"z", 0b100)
      .add(L"xy", 0b011)                   // union of x and y
      .add_union(L"yz", {L"y", L"z"})      // union of y and z, explicit
      .add_union(L"xyz", {L"xy", L"yz"});  // union of x, y, and z

  assert(isr.unIon(L"x", L"y") == isr.retrieve(L"xy"));
  assert(isr.intersection(L"xyz", L"y") == isr.retrieve(L"y"));

  // use the registry in global context to streamline composition
  set_default_context({.index_space_registry = std::move(isr)});
  Index xy1(L"xy_1");  // now can use space labels to define indices
}

void test2() {
  using namespace sequant;
  using namespace sequant::mbpt;
  // makes 2 base spaces, i and a, and their union
  auto isr = make_min_sr_spaces();
  set_default_context({.index_space_registry_shared_ptr = isr});

  // set theoretic operations on spaces
  auto i1 = Index(L"i_1");
  auto a1 = Index(L"a_1");
  assert(i1.space().attr().intersection(a1.space().attr()).type() ==
         IndexSpace::Type::null);
  assert(i1.space().attr().intersection(a1.space().attr()).qns() == Spin::any);
}

int main() {
  test0();
  test1();
  test2();

  return 0;
}
