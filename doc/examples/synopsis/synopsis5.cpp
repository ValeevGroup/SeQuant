//
// Created by Conner Masteran on 5/1/24.
//

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

void test0() {
  using namespace sequant;
  IndexBasisRegistry ibr;

  // base spaces
  ibr.add("i", 0b01).add("a", 0b10);
  // union of 2 base spaces
  // can create manually, as ibr.add("p", 0b11) , or explicitly ...
  ibr.add_union("p", {"i", "a"});  // union of i and a

  // can access unions and intersections of base and composite spaces
  SEQUANT_ASSERT(ibr.unIon("i", "a") == ibr.retrieve_space("p"));
  SEQUANT_ASSERT(ibr.intersection("p", "i") == ibr.retrieve_space("i"));

  // to use the vocabulary defined by ibr use it to make a Context object and
  // make it the default
  set_default_context({.index_basis_registry = std::move(ibr)});

  // now can use space labels to construct Index objects representing said
  // spaces
  Index i1("i_1");
  Index a1("a_1");
  Index p1("p_1");

  // set theoretic operations on spaces
  SEQUANT_ASSERT(i1.space().attr().intersection(a1.space().attr()) ==
                 IndexSpace::Attr::null);
}

void test1() {
  using namespace sequant;
  IndexBasisRegistry ibr;

  ibr.add("x", 0b001)
      .add("y", 0b010)
      .add("z", 0b100)
      .add("xy", 0b011)                 // union of x and y
      .add_union("yz", {"y", "z"})      // union of y and z, explicit
      .add_union("xyz", {"xy", "yz"});  // union of x, y, and z

  SEQUANT_ASSERT(ibr.unIon("x", "y") == ibr.retrieve_space("xy"));
  SEQUANT_ASSERT(ibr.intersection("xyz", "y") == ibr.retrieve_space("y"));

  // use the registry in global context to streamline composition
  set_default_context({.index_basis_registry = std::move(ibr)});
  Index xy1("xy_1");  // now can use space labels to define indices
}

void test2() {
  using namespace sequant;
  using namespace sequant::mbpt;
  // makes 2 base spaces, i and a, and their union
  auto ibr = make_min_sr_spaces();
  set_default_context({.index_basis_registry_shared_ptr = ibr});

  // set theoretic operations on spaces
  auto i1 = Index("i_1");
  auto a1 = Index("a_1");
  SEQUANT_ASSERT(i1.space().attr().intersection(a1.space().attr()).type() ==
                 IndexSpace::Type::null);
  SEQUANT_ASSERT(i1.space().attr().intersection(a1.space().attr()).qns() ==
                 Spin::any);
}

int main() {
  test0();
  test1();
  test2();

  return 0;
}
