#include <SeQuant/core/context.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/space_qns.hpp>

void v1() {
  // start-snippet-1
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
  set_default_context(Context(std::move(isr)));

  // now can use space labels to construct Index objects representing said
  // spaces
  Index i1(L"i_1");
  Index a1(L"a_1");
  Index p1(L"p_1");

  // set theoretic operations on spaces
  assert(i1.space().type().includes(a1.space().type()) == false);
  // end-snippet-1
}

void v2() {
  // start-snippet-2
  using namespace sequant;
  using namespace sequant::mbpt;

  // makes 2 base spaces, i and a, and their union
  auto isr = make_min_sr_spaces();
  set_default_context(Context(isr));

  // set theoretic operations on spaces
  auto i1 = Index(L"i_1");
  auto a1 = Index(L"a_1");
  assert(i1.space().attr().intersection(a1.space().attr()).type() ==
         IndexSpace::Type::null);
  assert(i1.space().attr().intersection(a1.space().attr()).qns() ==
         mbpt::Spin::any);
  // end-snippet-2
}

int main() {
  v1();
  v2();
  return 0;
}
