#include <SeQuant/core/context.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/space_qns.hpp>

int main() {
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
}
