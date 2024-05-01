#include <SeQuant/core/context.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/latex.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

int main() {
  using namespace sequant;

  set_default_context(
      Context{Vacuum::SingleProduct,
                              sequant::mbpt::make_F12_sr_subspaces(),
              IndexSpaceMetric::Unit, BraKetSymmetry::symm});

  get_default_context().index_space_registry()->relabel(L"i", L"y");
  get_default_context().index_space_registry()->relabel(L"α", L"z");

  IndexSpace sp = get_default_context().index_space_registry()->retrieve(L"p");
  Index p1(L"p_1", sp), p2(L"p_2", sp), p3(L"p_3", sp), p4(L"p_4", sp);

  auto cp1 = fcrex(p1), cp2 = fcrex(p2);
  auto ap3 = fannx(p3), ap4 = fannx(p4);

  std::wcout << to_latex(ap3 * cp1 * ap4 * cp2) << " = "
             << to_latex(FWickTheorem{ap3 * cp1 * ap4 * cp2}
                             .set_external_indices(std::array{p1, p2, p3, p4})
                             .full_contractions(false)
                             .compute())
             << std::endl;

  Index y21(L"y_21");  // <- represents IndexSpace::occupied
  Index z1(L"z_1");    // <- represents IndexSpace::complete_unoccupied
  auto op_oo_oo =
      ex<FNOperator>(WstrList{L"y_1", L"y_2"}, WstrList{L"y_3", L"y_4"});

  return 0;
}
