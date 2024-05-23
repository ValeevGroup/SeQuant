#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/latex.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/wick.hpp>

int main() {
  using namespace sequant;

  Index p1(L"p_1"), p2(L"p_2"), p3(L"p_3"), p4(L"p_4");

  auto cp1 = fcrex(p1), cp2 = fcrex(p2);
  auto ap3 = fannx(p3), ap4 = fannx(p4);

  std::wcout << to_latex(ap3 * ap4 * cp1 * cp2) << " = "
             << to_latex(FWickTheorem{ap3 * ap4 * cp1 * cp2}
                             .set_external_indices(std::array{p1, p2, p3, p4})
                             .full_contractions(false)
                             .compute())
             << std::endl;

  assert(FWickTheorem{ap3 * ap4 * cp1 * cp2}
             .set_external_indices(std::array{p1, p2, p3, p4})
             .full_contractions(false)
             .compute()
             ->size() == 7);

  auto nop1 = ex<FNOperator>(WstrList{L"p_1", L"p_2"}, std::array{p3, p4});
  auto nop2 = ex<FNOperator>(WstrList{L"p_5"}, WstrList{L"p_6", L"p_7"});

  std::wcout << to_latex(nop1 * nop2) << " = "
             << to_latex(FWickTheorem{nop1 * nop2}
                             .set_external_indices(
                                 WstrList{L"p_1", L"p_2", L"p_3", L"p_4",
                                          L"p_5", L"p_6", L"p_7"})
                             .full_contractions(false)
                             .compute())
             << std::endl;

  assert(FWickTheorem{nop1 * nop2}
             .set_external_indices(WstrList{L"p_1", L"p_2", L"p_3", L"p_4",
                                            L"p_5", L"p_6", L"p_7"})
             .full_contractions(false)
             .compute()
             ->size() == 3);

  auto nop3 = ex<BNOperator>(std::array{p1, p2}, std::array{p3, p4});
  auto nop4 = ex<BNOperator>(WstrList{L"p_5", L"p_6"}, WstrList{L"p_7"});

  std::wcout << to_latex(nop3 * nop4) << " = "
             << to_latex(BWickTheorem{nop3 * nop4}
                             .set_external_indices(
                                 WstrList{L"p_1", L"p_2", L"p_3", L"p_4",
                                          L"p_5", L"p_6", L"p_7"})
                             .full_contractions(false)
                             .compute())
             << std::endl;

  assert(BWickTheorem{nop3 * nop4}
             .set_external_indices(WstrList{L"p_1", L"p_2", L"p_3", L"p_4",
                                            L"p_5", L"p_6", L"p_7"})
             .full_contractions(false)
             .compute()
             ->size() == 7);

  return 0;
}
