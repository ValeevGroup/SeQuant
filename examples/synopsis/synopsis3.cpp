#include <SeQuant/core/context.hpp>
#include <SeQuant/core/wick.hpp>

int main() {
  using namespace sequant;

  set_default_context(Context{IndexSpaceRegistry{}
                                  .add(L"y", 0b01, is_vacuum_occupied)
                                  .add(L"z", 0b10)
                                  .add(L"p", 0b11, is_complete),
                              Vacuum::SingleProduct});

  auto cp1 = fcrex(L"p_1"), cp2 = fcrex(L"p_2");
  auto ap3 = fannx(L"p_3"), ap4 = fannx(L"p_4");

  std::wcout << to_latex(ap3 * cp1 * ap4 * cp2) << " = "
             << to_latex(FWickTheorem{ap3 * cp1 * ap4 * cp2}
                             .full_contractions(false)
                             .compute())
             << std::endl;

  return 0;
}
