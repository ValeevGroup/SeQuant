#include <SeQuant/core/context.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/wick.hpp>

int main() {
  using namespace sequant;

  set_default_context(
      {.index_basis_registry = IndexBasisRegistry{}
                                   .add("y", 0b01, {.vacuum_occupied = true})
                                   .add("z", 0b10)
                                   .add("p", 0b11, {.complete = true}),
       .vacuum = Vacuum::SingleProduct});

  auto cp1 = fcrex("p_1"), cp2 = fcrex("p_2");
  auto ap3 = fannx("p_3"), ap4 = fannx("p_4");

  std::cout << to_latex(ap3 * cp1 * ap4 * cp2) << " = "
            << to_latex(FWickTheorem{ap3 * cp1 * ap4 * cp2}
                            .full_contractions(false)
                            .compute())
            << std::endl;

  return 0;
}
