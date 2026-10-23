#include <SeQuant/core/index.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

int main() {
  using namespace sequant;
  using namespace sequant::mbpt;
  load(Convention::Minimal);

  Index i1("i_1");  // (active) occupied/hole
  Index a1("a_1");  // (active) unoccupied/particle
  Index p1("p_1");  // any state in computational basis
  // etc.
}
