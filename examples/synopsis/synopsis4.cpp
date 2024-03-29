#include <SeQuant/core/index.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

int main() {
  using namespace sequant;

  // load the QCiFS convention
  mbpt::
      set_default_convention();  // =mbpt::set_default_convention(mbpt::Convention::QCiFS)

  Index i1(L"i_1");  // active occupied
  Index a1(L"a_1");  // active unoccupied
  Index p1(L"p_1");  // any state in computational basis
  // etc.
}
