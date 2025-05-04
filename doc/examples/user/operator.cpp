//
// Created by Ajay Melekamburath on 5/3/25.
//

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

int main() {
  // start-snippet-0
  using namespace sequant;
  using namespace sequant::mbpt;
  set_default_context(Context(make_min_sr_spaces(), Vacuum::Physical));
  // end-snippet-0

  // start-snippet-1
  using sequant::mbpt::QuantumNumberChange;

  QuantumNumberChange qnc;
  // For a physical vacuum, we track two spaces
  qnc[0] = {0, 2};
  qnc[1] = {0, 2};

  // check if a given state is a vacuum
  auto is_vacuum = mbpt::is_vacuum(qnc);
  assert(!is_vacuum);

  // end-snippet-1

  return 0;
}
