#include <SeQuant/core/context.hpp>

#include <cassert>

int main() {
  // start-snippet-1
  using namespace sequant;

  // the default is to use genuine vacuum
  assert(get_default_context().vacuum() == Vacuum::Physical);
  // make default IndexSpaceRegistry
  IndexSpaceRegistry isr;
  // now set the context to a single product of SP states
  set_default_context(Context{isr, Vacuum::SingleProduct,
                              IndexSpaceMetric::Unit, BraKetSymmetry::symm});
  assert(get_default_context().vacuum() == Vacuum::SingleProduct);
  // reset the context back to the default
  reset_default_context();
  assert(get_default_context().vacuum() == Vacuum::Physical);
  // end-snippet-1

  // start-snippet-2
  isr.add(L"y", 0b01).vacuum_occupied_space(L"i");
  // end-snippet-2

  // start-snippet-3
  isr.add(L"y", 0b01, is_vacuum_occupied);
  // end-snippet-3

  // start-snippet-4
  isr.add(L"y", 0b01, is_vacuum_occupied)
      .add(L"z", 0b10)
      .add(L"p", 0b11, is_complete);
  // end-snippet-4

  return 0;
}
