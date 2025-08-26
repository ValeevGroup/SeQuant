#include <SeQuant/core/context.hpp>

#include <cassert>

int main() {
  // start-snippet-1
  using namespace sequant;

  // the default is to use genuine vacuum
  assert(get_default_context().vacuum() == Vacuum::Physical);
  // make default IndexSpaceRegistry
  auto isr = std::make_shared<IndexSpaceRegistry>();
  // now set the context to a single product of SP states
  set_default_context({.index_space_registry_shared_ptr = isr,
                       .vacuum = Vacuum::SingleProduct});
  assert(get_default_context().vacuum() == Vacuum::SingleProduct);
  // reset the context back to the default
  reset_default_context();
  assert(get_default_context().vacuum() == Vacuum::Physical);
  // end-snippet-1

  // start-snippet-2
  isr->add(L"y", 0b01, is_vacuum_occupied);
  // end-snippet-2

  // reset isr for the next example to work
  isr->clear();

  // start-snippet-3
  isr->add(L"y", 0b01, is_vacuum_occupied)
      .add(L"z", 0b10)
      .add(L"p", 0b11, is_complete);
  // end-snippet-3

  return 0;
}
