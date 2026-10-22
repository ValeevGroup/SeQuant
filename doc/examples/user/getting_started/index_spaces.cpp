#include <SeQuant/core/context.hpp>
#include <SeQuant/core/utility/macros.hpp>

int main() {
  // start-snippet-1
  using namespace sequant;

  // the default is to use genuine vacuum
  SEQUANT_ASSERT(get_default_context().vacuum() == Vacuum::Physical);
  // now set the context to a single product of SP states
  set_default_context({.index_basis_registry = IndexBasisRegistry{},
                       .vacuum = Vacuum::SingleProduct});
  SEQUANT_ASSERT(get_default_context().vacuum() == Vacuum::SingleProduct);
  // reset the context back to the default
  reset_default_context();
  SEQUANT_ASSERT(get_default_context().vacuum() == Vacuum::Physical);
  // end-snippet-1

  // start-snippet-2
  IndexBasisRegistry ibr;
  ibr.add("y", 0b01, {.vacuum_occupied = true});
  // end-snippet-2

  // reset ibr for the next example to work
  ibr.clear();

  // start-snippet-3
  ibr.add("y", 0b01, {.vacuum_occupied = true})
      .add("z", 0b10)
      .add("p", 0b11, {.complete = true});
  // end-snippet-3

  return 0;
}
