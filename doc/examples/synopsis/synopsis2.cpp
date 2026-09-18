#include <SeQuant/core/context.hpp>
#include <SeQuant/core/utility/macros.hpp>

int main() {
  using namespace sequant;

  // the default is to use genuine vacuum
  SEQUANT_ASSERT(get_default_context().vacuum() == Vacuum::Physical);
  // create a context with default IndexSpaceRegistry, using a vacuum of a
  // single product of SP states
  set_default_context({.index_space_registry = IndexSpaceRegistry{},
                       .vacuum = Vacuum::SingleProduct});
  SEQUANT_ASSERT(get_default_context().vacuum() == Vacuum::SingleProduct);
  // reset the context back to the default
  reset_default_context();
  SEQUANT_ASSERT(get_default_context().vacuum() == Vacuum::Physical);

  return 0;
}
