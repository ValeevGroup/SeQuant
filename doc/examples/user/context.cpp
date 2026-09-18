#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

int main() {
  using namespace sequant;

  // start-snippet-1
  // the one-line way to configure SeQuant for standard single-reference
  // quantum chemistry: populates the default IndexSpaceRegistry with the
  // usual occupied/virtual partitioning and switches the default Context to
  // the single-product (quasiparticle) vacuum
  mbpt::load(mbpt::Convention::SR, mbpt::SpinConvention::None);

  const Context& ctx = get_default_context();
  assert(ctx.vacuum() == Vacuum::SingleProduct);
  assert(ctx.index_space_registry() != nullptr);
  // end-snippet-1

  // start-snippet-2
  // most of SeQuant/domain/mbpt additionally needs an OpRegistry, which maps
  // operator labels ("t", "f", "g", ...) to their OpClass (excitation,
  // de-excitation, or general); this lives in a separate, mbpt-specific
  // context that must be configured on top of the core one
  mbpt::set_default_mbpt_context(
      {.op_registry_ptr = mbpt::make_legacy_registry()});

  assert(mbpt::to_op_class(L"t") == mbpt::OpClass::Ex);
  // end-snippet-2

  // start-snippet-3
  // temporarily switching context (e.g. to a spin-free basis for a single
  // calculation) without disturbing the enclosing default: the RAII
  // resetter restores the previous context when it goes out of scope
  {
    auto resetter = set_scoped_default_context({.spbasis = SPBasis::Spinfree});
    assert(get_default_context().spbasis() == SPBasis::Spinfree);
  }
  assert(get_default_context().spbasis() == SPBasis::Spinor);
  // end-snippet-3

  (void)ctx;

  return 0;
}
