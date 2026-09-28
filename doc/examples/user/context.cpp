#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/utility/macros.hpp>
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
  SEQUANT_ASSERT(ctx.vacuum() == Vacuum::SingleProduct);
  SEQUANT_ASSERT(ctx.index_space_registry() != nullptr);
  // end-snippet-1

  // start-snippet-2
  // most of SeQuant/domain/mbpt additionally needs an OpRegistry, which maps
  // operator labels ("t", "f", "g", ...) to their OpClass (excitation,
  // de-excitation, or general); this lives in a separate, mbpt-specific
  // context that must be configured on top of the core one
  mbpt::set_default_mbpt_context(
      {.op_registry_ptr = mbpt::make_legacy_registry()});

  SEQUANT_ASSERT(mbpt::to_op_class(L"t") == mbpt::OpClass::Ex);
  // end-snippet-2

  // start-snippet-4
  // CSV controls whether OpMaker builds excitation/de-excitation operators
  // (e.g. cluster amplitudes "t") with cluster-specific (index-dependent)
  // virtuals or plain, independent ones; it defaults to CSV::No
  auto csv_resetter = mbpt::set_scoped_default_mbpt_context(
      {.csv = mbpt::CSV::Yes, .op_registry_ptr = mbpt::make_legacy_registry()});
  SEQUANT_ASSERT(mbpt::get_default_mbpt_context().csv() == mbpt::CSV::Yes);
  // end-snippet-4

  // start-snippet-3
  // temporarily switching context (e.g. to a spin-free basis for a single
  // calculation) without disturbing the enclosing default: the RAII
  // resetter restores the previous context when it goes out of scope
  {
    auto resetter = set_scoped_default_context({.spbasis = SPBasis::Spinfree});
    SEQUANT_ASSERT(get_default_context().spbasis() == SPBasis::Spinfree);
  }
  SEQUANT_ASSERT(get_default_context().spbasis() == SPBasis::Spinor);
  // end-snippet-3

  (void)ctx;

  return 0;
}
