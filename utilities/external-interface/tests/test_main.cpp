#define CATCH_CONFIG_RUNNER
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

#include <limits>

int main(int argc, char* argv[]) {
  using namespace sequant;

  Catch::Session session;

  std::wcout.precision(std::numeric_limits<double>::max_digits10);
  std::wcerr.precision(std::numeric_limits<double>::max_digits10);
  sequant::set_locale();

  // Mirror the standard unit-test harness context so pipeline steps that need
  // the standard single-reference spaces (spintracing, projection, ...) behave
  // as they do in the external_interface tool; export-focused cases install
  // their own scoped registry via to_export_context() (test_export.hpp).
  sequant::set_default_context(
      {.index_basis_registry_shared_ptr = sequant::mbpt::make_sr_spaces(),
       .vacuum = Vacuum::SingleProduct,
       .spbasis = SPBasis::Spinor,
       .first_dummy_index_ordinal = 100,
       .deserialization_column_symmetry = ColumnSymmetry::Symm,
       .cardinal_tensor_labels = sequant::mbpt::cardinal_tensor_labels()});

  auto mbpt_ctx = mbpt::Context(
      {.csv = mbpt::CSV::No, .op_registry_ptr = mbpt::make_legacy_registry()});
  sequant::mbpt::set_default_mbpt_context(mbpt_ctx);

  return session.run(argc, argv);
}
