#include <benchmark/benchmark.h>

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

using namespace sequant;

int main(int argc, char *argv[]) {
  benchmark::MaybeReenterWithoutASLR(argc, argv);

  // Disable multithreading
  set_num_threads(1);
  set_locale();
  disable_thousands_separator();
  auto idxreg = mbpt::make_sr_spaces();
  Context fermi_ctx = Context({.index_space_registry_shared_ptr = idxreg,
                               .vacuum = Vacuum::SingleProduct,
                               .braket_symmetry = BraKetSymmetry::Nonsymm});
  set_default_context(fermi_ctx);

  Context bose_einstein_ctx =
      Context({.index_space_registry_shared_ptr = idxreg,
               .vacuum = Vacuum::Physical,
               .braket_symmetry = BraKetSymmetry::Nonsymm});

  set_default_context(bose_einstein_ctx, Statistics::BoseEinstein);

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) {
    return 1;
  }

  ::benchmark::RunSpecifiedBenchmarks();
  ::benchmark::Shutdown();

  return 0;
}
