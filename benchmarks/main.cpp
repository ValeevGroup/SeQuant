#include <benchmark/benchmark.h>

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

using namespace sequant;

int main(int argc, char *argv[]) {
  // Disable multithreading
  set_num_threads(1);
  set_locale();
  set_default_context(Context(mbpt::make_sr_spaces(), Vacuum::SingleProduct,
                              IndexSpaceMetric::Unit, BraKetSymmetry::nonsymm,
                              SPBasis::spinorbital));

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) {
    return 1;
  }

  ::benchmark::RunSpecifiedBenchmarks();
  ::benchmark::Shutdown();

  return 0;
}
