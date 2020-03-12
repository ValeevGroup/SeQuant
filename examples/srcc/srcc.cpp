#include "../sequant_setup.hpp"

using namespace  sequant;

int main(int argc, char* argv[]) {
  std::setlocale(LC_ALL, "en_US.UTF-8");
  std::wcout.precision(std::numeric_limits<double>::max_digits10);
  std::wcerr.precision(std::numeric_limits<double>::max_digits10);
  std::wcout.sync_with_stdio(false);
  std::wcerr.sync_with_stdio(false);
  std::wcout.imbue(std::locale("en_US.UTF-8"));
  std::wcerr.imbue(std::locale("en_US.UTF-8"));
  std::wcout.sync_with_stdio(true);
  std::wcerr.sync_with_stdio(true);
  sequant::detail::OpIdRegistrar op_id_registrar;

  mbpt::set_default_convention();

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
  //set_num_threads(1);

#ifndef NDEBUG
  const size_t DEFAULT_NMAX = 3;
#else
  const size_t DEFAULT_NMAX = 4;
#endif
  const size_t NMAX = argc > 1 ? std::atoi(argv[1]) : DEFAULT_NMAX;
  // change to true to print out the resulting equations
  constexpr bool print = true;
  // change to true to print stats
  Logger::get_instance().wick_stats = false;

  ranges::for_each(std::array<bool, 2>{false, true}, [=](const bool screen) {
    ranges::for_each(
        std::array<bool, 2>{false, true}, [=](const bool use_topology) {
          ranges::for_each(std::array<bool, 2>{false, true},
                           [=](const bool canonical_only) {
                             tpool.clear();
                             // comment out to run all possible combinations
                             if (screen && use_topology && canonical_only)
                               compute_all{NMAX}(print, screen, use_topology,
                                                 true, canonical_only);
                           });
        });
  });

  return 0;
}
