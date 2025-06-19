#include <SeQuant/core/math.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/timer.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>
#include <clocale>
#include <fstream>

using namespace sequant;

#define runtime_assert(tf)                                         \
  if (!(tf)) {                                                     \
    std::ostringstream oss;                                        \
    oss << "failed assert at line " << __LINE__                    \
        << " in closed-shell spin-traced coupled cluster example"; \
    throw std::runtime_error(oss.str().c_str());                   \
  }

int main(int argc, char* argv[]) {
  // std::setlocale(LC_ALL, "en_US.UTF-8");
  // std::wcout.precision(std::numeric_limits<double>::max_digits10);
  // std::wcerr.precision(std::numeric_limits<double>::max_digits10);
  // std::wcout.sync_with_stdio(false);
  // std::wcerr.sync_with_stdio(false);
  // std::wcout.imbue(std::locale("en_US.UTF-8"));
  // std::wcerr.imbue(std::locale("en_US.UTF-8"));
  // std::wcout.sync_with_stdio(true);
  // std::wcerr.sync_with_stdio(true);

  // // set global locale, keep precision as what it was, but global locale
  std::wcout.precision(std::numeric_limits<double>::max_digits10);
  std::wcerr.precision(std::numeric_limits<double>::max_digits10);
  sequant::set_locale();

  sequant::set_default_context(Context(
      mbpt::make_min_sr_spaces(), Vacuum::SingleProduct, IndexSpaceMetric::Unit,
      BraKetSymmetry::conjugate, SPBasis::spinorbital));
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

#ifndef NDEBUG
  const size_t DEFAULT_NMAX = 3;
#else
  const size_t DEFAULT_NMAX = 4;
#endif
  const size_t NMAX = argc > 1 ? std::atoi(argv[1]) : DEFAULT_NMAX;

  // Spin-orbital coupled cluster
  auto cc_r = sequant::mbpt::CC{NMAX}.t();
  for (auto i = 1; i < cc_r.size(); ++i) {
    std::cout << "Spin-orbital CC R" << i << " size: " << cc_r[i]->size()
              << "\n";
  }

  //
  // Closed-shell spintrace (fast)
  std::cout << "\nClosed-shell coupled cluster spintrace with biorthogonal "
               "transformation:\n";
  std::vector<ExprPtr> cc_st_r(cc_r.size());
  for (auto i = 1; i < cc_r.size(); ++i) {
    const auto tstart = std::chrono::high_resolution_clock::now();
    cc_st_r[i] = sequant::closed_shell_CC_spintrace(cc_r[i]);
    // cc_st_r[i] = sequant::closed_shell_CC_spintrace_compact_set(cc_r[i]);

    auto tstop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_elapsed = tstop - tstart;
    printf("CC R%d size: %lu time: %5.3f sec.\n", i, cc_st_r[i]->size(),
           time_elapsed.count());
  }

  if (NMAX == 4) {
    runtime_assert(cc_st_r.size() == 5)
        runtime_assert(cc_st_r.at(1)->size() == 30)    // T1
        runtime_assert(cc_st_r.at(2)->size() == 78)    // T2
        runtime_assert(cc_st_r.at(3)->size() == 567)   // T3
        runtime_assert(cc_st_r.at(4)->size() == 2150)  // T4
  } else if (NMAX == 3) {
    runtime_assert(cc_st_r.size() == 4)
        runtime_assert(cc_st_r.at(1)->size() == 30)   // T1
        runtime_assert(cc_st_r.at(2)->size() == 73)   // T2
        runtime_assert(cc_st_r.at(3)->size() == 490)  // T3
  }

  return 0;
}
