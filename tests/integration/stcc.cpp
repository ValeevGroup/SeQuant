#include <SeQuant/version.hpp>

#include <SeQuant/core/math.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/timer.hpp>
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
  std::wcout.precision(std::numeric_limits<double>::max_digits10);
  std::wcerr.precision(std::numeric_limits<double>::max_digits10);
  sequant::set_locale();

  std::cout << "SeQuant revision: " << sequant::git_revision() << "\n";
  std::cout << "Number of threads: " << sequant::num_threads() << "\n\n";

  sequant::set_default_context(
      {.index_space_registry_shared_ptr = mbpt::make_min_sr_spaces(),
       .vacuum = Vacuum::SingleProduct,
       .canonicalization_options =
           CanonicalizeOptions::default_options().copy_and_set(
               CanonicalizationMethod::Complete)});
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

#ifndef NDEBUG
  const size_t DEFAULT_NMAX = 3;
#else
  const size_t DEFAULT_NMAX = 4;
#endif
  const size_t NMAX = argc > 1 ? std::atoi(argv[1]) : DEFAULT_NMAX;

  using mbpt::BiorthogonalizationMethod;
  const std::map<std::string, BiorthogonalizationMethod> bm_str2type{
      {"v1", BiorthogonalizationMethod::V1},
      {"v2", BiorthogonalizationMethod::V2}};
  const std::string biorthogonalization_method_str = argc > 2 ? argv[2] : "v2";
  const auto biorthogonalization_method =
      bm_str2type.at(biorthogonalization_method_str);

  enum class SpintraceMethod { Standard, Naive };
  const std::map<std::string, SpintraceMethod> stm_str2type{
      {"standard", SpintraceMethod::Standard},
      {"naive", SpintraceMethod::Naive}};
  const std::string spintrace_method_str = argc > 3 ? argv[3] : "standard";
  const auto naive_spintrace =
      stm_str2type.at(spintrace_method_str) == SpintraceMethod::Naive;

  // Spin-orbital coupled cluster
  auto cc_r = sequant::mbpt::CC{NMAX}.t();
  for (auto i = 1; i < cc_r.size(); ++i) {
    std::cout << "Spin-orbital CC R" << i << " size: " << cc_r[i]->size()
              << "\n";
  }

  //
  // Closed-shell spintrace (fast)
  //
  std::cout << "\nClosed-shell coupled cluster spintrace with biorthogonal "
            << biorthogonalization_method_str << " transformation:\n";
  std::vector<ExprPtr> cc_st_r(cc_r.size());
  for (auto i = 1; i < cc_r.size(); ++i) {
    const auto tstart = std::chrono::high_resolution_clock::now();
    cc_st_r[i] = mbpt::closed_shell_CC_spintrace(
        cc_r[i], {.method = biorthogonalization_method,
                  .naive_spintrace = naive_spintrace});

    auto tstop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_elapsed = tstop - tstart;
    printf("CC R%d size: %lu time: %5.3f sec.\n", i, cc_st_r[i]->size(),
           time_elapsed.count());
  }

  switch (biorthogonalization_method) {
    case BiorthogonalizationMethod::V1:
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
      } else if (NMAX == 2) {
        runtime_assert(cc_st_r.size() == 3)
            runtime_assert(cc_st_r.at(1)->size() == 26)  // T1
            runtime_assert(cc_st_r.at(2)->size() == 55)  // T2
      }
      break;

    case BiorthogonalizationMethod::V2:
      if (NMAX == 4) {
        runtime_assert(cc_st_r.size() == 5)
            runtime_assert(cc_st_r.at(1)->size() == 30)   // T1
            runtime_assert(cc_st_r.at(2)->size() == 78)   // T2
            runtime_assert(cc_st_r.at(3)->size() == 111)  // T3
            runtime_assert(cc_st_r.at(4)->size() == 149)  // T4
      } else if (NMAX == 3) {
        runtime_assert(cc_st_r.size() == 4)
            runtime_assert(cc_st_r.at(1)->size() == 30)  // T1
            runtime_assert(cc_st_r.at(2)->size() == 73)  // T2
            runtime_assert(cc_st_r.at(3)->size() == 93)  // T3
      } else if (NMAX == 2) {
        runtime_assert(cc_st_r.size() == 3)
            runtime_assert(cc_st_r.at(1)->size() == 26)  // T1
            runtime_assert(cc_st_r.at(2)->size() == 55)  // T2
      }
      break;
    default:
      SEQUANT_ASSERT(false && "unreachable code reached");
      abort();
  }

  return 0;
}
