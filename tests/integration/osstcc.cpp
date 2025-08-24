#include <SeQuant/core/op.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/timer.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <clocale>

using namespace sequant;

#define runtime_assert(tf)                                       \
  if (!(tf)) {                                                   \
    std::ostringstream oss;                                      \
    oss << "failed assert at line " << __LINE__                  \
        << " in open-shell spin-traced coupled cluster example"; \
    throw std::runtime_error(oss.str().c_str());                 \
  }

int main(int argc, char* argv[]) {
  std::wcout.precision(std::numeric_limits<double>::max_digits10);
  std::wcerr.precision(std::numeric_limits<double>::max_digits10);
  sequant::set_locale();

  sequant::set_default_context(Context(
      mbpt::make_min_sr_spaces(), Vacuum::SingleProduct, IndexSpaceMetric::Unit,
      BraKetSymmetry::conjugate, SPBasis::spinor));
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

#ifndef NDEBUG
  const size_t DEFAULT_NMAX = 3;
#else
  const size_t DEFAULT_NMAX = 3;
#endif
  const size_t NMAX = argc > 1 ? std::atoi(argv[1]) : DEFAULT_NMAX;
  if (NMAX >= 4)
    throw std::logic_error(
        "spin-tracing rank-4 and "
        "higher-body equations is not currently supported");

  // Spin-orbital coupled cluster
  auto cc_r = sequant::mbpt::CC{NMAX}.t();
  for (auto i = 1; i < cc_r.size(); ++i) {
    std::cout << "Spin-orbital CC R" << i << " size: " << cc_r[i]->size()
              << "\n";
  }

  //
  // Open-shell spintrace
  //
  std::cout << "\nOpen-shell coupled cluster: nterms per spin blocks: "
            << std::endl;
  std::vector<std::vector<ExprPtr>> os_cc_st_r(cc_r.size());
  for (auto i = 1; i < cc_r.size(); ++i) {
    os_cc_st_r[i] = open_shell_CC_spintrace(cc_r[i]);
    std::cout << "CCK-" << NMAX << " R[" << i << "]  ";
    for (auto const& x : os_cc_st_r[i]) std::cout << x->size() << "  ";
    std::cout << std::endl;
  }

  if (NMAX == 4) {
    runtime_assert(os_cc_st_r.size() == 5)
        runtime_assert(os_cc_st_r.at(1).at(0)->size() == 30)   // T1a
        runtime_assert(os_cc_st_r.at(2).at(1)->size() == 130)  // T2ab
        runtime_assert(os_cc_st_r.at(2).at(2)->size() == 74)   // T2bb
        runtime_assert(os_cc_st_r.at(3).at(1)->size() == 249)  // T3aab
        runtime_assert(os_cc_st_r.at(3).at(3)->size() == 124)  // T3bbb
        runtime_assert(os_cc_st_r.at(4).at(1)->size() == 356)  // T4aaab
        runtime_assert(os_cc_st_r.at(4).at(2)->size() == 384)  // T4aabb
        runtime_assert(os_cc_st_r.at(4).at(4)->size() == 156)  // T4bbbb
  } else if (NMAX == 3) {
    runtime_assert(os_cc_st_r.size() == 4)
        runtime_assert(os_cc_st_r.at(1).at(0)->size() == 30)   // T1a
        runtime_assert(os_cc_st_r.at(2).at(0)->size() == 65)   // T2aa
        runtime_assert(os_cc_st_r.at(2).at(1)->size() == 122)  // T2ab
        runtime_assert(os_cc_st_r.at(3).at(2)->size() == 209)  // T3abb
        runtime_assert(os_cc_st_r.at(3).at(3)->size() == 75)   // T3bbb
  }

  return 0;
}
