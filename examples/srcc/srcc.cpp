#include <SeQuant/core/op.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/timer.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>

#include <clocale>

using namespace sequant;
using namespace sequant::mbpt::sr;

namespace {

#define runtime_assert(tf)                                         \
  if (!(tf)) {                                                     \
    std::ostringstream oss;                                        \
    oss << "failed assert at line " << __LINE__ << " in function " \
        << __func__;                                               \
    throw std::runtime_error(oss.str().c_str());                   \
  }

TimerPool<32> tpool;

/// types of CC equations to solve
enum class EqnType { t, λ };

/// maps equation type to string
inline const std::map<EqnType, std::wstring> type2str = {
    {EqnType::t, L"t"}, {EqnType::λ, L"lambda"}};

/// maps equation type string to enum
inline const std::map<std::string, EqnType> str2type = {{"t", EqnType::t},
                                                        {"lambda", EqnType::λ}};

/// maps unoccupied basis type string to enum
inline const std::map<std::string, mbpt::Context::CSV> str2uocc = {
    {"std", mbpt::Context::CSV::No}, {"csv", mbpt::Context::CSV::Yes}};

// profiles evaluation of all CC equations for a given ex rank N with projection
// ex rank PMIN .. P
class compute_cceqvec {
  size_t P, PMIN, N;
  EqnType type;

 public:
  compute_cceqvec(size_t p, size_t pmin, size_t n, EqnType t = EqnType::t)
      : P(p), PMIN(pmin), N(n), type(t) {}

  void operator()(bool print, bool screen, bool use_topology,
                  bool use_connectivity, bool canonical_only) {
    tpool.start(N);
    std::vector<ExprPtr> eqvec;
    switch (type) {
      case EqnType::t:
        eqvec = cceqs{N, P, PMIN}.t(screen, use_topology, use_connectivity,
                                    canonical_only);
        break;
      case EqnType::λ:
        eqvec = cceqs{N, P, PMIN}.λ(screen, use_topology, use_connectivity,
                                    canonical_only);
        break;
    }
    tpool.stop(N);
    std::wcout << std::boolalpha << "CC equations [type=" << type2str.at(type)
               << ",rank=" << N << ",screen=" << screen
               << ",use_topology=" << use_topology
               << ",use_connectivity=" << use_connectivity
               << ",canonical_only=" << canonical_only << "] computed in "
               << tpool.read(N) << " seconds" << std::endl;
    for (size_t R = PMIN; R <= P; ++R) {
      std::wcout << "R" << R << "(expS" << N << ") has " << eqvec[R]->size()
                 << " terms:" << std::endl;
      if (print) std::wcout << to_latex_align(eqvec[R], 20, 3) << std::endl;

      // validate known sizes of some CC residuals
      // N.B. # of equations depends on whether we use symmetric or
      // antisymmetric amplitudes
      if (mbpt::get_default_formalism().nbody_interaction_tensor_symm() ==
          mbpt::Context::NBodyInteractionTensorSymm::Yes) {
        if (type == EqnType::t) {
          if (R == 1 && N == 1) runtime_assert(eqvec[R]->size() == 8);
          if (R == 1 && N == 2) runtime_assert(eqvec[R]->size() == 14);
          if (R == 2 && N == 2) runtime_assert(eqvec[R]->size() == 31);
          if (R == 1 && N == 3) runtime_assert(eqvec[R]->size() == 15);
          if (R == 2 && N == 3) runtime_assert(eqvec[R]->size() == 37);
          if (R == 3 && N == 3) runtime_assert(eqvec[R]->size() == 47);
          if (R == 4 && N == 4) runtime_assert(eqvec[R]->size() == 74);
          if (R == 5 && N == 5) runtime_assert(eqvec[R]->size() == 99);
        }
      } else {
        if (type == EqnType::t) {
          if (R == 1 && N == 2) runtime_assert(eqvec[R]->size() == 26);
          if (R == 2 && N == 2) runtime_assert(eqvec[R]->size() == 55);
          if (R == 1 && N == 3) runtime_assert(eqvec[R]->size() == 30);
          if (R == 2 && N == 3) runtime_assert(eqvec[R]->size() == 73);
          if (R == 3 && N == 3) runtime_assert(eqvec[R]->size() == 93);
          if (R == 4 && N == 4) runtime_assert(eqvec[R]->size() == 149);
        }
      }
    }
  }

};  // class compute_cceqvec

// profiles evaluation of all CC equations with ex rank 2 .. N
class compute_all {
  size_t NMAX;
  EqnType type;

 public:
  compute_all(size_t nmax, EqnType t = EqnType::t) : NMAX(nmax), type(t) {}

  void operator()(bool print = true, bool screen = true,
                  bool use_topology = true, bool use_connectivity = true,
                  bool canonical_only = true) {
    for (size_t N = 1; N <= NMAX; ++N)
      compute_cceqvec{N, 1, N, type}(print, screen, use_topology,
                                     use_connectivity, canonical_only);
  }
};  // class compute_all

}  // namespace

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
  sequant::set_default_context(
      Context(Vacuum::SingleProduct, IndexSpaceMetric::Unit,
              BraKetSymmetry::conjugate, SPBasis::spinorbital));
  mbpt::set_default_convention();

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
  // set_num_threads(1);

#ifndef NDEBUG
  const size_t DEFAULT_NMAX = 3;
#else
  const size_t DEFAULT_NMAX = 4;
#endif
  const size_t NMAX = argc > 1 ? std::atoi(argv[1]) : DEFAULT_NMAX;
  const std::string eqn_type_str = argc > 2 ? argv[2] : "t";
  const EqnType eqn_type = str2type.at(eqn_type_str);
  const std::string uocc_type_str = argc > 3 ? argv[3] : "std";
  const mbpt::Context::CSV uocc_type = str2uocc.at(uocc_type_str);
  auto resetter = set_scoped_default_formalism(mbpt::Context(uocc_type));

  // change to true to print out the resulting equations
  constexpr bool print = false;
  // change to true to print stats
  Logger::get_instance().wick_stats = false;

  tpool.clear();
  // comment out to run all possible combinations
  compute_all{NMAX, eqn_type}(print);
}
