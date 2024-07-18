#include <SeQuant/core/op.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/timer.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <clocale>

using namespace sequant;
using namespace sequant::mbpt;

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
inline const std::map<std::string, mbpt::CSV> str2uocc = {
    {"std", mbpt::CSV::No}, {"csv", mbpt::CSV::Yes}};

/// maps SPBasis type string to enum
inline const std::map<std::string, SPBasis> str2spbasis = {
    {"so", SPBasis::spinorbital}, {"sf", SPBasis::spinfree}};

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
        eqvec = CC{N}.t(4, P, PMIN);
        break;
      case EqnType::λ:
        eqvec = CC{N}.λ(4);
        break;
    }
    tpool.stop(N);
    const bool spinfree = get_default_context().spbasis() == SPBasis::spinfree;
    std::wcout << std::boolalpha << "CC equations [type=" << type2str.at(type)
               << ",rank=" << N << ",spinfree=" << spinfree
               << ",screen=" << screen << ",use_topology=" << use_topology
               << ",use_connectivity=" << use_connectivity
               << ",canonical_only=" << canonical_only << "] computed in "
               << tpool.read(N) << " seconds" << std::endl;

    // validate spin-free equations against spin-traced spin-orbital equations
    std::vector<ExprPtr> eqvec_sf_ref;
    if (get_default_context().spbasis() == SPBasis::spinfree) {
      auto context_resetter = sequant::set_scoped_default_context(
          sequant::Context(make_min_sr_spaces(), Vacuum::SingleProduct,
                           IndexSpaceMetric::Unit, BraKetSymmetry::conjugate,
                           SPBasis::spinorbital));
      std::vector<ExprPtr> eqvec_so;
      switch (type) {
        case EqnType::t:
          eqvec_so = CC{N}.t(4, P, PMIN);
          break;
        case EqnType::λ:
          eqvec_so = CC{N}.λ();
          break;
      }

      eqvec_sf_ref.resize(eqvec_so.size());
      for (size_t R = PMIN; R <= P; ++R) {
        auto const ext_idxs =
            external_indices(eqvec_so[R]->at(0)->at(0)->as<Tensor>());
        eqvec_sf_ref[R] = closed_shell_spintrace(eqvec_so[R], ext_idxs);
        if (R == 1) {  // closed_shell_spintrace omits 1-body S
          using ranges::views::transform;
          auto bixs = ext_idxs | transform([](auto&& vec) { return vec[0]; });
          auto kixs = ext_idxs | transform([](auto&& vec) { return vec[1]; });
          auto s_tensor = ex<Tensor>(Tensor{L"S", kixs, bixs});
          eqvec_sf_ref[R] = s_tensor * eqvec_sf_ref[R];
          expand(eqvec_sf_ref[R]);
        }
      }
    }

    for (size_t R = PMIN; R <= P; ++R) {
      std::wcout << "R" << R << "(expS" << N << ") has " << eqvec[R]->size()
                 << " terms:" << std::endl;
      if (print) std::wcout << to_latex_align(eqvec[R], 20, 1) << std::endl;

      // validate known sizes of some CC residuals
      // N.B. # of equations depends on whether we use symmetric or
      // antisymmetric amplitudes
      if (get_default_context().spbasis() == SPBasis::spinorbital) {
        if (type == EqnType::t) {
          if (R == 1 && N == 1) runtime_assert(eqvec[R]->size() == 8);
          if (R == 1 && N == 2) runtime_assert(eqvec[R]->size() == 14);
          if (R == 2 && N == 2) runtime_assert(eqvec[R]->size() == 31);
          if (R == 1 && N == 3) runtime_assert(eqvec[R]->size() == 15);
          if (R == 2 && N == 3) runtime_assert(eqvec[R]->size() == 37);
          if (R == 3 && N == 3) runtime_assert(eqvec[R]->size() == 47);
          if (R == 4 && N == 4) runtime_assert(eqvec[R]->size() == 74);
          if (R == 5 && N == 5) runtime_assert(eqvec[R]->size() == 99);
        } else if (type == EqnType::λ) {
          if (R == 1 && N == 1) runtime_assert(eqvec[R]->size() == 14);
          if (R == 1 && N == 2) runtime_assert(eqvec[R]->size() == 45);
          if (R == 2 && N == 2) runtime_assert(eqvec[R]->size() == 32);
          if (R == 1 && N == 3) runtime_assert(eqvec[R]->size() == 83);
          if (R == 2 && N == 3) runtime_assert(eqvec[R]->size() == 71);
          if (R == 3 && N == 3) runtime_assert(eqvec[R]->size() == 32);
          if (R == 1 && N == 4) runtime_assert(eqvec[R]->size() == 134);
        }
      } else {  // spin-free
        if (type == EqnType::t) {
          if (R == 1 && N == 2) runtime_assert(eqvec[R]->size() == 26);
          if (R == 2 && N == 2) runtime_assert(eqvec[R]->size() == 110);
          if (R == 1 && N == 3) runtime_assert(eqvec[R]->size() == 30);
          if (R == 2 && N == 3) runtime_assert(eqvec[R]->size() == 146);
          if (R == 3 && N == 3) runtime_assert(eqvec[R]->size() == 490);
          if (R == 4 && N == 4) runtime_assert(eqvec[R]->size() == 2150);
        }

        // validate spin-free equations by spin-tracing spin-orbital equations
        const auto should_be_zero = simplify(eqvec_sf_ref[R] - eqvec[R]);
        if (should_be_zero != ex<Constant>(0))
          std::wcout << "Spin-free equations do not match spin-traced "
                        "spin-orbital equations: N="
                     << N << " R=" << R << ":\n"
                     << "spintraced-spinfree = "
                     << to_latex_align(should_be_zero, 0, 1) << std::endl;
        else
          std::wcout << "Spin-free equations match spin-traced "
                        "spin-orbital equations"
                     << std::endl;
        runtime_assert(should_be_zero == ex<Constant>(0));

        // validate sizes of spin-free t equations after biorthogonal transform
        if (type == EqnType::t) {
          auto const ext_idxs =
              external_indices(eqvec[R]->at(0)->at(0)->as<Tensor>());

          // Remove S operator
          for (auto& term : eqvec[R]->expr()) {
            if (term->is<Product>())
              term = remove_tensor(term->as<Product>(), L"S");
          }

          // Biorthogonal transformation
          eqvec[R] = biorthogonal_transform(eqvec[R], ext_idxs);

          // restore the particle symmetrizer
          auto bixs = ext_idxs | ranges::views::transform(
                                     [](auto&& vec) { return vec[0]; });
          auto kixs = ext_idxs | ranges::views::transform(
                                     [](auto&& vec) { return vec[1]; });
          // N.B. external_indices(expr) confuses bra and ket
          eqvec[R] = ex<Tensor>(Tensor{L"S", kixs, bixs}) * eqvec[R];
          eqvec[R] = expand(eqvec[R]);
          simplify(eqvec[R]);

          std::wcout << "biorthogonal spin-free R" << R << "(expS" << N
                     << ") has " << eqvec[R]->size() << " terms:" << std::endl;
          if (print) std::wcout << to_latex_align(eqvec[R], 20, 1) << std::endl;

          if (R == 1 && N == 2) runtime_assert(eqvec[R]->size() == 26);
          if (R == 2 && N == 2) runtime_assert(eqvec[R]->size() == 55);
          if (R == 1 && N == 3) runtime_assert(eqvec[R]->size() == 30);
          if (R == 2 && N == 3) runtime_assert(eqvec[R]->size() == 73);
          if (R == 3 && N == 3) runtime_assert(eqvec[R]->size() == 490);
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
  const mbpt::CSV uocc_type = str2uocc.at(uocc_type_str);
  auto mbpt_ctx = set_scoped_default_mbpt_context(mbpt::Context(uocc_type));

  const std::string spbasis_str = argc > 4 ? argv[4] : "so";
  const SPBasis spbasis = str2spbasis.at(spbasis_str);

  const std::string print_str = argc > 5 ? argv[5] : "noprint";
  const bool print = print_str == "print";

  sequant::detail::OpIdRegistrar op_id_registrar;
  sequant::set_default_context(sequant::Context(
      make_min_sr_spaces(), Vacuum::SingleProduct, IndexSpaceMetric::Unit,
      BraKetSymmetry::conjugate, spbasis));
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  // change to true to print stats
  Logger::instance().wick_stats = false;

  tpool.clear();
  // comment out to run all possible combinations
  compute_all{NMAX, eqn_type}(print);
}
