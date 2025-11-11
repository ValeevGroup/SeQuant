#include <SeQuant/version.hpp>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/timer.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/domain/mbpt/biorthogonalization.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <format>

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
    {"so", SPBasis::Spinor}, {"sf", SPBasis::Spinfree}};

// profiles evaluation of all CC equations for a given ex rank N with projection
// ex rank PMIN .. P
class compute_cceqvec {
  size_t P, PMIN, N;
  EqnType type;

 public:
  compute_cceqvec(size_t p, size_t pmin, size_t n, EqnType t = EqnType::t)
      : P(nₚ(p)), PMIN(pmin), N(n), type(t) {}

  void operator()(bool print, bool screen, bool use_topology,
                  bool use_connectivity) {
    tpool.start(N);
    std::vector<ExprPtr> eqvec;
    switch (type) {
      case EqnType::t:
        eqvec = CC{N, CC::Ansatz::T, screen, use_topology, use_connectivity}.t(
            4, P, PMIN);
        break;
      case EqnType::λ:
        eqvec =
            CC{N, CC::Ansatz::T, screen, use_topology, use_connectivity}.λ(4);
        break;
    }
    tpool.stop(N);
    const bool spinfree = get_default_context().spbasis() == SPBasis::Spinfree;
    std::wcout << std::format(
        L"\nCC equations [type={}, rank={}, spinfree={}, screen={}, "
        L"use_topology={}, use_connectivity={}] "
        L"computed in {} seconds",
        type2str.at(type), N, spinfree, screen, use_topology, use_connectivity,
        tpool.read(N));

    // validate spin-free equations against spin-traced spin-orbital equations
    std::vector<ExprPtr> eqvec_sf_ref;
    if (get_default_context().spbasis() == SPBasis::Spinfree) {
      auto context_resetter = sequant::set_scoped_default_context(
          {.index_space_registry_shared_ptr = make_min_sr_spaces(),
           .vacuum = Vacuum::SingleProduct});
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
          auto s_tensor = ex<Tensor>(Tensor{L"S", bra(kixs), ket(bixs)});
          eqvec_sf_ref[R] = s_tensor * eqvec_sf_ref[R];
          expand(eqvec_sf_ref[R]);
        }
      }
    }

    for (size_t R = PMIN; R <= P; ++R) {
      std::wcout << std::format(L"\nR{}(expS{}) has {} terms", R, N,
                                eqvec[R]->size());
      if (print) std::wcout << to_latex_align(eqvec[R], 20, 1) << std::endl;

      // validate known sizes of some CC residuals
      // N.B. # of equations depends on whether we use symmetric or
      // antisymmetric amplitudes
      if (get_default_context().spbasis() == SPBasis::Spinor) {
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
          std::wcout << std::format(
              L"Spin-free equations do not match spin-traced "
              L"spin-orbital equations: N={} R={}:\n"
              L"spintraced-spinfree = {}\n",
              N, R, to_latex_align(should_be_zero, 0, 1));
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
              term = remove_tensor(term.as_shared_ptr<Product>(), L"S");
          }

          // Biorthogonal transformation
          eqvec[R] = biorthogonal_transform(eqvec[R], ext_idxs);

          // restore the particle symmetrizer to then expand it in order to get
          // all the raw equations
          auto bixs = ext_idxs | ranges::views::transform(
                                     [](auto&& vec) { return vec[0]; });
          auto kixs = ext_idxs | ranges::views::transform(
                                     [](auto&& vec) { return vec[1]; });
          // N.B. external_indices(expr) confuses bra and ket
          if (bixs.size() > 1) {
            eqvec[R] =
                ex<Tensor>(Tensor{L"S", bra(kixs), ket(bixs)}) * eqvec[R];
          }
          simplify(eqvec[R]);

          // expand the particle symmetrizer to get all the raw equations
          eqvec[R] = S_maps(eqvec[R]);
          canonicalize(eqvec[R]);

          // apply WK_biorthogonalization_filter to get only terms with large
          // coefficients
          eqvec[R] = WK_biorthogonalization_filter(eqvec[R], ext_idxs);

          // restore the particle symmetrizer again to get the most compact set
          // of equations
          eqvec[R] = ex<Tensor>(Tensor{L"S", bra(kixs), ket(bixs)}) * eqvec[R];
          eqvec[R] = expand(eqvec[R]);

          // apply normalization and rescaling factors
          rational combined_factor;
          if (ext_idxs.size() <= 2) {
            combined_factor = rational(1, factorial(ext_idxs.size()));
          } else {
            auto fact_n = factorial(ext_idxs.size());
            combined_factor = rational(
                1, fact_n - 1);  // this is (1/fact_n) * (fact_n/(fact_n-1))
          }
          eqvec[R] = ex<Constant>(combined_factor) * eqvec[R];
          simplify(eqvec[R]);

          // WK_biorthogonalization_filter method removes the redundancy caused
          // by biorthogonal transformation and gives the most compact set of
          // equations. However, we need to restore the effects of those deleted
          // terms. So, after evaluate_symm call in sequant evaluation scope, we
          // need to call evaluate_biorthogonal_nns_project.

          std::wcout << std::format(
              L"biorthogonal spin-free R{}(expS{}) has {} terms:\n", R, N,
              eqvec[R]->size());
          if (print) std::wcout << to_latex_align(eqvec[R], 20, 1) << std::endl;

          if (R == 1 && N == 2) runtime_assert(eqvec[R]->size() == 26);
          if (R == 2 && N == 2) runtime_assert(eqvec[R]->size() == 55);
          if (R == 1 && N == 3) runtime_assert(eqvec[R]->size() == 30);
          if (R == 2 && N == 3) runtime_assert(eqvec[R]->size() == 73);
          if (R == 3 && N == 3) runtime_assert(eqvec[R]->size() == 93);
          if (R == 3 && N == 4) runtime_assert(eqvec[R]->size() == 111);
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
                  bool use_topology = true, bool use_connectivity = true) {
    for (size_t N = 1; N <= NMAX; ++N)
      compute_cceqvec{N, 1, N, type}(print, screen, use_topology,
                                     use_connectivity);
  }
};  // class compute_all

}  // namespace

int main(int argc, char* argv[]) {
  std::wcout.precision(std::numeric_limits<double>::max_digits10);
  std::wcerr.precision(std::numeric_limits<double>::max_digits10);
  sequant::set_locale();

  std::cout << "SeQuant revision: " << sequant::git_revision() << "\n";
  std::cout << "Number of threads: " << sequant::num_threads() << "\n\n";

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
  sequant::set_default_context(
      sequant::Context({.index_space_registry_shared_ptr =
                            make_min_sr_spaces(SpinConvention::None),
                        .vacuum = Vacuum::SingleProduct,
                        .spbasis = spbasis}));
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  // change to true to print stats
  Logger::instance().wick_stats = false;

  tpool.clear();
  // comment out to run all possible combinations
  compute_all{NMAX, eqn_type}(print, /*screen*/ true, /*use_topology*/ true,
                              /*use_connectivity*/ true);
}
