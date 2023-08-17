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
        eqvec = cceqs{N, P, PMIN}.t(screen, use_topology, use_connectivity,
                                    canonical_only);
        break;
      case EqnType::λ:
        eqvec = cceqs{N, P, PMIN}.λ(screen, use_topology, use_connectivity,
                                    canonical_only);
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
          Context(Vacuum::SingleProduct, IndexSpaceMetric::Unit,
                  BraKetSymmetry::conjugate, SPBasis::spinorbital));
      std::vector<ExprPtr> eqvec_so;
      switch (type) {
        case EqnType::t:
          eqvec_so = cceqs{N, P, PMIN}.t(screen, use_topology, use_connectivity,
                                         canonical_only);
          break;
        case EqnType::λ:
          eqvec_so = cceqs{N, P, PMIN}.λ(screen, use_topology, use_connectivity,
                                         canonical_only);
          break;
      }

      eqvec_sf_ref.resize(eqvec_so.size());
      for (size_t R = PMIN; R <= P; ++R) {
        auto const ext_idxs = external_indices(eqvec_so[R]);
        eqvec_sf_ref[R] = closed_shell_spintrace(eqvec_so[R], ext_idxs);
      }
    }

    for (size_t R = PMIN; R <= P; ++R) {
      std::wcout << "R" << R << "(expS" << N << ") has " << eqvec[R]->size()
                 << " terms:" << std::endl;
      if (print) std::wcout << to_latex_align(eqvec[R], 20, 3) << std::endl;

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
          // Biorthogonal transformation
          auto eq_biorth = biorthogonalize(
              eqvec[R], BiorthogonalizationMethod::Pseudoinverse);

          std::wcout << "biorthogonal(PI) spin-free R" << R << "(expS" << N
                     << ") has " << eq_biorth->size() << " terms:" << std::endl;
          if (print)
            std::wcout << to_latex_align(eq_biorth, 20, 3) << std::endl;

          // check F*T->R terms
          if (print) {
            auto extract_F_TR_terms = [&](const auto& expr) {
              assert(expr.template is<Sum>());
              auto F_TR_terms =
                  *expr | ranges::views::filter([&](const auto& term) {
                    return term.template is<Product>() &&
                           ranges::any_of(
                               *term,
                               [](const auto& factor) {
                                 return factor.template is<Tensor>() &&
                                        factor.template as<Tensor>().label() ==
                                            L"f";
                               }) &&
                           ranges::any_of(
                               *term,
                               [&](const auto& factor) {
                                 return factor.template is<Tensor>() &&
                                        factor.template as<Tensor>().label() ==
                                            L"t" &&
                                        factor.template as<Tensor>().rank() ==
                                            R;
                               }) &&
                           (ranges::count_if(*term, [&](const auto& factor) {
                              return factor.template is<Tensor>() &&
                                     factor.template as<Tensor>().label() ==
                                         L"t";
                            }) == 1);
                  }) |
                  ranges::to_vector;
              return ex<Sum>(F_TR_terms);
            };
            std::wcout << "F_TR terms (from spin-free equations): "
                       << extract_F_TR_terms(eqvec[R]).to_latex() << std::endl;
            std::wcout << "F_TR terms (from spin-integrated equations): "
                       << extract_F_TR_terms(eqvec_sf_ref[R]).to_latex()
                       << std::endl;
            std::wcout << "F_TR terms (from biorthogonal spin-free equations): "
                       << extract_F_TR_terms(eq_biorth).to_latex() << std::endl;

            // check biorthogonal transformation variants
            auto eq_biorth_qr =
                biorthogonalize(eqvec[R], BiorthogonalizationMethod::QR);
            auto F_TR_pinv_minus_qr =
                simplify(extract_F_TR_terms(eq_biorth) -
                         extract_F_TR_terms(eq_biorth_qr));
            if (F_TR_pinv_minus_qr->size() > 0) {
              std::wcout
                  << "F_TR terms biorthogonal terms (PI minus QR) spin-free R"
                  << R << "(expS" << N << ") has " << F_TR_pinv_minus_qr->size()
                  << " terms:" << std::endl;
              std::wcout << to_latex_align(F_TR_pinv_minus_qr, 20, 3)
                         << std::endl;
            }
          }

          if (R == 1 && N == 2) runtime_assert(eq_biorth->size() == 26);
          if (R == 2 && N == 2) runtime_assert(eq_biorth->size() == 55);
          if (R == 1 && N == 3) runtime_assert(eq_biorth->size() == 30);
          if (R == 2 && N == 3) runtime_assert(eq_biorth->size() == 73);
          if (R == 3 && N == 3) runtime_assert(eq_biorth->size() == 490);
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
  auto resetter = set_scoped_default_formalism(mbpt::Context(uocc_type));

  const std::string spbasis_str = argc > 4 ? argv[4] : "so";
  const SPBasis spbasis = str2spbasis.at(spbasis_str);

  sequant::detail::OpIdRegistrar op_id_registrar;
  sequant::set_default_context(Context(Vacuum::SingleProduct,
                                       IndexSpaceMetric::Unit,
                                       BraKetSymmetry::conjugate, spbasis));
  mbpt::set_default_convention();
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  // change to true to print out the resulting equations
  constexpr bool print = false;
  // change to true to print stats
  Logger::get_instance().wick_stats = false;

  tpool.clear();
  // comment out to run all possible combinations
  compute_all{NMAX, eqn_type}(print);
}
