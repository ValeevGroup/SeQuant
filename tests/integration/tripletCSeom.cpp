#include <SeQuant/version.hpp>

#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/timer.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <range/v3/algorithm/transform.hpp>

#include <chrono>
#include <sstream>
#include <vector>

using namespace sequant;
using namespace sequant::mbpt;

namespace {
#define runtime_assert(tf)                                           \
  if (!(tf)) {                                                       \
    std::ostringstream oss;                                          \
    oss << "failed assert at line " << __LINE__                      \
        << " in closed-shell triplet equation-of-motion CC example"; \
    throw Exception(oss.str());                                      \
  }

TimerPool<32> timer_pool;

std::pair<size_t, size_t> parse_excitation_manifold(std::string& str) {
  std::pair<size_t, size_t> result;

  ranges::transform(str, str.begin(), ::tolower);
  const auto h_pos = str.find('h');
  const auto p_pos = str.find('p');

  if (h_pos == std::string::npos && p_pos == std::string::npos) {
    throw Exception(
        "Invalid excitation manifold string: must contain 'h' or 'p'");
  }

  if (h_pos != std::string::npos && p_pos == std::string::npos) {
    result.first = std::stoi(str.substr(0, h_pos));
    result.second = 0;
  } else if (p_pos != std::string::npos && h_pos == std::string::npos) {
    result.first = 0;
    result.second = std::stoi(str.substr(0, p_pos));
  } else {
    if (h_pos < p_pos) {
      result.first = std::stoi(str.substr(0, h_pos));
      result.second = std::stoi(str.substr(h_pos + 1, p_pos - h_pos - 1));
    } else {
      result.first = std::stoi(str.substr(p_pos + 1, h_pos - p_pos - 1));
      result.second = std::stoi(str.substr(0, p_pos));
    }
  }

  if (result.first == 0 && result.second == 0)
    throw Exception(
        "Invalid excitation manifold: both particle and hole ranks cannot be "
        "zero");

  return result;
}

enum class EqnType { left, right };

inline const container::map<std::string, EqnType> str2type = {
    {"L", EqnType::left}, {"R", EqnType::right}};

inline const container::map<EqnType, std::wstring> type2wstr = {
    {EqnType::left, L"L"}, {EqnType::right, L"R"}};

class compute_eomcc_closedshell_triplet {
  size_t N, np, nh;
  std::string manifold;
  EqnType type;

 public:
  compute_eomcc_closedshell_triplet(size_t n, const std::string& exc_manifold,
                                    EqnType t = EqnType::right)
      : N(n), manifold(exc_manifold), type(t) {
    std::tie(nh, np) = parse_excitation_manifold(manifold);

    // triplet spintrace currently supports particle-conserving EE only
    if (nh != np) {
      throw Exception(
          "Closed-shell triplet EOM spintrace only supports particle-"
          "conserving (EE) manifolds; got " +
          manifold);
    }
  }

  void operator()(bool print) {
    SEQUANT_ASSERT(get_default_context().spbasis() == SPBasis::Spinor);

    // generate spin-orbital EOM eqs first
    timer_pool.start(N);
    std::vector<ExprPtr> eqvec;
    switch (type) {
      case EqnType::right:
        eqvec = CC{N}.eom_r(nₚ(np), nₕ(nh));
        break;
      case EqnType::left:
        eqvec = CC{N}.eom_l(nₚ(np), nₕ(nh));
        break;
    }
    timer_pool.stop(N);

    std::wcout << std::boolalpha
               << "EOM-CC Equations [type=" << type2wstr.at(type)
               << ", CC rank=" << N
               << ", manifold=" << sequant::toUtf16(manifold) << "]"
               << " computed in " << timer_pool.read(N) << " s\n";

    for (size_t i = 0; i < eqvec.size(); ++i) {
      if (eqvec[i] == nullptr) continue;
      std::wcout << "Spin-orbital R[" << i << "] size: " << eqvec[i]->size()
                 << "\n";
    }

    std::wcout << "\nClosed-shell triplet EOM-CC spintrace:\n";

    auto term_count = [](const ExprPtr& e) -> size_t {
      if (e->is<Constant>()) return e->as<Constant>().value() == 0 ? 0 : 1;
      if (e->is<Sum>()) return e->size();
      return 1;
    };

    timer_pool.start(N + 16);
    for (size_t i = 0; i < eqvec.size(); ++i) {
      if (eqvec[i] == nullptr) continue;

      auto tstart = std::chrono::high_resolution_clock::now();
      const auto st = closed_shell_CC_triplet_spintrace(
          eqvec[i],
          {.compact = false, .residual = TripletResidualKind::Combined});
      auto tstop = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> dt = tstop - tstart;
      std::wcout << "R[" << i << "] size: " << term_count(st)
                 << " time: " << dt.count() << " s\n";

      // validated term counts of the full triplet residual
      const auto n_st = term_count(st);
      if (N == 2 && type == EqnType::right && np == 2 && nh == 2) {
        if (i == 1) runtime_assert(n_st == 42);
        if (i == 2) runtime_assert(n_st == 540);
      }
      if (N == 3 && type == EqnType::right && np == 3 && nh == 3) {
        if (i == 1) runtime_assert(n_st == 54);
        if (i == 2) runtime_assert(n_st == 920);
        if (i == 3) runtime_assert(n_st == 11592);
      }

      // compact residual (the production form)
      tstart = std::chrono::high_resolution_clock::now();
      const auto compact = closed_shell_CC_triplet_spintrace(
          eqvec[i],
          {.compact = true, .residual = TripletResidualKind::Combined});
      tstop = std::chrono::high_resolution_clock::now();
      dt = tstop - tstart;
      std::wcout << "R[" << i << "] compact size: " << term_count(compact)
                 << " time: " << dt.count() << " s\n";

      // validated term counts of the compact triplet eom residual
      const auto n_compact = term_count(compact);
      if (N == 2 && type == EqnType::right && np == 2 && nh == 2) {
        if (i == 1) runtime_assert(n_compact == 42);
        if (i == 2) runtime_assert(n_compact == 135);
      }
      if (N == 3 && type == EqnType::right && np == 3 && nh == 3) {
        if (i == 1) runtime_assert(n_compact == 54);
        if (i == 2) runtime_assert(n_compact == 230);
        if (i == 3) runtime_assert(n_compact == 429);
      }

      // bare-TE residual (doubles only): rebuilding the full residual from it
      // via Omega = te + (1/4)(bra_swap(te) + ket_swap(te)) must be exact
      const auto ext_idxs = external_indices(eqvec[i]);
      if (ext_idxs.size() == 2 && N <= 2) {
        auto te = closed_shell_CC_triplet_spintrace(
            eqvec[i],
            {.compact = false, .residual = TripletResidualKind::BareTE});
        simplify(te);

        const Index b0 = get_bra_idx(ext_idxs.at(0));
        const Index b1 = get_bra_idx(ext_idxs.at(1));
        const Index k0 = get_ket_idx(ext_idxs.at(0));
        const Index k1 = get_ket_idx(ext_idxs.at(1));
        const container::map<Index, Index> bra_swap{{b0, b1}, {b1, b0}};
        const container::map<Index, Index> ket_swap{{k0, k1}, {k1, k0}};

        ExprPtr diff =
            st->clone() - (te->clone() + ex<Constant>(ratio(1, 4)) *
                                             (transform_expr(te, bra_swap) +
                                              transform_expr(te, ket_swap)));
        canonicalize(diff);
        simplify(diff);
        std::wcout << "R[" << i
                   << "] bare-TE reconstruction - full: " << term_count(diff)
                   << " terms (expect 0)\n";
        runtime_assert(term_count(diff) == 0);
      }

      if (print) {
        std::wcout << "\n R[" << i << "] equations:\n"
                   << to_latex_align(st, 20, 1) << "\n";
      }
    }
    timer_pool.stop(N + 16);

    std::wcout << "\nClosed-shell triplet spintracing completed in "
               << timer_pool.read(N + 16) << " s\n";
  }
};
}  // namespace

int main(int argc, char* argv[]) {
  std::wcout.precision(std::numeric_limits<double>::max_digits10);
  std::wcerr.precision(std::numeric_limits<double>::max_digits10);
  sequant::set_locale();

  std::cout << "SeQuant revision: " << sequant::git_revision() << "\n";
  std::cout << "Number of threads: " << sequant::num_threads() << "\n\n";

#ifndef NDEBUG
  constexpr size_t DEFAULT_NMAX = 2;
#else
  constexpr size_t DEFAULT_NMAX = 3;
#endif

  // command line arguments:
  //   argv[1]: NMAX (CC rank)
  //   argv[2]: excitation manifold (e.g. "1h1p", "2h2p"). Must be EE
  //   argv[3]: equation type ("R" or "L")
  //   argv[4]: "print" or "noprint"
  const size_t NMAX = argc > 1 ? std::stoi(argv[1]) : DEFAULT_NMAX;
  SEQUANT_ASSERT(NMAX > 0 && "Invalid NMAX");
  const std::string exc_manifold =
      argc > 2 ? argv[2]
               : (std::to_string(NMAX) + "h" + std::to_string(NMAX) + "p");
  SEQUANT_ASSERT(!exc_manifold.empty() && "Invalid excitation manifold");
  const std::string eqn_type = argc > 3 ? argv[3] : "R";
  const std::string print_str = argc > 4 ? argv[4] : "noprint";
  const bool print = print_str == "print";

  sequant::set_default_context(sequant::Context(
      {.index_space_registry_shared_ptr = make_min_sr_spaces(),
       .vacuum = Vacuum::SingleProduct,
       .canonicalization_options = CanonicalizeOptions().copy_and_set(
           CanonicalizationMethod::Complete)}));
  mbpt::set_default_mbpt_context(
      {.op_registry_ptr = mbpt::make_minimal_registry()});

  Logger::instance().wick_stats = false;

  compute_eomcc_closedshell_triplet{NMAX, exc_manifold,
                                    str2type.at(eqn_type)}(print);
}
