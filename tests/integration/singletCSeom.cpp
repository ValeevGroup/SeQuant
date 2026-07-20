#include <SeQuant/version.hpp>

#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/timer.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <range/v3/algorithm/transform.hpp>

#include <chrono>

using namespace sequant;
using namespace sequant::mbpt;

namespace {
#define runtime_assert(tf)                                   \
  if (!(tf)) {                                               \
    std::ostringstream oss;                                  \
    oss << "failed assert at line " << __LINE__              \
        << " in equation-of-motion coupled cluster example"; \
    throw std::runtime_error(oss.str().c_str());             \
  }

TimerPool<32> timer_pool;

std::pair<size_t, size_t> parse_excitation_manifold(std::string& str) {
  std::pair<size_t, size_t> result;

  ranges::transform(str, str.begin(), ::tolower);
  const auto h_pos = str.find('h');
  const auto p_pos = str.find('p');

  if (h_pos == std::string::npos && p_pos == std::string::npos) {
    throw std::runtime_error(
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
    throw std::runtime_error(
        "Invalid excitation manifold: both particle and hole ranks cannot be "
        "zero");

  return result;
}

enum class EqnType { left, right };

inline const container::map<std::string, EqnType> str2type = {
    {"L", EqnType::left}, {"R", EqnType::right}};

inline const container::map<EqnType, std::wstring> type2wstr = {
    {EqnType::left, L"L"}, {EqnType::right, L"R"}};

class compute_eomcc_closedshell {
  size_t N, np, nh;
  std::string manifold;
  EqnType type;
  BiorthogonalizationMethod biorth_method;

 public:
  compute_eomcc_closedshell(
      size_t n, const std::string& exc_manifold, EqnType t = EqnType::right,
      BiorthogonalizationMethod bm = BiorthogonalizationMethod::V2)
      : N(n), manifold(exc_manifold), type(t), biorth_method(bm) {
    std::tie(nh, np) = parse_excitation_manifold(manifold);
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

    // report spin-orbital sizes (analogous to the ground-state CC test)
    for (size_t i = 0; i < eqvec.size(); ++i) {
      if (eqvec[i] == nullptr) continue;
      std::wcout << "Spin-orbital R[" << i << "] size: " << eqvec[i]->size()
                 << "\n";
    }

    // closed-shell spintrace: one ExprPtr per equation (a Sum of terms),
    // mirroring the ground-state CC pattern
    const std::wstring biorth_label =
        (biorth_method == BiorthogonalizationMethod::V1) ? L"v1" : L"v2";
    std::wcout << "\nClosed-shell EOM-CC spintrace with biorthogonal "
               << biorth_label << " transformation:\n";

    std::vector<ExprPtr> cs_st_eom(eqvec.size());

    timer_pool.start(N + 16);  // distinct slot for spintracing timer
    for (size_t i = 0; i < eqvec.size(); ++i) {
      if (eqvec[i] == nullptr) {
        cs_st_eom[i] = nullptr;
        continue;
      }

      const auto tstart = std::chrono::high_resolution_clock::now();
      auto st = closed_shell_CC_spintrace(eqvec[i], {.method = biorth_method});
      // if (i > 1) {
      //   st = S_maps(st);
      //   simplify(st);
      // }
      for (auto& term : *st) {
        if (term->is<Product>())
          term = remove_tensor(term.as_shared_ptr<Product>(),
                               reserved::symm_label());
      }
      const auto tstop = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> dt = tstop - tstart;

      std::wcout << "R[" << i << "] size: " << st->size()
                 << " time: " << dt.count() << " s\n";

      if (print) {
        std::wcout << "\n R[" << i << "] equations:\n"
                   << to_latex_align(st, 0, 4) << "\n";
      }

      cs_st_eom[i] = std::move(st);
    }
    timer_pool.stop(N + 16);

    std::wcout << "\nClosed-shell spintracing completed in "
               << timer_pool.read(N + 16) << " s\n";
  }
};

class compute_all_closedshell {
  size_t NMAX;
  std::string manifold;
  EqnType type;
  BiorthogonalizationMethod biorth_method;

 public:
  compute_all_closedshell(
      size_t nmax, const std::string manifold, EqnType t = EqnType::right,
      BiorthogonalizationMethod bm = BiorthogonalizationMethod::V2)
      : NMAX(nmax), manifold(manifold), type(t), biorth_method(bm) {}

  void operator()(bool print = false) {
    for (size_t N = 1; N <= NMAX; ++N) {
      std::vector<std::string> manifold_vec;
      auto [Nh, Np] = parse_excitation_manifold(manifold);
      while (Nh > 0 || Np > 0) {
        if (Nh == 0 && Np == 0) break;
        manifold_vec.push_back(std::to_string(Nh) + "h" + std::to_string(Np) +
                               "p");
        if (Nh == 0 || Np == 0) break;
        Nh--;
        Np--;
      }
      for (auto it = manifold_vec.rbegin(); it != manifold_vec.rend(); ++it) {
        compute_eomcc_closedshell{N, *it, type, biorth_method}(print);
      }
    }
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
  constexpr size_t DEFAULT_NMAX = 3;
#else
  constexpr size_t DEFAULT_NMAX = 4;
#endif

  // command line arguments:
  //   argv[1]: NMAX (CC rank)
  //   argv[2]: excitation manifold (e.g. "1h1p", "2h1p")
  //   argv[3]: equation type ("R" or "L")
  //   argv[4]: biorthogonalization method ("V1" or "V2")
  //   argv[5]: "print" or "noprint"
  const size_t NMAX = argc > 1 ? std::stoi(argv[1]) : DEFAULT_NMAX;
  SEQUANT_ASSERT(NMAX > 0 && "Invalid NMAX");
  const std::string exc_manifold =
      argc > 2 ? argv[2]
               : (std::to_string(NMAX) + "h" + std::to_string(NMAX) + "p");
  SEQUANT_ASSERT(!exc_manifold.empty() && "Invalid excitation manifold");
  const std::string eqn_type = argc > 3 ? argv[3] : "R";
  const std::string biorth_str = argc > 4 ? argv[4] : "V2";
  const std::string print_str = argc > 5 ? argv[5] : "noprint";
  const bool print = print_str == "print";

  BiorthogonalizationMethod biorth_method;
  if (biorth_str == "V1") {
    biorth_method = BiorthogonalizationMethod::V1;
  } else if (biorth_str == "V2") {
    biorth_method = BiorthogonalizationMethod::V2;
  } else {
    throw std::runtime_error("Invalid biorthogonalization method: " +
                             biorth_str + " (expected V1 or V2)");
  }

  sequant::detail::OpIdRegistrar op_id_registrar;

  sequant::set_default_context(sequant::Context(
      {.index_space_registry_shared_ptr = make_min_sr_spaces(),
       .vacuum = Vacuum::SingleProduct,
       .canonicalization_options = CanonicalizeOptions().copy_and_set(
           CanonicalizationMethod::Complete)}));
  mbpt::set_default_mbpt_context(
      {.op_registry_ptr = mbpt::make_minimal_registry()});
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  Logger::instance().wick_stats = false;

  // compute_all_closedshell{NMAX, exc_manifold, str2type.at(eqn_type),
  //                         biorth_method}(print);
  compute_eomcc_closedshell{NMAX, exc_manifold, str2type.at(eqn_type),
                            biorth_method}(print);
}
