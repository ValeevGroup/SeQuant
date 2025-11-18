#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/permutation.hpp>
#include <SeQuant/core/utility/timer.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <algorithm>
#include <functional>
#include <iterator>
#include <stdexcept>
#include <string_view>
#include <unordered_map>

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
    // holes only
    result.first = std::stoi(str.substr(0, h_pos));
    result.second = 0;
  } else if (p_pos != std::string::npos && h_pos == std::string::npos) {
    // particles only
    result.first = 0;
    result.second = std::stoi(str.substr(0, p_pos));
  } else {
    // hp/ph cases
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

class compute_eomcc_openshell {
  size_t N, np, nh;
  std::string manifold;
  EqnType type;

 public:
  compute_eomcc_openshell(size_t n, const std::string& exc_manifold,
                          EqnType t = EqnType::right)
      : N(n), manifold(exc_manifold), type(t) {
    std::tie(nh, np) = parse_excitation_manifold(manifold);
  }

  void operator()(bool print) {
    SEQUANT_ASSERT(get_default_context().spbasis() == SPBasis::Spinor);
    // generate so EOM eqs like how Ajay constructed them
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
               << ", manifold=" << sequant::to_wstring(manifold) << "]"
               << " computed in " << timer_pool.read(N) << " s\n";

    if (print) std::wcout << "\n";
    // to_latex_align(eqvec[i], 20, 1) << "\n";

    // open-shell spin-tracing
    std::vector<std::vector<ExprPtr>> spintraced_results;

    for (size_t i = 0; i < eqvec.size(); ++i) {
      if (eqvec[i] == nullptr) continue;
      //     // spintraced_results.push_back({});
      std::wcout << "R[" << i << "] has " << eqvec[i].size() << " terms\n";

      auto expr = eqvec[i];
      Tensor A = expr->at(0)->at(0)->as<Tensor>();
      SEQUANT_ASSERT(A.label() == L"A");
      size_t rank = A.bra_rank();
      runtime_assert(rank == A.ket_rank());
      auto ext_idxs = external_indices(A);

      auto P_vec = open_shell_P_op_vector(A);
      SEQUANT_ASSERT(P_vec.size() == rank + 1);

      auto A_vec = open_shell_A_op(A);
      SEQUANT_ASSERT(A_vec.size() == rank + 1);
      std::vector<Sum> concat_terms(rank + 1);

      std::wcout << "processing " << expr->size() << " terms:\n";
      for (auto& product_term : *expr) {
        // remove A from this term, like how we did for S in v2
        auto term = remove_tensor(product_term.as_shared_ptr<Product>(), L"A");

        std::vector<ExprPtr> os_st(rank + 1);

        for (size_t s = 0; s <= rank; ++s) {
          os_st[s] = P_vec[s] * term;
          expand(os_st[s]);
          os_st[s] = expand_P_op(os_st[s]);

          // os_st[s] = open_shell_spintrace(os_st[s], ext_idxs, s)[0];
          auto spintraced = open_shell_spintrace(os_st[s], ext_idxs, s);
          if (spintraced.empty()) {
            throw std::runtime_error(
                "open_shell_spintrace returned empty vector");
          }
          os_st[s] = spintraced[0];

          if (rank > 2) {
            os_st[s] = A_vec[s] * os_st[s];
            simplify(os_st[s]);
            os_st[s] = remove_tensor(os_st[s], L"A");
          }
        }

        for (size_t s = 0; s <= rank; ++s) {
          concat_terms[s].append(os_st[s]);
        }
      }

      std::vector<ExprPtr> level_spin_cases;
      for (auto& spin_case : concat_terms) {
        auto ptr = ex<Sum>(spin_case);
        simplify(ptr);
        canonicalize(ptr);
        level_spin_cases.push_back(ptr);
      }
      spintraced_results.push_back(level_spin_cases);
    }

    for (size_t i = 0; i < spintraced_results.size(); ++i) {
      const auto& spin_cases = spintraced_results[i];

      if (spin_cases.empty()) {
        std::wcout << type2wstr.at(type) << i << ": empty\n";
        continue;
      }

      if (eqvec[i] != nullptr) {
        std::wcout << "original (spin-orbital): " << eqvec[i]->size()
                   << " terms\n";
      }
      std::wcout << "spin-traced cases: " << spin_cases.size() << "\n";

      for (size_t sc = 0; sc < spin_cases.size(); ++sc) {
        if (spin_cases[sc] == nullptr) {
          std::wcout << "case " << sc << ": null\n";
          continue;
        }
        std::wcout << "case " << sc << " : " << spin_cases[sc]->size()
                   << " terms\n";
        // std::wcout << "check spintraced eqs: "
        // <<to_latex_align(spin_cases[sc], 20, 1) << "\n";
      }
    }
  }
};

class compute_all_openshell {
  size_t NMAX;
  std::string manifold;
  EqnType type;

 public:
  compute_all_openshell(size_t nmax, const std::string& manifold,
                        EqnType t = EqnType::right)
      : NMAX(nmax), manifold(manifold), type(t) {}

  void operator()(bool print = false) {
    for (size_t N = 1; N <= NMAX; ++N) {
      std::vector<std::string> manifold_vec;
      auto [Nh, Np] = parse_excitation_manifold(manifold);
      // generate all possible manifolds
      while (Nh > 0 || Np > 0) {
        if (Nh == 0 && Np == 0) break;
        manifold_vec.push_back(std::to_string(Nh) + "h" + std::to_string(Np) +
                               "p");
        if (Nh == 0 || Np == 0) break;
        Nh--;
        Np--;
      }
      for (auto it = manifold_vec.rbegin(); it != manifold_vec.rend(); ++it) {
        compute_eomcc_openshell{N, *it, type}(print);
      }
    }
  }
};
}  // namespace

int main(int argc, char* argv[]) {
  std::wcout.precision(std::numeric_limits<double>::max_digits10);
  std::wcerr.precision(std::numeric_limits<double>::max_digits10);
  sequant::set_locale();

#ifndef NDEBUG
  constexpr size_t DEFAULT_NMAX = 3;
#else
  constexpr size_t DEFAULT_NMAX = 4;
#endif

  // read command line arguments
  const size_t NMAX = argc > 1 ? std::stoi(argv[1]) : DEFAULT_NMAX;
  SEQUANT_ASSERT(NMAX > 0 && "Invalid NMAX");
  const std::string exc_manifold =
      argc > 2 ? argv[2]
               : (std::to_string(NMAX) + "h" + std::to_string(NMAX) + "p");
  SEQUANT_ASSERT(!exc_manifold.empty() && "Invalid excitation manifold");
  const std::string eqn_type = argc > 3 ? argv[3] : "R";
  const std::string print_str = argc > 4 ? argv[4] : "noprint";
  const bool print = print_str == "print";

  sequant::detail::OpIdRegistrar op_id_registrar;

  sequant::set_default_context(sequant::Context(
      {.index_space_registry_shared_ptr = make_min_sr_spaces(),
       .vacuum = Vacuum::SingleProduct,
       .canonicalization_options = CanonicalizeOptions().copy_and_set(
           CanonicalizationMethod::Complete)}));

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  Logger::instance().wick_stats = false;

  // call the compute_all function here
  compute_all_openshell{NMAX, exc_manifold, str2type.at(eqn_type)}(print);
}
