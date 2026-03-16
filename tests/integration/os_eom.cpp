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
    // generate so EOM eqs
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

    if (print) std::wcout << "\n";
    // to_latex_align(eqvec[i], 20, 1) << "\n";

    std::vector<std::vector<ExprPtr>> os_st_eom;
    std::vector<container::svector<container::svector<SlottedIndex>>>
        ext_index_groups_from_A;
    std::vector<size_t> n_paired_from_A;

    for (size_t i = 0; i < eqvec.size(); ++i) {
      if (eqvec[i] == nullptr) continue;

      Tensor A = eqvec[i].is<Sum>() ? eqvec[i]->at(0)->at(0)->as<Tensor>()
                                    : eqvec[i]->at(0)->as<Tensor>();
      std::wcout << "\n R[" << i << "]: A tensor (first tensor): \n"
                 << to_latex_align(ex<Tensor>(A)) << "\n";
      // std::wcout << "before spintracing: " << to_latex_align(eqvec[i], 20, 1)
      // << "\n";

      ext_index_groups_from_A.push_back(external_indices(A));
      n_paired_from_A.push_back(std::min(A.bra_rank(), A.ket_rank()));

      std::wcout << "R[" << i << "] has " << eqvec[i].size() << " terms\n";
      os_st_eom.push_back(open_shell_CC_spintrace(eqvec[i]));
      // std::wcout << to_latex_align(os_st_eom[i][0], 20, 1) << "\n";
      // std::wcout << "how many elements does it have? " << os_st_eom[i].size()
      // << "\n";
    }

    const auto alpha_qns = IndexSpace::QuantumNumbers{1};
    for (size_t i = 0; i < os_st_eom.size(); ++i) {
      const auto& spin_cases = os_st_eom[i];
      const auto& ext_from_A = ext_index_groups_from_A[i];
      const auto n_groups = ext_from_A.size();
      const size_t n_paired = n_paired_from_A[i];
      const size_t n_half = n_groups - n_paired;
      const bool orbit_encoding =
          (spin_cases.size() == (n_paired + 1) * (n_half + 1));

      for (size_t sc = 0; sc < spin_cases.size(); ++sc) {
        container::svector<int> spins(n_groups, 0);
        if (orbit_encoding) {
          const size_t i_f = sc % (n_paired + 1);
          const size_t i_h = sc / (n_paired + 1);
          if (i_f > 0)
            std::fill(spins.begin() + (n_paired - i_f),
                      spins.begin() + n_paired, 1);
          if (i_h > 0)
            std::fill(spins.begin() + (n_groups - i_h), spins.end(), 1);
        } else if (sc > 0) {
          std::fill(spins.end() - sc, spins.end(), 1);
        }

        std::wstring spin_label;
        for (size_t g = 0; g < ext_from_A.size(); ++g)
          for (const auto& slotted : ext_from_A[g]) {
            const auto spin_idx = spins[g] == 0
                                      ? make_spinalpha(slotted.index())
                                      : make_spinbeta(slotted.index());
            spin_label += (spin_idx.space().qns() == alpha_qns) ? L"α" : L"β";
          }

        std::wcout << "\n R[" << i << "] case [" << sc << "] spin=("
                   << spin_label << ") : " << spin_cases[sc]->size()
                   << " terms\n";
        // std::wcout << to_latex_align(spin_cases[sc]) << std::endl;

        // std::wcout << "  ext_index_groups (" << ext_from_A.size()
        //            << " groups):\n";
        // for (size_t g = 0; g < ext_from_A.size(); ++g) {
        //   std::wcout << "    group[" << g << "]: ";
        //   for (const auto& slotted : ext_from_A[g]) {
        //     const auto spin_idx = spins[g] == 0
        //                               ? make_spinalpha(slotted.index())
        //                               : make_spinbeta(slotted.index());
        //     std::wcout << spin_idx.to_latex() << "("
        //                << ((spin_idx.space().qns() == alpha_qns) ? L"α" :
        //                L"β")
        //                << ") ";
        //   }
        //   std::wcout << "\n";
        // }
      }
    }
  }
};

class compute_all_openshell {
  size_t NMAX;
  std::string manifold;
  EqnType type;

 public:
  compute_all_openshell(size_t nmax, const std::string manifold,
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

  std::cout << "SeQuant revision: " << sequant::git_revision() << "\n";
  std::cout << "Number of threads: " << sequant::num_threads() << "\n\n";

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
  mbpt::set_default_mbpt_context(
      {.op_registry_ptr = mbpt::make_minimal_registry()});
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  // change to true to print stats
  Logger::instance().wick_stats = false;

  // call the compute_all function here
  // compute_all_openshell{NMAX, exc_manifold, str2type.at(eqn_type)}(print);
  compute_eomcc_openshell{NMAX, exc_manifold, str2type.at(eqn_type)}(print);
}
