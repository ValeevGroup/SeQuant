#include <SeQuant/version.hpp>

#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/timer.hpp>
#include <SeQuant/domain/mbpt/biorthogonalization.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <chrono>

using namespace sequant;
using namespace sequant::mbpt;

namespace {
#define runtime_assert(tf)                                           \
  if (!(tf)) {                                                       \
    std::ostringstream oss;                                          \
    oss << "failed assert at line " << __LINE__                      \
        << " in closed-shell triplet equation-of-motion CC example"; \
    throw std::runtime_error(oss.str().c_str());                     \
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

std::size_t count_distinct_hashes(ExprPtr expr) {
  if (expr->is<Sum>()) {
    for (auto& term : *expr) {
      if (term->is<Product>())
        term = remove_tensor(term.as_shared_ptr<Product>(),
                             reserved::symm_label());
    }
  }
  canonicalize(expr);
  simplify(expr);

  if (!expr->is<Sum>()) return expr->is<Constant>() ? 0 : 1;

  container::set<std::size_t> hashes;
  for (const auto& term : *expr) {
    if (!term->is<Product>()) continue;
    auto product = term.as_shared_ptr<Product>();
    sequant::TensorNetwork tn(*product);
    auto hash =
        tn.canonicalize_slots(TensorCanonicalizer::cardinal_tensor_labels())
            .hash_value();
    hashes.insert(hash);
  }
  return hashes.size();
}

container::svector<container::svector<Index>> unwrap_ext_groups(
    const container::svector<container::svector<SlottedIndex>>& ext_idxs) {
  container::svector<container::svector<Index>> ext_groups;
  for (const auto& g : ext_idxs) {
    container::svector<Index> grp;
    grp.reserve(g.size());
    for (const auto& s : g) grp.push_back(s.index());
    ext_groups.push_back(std::move(grp));
  }
  return ext_groups;
}

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
      throw std::runtime_error(
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

    // Per-equation: os_eom.cpp sector/singlet checks, then closed-shell
    // triplet.
    std::wcout
        << "\nClosed-shell triplet EOM-CC (os_eom-aligned diagnostics):\n";

    std::vector<ExprPtr> cs_st_eom(eqvec.size());

    timer_pool.start(N + 16);
    for (size_t i = 0; i < eqvec.size(); ++i) {
      if (eqvec[i] == nullptr) {
        cs_st_eom[i] = nullptr;
        continue;
      }
      if (eqvec[i]->is<Constant>() || eqvec[i]->is<Variable>()) {
        cs_st_eom[i] = nullptr;
        continue;
      }

      Tensor A = eqvec[i]->is<Sum>() ? eqvec[i]->at(0)->at(0)->as<Tensor>()
                                     : eqvec[i]->at(0)->as<Tensor>();
      std::wcout << "\n R[" << i << "]: A tensor (first tensor): \n"
                 << to_latex_align(ex<Tensor>(A)) << "\n";
      std::wcout << "R[" << i << "] has " << eqvec[i]->size() << " terms\n";

      const auto ext_idxs = external_indices(eqvec[i]);
      const auto ext_groups = unwrap_ext_groups(ext_idxs);

      // ----- os_eom main path: open_shell_CC_spintrace_by_sector + sum -----
      auto os_sectors = open_shell_CC_spintrace_by_sector(eqvec[i]);
      const size_t n_cases = os_sectors.size();
      std::wcout << "number of spin cases " << n_cases << "\n";
      SEQUANT_ASSERT(n_cases >= 2);

      auto summed_spinfree = std::make_shared<Sum>();
      for (size_t sc = 0; sc < n_cases; ++sc) {
        ExprPtr stripped = os_sectors[sc]->clone();
        expand(stripped);
        std::wcout << "sc" << sc << "\n";
        std::wcout << "spin-free sector R[" << i << "] has " << stripped->size()
                   << " terms\n";
        canonicalize(stripped);
        summed_spinfree->append(stripped);
      }

      ExprPtr summed = summed_spinfree;
      simplify(summed);
      std::wcout << "R[" << i
                 << "] open-shell sum (Ŝ NOT expanded): " << summed->size()
                 << " terms\n";

      auto singlet_ref = closed_shell_CC_spintrace(
          eqvec[i], {.method = BiorthogonalizationMethod::V2});
      simplify(singlet_ref);
      std::wcout << "R[" << i << "] reference closed-shell (CC spintrace): "
                 << singlet_ref->size() << " terms\n";

      ExprPtr os_singlet =
          remove_tensor(summed->clone(), reserved::symm_label());
      simplify(os_singlet);
      ExprPtr singlet_diff = os_singlet - singlet_ref;
      canonicalize(singlet_diff);
      simplify(singlet_diff);
      std::wcout << "R[" << i << "] (open-shell singlet) - (closed-shell ref): "
                 << singlet_diff->size() << " terms\n";

      // ----- os_eom independent path: spintrace_by_sector -----------------
      {
        auto sectors = spintrace_by_sector(eqvec[i], ext_groups);
        std::wcout << L"\n========== R[" << i << L"] per-sector spin trace ("
                   << sectors.size() << L" sectors) ==========\n";
        for (auto& [label, sec] : sectors) {
          std::wcout << L"  sector " << label << L": " << sec->size()
                     << L" terms, " << count_distinct_hashes(sec->clone())
                     << L" distinct hashes\n";
        }

        auto sector_total = std::make_shared<Sum>();
        for (auto& [label, sec] : sectors) sector_total->append(sec->clone());

        ExprPtr so = eqvec[i]->clone();
        so->visit(
            [](ExprPtr& n) {
              if (n->is<Tensor>()) n->as<Tensor>().reset_tags();
            },
            /*atoms_only=*/true);
        ExprPtr generic =
            spintrace(so, ext_groups, /*spinfree_index_spaces=*/false);
        canonicalize(generic);
        simplify(generic);
        generic = remove_spin_with_relabel(generic);
        canonicalize(generic);
        simplify(generic);

        ExprPtr sector_generic_diff = ExprPtr(sector_total) - generic;
        canonicalize(sector_generic_diff);
        simplify(sector_generic_diff);
        std::wcout << L"  R[" << i << L"] (Σ sectors) - (generic) : "
                   << sector_generic_diff->size() << L" terms\n";
      }

      // ----- triplet: first − last with triplet_R (matches production) ----
      auto triplet_sectors =
          spintrace_by_sector(eqvec[i], ext_groups, /*triplet_R=*/true);
      ExprPtr T_diag = triplet_sectors.front().second->clone() -
                       triplet_sectors.back().second->clone();
      canonicalize(T_diag);
      simplify(T_diag);
      std::wcout << "R[" << i
                 << "] triplet (triplet_R, first - last): " << T_diag->size()
                 << " terms\n";

      try {
        const auto tstart = std::chrono::high_resolution_clock::now();
        auto st = closed_shell_EOM_triplet_spintrace(
            eqvec[i], {.method = BiorthogonalizationMethod::V2});
        const auto tstop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> dt = tstop - tstart;

        ExprPtr T_diff = st - T_diag;
        canonicalize(T_diff);
        simplify(T_diff);
        std::wcout << "R[" << i
                   << "] closed_shell_EOM_triplet_spintrace: " << st->size()
                   << " terms, time: " << dt.count() << " s\n";
        std::wcout << "R[" << i
                   << "] (production triplet) - (diag first-last): "
                   << T_diff->size() << " terms\n";

        if (print) {
          std::wcout << "\n open-shell singlet sum:\n"
                     << to_latex_align(summed, 20, 1) << "\n";
          std::wcout << "\n triplet (production):\n"
                     << to_latex_align(st, 20, 1) << "\n";
        }

        cs_st_eom[i] = std::move(st);
      } catch (const std::exception& e) {
        std::wcerr << "\n*** FAILED on R[" << i << "] (CC rank=" << N
                   << ", manifold=" << sequant::toUtf16(manifold) << ") ***\n";
        std::wcerr << "Exception: " << e.what() << "\n";
        std::wcerr << "Input spin-orbital expression:\n"
                   << to_latex_align(eqvec[i], 20, 1) << "\n";
        std::wcerr.flush();
        throw;
      }
    }
    timer_pool.stop(N + 16);

    std::wcout << "\nClosed-shell triplet spintracing completed in "
               << timer_pool.read(N + 16) << " s\n";
  }
};

class compute_all_closedshell_triplet {
  size_t NMAX;
  std::string manifold;
  EqnType type;

 public:
  compute_all_closedshell_triplet(size_t nmax, const std::string manifold,
                                  EqnType t = EqnType::right)
      : NMAX(nmax), manifold(manifold), type(t) {}

  void operator()(bool print = false) {
    for (size_t N = 1; N <= NMAX; ++N) {
      std::vector<std::string> manifold_vec;
      auto [Nh, Np] = parse_excitation_manifold(manifold);
      // triplet only supports particle-conserving manifolds; collapse to the
      // diagonal sequence (NhNp, ..., 1h1p)
      while (Nh > 0 && Np > 0) {
        manifold_vec.push_back(std::to_string(Nh) + "h" + std::to_string(Np) +
                               "p");
        Nh--;
        Np--;
      }
      for (auto it = manifold_vec.rbegin(); it != manifold_vec.rend(); ++it) {
        compute_eomcc_closedshell_triplet{N, *it, type}(print);
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

  // compute_all_closedshell_triplet{NMAX, exc_manifold,
  //                                 str2type.at(eqn_type)}(print);
  compute_eomcc_closedshell_triplet{NMAX, exc_manifold,
                                    str2type.at(eqn_type)}(print);
}
