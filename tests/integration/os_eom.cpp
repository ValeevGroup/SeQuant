#include <SeQuant/version.hpp>

#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/timer.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <range/v3/algorithm/transform.hpp>
#include "SeQuant/domain/mbpt/biorthogonalization.hpp"

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

std::size_t count_distinct_hashes(ExprPtr expr) {
  // Strip symmetrizer so hashing matches the closed-shell / WK filter path
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

// ===========================================================================
//  Diagnostic harness: brute-force ground-truth comparison
// ---------------------------------------------------------------------------
//  Implements the approach your PI described (transcript p.8): run the generic
//  spintrace() (loops over all 2^n spin cases, drops non-conserving terms via
//  can_expand(), assigns ONE spin per index per bitstring), then manually
//  remove all spins with plain remove_spin() -- NOT remove_spin_with_relabel().
//
//  Because every index in a product already has a distinct base label, plain
//  remove_spin() never hits the alpha/beta -> spinfree collision that
//  remove_spin_with_relabel exists to patch, so this path never relabels
//  internals. It is therefore the clean ground truth.
//
//  This function only INSTRUMENTS and PRINTS; it does not try to "fix" the
//  discrepancy. Read the differing terms by hand.
//
//  CAVEAT: generic spintrace() also sums over the EXTERNAL-index spin cases,
//  whereas closed_shell_(CC_)spintrace fixes external spins and cycle-counts a
//  single representative. So (generic) - (closed-shell) may differ by an
//  overall scalar even when structurally identical -- read that comparison for
//  structure (distinct-hash count, which monomials appear), not an exact zero.
//  The decisive test is (generic) - (open-shell sum): both sum over spin
//  cases, so any nonzero diff localizes the issue to your open-shell path.
// ===========================================================================
void diagnose_against_generic_spintrace(const ExprPtr& spinorbital_rhs,
                                        size_t i,
                                        const ExprPtr& closed_shell_ref,
                                        const ExprPtr& open_shell_sum/*,
                                        bool dump_full = true*/) {
  if (!spinorbital_rhs) return;

  auto banner = [&](const std::wstring& s) {
    std::wcout << L"\n=== R[" << i << L"] " << s << L" ==========\n";
  };

  // External index groups. Use the expression-level call so it matches what
  // spintrace() determines internally (its internal assert compares against
  // external_indices(expression)).
  // Work on a clone, then clear all index tags. Earlier passes (open/closed-
  // shell spintrace, external_indices, canonicalize) leave int tags on the
  // shared expression; generic spintrace()'s internal phase tagging would then
  // trip SEQUANT_ASSERT(!tag_.has_value()).
  ExprPtr so = spinorbital_rhs->clone();
  const auto ext_groups = external_indices(so);
  so->visit(
      [](ExprPtr& node) {
        if (node->is<Tensor>()) node->as<Tensor>().reset_tags();
      },
      /*atoms_only=*/true);

  // -- (1) Generic brute-force spintrace -----------------------------------
  // Keep spin labels (spinfree_index_spaces = false) so we remove them
  // ourselves
  banner(L"generic spintrace (brute force, 2^n spin cases)");
  ExprPtr generic = spintrace(so, ext_groups,
                              /*spinfree_index_spaces=*/false);
  canonicalize(generic);
  simplify(generic);
  std::wcout << L"  with spin labels : " << generic->size() << L" terms\n";

  // Manual spin removal -- plain remove_spin(), NOT remove_spin_with_relabel().
  generic = remove_spin_with_relabel(generic);
  canonicalize(generic);
  simplify(generic);
  generic = biorthogonal_transform_pre_nnsproject(generic, ext_groups);
  std::wcout << L"  after remove_spin: " << generic->size() << L" terms, "
             << count_distinct_hashes(generic->clone())
             << L" distinct hashes\n";
  // if (dump_full)
  //   std::wcout << L"  --- brute-force result ---\n"
  //              << to_latex_align(generic, 0, 4) << L"\n";

  // Strip any leftover symmetrizer so all three expressions are bare sums.
  auto strip = [](const ExprPtr& e) {
    ExprPtr r = e->clone();
    canonicalize(r);
    simplify(r);
    return r;
  };

  auto report_diff = [&](const std::wstring& label, const ExprPtr& other) {
    if (!other) return;
    ExprPtr o = strip(other);
    std::wcout << L"  " << label << L" : " << o->size() << L" terms, "
               << count_distinct_hashes(o->clone()) << L" distinct hashes\n";
    ExprPtr d = generic - o;
    canonicalize(d);
    simplify(d);
    std::wcout << L"  (generic) - (" << label << L") : " << d->size()
               << L" terms\n";
    // if (d->size() != 0) std::wcout << to_latex_align(d, 0, 4) << L"\n";
  };

  // -- (2) Diff vs closed-shell path (structure only; see CAVEAT) ----------
  banner(L"diff vs closed-shell path");
  report_diff(L"closed-shell", closed_shell_ref);

  // -- (3) Diff vs open-shell sum path (THE decisive comparison) -----------
  banner(L"diff vs open-shell sum path (decisive)");
  report_diff(L"open-shell sum", open_shell_sum);
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
    std::vector<ExprPtr> os_spinfree_summed_eom;
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

      std::wcout << "R[" << i << "] has " << eqvec[i].size() << " terms\n";

      // Singlet/triplet sector construction (and all the comparison/diagnostic
      // blocks below) only make sense for particle-conserving manifolds where
      // the resulting state has well-defined total spin (S=0/S=1).
      if (np != nh) continue;

      ext_index_groups_from_A.push_back(external_indices(A));
      n_paired_from_A.push_back(std::min(A.bra_rank(), A.ket_rank()));

      os_st_eom.push_back(open_shell_CC_spintrace_by_sector(eqvec[i]));
      auto summed_spinfree = std::make_shared<Sum>();

      const size_t n_cases = os_st_eom.back().size();
      std::wcout << "number of spin cases " << n_cases << "\n";

      for (size_t sc = 0; sc < n_cases; ++sc) {
        auto& spin_case_expr = os_st_eom.back()[sc];

        // open_shell_CC_spintrace_by_sector: already spin-free per sector.
        ExprPtr stripped = spin_case_expr->clone();
        expand(stripped);
        std::wcout << "sc" << sc << "\n";
        std::wcout << "spin-free sector R[" << i << "] has " << stripped.size()
                   << " terms\n";
        canonicalize(stripped);

        // const bool is_mixed = (n_cases == 3 && sc == 1);
        // if (is_mixed) {
        //   const auto nf = ex<Constant> (factorial(2));
        //   stripped = nf* stripped;
        //   simplify(stripped);
        // }

        // if (!is_mixed) {
        // ExprPtr swapped = swap_spin(stripped);
        // canonicalize(swapped);
        // simplify(swapped);
        summed_spinfree->append(stripped);  // αβ
        // summed_spinfree->append(swapped);          // βα
        // continue;  // skip the normal append below
        // }
        // std::wcout << "sub set size  R[" << i << "] has " << stripped.size()
        //            << " terms\n";
        // if (!is_mixed) {
        // summed_spinfree->append(stripped);
        // }
      }
      // ExprPtr sc0_result, sc2_result;
      //
      // for (size_t sc = 0; sc < n_cases; ++sc) {
      //   ExprPtr spin_case_expr = os_st_eom.back()[sc];
      //   canonicalize(spin_case_expr);
      //   simplify(spin_case_expr);
      //
      //   ExprPtr stripped = expand_antisymm(spin_case_expr);
      //   expand(stripped);
      //   canonicalize(stripped);
      //   simplify(stripped);
      //   stripped = remove_spin_with_relabel(stripped);
      //   canonicalize(stripped);
      //   simplify(stripped);
      //
      //   if (sc == 0) sc0_result = stripped;
      //   if (sc == n_cases - 1) sc2_result = stripped;
      //
      //   summed_spinfree->append(stripped);
      // }
      //
      // // subtract sc0 and sc2 to check consistency
      // if (sc0_result && sc2_result) {
      //   ExprPtr diff = sc0_result + ex<Constant>(-1) * sc2_result;
      //   simplify(diff);
      //   std::wcout << "sc0 - sc_last difference: " << diff->size()
      //              << " terms\n";
      //   std::wcout << to_latex_align(diff, 0, 4) << "\n";
      // }

      // os_spinfree_summed_eom.push_back(std::move(summed_spinfree));
      // expand_antisymm(os_spinfree_summed_eom.back());
      // // simplify(os_spinfree_summed_eom.back());

      // //ExprPtr summed = expand_antisymm(summed_spinfree);

      // The fixed open_shell_CC_spintrace attaches S internally (via
      // symmetrize_expr), so summed_spinfree already carries one S per term
      // for rank >= 2. No outer S multiplication is needed; we only call
      // S_maps to expand the existing S into particle permutations.
      ExprPtr summed = summed_spinfree;
      using ranges::views::transform;
      auto const ext_idxs = external_indices(summed);
      simplify(summed);
      summed = biorthogonal_transform_pre_nnsproject(summed, ext_idxs);

      // ExprPtr summed_smaps = summed->clone();
      // if (has_tensor(summed_smaps, reserved::symm_label())) {
      //   summed_smaps = S_maps(summed_smaps);
      //   simplify(summed_smaps);
      // }
      std::wcout << "R[" << i
                 << "] open-shell sum (Ŝ NOT expanded): " << summed->size()
                 << " terms\n";
      // std::wcout << "R[" << i << "] open-shell sum (Ŝ expanded via S_maps): "
      //            << summed_smaps->size() << " terms\n";

      os_spinfree_summed_eom.push_back(std::move(summed));

      auto singlet_ref = closed_shell_CC_spintrace(
          eqvec[i], {.method = BiorthogonalizationMethod::V2});
      // singlet_ref = remove_tensor(singlet_ref, reserved::symm_label());
      simplify(singlet_ref);
      std::wcout << "R[" << i << "] reference closed-shell (CC spintrace): "
                 << singlet_ref->size() << " terms\n";

      ExprPtr os_singlet = os_spinfree_summed_eom.back();
      simplify(os_singlet);

      // subtract
      ExprPtr diff = os_singlet - singlet_ref;
      canonicalize(diff);
      simplify(diff);

      std::wcout << "R[" << i << "] (open-shell singlet) - (closed-shell ref): "
                 << diff->size() << " terms\n";
      // if (diff->size() != 0)
      // std::wcout << to_latex_align(os_singlet, 0, 4) << "\n";

      // ----- PI-suggested brute-force ground-truth comparison -------------
      // Generic spintrace + manual remove_spin, diffed against both the
      // closed-shell path and your open-shell sum path. Read the (generic) -
      // (open-shell sum) block: it is the decisive comparison.
      diagnose_against_generic_spintrace(eqvec[i], i, singlet_ref,
                                         os_spinfree_summed_eom.back());

      // ----- INDEPENDENT PATH: per-sector spin trace ----------------------
      // Self-contained; does not touch os_st_eom / os_spinfree_summed_eom or
      // any open-shell state. Safe to delete without affecting the rest.
      {
        // external_indices() yields SlottedIndex groups; unwrap to Index.
        container::svector<container::svector<Index>> ext_groups;
        for (const auto& g : external_indices(eqvec[i])) {
          container::svector<Index> grp;
          grp.reserve(g.size());
          for (const auto& s : g) grp.push_back(s.index());
          ext_groups.push_back(std::move(grp));
        }

        auto sectors = spintrace_by_sector(eqvec[i], ext_groups);

        std::wcout << L"\n========== R[" << i << L"] per-sector spin trace ("
                   << sectors.size() << L" sectors) ==========\n";
        for (auto& [label, sec] : sectors) {
          std::wcout << L"  sector " << label << L": " << sec->size()
                     << L" terms, " << count_distinct_hashes(sec->clone())
                     << L" distinct hashes\n";
          // std::wcout << to_latex_align(sec, 0, 4) << L"\n";
        }

        // Consistency check: Σ sectors must equal generic spintrace.
        auto sector_sum = std::make_shared<Sum>();
        for (auto& [label, sec] : sectors) sector_sum->append(sec->clone());

        ExprPtr sector_total = sector_sum;
        sector_total = biorthogonal_transform_pre_nnsproject(
            sector_total, external_indices(eqvec[i]));

        ExprPtr so = eqvec[i]->clone();
        so->visit(
            [](ExprPtr& n) {
              if (n->is<Tensor>()) n->as<Tensor>().reset_tags();
            },
            /*atoms_only=*/true);
        // Match diagnose_against_generic_spintrace: spintrace keeps labels,
        // then remove_spin_with_relabel (handles α/β collisions after summing
        // spin configurations).
        ExprPtr generic =
            spintrace(so, ext_groups, /*spinfree_index_spaces=*/false);
        canonicalize(generic);
        simplify(generic);
        generic = remove_spin_with_relabel(generic);
        canonicalize(generic);
        simplify(generic);
        generic = biorthogonal_transform_pre_nnsproject(generic, ext_idxs);

        ExprPtr diff = ExprPtr(sector_total) - generic;
        canonicalize(diff);
        simplify(diff);
        std::wcout << L"  R[" << i << L"] (Σ sectors) - (generic) : "
                   << diff->size() << L" terms\n";
        // if (diff->size() != 0) std::wcout << to_latex_align(diff, 0, 4) <<
        // L"\n";
      }
    }

    const auto alpha_qns = IndexSpace::QuantumNumbers{1};
    for (size_t i = 0; i < os_st_eom.size(); ++i) {
      const auto& spin_cases = os_st_eom[i];
      const auto& ext_from_A = ext_index_groups_from_A[i];
      const auto n_groups = ext_from_A.size();
      const size_t n_paired = n_paired_from_A[i];
      const size_t n_half = n_groups - n_paired;
      const bool full_external_spin_encoding =
          (spin_cases.size() == (size_t{1} << n_groups));
      const bool orbit_encoding =
          !full_external_spin_encoding &&
          (spin_cases.size() == (n_paired + 1) * (n_half + 1));

      for (size_t sc = 0; sc < spin_cases.size(); ++sc) {
        container::svector<int> spins(n_groups, 0);
        if (full_external_spin_encoding) {
          for (size_t g = 0; g < n_groups; ++g)
            spins[g] = static_cast<int>((sc >> g) & 1u);
        } else if (orbit_encoding) {
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

      std::wcout << "\n R[" << i << "] spin-free summed (all spin cases): "
                 << os_spinfree_summed_eom[i]->size() << " terms\n";
      // std::wcout << to_latex_align(os_spinfree_summed_eom[i], 0, 4) << "\n";
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
       // .assert_strict_braket_symmetry = false,
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
