#include <SeQuant/version.hpp>

#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/timer.hpp>
#include <SeQuant/domain/mbpt/biorthogonalization.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <range/v3/algorithm/transform.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <map>
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

void prepare_expr_for_hashing(ExprPtr& expr) {
  if (expr->is<Sum>()) {
    for (auto& term : *expr) {
      if (term->is<Product>())
        term = remove_tensor(term.as_shared_ptr<Product>(),
                             reserved::symm_label());
    }
  }
  canonicalize(expr);
  simplify(expr);
}

std::size_t term_network_hash(const ExprPtr& term) {
  SEQUANT_ASSERT(term->is<Product>());
  sequant::TensorNetwork tn(*term.as_shared_ptr<Product>());
  return tn.canonicalize_slots(TensorCanonicalizer::cardinal_tensor_labels())
      .hash_value();
}

std::size_t count_distinct_hashes(ExprPtr expr) {
  prepare_expr_for_hashing(expr);

  if (!expr->is<Sum>()) return expr->is<Constant>() ? 0 : 1;

  container::set<std::size_t> hashes;
  for (const auto& term : *expr) {
    if (!term->is<Product>()) continue;
    hashes.insert(term_network_hash(term));
  }
  return hashes.size();
}

struct HashGroupStats {
  std::size_t n_terms = 0;
  std::size_t n_groups = 0;
  container::map<std::size_t, std::size_t> size_histogram;
};

HashGroupStats hash_group_stats(const ExprPtr& expr_in) {
  HashGroupStats stats;
  ExprPtr expr = expr_in->clone();
  prepare_expr_for_hashing(expr);

  if (!expr->is<Sum>()) {
    stats.n_terms = expr->is<Constant>() ? 0 : 1;
    stats.n_groups = stats.n_terms;
    if (stats.n_terms) stats.size_histogram[1] = 1;
    return stats;
  }

  container::map<std::size_t, container::vector<ExprPtr>> groups;
  for (const auto& term : *expr) {
    if (!term->is<Product>()) continue;
    ++stats.n_terms;
    groups[term_network_hash(term)].push_back(term);
  }
  stats.n_groups = groups.size();
  for (const auto& [_, terms] : groups) ++stats.size_histogram[terms.size()];

  return stats;
}

container::map<std::size_t, container::vector<ExprPtr>> group_by_hash(
    ExprPtr expr_in) {
  ExprPtr expr = expr_in->clone();
  prepare_expr_for_hashing(expr);

  container::map<std::size_t, container::vector<ExprPtr>> groups;
  if (!expr->is<Sum>()) {
    if (expr->is<Product>()) groups[term_network_hash(expr)].push_back(expr);
    return groups;
  }

  for (const auto& term : *expr) {
    if (!term->is<Product>()) continue;
    groups[term_network_hash(term)].push_back(term);
  }
  return groups;
}

std::wstring coeff_to_wstring(const Product::scalar_type& c) {
  return to_latex_align(ex<Constant>(c), 1, 1);
}

// check the compact/reconstruct. If the compactness is not correct
// report which tensor-network groups disagree
void report_reconstruction_mismatch(const ExprPtr& diff, const ExprPtr& full,
                                    const ExprPtr& compact) {
  if (!diff->is<Sum>()) return;
  container::set<std::size_t> bad_hashes;
  for (const auto& t : *diff) {
    if (!t->is<Product>()) continue;
    const auto h = term_network_hash(t);
    bad_hashes.insert(h);
    std::wcout << L"    DIFF hash=0x" << std::hex << h << std::dec << L" coeff="
               << coeff_to_wstring(t->as<Product>().scalar()) << L"\n";
  }
  const auto full_groups = group_by_hash(full->clone());
  const auto compact_groups = group_by_hash(compact->clone());
  for (const auto h : bad_hashes) {
    std::wcout << L"    BAD group 0x" << std::hex << h << std::dec << L"\n";
    if (const auto it = full_groups.find(h); it != full_groups.end())
      for (const auto& t : it->second)
        std::wcout << L"      full: " << t->to_latex() << L"\n";
    if (const auto it = compact_groups.find(h); it != compact_groups.end())
      for (const auto& t : it->second)
        std::wcout << L"      kept: " << t->to_latex() << L"\n";
  }
}

std::size_t expr_term_count(const ExprPtr& e) {
  if (e->is<Constant>())
    return e->as<Constant>().value() == Constant::scalar_type{0} ? 0 : 1;
  if (e->is<Sum>()) return e->size();
  return 1;
}

void print_hash_histogram(const std::wstring& stage_label,
                          const ExprPtr& expr) {
  const auto stats = hash_group_stats(expr);
  std::wcout << L"  [" << stage_label << L"] " << stats.n_terms << L" terms, "
             << stats.n_groups << L" distinct hashes\n";
  std::wcout << L"    group-size histogram: ";
  for (const auto& [sz, cnt] : stats.size_histogram) {
    std::wcout << L"size=" << sz << L" x" << cnt << L"  ";
  }
  std::wcout << L"\n";
}

enum class ExternalSwapKind { Identity, BraSwap, KetSwap, PairSwap, Other };

ExternalSwapKind classify_external_swap(
    const container::svector<container::svector<SlottedIndex>>& ext_idxs,
    const container::svector<Index>& ref_bra,
    const container::svector<Index>& ref_ket,
    const container::svector<Index>& bra,
    const container::svector<Index>& ket) {
  if (ext_idxs.size() != 2 || ref_bra.size() != 2 || ref_ket.size() != 2 ||
      bra.size() != 2 || ket.size() != 2) {
    return ExternalSwapKind::Other;
  }

  const Index b0 = get_bra_idx(ext_idxs.at(0));
  const Index b1 = get_bra_idx(ext_idxs.at(1));
  const Index k0 = get_ket_idx(ext_idxs.at(0));
  const Index k1 = get_ket_idx(ext_idxs.at(1));

  auto matches = [&](const container::svector<Index>& b,
                     const container::svector<Index>& k) {
    return b.at(0) == bra.at(0) && b.at(1) == bra.at(1) &&
           k.at(0) == ket.at(0) && k.at(1) == ket.at(1);
  };

  if (matches(ref_bra, ref_ket)) return ExternalSwapKind::Identity;
  if (matches(container::svector<Index>{b1, b0}, ref_ket))
    return ExternalSwapKind::BraSwap;
  if (matches(ref_bra, container::svector<Index>{k1, k0}))
    return ExternalSwapKind::KetSwap;
  if (matches(container::svector<Index>{b1, b0},
              container::svector<Index>{k1, k0}))
    return ExternalSwapKind::PairSwap;
  return ExternalSwapKind::Other;
}

std::wstring swap_kind_label(ExternalSwapKind k) {
  switch (k) {
    case ExternalSwapKind::Identity:
      return L"id";
    case ExternalSwapKind::BraSwap:
      return L"bra";
    case ExternalSwapKind::KetSwap:
      return L"ket";
    case ExternalSwapKind::PairSwap:
      return L"pair";
    case ExternalSwapKind::Other:
      return L"?";
  }
  return L"?";
}

// Collect bra/ket index labels from the first EOM amplitude tensor (R/L) in a
// product, if present; otherwise from the first tensor with rank >= 2.
std::pair<container::svector<Index>, container::svector<Index>>
external_slot_labels(const Product& p) {
  for (const auto& f : p) {
    if (!f->is<Tensor>()) continue;
    const auto& t = f->as<Tensor>();
    if (t.label() == L"R" || t.label() == L"L") {
      return {container::svector<Index>(t.bra().begin(), t.bra().end()),
              container::svector<Index>(t.ket().begin(), t.ket().end())};
    }
  }
  for (const auto& f : p) {
    if (!f->is<Tensor>()) continue;
    const auto& t = f->as<Tensor>();
    if (t.rank() >= 2) {
      return {container::svector<Index>(t.bra().begin(), t.bra().end()),
              container::svector<Index>(t.ket().begin(), t.ket().end())};
    }
  }
  return {};
}

void dump_hash_groups(
    const std::wstring& stage_label, const ExprPtr& expr,
    const container::svector<container::svector<SlottedIndex>>& ext_idxs,
    std::size_t max_groups, std::size_t max_terms_per_group) {
  print_hash_histogram(stage_label, expr);

  const auto groups = group_by_hash(expr->clone());
  std::size_t printed = 0;
  for (const auto& [hash, terms] : groups) {
    if (printed >= max_groups) break;
    ++printed;

    std::wcout << L"\n    --- hash group 0x" << std::hex << hash << std::dec
               << L" (" << terms.size() << L" terms) ---\n";

    container::svector<Index> ref_bra, ref_ket;
    if (!terms.empty() && terms.front()->is<Product>()) {
      std::tie(ref_bra, ref_ket) =
          external_slot_labels(terms.front()->as<Product>());
    }

    std::size_t tidx = 0;
    for (const auto& term : terms) {
      if (tidx >= max_terms_per_group) {
        std::wcout << L"      ... (" << (terms.size() - max_terms_per_group)
                   << L" more terms)\n";
        break;
      }
      const auto& p = term->as<Product>();
      const auto coeff = p.scalar();
      const auto [bra, ket] = external_slot_labels(p);
      const auto swap =
          classify_external_swap(ext_idxs, ref_bra, ref_ket, bra, ket);

      std::wcout << L"      [" << tidx << L"] coeff=" << coeff_to_wstring(coeff)
                 << L" swap=" << swap_kind_label(swap) << L" bra=(";
      for (std::size_t i = 0; i < bra.size(); ++i) {
        if (i) std::wcout << L",";
        std::wcout << bra.at(i).full_label();
      }
      std::wcout << L") ket=(";
      for (std::size_t i = 0; i < ket.size(); ++i) {
        if (i) std::wcout << L",";
        std::wcout << ket.at(i).full_label();
      }
      std::wcout << L")\n        " << to_latex_align(term, 1, 1) << L"\n";
      ++tidx;
    }
  }
  if (groups.size() > max_groups) {
    std::wcout << L"    ... (" << (groups.size() - max_groups)
               << L" more hash groups not shown)\n";
  }
}

// For each group, try to express every term as an external-index transform of
// the first term. Report groups where all members are {id,bra,ket,pair} swaps.
void analyze_group_recovery(
    const std::wstring& stage_label, const ExprPtr& expr,
    const container::svector<container::svector<SlottedIndex>>& ext_idxs,
    std::size_t target_group_size) {
  const auto groups = group_by_hash(expr->clone());

  std::size_t n_target = 0;
  std::size_t n_swap_classified = 0;
  std::size_t n_permutation_recovery = 0;
  std::size_t n_coeff_sum_zero = 0;
  std::size_t n_symbolic_group_sum_zero = 0;

  for (const auto& [_, terms] : groups) {
    if (terms.size() != target_group_size) continue;
    ++n_target;

    sequant::rational coeff_sum{0};
    bool coeff_sum_exact = true;
    for (const auto& t : terms) {
      if (!t->is<Product>()) {
        coeff_sum_exact = false;
        break;
      }
      const auto s = t->as<Product>().scalar();
      if (s.imag() != 0) {
        coeff_sum_exact = false;
        break;
      }
      coeff_sum += s.real();
    }
    if (coeff_sum_exact && coeff_sum == sequant::rational{0})
      ++n_coeff_sum_zero;

    ExprPtr group_sum = ex<Constant>(sequant::rational{0});
    for (const auto& t : terms) group_sum = group_sum + t;
    canonicalize(group_sum);
    simplify(group_sum);
    if (expr_term_count(group_sum) == 0) ++n_symbolic_group_sum_zero;

    if (target_group_size != 4 || !terms.front()->is<Product>()) continue;

    const auto [ref_bra, ref_ket] =
        external_slot_labels(terms.front()->as<Product>());

    bool all_classified = true;
    container::set<ExternalSwapKind> kinds;
    for (const auto& t : terms) {
      if (!t->is<Product>()) {
        all_classified = false;
        break;
      }
      const auto [bra, ket] = external_slot_labels(t->as<Product>());
      const auto kind =
          classify_external_swap(ext_idxs, ref_bra, ref_ket, bra, ket);
      kinds.insert(kind);
      if (kind == ExternalSwapKind::Other) all_classified = false;
    }
    if (all_classified && kinds.size() == 4) ++n_swap_classified;

    // Permutation recovery: term_j == P(term_0) for one of four swaps
    if (ext_idxs.size() == 2) {
      const Index b0 = get_bra_idx(ext_idxs.at(0));
      const Index b1 = get_bra_idx(ext_idxs.at(1));
      const Index k0 = get_ket_idx(ext_idxs.at(0));
      const Index k1 = get_ket_idx(ext_idxs.at(1));
      const container::map<Index, Index> swap_bra{{b0, b1}, {b1, b0}};
      const container::map<Index, Index> swap_ket{{k0, k1}, {k1, k0}};
      container::map<Index, Index> swap_pair = swap_bra;
      for (const auto& [a, b] : swap_ket) swap_pair.emplace(a, b);

      const auto ref = terms.front()->clone();
      const container::map<Index, Index> id_map;
      const container::map<Index, Index>* maps[] = {&id_map, &swap_bra,
                                                    &swap_ket, &swap_pair};

      bool group_ok = true;
      for (std::size_t j = 1; j < terms.size(); ++j) {
        bool matched = false;
        for (const auto* m : maps) {
          ExprPtr transformed = transform_expr(ref, *m);
          canonicalize(transformed);
          simplify(transformed);
          ExprPtr diff = transformed - terms.at(j);
          canonicalize(diff);
          simplify(diff);
          if (diff->is<Constant>() &&
              diff->as<Constant>().value() == sequant::rational{0}) {
            matched = true;
            break;
          }
        }
        if (!matched) {
          group_ok = false;
          break;
        }
      }
      if (group_ok) ++n_permutation_recovery;
    }
  }

  std::wcout << L"  [" << stage_label << L"] recovery summary (size-"
             << target_group_size << L" groups): " << n_target << L" groups\n";
  std::wcout << L"    coeff sums zero: " << n_coeff_sum_zero << L" / "
             << n_target << L"\n";
  std::wcout << L"    symbolic group sums zero: " << n_symbolic_group_sum_zero
             << L" / " << n_target << L"\n";
  if (target_group_size == 4) {
    std::wcout << L"    all 4 swap kinds {id,bra,ket,pair}: "
               << n_swap_classified << L" / " << n_target << L"\n";
    std::wcout << L"    permutation recovery from 1st term: "
               << n_permutation_recovery << L" / " << n_target << L"\n";
  }
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
  bool hashgroups_;

 public:
  compute_eomcc_closedshell_triplet(size_t n, const std::string& exc_manifold,
                                    EqnType t = EqnType::right,
                                    bool hashgroups = false)
      : N(n), manifold(exc_manifold), type(t), hashgroups_(hashgroups) {
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

      // triplet: explicitly spin-coupled basis (Hattig/Kohn/Hald + Faber
      // paper), with a different approach
      if (N > 3 || ext_groups.size() > 3) {
        std::wcout << "R[" << i
                   << "] triplet: skipped (explicitly spin-coupled triplet "
                      "manifold implemented through triples only)\n";
        cs_st_eom[i] = nullptr;
        continue;
      }

      auto term_count = [](const ExprPtr& e) -> size_t {
        if (e->is<Constant>()) return e->as<Constant>().value() == 0 ? 0 : 1;
        if (e->is<Sum>()) return e->size();
        return 1;
      };

      try {
        const auto tstart = std::chrono::high_resolution_clock::now();
        auto st = closed_shell_EOM_triplet_spintrace(
            eqvec[i], {.method = BiorthogonalizationMethod::V2});
        const auto tstop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> dt = tstop - tstart;
        std::wcout << "R[" << i
                   << "] closed_shell_EOM_triplet_spintrace: " << st->size()
                   << " terms, and " << count_distinct_hashes(st->clone())
                   << " distinct hashes, and time: " << dt.count() << " s\n";

        // validated term counts of the production triplet residual
        const auto n_st = term_count(st);
        if (N == 2 && type == EqnType::right) {
          if (np == 2 && nh == 2) {  // triplet EOM-CCSD(2h2p)
            if (i == 1) runtime_assert(n_st == 42);
            if (i == 2) runtime_assert(n_st == 540);
          }
        }
        if (N == 3 && type == EqnType::right) {
          if (np == 3 && nh == 3) {  // triplet EOM-CCSDT(3h3p)
            if (i == 1) runtime_assert(n_st == 54);
            if (i == 2) runtime_assert(n_st == 920);
            if (i == 3) runtime_assert(n_st == 11592);
          }
        }

        // compact residual + symbolic reconstruction
        // defined for the rank-2 and rank-3 manifolds (doubles:
        // the -3c member of each {c,c,c,-3c} group; triples: one
        // scaled member per 36 slot perms). The singles residual
        // has nothing to compact (one term per hash group, with factor 1/2).
        ExprPtr compact;
        if (ext_groups.size() == 2 || ext_groups.size() == 3) {
          const auto cstart = std::chrono::high_resolution_clock::now();
          compact = closed_shell_EOM_triplet_spintrace(
              eqvec[i], {.method = BiorthogonalizationMethod::V2,
                         .triplet_doubles_compact = true});
          const auto cstop = std::chrono::high_resolution_clock::now();
          std::chrono::duration<double> cdt = cstop - cstart;
          std::wcout << "R[" << i
                     << "] triplet compact (WK factor): " << term_count(compact)
                     << " terms, and "
                     << count_distinct_hashes(compact->clone())
                     << " distinct hashes, and time: " << cdt.count() << " s\n";

          const ExprPtr recon =
              triplet_symbolic_reconstruct(compact, ext_groups);
          ExprPtr recon_diff = recon->clone() - st->clone();
          canonicalize(recon_diff);
          simplify(recon_diff);
          std::wcout << "R[" << i
                     << "] reconstruct(compact): " << term_count(recon)
                     << " terms; reconstruct - full: " << term_count(recon_diff)
                     << " terms (expect 0)\n";
          if (term_count(recon_diff) > 0)
            report_reconstruction_mismatch(recon_diff, st, compact);
          runtime_assert(term_count(recon_diff) == 0);
        }

        // ===== EFV experiment: bare-TE residual (CCSD doubles only for now)
        // ===== te_only is a doubles-only method for now; in an N >= 3 theory
        // even the doubles residual contains rank-3 R amplitudes it does not
        // support.
        if (ext_groups.size() == 2 && N <= 2) {
          // te_only: drop the external pair-swap TE_ps (or we can say ET)->
          // residual = TE/4.
          auto te_a = closed_shell_EOM_triplet_spintrace(
              eqvec[i], {.method = BiorthogonalizationMethod::V2,
                         .triplet_te_only = true});
          simplify(te_a);

          std::wcout << "\n----- EFV experiment (TE-only) comparison R[" << i
                     << "] -----\n";
          std::wcout << "  production Omega = (3TE-TE_ps)/16 : "
                     << term_count(st) << " terms, "
                     << count_distinct_hashes(st->clone())
                     << " distinct hashes\n";
          std::wcout << "  compact (-3c representative)      : "
                     << term_count(compact) << " terms, "
                     << count_distinct_hashes(compact->clone())
                     << " distinct hashes\n";
          std::wcout << "  TE-only                 : " << term_count(te_a)
                     << " terms, " << count_distinct_hashes(te_a->clone())
                     << " distinct hashes\n";

          if (hashgroups_) {
            dump_hash_groups(L"te_a (TE-only)", te_a, ext_idxs, 5, 4);
            analyze_group_recovery(L"te_a (TE-only)", te_a, ext_idxs, 4);
          }

          // external single/pair swap maps on the projector (mu) indices
          const auto& eg0 = ext_idxs.at(0);
          const auto& eg1 = ext_idxs.at(1);
          const Index b0 = get_bra_idx(eg0);
          const Index b1 = get_bra_idx(eg1);
          const Index k0 = get_ket_idx(eg0);
          const Index k1 = get_ket_idx(eg1);
          const container::map<Index, Index> e_swap_bra{{b0, b1}, {b1, b0}};
          const container::map<Index, Index> e_swap_ket{{k0, k1}, {k1, k0}};
          const container::map<Index, Index> e_swap_pair{
              {b0, b1}, {b1, b0}, {k0, k1}, {k1, k0}};

          // POSTPROCESSING test: TE-only lives in the same 135
          // hash groups as Omega, so a postprocessing should rebuild
          // Omega via the null-space identity TE + TE_ps + TE_bs + TE_ks = 0:
          //   Omega == te_a + (1/4)( bra_swap(te_a) + ket_swap(te_a) ).
          ExprPtr te_recon =
              te_a->clone() +
              ex<Constant>(ratio(1, 4)) * (transform_expr(te_a, e_swap_bra) +
                                           transform_expr(te_a, e_swap_ket));
          ExprPtr te_recon_diff = st->clone() - te_recon;
          canonicalize(te_recon_diff);
          simplify(te_recon_diff);
          std::wcout << "  postproc: Omega - [te_a + 1/4(bs+ks)te_a] = "
                     << term_count(te_recon_diff) << " terms (expect 0)\n";
          runtime_assert(term_count(te_recon_diff) == 0);

          // What TE-only drops (raw, expected nonzero): D = Omega - TE/4.
          ExprPtr D = st->clone() - te_a->clone();
          canonicalize(D);
          simplify(D);
          std::wcout << "  raw difference Omega - TE-only   : " << term_count(D)
                     << " terms (nonzero => not a free identity)\n";

          // Physical-subspace test: project the difference onto the
          // Non-Null-Space, P_nns = 1 - 1/4 (1+Pij)(1+Pab), acting on the
          // external mu indices (Pij = bra swap, Pab = ket swap). If P_nns D =
          // 0 the two residuals agree on the physical triplet bra space => the
          // EFV experiment identity holds; nonzero => the single-cross
          // symmetrization is load-bearing.
          auto project_nns = [&](const ExprPtr& e) -> ExprPtr {
            ExprPtr avg = e->clone() + transform_expr(e, e_swap_bra) +
                          transform_expr(e, e_swap_ket) +
                          transform_expr(e, e_swap_pair);
            ExprPtr res = e->clone() - ex<Constant>(ratio(1, 4)) * avg;
            canonicalize(res);
            simplify(res);
            return res;
          };
          ExprPtr nns_D = project_nns(D);
          std::wcout << "  P_nns (Omega - TE-only)          : "
                     << term_count(nns_D)
                     << " terms (0 => TE-only physically equivalent)\n";
          std::wcout << "-----------------------------------------------\n";
        }

        if (hashgroups_ && ext_groups.size() <= 3) {
          dump_hash_groups(L"st (production spintrace)", st, ext_idxs, 5, 4);
          analyze_group_recovery(L"st (production spintrace)", st, ext_idxs, 4);
        }

        if (print) {
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
  const bool hashgroups = (argc > 5 && std::string(argv[5]) == "hashgroups") ||
                          print_str == "hashgroups";

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

  compute_eomcc_closedshell_triplet{NMAX, exc_manifold, str2type.at(eqn_type),
                                    hashgroups}(print);
}
