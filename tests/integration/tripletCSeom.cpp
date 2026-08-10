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

std::size_t expr_term_count(const ExprPtr& e) {
  if (e->is<Constant>())
    return e->as<Constant>().value() == Constant::scalar_type{0} ? 0 : 1;
  if (e->is<Sum>()) return e->size();
  return 1;
}

// Map each (scalar-stripped) product term to its real coefficient. Key is the
// concatenated canonical latex of the non-scalar factors, so structurally
// identical terms across expressions share a key.
std::map<std::wstring, double> term_coeff_map(ExprPtr e) {
  prepare_expr_for_hashing(e);
  std::map<std::wstring, double> m;
  auto add = [&m](const Product& p) {
    std::wstring key;
    for (const auto& f : p.nonscalar_factors()) key += f->to_latex();
    m[key] += Constant(p.scalar()).value<double>();
  };
  if (e->is<Product>())
    add(e->as<Product>());
  else if (e->is<Sum>())
    for (const auto& t : *e)
      if (t->is<Product>()) add(t->as<Product>());
  return m;
}

// Least-squares fit of `target` by the columns `basis` over the union of term
// keys (solved via 4x4 normal equations, Gauss elimination). Returns the
// coefficients and the max |residual| over all keys. residual ~ 0 means
// `target` lies exactly in the span of `basis` (i.e. a reconstruction exists).
struct FitResult {
  std::vector<double> weights;
  double max_residual = 0;
};
FitResult fit_in_span(const std::vector<ExprPtr>& basis,
                      const ExprPtr& target) {
  const std::size_t n = basis.size();
  std::vector<std::map<std::wstring, double>> B;
  B.reserve(n);
  for (const auto& b : basis) B.push_back(term_coeff_map(b->clone()));
  const auto T = term_coeff_map(target->clone());

  container::set<std::wstring> keys;
  for (const auto& bm : B)
    for (const auto& [k, _] : bm) keys.insert(k);
  for (const auto& [k, _] : T) keys.insert(k);

  // normal equations  (BᵀB) x = Bᵀ t
  std::vector<std::vector<double>> A(n, std::vector<double>(n + 1, 0.0));
  for (const auto& k : keys) {
    std::vector<double> row(n);
    for (std::size_t j = 0; j < n; ++j) {
      auto it = B[j].find(k);
      row[j] = it == B[j].end() ? 0.0 : it->second;
    }
    double tv = 0.0;
    if (auto it = T.find(k); it != T.end()) tv = it->second;
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j) A[i][j] += row[i] * row[j];
      A[i][n] += row[i] * tv;
    }
  }
  // Gauss elimination with partial pivoting
  for (std::size_t c = 0; c < n; ++c) {
    std::size_t piv = c;
    for (std::size_t r = c + 1; r < n; ++r)
      if (std::abs(A[r][c]) > std::abs(A[piv][c])) piv = r;
    std::swap(A[c], A[piv]);
    if (std::abs(A[c][c]) < 1e-14) continue;
    for (std::size_t r = 0; r < n; ++r) {
      if (r == c) continue;
      const double f = A[r][c] / A[c][c];
      for (std::size_t j = c; j <= n; ++j) A[r][j] -= f * A[c][j];
    }
  }
  FitResult res;
  res.weights.resize(n, 0.0);
  for (std::size_t i = 0; i < n; ++i)
    res.weights[i] = std::abs(A[i][i]) < 1e-14 ? 0.0 : A[i][n] / A[i][i];

  for (const auto& k : keys) {
    double v = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
      auto it = B[j].find(k);
      if (it != B[j].end()) v += res.weights[j] * it->second;
    }
    double tv = 0.0;
    if (auto it = T.find(k); it != T.end()) tv = it->second;
    res.max_residual = std::max(res.max_residual, std::abs(v - tv));
  }
  return res;
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

ExprPtr unit_product(const ExprPtr& term) {
  auto out = term->clone();
  auto& p = out->as<Product>();
  p.scale(Product::scalar_type{1} / p.scalar());
  canonicalize(out);
  simplify(out);
  return out;
}

// Match a term against a pool up to scalar factor (exact product structure).
bool term_matches_pool(const ExprPtr& term,
                       const container::vector<ExprPtr>& pool) {
  if (!term->is<Product>()) return false;
  const auto term_unit = unit_product(term);
  for (const auto& cand : pool) {
    if (!cand->is<Product>()) continue;
    ExprPtr diff = term_unit - unit_product(cand);
    canonicalize(diff);
    simplify(diff);
    if (expr_term_count(diff) == 0) return true;
  }
  return false;
}

// Trace how paper combine (3V - V_ps)/16 maps V size-3 hash groups to Omega
// size-4 groups.
void analyze_v_to_omega_transition(
    const ExprPtr& V, const ExprPtr& T_ref,
    const container::map<Index, Index>& pair_swap) {
  const auto v_groups = group_by_hash(V->clone());
  ExprPtr V_ps = transform_expr(V, pair_swap);
  canonicalize(V_ps);
  simplify(V_ps);
  const auto vps_groups = group_by_hash(V_ps->clone());
  const auto t_groups = group_by_hash(T_ref->clone());

  std::size_t n_v3_t4 = 0;
  std::size_t n_from_v_only = 0;
  std::size_t n_from_vps_only = 0;
  std::size_t n_from_both = 0;

  container::vector<ExprPtr> all_v_terms;
  for (const auto& [_, terms] : v_groups)
    for (const auto& t : terms) all_v_terms.push_back(t);
  container::vector<ExprPtr> all_vps_terms;
  for (const auto& [_, terms] : vps_groups)
    for (const auto& t : terms) all_vps_terms.push_back(t);

  for (const auto& [hash, t_terms] : t_groups) {
    const auto v_it = v_groups.find(hash);
    if (v_it != v_groups.end() && v_it->second.size() == 3 &&
        t_terms.size() == 4)
      ++n_v3_t4;

    std::size_t n_v = 0;
    std::size_t n_vps = 0;
    const auto v_pool = v_it != v_groups.end() ? v_it->second : all_v_terms;
    const auto vps_pool = [&]() -> const container::vector<ExprPtr>& {
      const auto it = vps_groups.find(hash);
      return it != vps_groups.end() ? it->second : all_vps_terms;
    }();
    for (const auto& t : t_terms) {
      if (term_matches_pool(t, v_pool)) ++n_v;
      if (term_matches_pool(t, vps_pool)) ++n_vps;
    }
    if (n_v == 4 && n_vps == 0)
      ++n_from_v_only;
    else if (n_v == 0 && n_vps == 4)
      ++n_from_vps_only;
    else if (n_v > 0 && n_vps > 0)
      ++n_from_both;
  }

  std::wcout << L"  [V -> Omega] paper-combine transition:\n";
  std::wcout << L"    hash groups with V size=3 and Omega size=4: " << n_v3_t4
             << L" / " << t_groups.size() << L"\n";
  std::wcout << L"    Omega terms traceable to V pool only: " << n_from_v_only
             << L" groups\n";
  std::wcout << L"    Omega terms traceable to V_ps pool only: "
             << n_from_vps_only << L" groups\n";
  std::wcout << L"    Omega groups with terms from both V and V_ps: "
             << n_from_both << L" / " << t_groups.size() << L"\n";
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

      // Just test ----- os_eom main path: open_shell_CC_spintrace_by_sector +
      // sum ----- auto os_sectors =
      // open_shell_CC_spintrace_by_sector(eqvec[i]); const size_t n_cases =
      // os_sectors.size(); std::wcout << "number of spin cases " << n_cases <<
      // "\n"; SEQUANT_ASSERT(n_cases >= 2);
      //
      // auto summed_spinfree = std::make_shared<Sum>();
      // for (size_t sc = 0; sc < n_cases; ++sc) {
      //   ExprPtr stripped = os_sectors[sc]->clone();
      //   expand(stripped);
      //   std::wcout << "sc" << sc << "\n";
      //   std::wcout << "spin-free sector R[" << i << "] has " <<
      //   stripped->size()
      //              << " terms\n";
      //   canonicalize(stripped);
      //   summed_spinfree->append(stripped);
      // }
      //
      // ExprPtr summed = summed_spinfree;
      // simplify(summed);
      // summed = biorthogonal_transform_pre_nnsproject(summed, ext_idxs);
      // std::wcout << "R[" << i
      //            << "] open-shell sum (Ŝ NOT expanded): " << summed->size()
      //            << " terms\n";
      //
      // auto singlet_ref = closed_shell_CC_spintrace(
      //     eqvec[i], {.method = BiorthogonalizationMethod::V2});
      // simplify(singlet_ref);
      // std::wcout << "R[" << i << "] reference closed-shell (CC spintrace): "
      //            << singlet_ref->size() << " terms\n";
      //
      // ExprPtr os_singlet = summed->clone();
      // simplify(os_singlet);
      // ExprPtr singlet_diff = os_singlet - singlet_ref;
      // canonicalize(singlet_diff);
      // simplify(singlet_diff);
      // std::wcout << "R[" << i << "] (open-shell singlet) - (closed-shell
      // ref): "
      //            << singlet_diff->size() << " terms\n";

      // ----- os_eom independent path: spintrace_by_sector -----------------
      // {
      //   auto sectors = spintrace_by_sector(eqvec[i], ext_groups);
      //   std::wcout << L"\n========== R[" << i << L"] per-sector spin trace ("
      //              << sectors.size() << L" sectors) ==========\n";
      //   for (auto& [label, sec] : sectors) {
      //     std::wcout << L"  sector " << label << L": " << sec->size()
      //                << L" terms, " << count_distinct_hashes(sec->clone())
      //                << L" distinct hashes\n";
      //   }
      //
      //   auto sector_total = std::make_shared<Sum>();
      //   for (auto& [label, sec] : sectors)
      //   sector_total->append(sec->clone());
      //
      //   ExprPtr so = eqvec[i]->clone();
      //   so->visit(
      //       [](ExprPtr& n) {
      //         if (n->is<Tensor>()) n->as<Tensor>().reset_tags();
      //       },
      //       /*atoms_only=*/true);
      //   ExprPtr generic =
      //       spintrace(so, ext_groups, /*spinfree_index_spaces=*/false);
      //   canonicalize(generic);
      //   simplify(generic);
      //   generic = remove_spin(generic, true);
      //   canonicalize(generic);
      //   simplify(generic);
      //
      //   ExprPtr sector_generic_diff = ExprPtr(sector_total) - generic;
      //   canonicalize(sector_generic_diff);
      //   simplify(sector_generic_diff);
      //   std::wcout << L"  R[" << i << L"] (Σ sectors) - (generic) : "
      //              << sector_generic_diff->size() << L" terms\n";
      // }

      // ----- triplet: explicitly spin-coupled basis (Hattig/Kohn/Hald) -----
      // doubles use the T (x) E coupling; triples the T (x) E (x) E coupling
      // (18-op orbit, tools/triplet_triples_check.py)
      if (N > 3 || ext_groups.size() > 3) {
        std::wcout << "R[" << i
                   << "] triplet: skipped (explicitly spin-coupled triplet "
                      "manifold implemented through triples only)\n";
        cs_st_eom[i] = nullptr;
        continue;
      }

      // V_mu = sum over external spin sectors weighted by the sign of the
      // spin of external group 0 (the line carrying T = E(alpha) - E(beta))
      auto triplet_sectors =
          spintrace_by_sector(eqvec[i], ext_groups, /*triplet_R=*/true);
      auto V_sum = std::make_shared<Sum>();
      for (size_t sc = 0; sc < triplet_sectors.size(); ++sc) {
        auto sector = triplet_sectors[sc].second->clone();
        if (sc & 1u) sector = ex<Constant>(-1) * sector;
        V_sum->append(sector);
      }
      ExprPtr V = V_sum;
      canonicalize(V);
      simplify(V);

      auto term_count = [](const ExprPtr& e) -> size_t {
        if (e->is<Constant>()) return e->as<Constant>().value() == 0 ? 0 : 1;
        if (e->is<Sum>()) return e->size();
        return 1;
      };

      ExprPtr T_ref;
      if (ext_groups.size() == 1) {
        // singles: the rank-1 biorthogonal coefficient (1/2) coincides with
        // the singlet one
        T_ref = biorthogonal_transform_pre_nnsproject(V, ext_idxs);
      } else if (ext_groups.size() == 3) {
        // triples: P3 = (1/80)[6 V - V_ps01 - V_ps02 + 2 V_ks12] with an
        // extra 1/2 (metric idempotency: 36 ordered labels cover the 18-op
        // orbit 2:1) -> 1/160 (tools/triplet_triples_check.py, checks 2+6)
        const auto& g0 = ext_idxs.at(0);
        const auto& g1 = ext_idxs.at(1);
        const auto& g2 = ext_idxs.at(2);
        const std::array<Index, 3> b{get_bra_idx(g0), get_bra_idx(g1),
                                     get_bra_idx(g2)};
        const std::array<Index, 3> k{get_ket_idx(g0), get_ket_idx(g1),
                                     get_ket_idx(g2)};
        auto whole_pair_swap = [&](int m, int n) {
          return container::map<Index, Index>{
              {b[m], b[n]}, {b[n], b[m]}, {k[m], k[n]}, {k[n], k[m]}};
        };
        const container::map<Index, Index> ket_swap_12{{k[1], k[2]},
                                                       {k[2], k[1]}};

        // null-space identity: the sum over the whole 18-op orbit (as the 36
        // independent bra x ket external permutations, a 2:1 cover) vanishes
        // -- the rank-3 analog of V + V_ps + V_bs + V_ks = 0
        {
          constexpr std::array<std::array<int, 3>, 6> s3{{{0, 1, 2},
                                                          {0, 2, 1},
                                                          {1, 0, 2},
                                                          {1, 2, 0},
                                                          {2, 0, 1},
                                                          {2, 1, 0}}};
          auto orbit_sum = std::make_shared<Sum>();
          for (const auto& pb : s3) {
            for (const auto& pk : s3) {
              container::map<Index, Index> m;
              for (int n = 0; n != 3; ++n) {
                if (pb[n] != n) m.emplace(b[n], b[pb[n]]);
                if (pk[n] != n) m.emplace(k[n], k[pk[n]]);
              }
              orbit_sum->append(m.empty() ? V->clone() : transform_expr(V, m));
            }
          }
          ExprPtr null_check = orbit_sum;
          canonicalize(null_check);
          simplify(null_check);
          std::wcout << "R[" << i
                     << "] triplet triples null-space identity (Σ 18-op "
                        "orbit): "
                     << term_count(null_check) << " terms (expect 0)\n";
          runtime_assert(term_count(null_check) == 0);
        }

        ExprPtr V_ps01 = transform_expr(V, whole_pair_swap(0, 1));
        ExprPtr V_ps02 = transform_expr(V, whole_pair_swap(0, 2));
        ExprPtr V_ks12 = transform_expr(V, ket_swap_12);
        T_ref =
            ex<Constant>(ratio(1, 160)) *
            (ex<Constant>(6) * V - V_ps01 - V_ps02 + ex<Constant>(2) * V_ks12);
      } else {
        // assemble the paper-native combined R2 residual (Eqs. 7-8, 1/8):
        //   V^{(1)} = (1/8)(V + V_{pair-swap}), V^{(2)} = (1/8)(V -
        //   V_{pair-swap}), Omega = (1/2) V^{(1)} + V^{(2)} = (3V -
        //   V_{pair-swap})/16
        const auto& g0 = ext_idxs.at(0);
        const auto& g1 = ext_idxs.at(1);
        const Index b0 = get_bra_idx(g0);
        const Index b1 = get_bra_idx(g1);
        const Index k0 = get_ket_idx(g0);
        const Index k1 = get_ket_idx(g1);
        const container::map<Index, Index> swap_pair{
            {b0, b1}, {b1, b0}, {k0, k1}, {k1, k0}};
        const container::map<Index, Index> swap_bra{{b0, b1}, {b1, b0}};
        const container::map<Index, Index> swap_ket{{k0, k1}, {k1, k0}};

        ExprPtr V_ps = transform_expr(V, swap_pair);

        ExprPtr null_check = V->clone() + V_ps->clone() +
                             transform_expr(V, swap_bra) +
                             transform_expr(V, swap_ket);
        canonicalize(null_check);
        simplify(null_check);
        std::wcout << "R[" << i
                   << "] triplet null-space identity (V + V_ps + V_bs + "
                      "V_ks): "
                   << term_count(null_check) << " terms (expect 0)\n";
        runtime_assert(term_count(null_check) == 0);

        ExprPtr V_ch1 = ex<Constant>(ratio(1, 8)) * (V + V_ps);
        ExprPtr V_ch2 = ex<Constant>(ratio(1, 8)) * (V - V_ps);
        T_ref = ex<Constant>(ratio(1, 2)) * V_ch1 + V_ch2;
      }
      simplify(T_ref);
      std::wcout << "R[" << i
                 << "] triplet (reference assembly): " << term_count(T_ref)
                 << " terms, and " << count_distinct_hashes(T_ref->clone())
                 << " distinct hashes\n";

      if (hashgroups_ && ext_groups.size() == 2) {
        const auto& g0 = ext_idxs.at(0);
        const auto& g1 = ext_idxs.at(1);
        const container::map<Index, Index> pair_swap{
            {get_bra_idx(g0), get_bra_idx(g1)},
            {get_bra_idx(g1), get_bra_idx(g0)},
            {get_ket_idx(g0), get_ket_idx(g1)},
            {get_ket_idx(g1), get_ket_idx(g0)}};

        std::wcout
            << L"\n========== R[" << i
            << L"] hash-group diagnostics (triplet doubles) ==========\n";
        dump_hash_groups(L"V (sector sum)", V, ext_idxs, 5, 4);
        analyze_group_recovery(L"V (sector sum)", V, ext_idxs, 3);

        dump_hash_groups(L"T_ref (paper Omega)", T_ref, ext_idxs, 5, 4);
        analyze_group_recovery(L"T_ref (paper Omega)", T_ref, ext_idxs, 4);
        analyze_v_to_omega_transition(V, T_ref, pair_swap);
      }

      try {
        const auto tstart = std::chrono::high_resolution_clock::now();
        auto st = closed_shell_EOM_triplet_spintrace(
            eqvec[i], {.method = BiorthogonalizationMethod::V2});
        const auto tstop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> dt = tstop - tstart;

        ExprPtr T_diff = st - T_ref;
        canonicalize(T_diff);
        simplify(T_diff);
        std::wcout << "R[" << i
                   << "] closed_shell_EOM_triplet_spintrace: " << st->size()
                   << " terms, and " << count_distinct_hashes(st->clone())
                   << " distinct hashes, and time: " << dt.count() << " s\n";
        std::wcout << "R[" << i
                   << "] (production triplet) - (reference assembly): "
                   << term_count(T_diff) << " terms (expect 0)\n";
        runtime_assert(term_count(T_diff) == 0);

        // validated term counts of the production triplet residual; pinning
        // `st` also pins T_ref, since the two were just asserted equal
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

        // compact + te/ter experiment knobs are defined for the CCSD doubles
        // study only; in an N >= 3 theory even the doubles residual contains
        // rank-3 R amplitudes which those knobs do not support
        if (ext_groups.size() == 2 && N <= 2) {
          auto compact = closed_shell_EOM_triplet_spintrace(
              eqvec[i], {.method = BiorthogonalizationMethod::V2,
                         .triplet_doubles_compact = true});
          const ExprPtr recon =
              triplet_symbolic_reconstruct(compact, ext_groups);
          runtime_assert(recon == st);
          if (recon == st) {
            std::wcout << "recon == st\n";
          }

          std::wcout << "R[" << i
                     << "] triplet compact (WK factor): " << term_count(compact)
                     << " terms\n";

          // ===== EFV experiment: bare-TE residual =================
          // te_only: drop the external pair-swap TE_ps -> residual = TE/4.
          // ter_only: also drop the column-swapped R amplitude partner.
          auto te_a = closed_shell_EOM_triplet_spintrace(
              eqvec[i], {.method = BiorthogonalizationMethod::V2,
                         .triplet_te_only = true});
          auto te_ab = closed_shell_EOM_triplet_spintrace(
              eqvec[i], {.method = BiorthogonalizationMethod::V2,
                         .triplet_te_only = true,
                         .triplet_amp_no_swap = true});
          simplify(te_a);
          simplify(te_ab);

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
          std::wcout << "  TE-only + no R_swap  : " << term_count(te_ab)
                     << " terms, " << count_distinct_hashes(te_ab->clone())
                     << " distinct hashes\n";
          ExprPtr ab_diff = te_a->clone() - te_ab->clone();
          canonicalize(ab_diff);
          simplify(ab_diff);
          std::wcout << "  te_only vs ter_only (te_a - te_ab)     : "
                     << term_count(ab_diff)
                     << " terms (0 => R_only is a no-op on the residual)\n";

          if (hashgroups_ && ext_groups.size() == 2) {
            std::wcout << "hashgroups for ter_only\n";
            dump_hash_groups(L"te_ab (TE-only, no R_swap)", te_ab, ext_idxs, 5,
                             4);
            analyze_group_recovery(L"te_ab (TE-only, no R_swap)", te_ab,
                                   ext_idxs, 4);
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

          // Sanity identity: Omega == TE/4 + (TE_bs + TE_ks)/16.
          // (Uses the null-space identity TE + TE_ps + TE_bs + TE_ks = 0.)
          ExprPtr V_bs = transform_expr(V, e_swap_bra);
          ExprPtr V_ks = transform_expr(V, e_swap_ket);
          ExprPtr identity_check =
              st->clone() - (ex<Constant>(ratio(1, 4)) * V->clone() +
                             ex<Constant>(ratio(1, 16)) * (V_bs + V_ks));
          canonicalize(identity_check);
          simplify(identity_check);
          std::wcout << "  sanity: Omega - [TE/4 + (TE_bs+TE_ks)/16] = "
                     << term_count(identity_check) << " terms (expect 0)\n";
          runtime_assert(term_count(identity_check) == 0);

          // POSTPROCESSING test (user's idea): TE-only lives in the same 135
          // hash groups as Omega, so a Klein-four postprocessing should rebuild
          // Omega. Since te_a == TE/4, the sanity identity becomes
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

          // ===== TER-only (te_ab, R-only) reconstructability ================
          // Does ANY Klein-four postprocessing rebuild Omega from te_ab?
          // Fit Omega in span{te_ab, bs te_ab, ks te_ab, ps te_ab}; residual ~0
          // => a kernel exists (weights printed), else te_ab is NOT
          // Klein-four-reconstructable (dropped R_swap is unrecoverable).
          std::wcout << "\n  --- TER-only (R-only) reconstructability ---\n";
          ExprPtr ter_naive = st->clone() - te_ab->clone();
          canonicalize(ter_naive);
          simplify(ter_naive);
          std::wcout << "  raw  : Omega - te_ab            = "
                     << term_count(ter_naive) << " terms\n";

          ExprPtr ter_kA =
              st->clone() -
              (te_ab->clone() +
               ex<Constant>(ratio(1, 4)) * (transform_expr(te_ab, e_swap_bra) +
                                            transform_expr(te_ab, e_swap_ket)));
          canonicalize(ter_kA);
          simplify(ter_kA);
          std::wcout << "  kА   : Omega - [te_ab+1/4(bs+ks)te_ab] = "
                     << term_count(ter_kA) << " terms\n";

          const std::vector<ExprPtr> ter_basis{
              te_ab->clone(), transform_expr(te_ab, e_swap_bra),
              transform_expr(te_ab, e_swap_ket),
              transform_expr(te_ab, e_swap_pair)};
          const auto fit = fit_in_span(ter_basis, st);
          std::wcout << "  fit  : Omega = a*te_ab + b*bs + c*ks + d*ps, "
                     << "weights {a,b,c,d} = {" << fit.weights.at(0) << ", "
                     << fit.weights.at(1) << ", " << fit.weights.at(2) << ", "
                     << fit.weights.at(3)
                     << "}, max|resid| = " << fit.max_residual << "\n";
          std::wcout << "  => TER-only "
                     << (fit.max_residual < 1e-9 ? "IS" : "is NOT")
                     << " Klein-four reconstructable\n";

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

        if (hashgroups_ && ext_groups.size() == 2) {
          dump_hash_groups(L"st (production spintrace)", st, ext_idxs, 5, 4);
          analyze_group_recovery(L"st (production spintrace)", st, ext_idxs, 4);
        }

        if (print) {
          // std::wcout << "\n open-shell singlet sum:\n"
          //            << to_latex_align(summed, 20, 1) << "\n";
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

  // compute_all_closedshell_triplet{NMAX, exc_manifold,
  //                                 str2type.at(eqn_type)}(print);
  compute_eomcc_closedshell_triplet{NMAX, exc_manifold, str2type.at(eqn_type),
                                    hashgroups}(print);
}
