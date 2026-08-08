//
// Tests for the genetic factorizer (SeQuant/core/optimize/ga). The parity
// sections validate the port against the Python prototype's exact reference
// numbers (see ga_reference.hpp) on the PAO-CCSD-DF f-only and g-only E+R1
// universes parsed from ga_pao_ccsd_df.md.
//

#include "catch2_sequant.hpp"
#include "ga_reference.hpp"

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node_compare.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/ga/cost.hpp>
#include <SeQuant/core/optimize/ga/emit.hpp>
#include <SeQuant/core/optimize/ga/fitness.hpp>
#include <SeQuant/core/optimize/ga/ga.hpp>
#include <SeQuant/core/optimize/ga/genome.hpp>
#include <SeQuant/core/optimize/ga/key_table.hpp>
#include <SeQuant/core/optimize/ga/optimize_ga.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <catch2/catch_test_macros.hpp>

#include <bit>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <functional>
#include <random>
#include <set>
#include <sstream>
#include <unordered_map>

namespace {

using namespace sequant;
namespace ga = sequant::opt::ga;

// ---------------------------------------------------------------- codec ----

std::size_t double_factorial(int n) {  // (2n-3)!!
  std::size_t r = 1;
  for (int k = 2 * n - 3; k > 1; k -= 2) r *= k;
  return r;
}

// ------------------------------------------------- universe from markdown --

struct MdTerm {
  std::wstring target;
  ExprPtr summand;  // Product of the data tensors, with the scalar prefactor
  ga::FaceSet ext;
  container::svector<std::wstring> labels;
};

void collect_factors(ExprPtr const& e, Constant::scalar_type& scalar,
                     container::svector<ExprPtr>& tensors) {
  if (e->is<Tensor>())
    tensors.push_back(e);
  else if (e->is<Constant>())
    scalar *= e->as<Constant>().value();
  else if (e->is<Product>()) {
    auto const& p = e->as<Product>();
    scalar *= p.scalar();
    for (auto const& f : p.factors()) collect_factors(f, scalar, tensors);
  } else
    FAIL("unexpected expression kind in term");
}

container::svector<MdTerm> load_md_terms(std::string const& path) {
  std::ifstream in(path);
  REQUIRE(in);
  container::svector<MdTerm> out;
  std::wstring target;
  std::string line;
  while (std::getline(in, line)) {
    if (line.rfind("## E", 0) == 0) {
      target = L"E";
      continue;
    }
    if (line.rfind("## R1", 0) == 0) {
      target = L"R1";
      continue;
    }
    if (line.rfind("## R2", 0) == 0) {
      target = L"R2";
      continue;
    }
    if (target.empty() || line.empty() || !std::isdigit(line.front()))
      continue;
    const auto lo = line.find('`'), hi = line.rfind('`');
    if (lo == std::string::npos || hi <= lo) continue;
    const std::string body = line.substr(lo + 1, hi - lo - 1);
    ExprPtr expr = deserialize(body);
    Constant::scalar_type scalar = 1;
    container::svector<ExprPtr> tensors;
    collect_factors(expr, scalar, tensors);
    MdTerm term;
    term.target = target;
    Product prod(scalar, ExprPtrList{});
    for (auto const& t : tensors) {
      if (t->as<Tensor>().label() == L"Ŝ") {
        for (auto const& ix : t->as<Tensor>().const_indices())
          term.ext.emplace(ix);
      } else {
        prod.append(1, t, Product::Flatten::No);
        term.labels.push_back(std::wstring(t->as<Tensor>().label()));
      }
    }
    term.summand = ex<Product>(std::move(prod));
    out.push_back(std::move(term));
  }
  return out;
}

// mode: 0 = f-only (f and no g, E+R1), 1 = g-only (g and no f, E+R1),
// 2 = everything (full CCSD: E+R1+R2), 3 = R1+R2 -- the 81-term universe the
// GATE-B benchmark and the mpqc4 CSV-CCSD job actually optimize
container::svector<ga::TargetInput> universe(
    container::svector<MdTerm> const& terms, int mode) {
  auto has = [](MdTerm const& t, std::wstring const& lbl) {
    return ranges::find(t.labels, lbl) != ranges::end(t.labels);
  };
  container::svector<ga::TargetInput> targets;
  for (auto const& t : terms) {
    if (mode == 3 && t.target == L"E") continue;
    if (mode < 2) {
      if (t.target != L"E" && t.target != L"R1") continue;
      if (mode == 0 ? (!has(t, L"f") || has(t, L"g"))
                    : (!has(t, L"g") || has(t, L"f")))
        continue;
    }
    auto it = std::find_if(targets.begin(), targets.end(),
                           [&](auto const& tg) { return tg.label == t.target; });
    if (it == targets.end()) {
      targets.push_back({t.target, {}, {}});
      it = targets.end() - 1;
    }
    it->summands.push_back(t.summand);
    it->ext.push_back(t.ext);
  }
  return targets;
}

// ------------------------------------------- cache-aware runtime coster ----
// Costs a binarized forest the way evaluate() + CacheManager do: every
// structurally repeated internal node (TreeNodeHasher/-EqualityComparator
// identity -- exactly what the runtime cache is keyed on) is computed once
// and is a cache hit thereafter. Same flop convention as
// CostModel::native(). Producer-substituted emission makes subtree identity
// coincide with the schedule's array identity, so the unique contraction
// nodes must reproduce Schedule::l1 exactly; sum-carrying structures are not
// keyed by the schedule (the GA charges them once per target occurrence), so
// structurally duplicate sums the cache shares anyway are the only -- and
// exactly enumerable -- difference.
struct RuntimeCoster {
  using Node = EvalNode<EvalExpr>;
  template <typename V>
  using NodeMap =
      std::unordered_map<Node const*, V, TreeNodeHasher<Node>,
                         TreeNodeEqualityComparator<Node>>;

  std::function<std::size_t(Index const&)> extent;
  NodeMap<std::size_t> uses;  ///< cache-semantics use counts (== max_life)
  NodeMap<bool> sum_below;
  double cache_total = 0;  ///< every unique node charged once
  double l1_unique = 0;    ///< unique contraction nodes with no Sum below
  double l2_unique = 0;    ///< unique sum-carrying nodes
  double l2_every = 0;     ///< sum-carrying nodes, every occurrence

  double volume(Index::index_vector const& a,
                Index::index_vector const& b = {}) {
    container::set<Index> closure;
    for (auto const* ixs : {&a, &b})
      for (auto const& ix : *ixs) {
        closure.emplace(ix);
        for (auto const& p : ix.proto_indices()) closure.emplace(p);
      }
    double v = 1;
    for (auto const& ix : closure) v *= static_cast<double>(extent(ix));
    return v;
  }

  double cost(Node const& n) {
    if (n.leaf()) return 0;
    if (*n->op_type() == EvalOp::Sum) return volume(n->canon_indices());
    if (n.left()->is_constant() || n.right()->is_constant()) return 0;
    const double v =
        volume(n.left()->canon_indices(), n.right()->canon_indices());
    return v == 1 ? 0 : v;
  }

  bool carries_sum(Node const& n) {
    if (n.leaf()) return false;
    if (auto it = sum_below.find(&n); it != sum_below.end()) return it->second;
    const bool r = *n->op_type() == EvalOp::Sum || carries_sum(n.left()) ||
                   carries_sum(n.right());
    sum_below.emplace(&n, r);
    return r;
  }

  // cache semantics: charge on first encounter, never descend into repeats
  void walk_cached(Node const& n) {
    if (n.leaf()) return;
    if (++uses[&n] > 1) return;
    const double c = cost(n);
    cache_total += c;
    (carries_sum(n) ? l2_unique : l1_unique) += c;
    walk_cached(n.left());
    walk_cached(n.right());
  }

  // the schedule's per-occurrence L2 accounting: sum-carrying nodes are
  // upward-closed, so pure contraction subtrees need not be entered
  void walk_l2(Node const& n) {
    if (n.leaf() || !carries_sum(n)) return;
    l2_every += cost(n);
    walk_l2(n.left());
    walk_l2(n.right());
  }
};

std::size_t python_extent(Index const& ix) {
  switch (ix.proto_indices().size()) {
    case 1:
      return 100;  // OSV
    case 2:
      return 50;  // PNO
    default:
      break;
  }
  auto const& base = ix.space().base_key();
  if (base == L"i") return 180;
  if (base == L"μ̃") return 1000;
  if (base == L"Κ") return 3000;
  FAIL("unexpected index space in the PAO-CCSD-DF universe");
  return 1;
}

/// The C4H10 / cc-pVDZ extents the GATE-B benchmark and the mpqc4 job use
/// (occ 17 / PAO 106 / aux 364 / CSV 21). Verbatim from
/// `GA_scratch/dag_batching/bench_ga.cpp`, so the 81-term universe below is
/// the one the fingerprint is taken on.
std::size_t dz_extent(Index const& ix) {
  if (ix.proto_indices().size() >= 1) return 21;
  auto const& base = ix.space().base_key();
  if (base == L"i") return 17;
  if (base == L"μ̃") return 106;
  if (base == L"Κ") return 364;
  return 1;
}

// cost of one term's tree under the prototype convention
double tree_cost(ga::KeyTable const& kt, ga::CostModel const& cm, int d,
                 ga::TreeCode const& code) {
  auto const& T = kt.terms[d];
  auto fam = ga::decode_tree(code, static_cast<int>(T.n()));
  double c = 0;
  for (ga::NodeMask S : fam) {
    if (std::popcount(S) < 2) continue;
    auto [c1, c2] = ga::tree_children(fam, S);
    c += cm.merge(T, c1, c2);
  }
  return c;
}

ga::Genome slice_genome(std::vector<int> const& g, std::vector<int> const& h) {
  ga::Genome out;
  out.g.assign(g.begin(), g.end());
  out.h.assign(h.begin(), h.end());
  return out;
}

// ------------------------------------------ cost path == explain (T-A5) ----

/// All genes zero: a valid genome for any layout (gene k of a block has range
/// 2k-3, so 0 is always legal) and a deterministic starting point that costs
/// nothing to build -- unlike `seed_genome`, whose per-term DP is 21 s on the
/// 81-term DZ universe.
ga::Genome zero_genome(ga::GenomeLayout const& lay) {
  ga::Genome out;
  out.g.assign(lay.g_ranges.size(), 0);
  out.h.assign(lay.h_ranges.size(), 0);
  return out;
}

/// `k` NNI moves on randomly chosen L1/L2 blocks -- ga.cpp's `nni_kick`,
/// restated here (it is file-local there) against a caller-owned RNG, so the
/// perturbed genomes are reproducible from a fixed seed.
void nni_kick_local(ga::KeyTable const& kt, ga::GenomeLayout const& lay,
                    ga::Genome& x, std::mt19937_64& rng, int k) {
  auto below = [&](std::size_t n) {
    return std::uniform_int_distribution<std::size_t>{0, n - 1}(rng);
  };
  auto kick_block = [&](container::svector<int>& code,
                        std::pair<int, int> slice, int n) {
    if (n < 3) return;
    ga::TreeCode block(code.begin() + slice.first, code.begin() + slice.second);
    auto const moves = ga::nni_moves(ga::decode_tree(block, n));
    if (moves.empty()) return;
    auto const next = ga::encode_tree(moves[below(moves.size())], n);
    std::copy(next.begin(), next.end(), code.begin() + slice.first);
  };
  for (int i = 0; i < k; ++i) {
    if (std::uniform_real_distribution<>{}(rng) < 0.85) {  // ga.cpp l1_move_prob
      const auto d = below(kt.terms.size());
      kick_block(x.g, lay.g_slice[d], static_cast<int>(kt.terms[d].n()));
    } else {
      const auto t = below(kt.targets.size());
      kick_block(x.h, lay.h_slice[t],
                 static_cast<int>(kt.targets[t].terms.size()));
    }
  }
}

/// Node census of an L2 value. The tripwire below is only meaningful if the
/// genomes it walks actually reach `build_node`'s extraction branch (i.e.
/// `find_beta` succeeds somewhere) and reach it with a multi-tensor shared
/// factor -- so it asserts on these counts rather than trusting the kicks.
struct L2Shape {
  std::size_t fx = 0;       ///< extraction nodes
  std::size_t fx_wide = 0;  ///< ... whose extracted cluster is itself a merge
  std::size_t sum = 0;      ///< plain additions
  std::size_t cl = 0;       ///< L1 arrays reaching L2
  std::size_t cl_wide = 0;  ///< ... that are merges (i.e. demanded of L1)
  void operator+=(L2Shape const& o) {
    fx += o.fx;
    fx_wide += o.fx_wide;
    sum += o.sum;
    cl += o.cl;
    cl_wide += o.cl_wide;
  }
};

L2Shape l2_shape(ga::ValPtr const& v) {
  L2Shape s;
  if (!v) return s;
  switch (v->kind) {
    case ga::Val::Cl:
      s.cl = 1;
      s.cl_wide = std::popcount(v->S) >= 2 ? 1 : 0;
      return s;
    case ga::Val::Fx:
      s += l2_shape(v->inner);
      ++s.fx;
      if (std::popcount(v->V) >= 2) ++s.fx_wide;
      return s;
    case ga::Val::Sum:
      s += l2_shape(v->s1.val);
      s += l2_shape(v->s2.val);
      ++s.sum;
      return s;
  }
  return s;
}

/// THE TRIPWIRE for T-A5's two evaluation paths. `Fitness::operator()` (the
/// cost-only recursion the ~164k-evaluation search runs) and
/// `Fitness::explain` (the rich `Schedule` mpqc4 emits from) must return the
/// SAME double, forever: the search optimizes the former and the emitter
/// realizes the latter, so any divergence silently emits a schedule that is
/// not the one the search chose -- and no other gate would attribute it here.
///
/// Exact double equality, never a tolerance: the two paths run the same
/// additions in the same order by construction, so anything but bit equality
/// is a real defect, not rounding.
void check_cost_path_agrees(ga::KeyTable const& kt, ga::Fitness const& F,
                            container::svector<ga::Genome> const& bases,
                            std::uint64_t rng_seed) {
  std::set<double> totals;
  L2Shape shape;
  auto one = [&](ga::Genome const& g) {
    auto const sch = F.explain(g);
    CHECK(F(g) == sch.total);  // <-- the tripwire
    totals.insert(sch.total);
    for (auto const& root : sch.roots) shape += l2_shape(root.val);
  };
  std::mt19937_64 rng(rng_seed);
  for (auto const& base : bases) {
    one(base);
    ga::Genome g = base;
    for (int i = 0; i < 200 / static_cast<int>(bases.size()); ++i) {
      nni_kick_local(kt, F.layout(), g, rng, 1 + i % 12);
      one(g);
    }
  }
  // Coverage, not decoration: without these the equalities above could be 201
  // checks of the same trivial schedule, and -- measured the hard way, by
  // mutating the cost path and watching this pass -- an `fx > 0` bound alone
  // does NOT imply the extracted cluster is ever a merge, which is the branch
  // that feeds `demanded` from L2.
  CHECK(totals.size() > 1);
  CHECK(shape.sum > 0);
  CHECK(shape.fx > 0);
  CHECK(shape.fx_wide > 0);
  CHECK(shape.cl_wide > 0);
}

/// The bit-stripping canon_less that genome.hpp's branchless version replaced;
/// kept here so the replacement stays pinned to it exhaustively (T-A1).
bool canon_less_ref(ga::NodeMask a, ga::NodeMask b) {
  const int pa = std::popcount(a), pb = std::popcount(b);
  if (pa != pb) return pa < pb;
  while (a && b) {
    const int la = std::countr_zero(a), lb = std::countr_zero(b);
    if (la != lb) return la < lb;
    a &= a - 1;
    b &= b - 1;
  }
  return false;
}

/// The pre-T-A2 `ins` body: rebuild every member, then a full canon_sort.
/// Kept here so the merge-based replacement stays pinned to it exhaustively.
void ins_ref(ga::Laminar& fam, ga::NodeMask T, ga::NodeMask l) {
  ga::Laminar out;
  out.reserve(fam.size() + 2);
  for (ga::NodeMask C : fam) {
    const bool superset = (C & T) == T;
    if (!(superset && C != T)) out.push_back(C);  // not a proper superset
    if (superset) out.push_back(C | l);
  }
  out.push_back(l);
  ga::canon_sort(out);
  fam = std::move(out);
}

/// decode_tree driven by ins_ref instead of ins.
ga::Laminar decode_tree_ref(ga::TreeCode const& code, int n) {
  ga::Laminar fam{ga::NodeMask{1}};
  for (int k = 2; k <= n; ++k)
    ins_ref(fam, fam[code[k - 2]], ga::NodeMask{1} << (k - 1));
  return fam;
}

/// The pre-T-A2b `nni_moves` body: copy, erase A, push_back the new cluster,
/// then a full canon_sort of the whole family. Kept here so the splicing
/// replacement stays pinned to it -- including the ORDER of the move list,
/// which `hill_climb` reads (it breaks ties first-wins on a strict `<`).
container::svector<ga::Laminar> nni_moves_ref(ga::Laminar const& fam) {
  container::svector<ga::Laminar> out;
  for (ga::NodeMask B : fam) {
    if (std::popcount(B) < 2) continue;
    const auto [b1, b2] = ga::tree_children(fam, B);
    for (int i = 0; i < 2; ++i) {
      const ga::NodeMask A = i ? b2 : b1, A3 = i ? b1 : b2;
      if (std::popcount(A) < 2) continue;
      const auto [A1, A2] = ga::tree_children(fam, A);
      for (ga::NodeMask keep : {A2, A1}) {  // swap the other child out
        ga::Laminar next = fam;
        next.erase(std::find(next.begin(), next.end(), A));
        next.push_back(keep | A3);
        ga::canon_sort(next);
        out.push_back(std::move(next));
      }
    }
  }
  return out;
}

/// The pre-T-A2b `encode_tree` body: rebuild the peeled family, canon_sort it,
/// then std::unique it. Kept here so the 2-way-merge replacement stays pinned.
ga::TreeCode encode_tree_ref(ga::Laminar fam, int n) {
  ga::TreeCode code(n - 1);
  for (int k = n; k >= 2; --k) {
    const ga::NodeMask sl = ga::NodeMask{1} << (k - 1);
    ga::NodeMask P = 0;
    for (ga::NodeMask C : fam)
      if ((C & sl) && C != sl && (!P || std::popcount(C) < std::popcount(P)))
        P = C;
    const ga::NodeMask T = P & ~sl;
    ga::Laminar next;
    for (ga::NodeMask C : fam)
      if (C != P && C != sl) next.push_back(C & ~sl);
    ga::canon_sort(next);
    next.erase(std::unique(next.begin(), next.end()), next.end());
    const auto it =
        std::lower_bound(next.begin(), next.end(), T, ga::canon_less);
    code[k - 2] = static_cast<int>(it - next.begin());
    fam = std::move(next);
  }
  return code;
}

}  // namespace

TEST_CASE("ga genome codec", "[ga]") {
  using namespace ga;
  SECTION("canon_less matches the bit-stripping reference, exhaustive 10 bit") {
    // 2^20 ordered pairs, reported as one aggregate CHECK so a regression does
    // not emit a million assertions.
    bool all_equal = true;
    NodeMask bad_a = 0, bad_b = 0;
    for (NodeMask a = 0; a < 1024 && all_equal; ++a)
      for (NodeMask b = 0; b < 1024; ++b)
        if (canon_less(a, b) != canon_less_ref(a, b)) {
          bad_a = a;
          bad_b = b;
          all_equal = false;
          break;
        }
    INFO("first mismatch at a = " << bad_a << ", b = " << bad_b);
    CHECK(all_equal);
  }
  SECTION("ins merges exactly like the re-sorting reference, exhaustive n <= 7") {
    // T-A2: `ins` replaced canon_sort by a 3-way merge of the three sorted runs
    // it produces. Pin it to the old body over EVERY code of EVERY length
    // n <= 7 -- 1 + 3 + 15 + 105 + 945 + 10395 = 11464 codes, i.e. every
    // (2n-3)!! binary tree reached through every insertion order -- plus a
    // handful of n = 55 codes, the size of the real L2 block (the exhaustive
    // part only reaches 13-mask families, so it cannot exercise deep chains).
    // Reported as aggregate CHECKs rather than 11464+ assertions.
    std::size_t n_codes = 0;
    bool all_equal = true;
    for (int n = 2; n <= 7 && all_equal; ++n) {
      container::svector<int> ranges;
      for (int k = 2; k <= n; ++k) ranges.push_back(2 * k - 3);
      TreeCode code(n - 1, 0);
      while (true) {
        ++n_codes;
        if (decode_tree(code, n) != decode_tree_ref(code, n)) {
          INFO("n = " << n << ", gene 0 = " << code[0]);
          all_equal = false;
          break;
        }
        int i = 0;
        for (; i < n - 1; ++i) {
          if (++code[i] < ranges[i]) break;
          code[i] = 0;
        }
        if (i == n - 1) break;
      }
    }
    CHECK(n_codes == 11464u);
    CHECK(all_equal);

    constexpr int big_n = 55;
    std::uint64_t rng = 0x9e3779b97f4a7c15ull;
    bool big_equal = true;
    for (int rep = 0; rep < 256 && big_equal; ++rep) {
      TreeCode code(big_n - 1);
      for (int k = 2; k <= big_n; ++k) {
        rng = rng * 6364136223846793005ull + 1442695040888963407ull;
        code[k - 2] = static_cast<int>((rng >> 33) % (2 * k - 3));
      }
      big_equal = decode_tree(code, big_n) == decode_tree_ref(code, big_n);
    }
    CHECK(big_equal);
  }
  SECTION("nni_moves/encode_tree splice like the sorting refs, exhaustive n<=7") {
    // T-A2b: `nni_moves` replaced its per-move canon_sort by a single splice
    // (erase A, insert keep|A3 at its lower_bound position), and `encode_tree`
    // replaced its per-leaf canon_sort + unique by a 2-way merge of the two
    // sorted, provably disjoint runs the peel produces. Both are pinned to the
    // old bodies over EVERY code of EVERY length n <= 7 -- the same 11464
    // codes the `ins` section walks, i.e. every (2n-3)!! tree through every
    // insertion order. The move list is compared ELEMENT FOR ELEMENT IN ORDER:
    // `hill_climb` breaks ties first-wins on a strict `<`, so a permuted move
    // list would silently select a different block. `encode_tree` is checked
    // both on the tree itself and on each of its NNI neighbours -- the latter
    // is exactly how hill_climb/nni_kick call it. Then a few n = 55 codes, the
    // real L2 block size, since n <= 7 only reaches 13-mask families and
    // cannot exercise deep ancestor chains. Aggregate CHECKs, not millions of
    // assertions.
    std::size_t n_codes = 0, n_moves = 0;
    bool moves_equal = true, codes_equal = true;
    for (int n = 2; n <= 7 && moves_equal && codes_equal; ++n) {
      container::svector<int> ranges;
      for (int k = 2; k <= n; ++k) ranges.push_back(2 * k - 3);
      TreeCode code(n - 1, 0);
      while (true) {
        ++n_codes;
        const Laminar fam = decode_tree(code, n);
        auto const mv = nni_moves(fam);
        if (mv != nni_moves_ref(fam)) moves_equal = false;
        if (encode_tree(fam, n) != encode_tree_ref(fam, n)) codes_equal = false;
        for (auto const& m : mv) {
          ++n_moves;
          if (encode_tree(m, n) != encode_tree_ref(m, n)) codes_equal = false;
        }
        if (!moves_equal || !codes_equal) {
          INFO("n = " << n << ", gene 0 = " << code[0]);
          break;
        }
        int i = 0;
        for (; i < n - 1; ++i) {
          if (++code[i] < ranges[i]) break;
          code[i] = 0;
        }
        if (i == n - 1) break;
      }
    }
    CHECK(n_codes == 11464u);
    // every tree on n leaves has exactly 2(n-2) NNI moves (one per internal
    // non-root node, times the two children it can swap out), so this is
    // 3*2 + 15*4 + 105*6 + 945*8 + 10395*10.
    CHECK(n_moves == 112206u);
    CHECK(moves_equal);
    CHECK(codes_equal);

    constexpr int big_n = 55;
    std::uint64_t rng = 0x9e3779b97f4a7c15ull;
    bool big_equal = true;
    std::size_t big_moves = 0;
    for (int rep = 0; rep < 16 && big_equal; ++rep) {
      TreeCode code(big_n - 1);
      for (int k = 2; k <= big_n; ++k) {
        rng = rng * 6364136223846793005ull + 1442695040888963407ull;
        code[k - 2] = static_cast<int>((rng >> 33) % (2 * k - 3));
      }
      const Laminar fam = decode_tree(code, big_n);
      auto const mv = nni_moves(fam);
      big_moves += mv.size();
      if (mv != nni_moves_ref(fam)) big_equal = false;
      if (encode_tree(fam, big_n) != encode_tree_ref(fam, big_n))
        big_equal = false;
      for (auto const& m : mv)
        if (encode_tree(m, big_n) != encode_tree_ref(m, big_n))
          big_equal = false;
    }
    CHECK(big_moves == 1696u);  // 16 x 2 x (53 internal non-root nodes)
    CHECK(big_equal);
  }
  SECTION("nni_move_count/nni_move_at == the move list, exhaustive n <= 7") {
    // T-A2c: `nni_kick` materialized every NNI neighbour of a block -- ~104
    // families of 109 masks for the 55-leaf L2 block -- and then used exactly
    // ONE, `moves[rng.below(moves.size())]`. It now counts the neighbours,
    // makes that same draw against the count, and builds only the drawn one.
    // That is bit-identical only if BOTH halves are exact:
    //
    //   * `nni_move_count(fam) == nni_moves(fam).size()`, or the draw is over
    //     a different range and every later decision in the search shifts --
    //     silently, since a genome with the wrong move applied is still a
    //     perfectly valid genome;
    //   * `nni_move_at(fam, k) == nni_moves(fam)[k]` for EVERY k, or the drawn
    //     index names a different neighbour.
    //
    // Both are pinned here against `nni_moves_ref` -- the pre-T-A2b body, a
    // full canon_sort per neighbour -- rather than against `nni_moves`, so the
    // chain runs all the way back to the original code. Driven over EVERY code
    // of EVERY length n <= 7, the same 11464 = 1+3+15+105+945+10395
    // mixed-radix codes the two sections above walk (every (2n-3)!! tree
    // through every insertion order), at every index of every one of them;
    // then a few n = 55 codes, the real L2 block size, since n <= 7 only
    // reaches 13-mask families. Aggregate CHECKs, not 100k assertions.
    std::size_t n_codes = 0, n_moves = 0;
    bool counts_equal = true, at_equal = true;
    for (int n = 2; n <= 7 && counts_equal && at_equal; ++n) {
      container::svector<int> ranges;
      for (int k = 2; k <= n; ++k) ranges.push_back(2 * k - 3);
      TreeCode code(n - 1, 0);
      while (true) {
        ++n_codes;
        const Laminar fam = decode_tree(code, n);
        auto const ref = nni_moves_ref(fam);
        if (nni_move_count(fam) != ref.size()) counts_equal = false;
        // the ChildTable-driven walk (what nni_kick/hill_climb run, with the
        // table from the decode memo) is pinned against the same reference as
        // the fam-only forms, over the same exhaustive code sweep
        ChildTable ch;
        build_children(fam, ch);
        if (nni_moves(fam, ch) != ref) at_equal = false;
        for (std::size_t k = 0; k < ref.size(); ++k) {
          ++n_moves;
          if (nni_move_at(fam, k) != ref[k]) at_equal = false;
          if (nni_move_at(fam, ch, k) != ref[k]) at_equal = false;
        }
        if (!counts_equal || !at_equal) {
          INFO("n = " << n << ", gene 0 = " << code[0]);
          break;
        }
        int i = 0;
        for (; i < n - 1; ++i) {
          if (++code[i] < ranges[i]) break;
          code[i] = 0;
        }
        if (i == n - 1) break;
      }
    }
    CHECK(n_codes == 11464u);
    // same tally as the section above: 3*2 + 15*4 + 105*6 + 945*8 + 10395*10
    CHECK(n_moves == 112206u);
    CHECK(counts_equal);
    CHECK(at_equal);

    constexpr int big_n = 55;
    std::uint64_t rng = 0x9e3779b97f4a7c15ull;
    bool big_equal = true;
    std::size_t big_moves = 0;
    for (int rep = 0; rep < 16 && big_equal; ++rep) {
      TreeCode code(big_n - 1);
      for (int k = 2; k <= big_n; ++k) {
        rng = rng * 6364136223846793005ull + 1442695040888963407ull;
        code[k - 2] = static_cast<int>((rng >> 33) % (2 * k - 3));
      }
      const Laminar fam = decode_tree(code, big_n);
      auto const ref = nni_moves_ref(fam);
      if (nni_move_count(fam) != ref.size()) big_equal = false;
      for (std::size_t k = 0; k < ref.size() && big_equal; ++k) {
        ++big_moves;
        if (nni_move_at(fam, k) != ref[k]) big_equal = false;
      }
    }
    CHECK(big_moves == 1696u);  // 16 x 2 x (53 internal non-root nodes)
    CHECK(big_equal);
  }
  SECTION("bijection onto binary trees, exhaustive for n <= 6") {
    for (int n = 2; n <= 6; ++n) {
      container::svector<int> ranges;
      for (int k = 2; k <= n; ++k) ranges.push_back(2 * k - 3);
      container::svector<Laminar> seen;
      TreeCode code(n - 1, 0);
      while (true) {
        Laminar fam = decode_tree(code, n);
        REQUIRE(fam.size() == std::size_t(2 * n - 1));
        REQUIRE(encode_tree(fam, n) == code);
        if (ranges::find(seen, fam) == ranges::end(seen))
          seen.push_back(fam);
        int i = 0;
        for (; i < n - 1; ++i) {
          if (++code[i] < ranges[i]) break;
          code[i] = 0;
        }
        if (i == n - 1) break;
      }
      REQUIRE(seen.size() == double_factorial(n));
    }
  }
  SECTION("NNI moves are valid trees and round-trip") {
    TreeCode code{0, 2, 3, 6};  // some tree on 5 leaves
    auto fam = decode_tree(code, 5);
    auto moves = nni_moves(fam);
    REQUIRE(!moves.empty());
    for (auto const& mv : moves) {
      REQUIRE(mv.size() == fam.size());
      REQUIRE(decode_tree(encode_tree(mv, 5), 5) == mv);
      REQUIRE(mv != fam);
    }
  }
  SECTION("decode memo reproduces a fresh decode, exhaustive n <= 7") {
    // T-A3: the memo is keyed on the FULL code slice (the 64-bit hash only
    // picks the probe chain; `n` and every gene are compared before a hit is
    // accepted), so a hit must be byte-identical to decode_tree +
    // build_children -- same masks, same ORDER, same children pairing. Drive
    // every code of every length n <= 7 (the same 11464 codes the `ins`
    // section uses) through three memos: a default one twice (pass 0 = all
    // misses, pass 1 = all hits), one capped at 64 entries so clear-on-full is
    // exercised, and one that is asked for `fam` only first so the lazy
    // children fill-in path runs. Aggregate CHECKs, not 11464 assertions.
    TreeMemo memo, tiny(64), lazy;
    bool ok = true;
    std::size_t n_codes = 0;
    for (int n = 2; n <= 7 && ok; ++n) {
      container::svector<int> ranges;
      for (int k = 2; k <= n; ++k) ranges.push_back(2 * k - 3);
      TreeCode code(n - 1, 0);
      while (true) {
        ++n_codes;
        const Laminar ref = decode_tree(code, n);
        ChildTable ref_ch;
        build_children(ref, ref_ch);
        Laminar fam;
        ChildTable ch;
        for (TreeMemo* m : {&memo, &tiny, &memo})
          for (int pass = 0; pass < 2; ++pass) {
            m->decode(code.data(), n, fam, ch);
            if (fam != ref || ch != ref_ch) ok = false;
          }
        lazy.decode(code.data(), n, fam);  // fam only: no children stored
        if (fam != ref) ok = false;
        lazy.decode(code.data(), n, fam, ch);  // ... now fill them in
        if (fam != ref || ch != ref_ch) ok = false;
        lazy.decode(code.data(), n, fam, ch);  // ... and read them back
        if (fam != ref || ch != ref_ch) ok = false;
        if (!ok) {
          INFO("n = " << n << ", gene 0 = " << code[0]);
          break;
        }
        int i = 0;
        for (; i < n - 1; ++i) {
          if (++code[i] < ranges[i]) break;
          code[i] = 0;
        }
        if (i == n - 1) break;
      }
    }
    CHECK(n_codes == 11464u);
    CHECK(ok);
    CHECK(memo.hits() > 0);
    CHECK(tiny.clears() > 0);  // the wholesale-eviction path really ran
  }
  SECTION("memo seeding is sound: encode(mv) decodes back to mv, n <= 6") {
    // T-A3 seeds the memo in nni_kick/hill_climb, where `encode_tree(mv, n)`
    // is computed with the family `mv` already in hand. That entry is only
    // sound if the codec round-trips ORDER-exactly -- otherwise the memo would
    // answer a later decode of that code with a DIFFERENT tree and silently
    // move the objective. Pin it over every NNI move of every tree on n <= 6
    // leaves, and check the seeded memo (whose children table is filled
    // lazily, since the seeder has none) matches a fresh decode.
    bool ok = true;
    std::size_t n_moves = 0;
    TreeMemo memo;
    for (int n = 3; n <= 6 && ok; ++n) {
      container::svector<int> ranges;
      for (int k = 2; k <= n; ++k) ranges.push_back(2 * k - 3);
      TreeCode code(n - 1, 0);
      while (true) {
        for (auto const& mv : nni_moves(decode_tree(code, n))) {
          ++n_moves;
          auto const cand = encode_tree(mv, n);
          if (decode_tree(cand, n) != mv) ok = false;
          memo.seed(cand.data(), n, mv);
          Laminar fam;
          ChildTable ch, ref_ch;
          build_children(mv, ref_ch);
          memo.decode(cand.data(), n, fam, ch);
          if (fam != mv || ch != ref_ch) ok = false;
          if (!ok) break;
        }
        if (!ok) {
          INFO("n = " << n << ", gene 0 = " << code[0]);
          break;
        }
        int i = 0;
        for (; i < n - 1; ++i) {
          if (++code[i] < ranges[i]) break;
          code[i] = 0;
        }
        if (i == n - 1) break;
      }
    }
    // 2(n-2) moves per tree (one per internal non-root node, two ways),
    // (2n-3)!! codes per n: 3*2 + 15*4 + 105*6 + 945*8.
    CHECK(n_moves == 8256u);
    CHECK(ok);
  }
}

TEST_CASE("ga python parity", "[ga]") {
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg);

  auto const terms = load_md_terms(SEQUANT_UNIT_TESTS_SOURCE_DIR
                                   "/ga_pao_ccsd_df.md");
  REQUIRE(terms.size() == 86);

  SECTION("f-only universe") {
    auto const& ref = ga_reference::fonly;
    auto targets = universe(terms, 0);
    auto kt = ga::build_key_table(targets, python_extent);
    REQUIRE(kt.terms.size() == ref.n_terms);
    std::size_t nodes = 0, pairs = 0;
    for (auto const& T : kt.terms) {
      nodes += T.n();
      pairs += (std::size_t{1} << T.n()) - 1;
    }
    CHECK(nodes == ref.n_nodes);
    CHECK(pairs == ref.key_pairs);
    CHECK(kt.n_keys == ref.n_keys);

    ga::Fitness F(kt, ga::CostModel::python_parity(),
                  ga::ProducerResolution::Greedy);
    ga::Fitness Fx(kt, ga::CostModel::python_parity(),
                   ga::ProducerResolution::Exact);

    // per-term optimal tree costs, on the prototype's own optimal trees
    for (std::size_t d = 0; d < kt.terms.size(); ++d) {
      auto const [lo, hi] = F.layout().g_slice[d];
      ga::TreeCode block(ref.g_opt.begin() + lo, ref.g_opt.begin() + hi);
      CHECK(tree_cost(kt, F.cost(), static_cast<int>(d), block) ==
            ref.per_term_opt[d]);
    }

    auto seed = slice_genome(ref.g_opt, ref.h0);
    CHECK(F(seed) == ref.seed_greedy);
    CHECK(Fx(seed) == ref.seed_exact);

    auto control = slice_genome(ref.g_ctrl, ref.h_ctrl);
    CHECK(Fx(control) == ref.control_exact);
    CHECK(F(control) == ref.control_greedy);

    // our own seed reproduces the same per-term optima (values, not trees),
    // TERM BY TERM. Not a sum: a sum lets one term's loss cancel another's
    // gain, and per-term optimality is precisely what `seed_genome`'s mask DP
    // (T-B1) has to preserve -- it replaced `single_term_opt`, whose per-subset
    // index sets come from a different closure convention.
    auto mine = ga::seed_genome(kt);
    for (std::size_t d = 0; d < kt.terms.size(); ++d) {
      auto const [lo, hi] = F.layout().g_slice[d];
      CHECK(tree_cost(kt, F.cost(), static_cast<int>(d),
                      {mine.g.begin() + lo, mine.g.begin() + hi}) ==
            ref.per_term_opt[d]);
    }
  }

  SECTION("g-only universe") {
    auto const& ref = ga_reference::gonly;
    auto targets = universe(terms, 1);
    auto kt = ga::build_key_table(targets, python_extent);
    REQUIRE(kt.terms.size() == ref.n_terms);
    std::size_t nodes = 0, pairs = 0;
    for (auto const& T : kt.terms) {
      nodes += T.n();
      pairs += (std::size_t{1} << T.n()) - 1;
    }
    CHECK(nodes == ref.n_nodes);
    CHECK(pairs == ref.key_pairs);
    CHECK(kt.n_keys == ref.n_keys);

    ga::Fitness F(kt, ga::CostModel::python_parity(),
                  ga::ProducerResolution::Greedy);

    for (std::size_t d = 0; d < kt.terms.size(); ++d) {
      auto const [lo, hi] = F.layout().g_slice[d];
      ga::TreeCode block(ref.g_opt.begin() + lo, ref.g_opt.begin() + hi);
      CHECK(tree_cost(kt, F.cost(), static_cast<int>(d), block) ==
            ref.per_term_opt[d]);
    }

    auto seed = slice_genome(ref.g_opt, ref.h0);
    CHECK(F(seed) == ref.seed_greedy);

    // our own seed reproduces the same per-term optima, term by term (see the
    // f-only section)
    auto mine = ga::seed_genome(kt);
    for (std::size_t d = 0; d < kt.terms.size(); ++d) {
      auto const [lo, hi] = F.layout().g_slice[d];
      CHECK(tree_cost(kt, F.cost(), static_cast<int>(d),
                      {mine.g.begin() + lo, mine.g.begin() + hi}) ==
            ref.per_term_opt[d]);
    }

    // the prototype's GA winner decodes to the same schedule and cost
    auto winner = slice_genome(ref.winner_g, ref.winner_h);
    auto sch = F.explain(winner);
    CHECK(sch.total == ref.winner_total);
    CHECK(sch.l1 == ref.winner_l1);
    CHECK(sch.l2 == ref.winner_l2);
    CHECK(sch.pick.size() == ref.winner_n_arrays);
  }
}

// T-A5 gave `Fitness::operator()` its own cost-only L2 recursion
// (`build_node_cost` + `cost_of_cost_val` over `CostVal`s) so the ~164k
// evaluations of a search stop materializing a `Schedule` they throw away.
// That leaves TWO evaluation paths in the file, and this test case is the
// standing guarantee that they never drift apart: it must be extended by any
// change to either path and must never be relaxed to a tolerance. See the
// comments on `check_cost_path_agrees` and on `Fitness::operator()`.
TEST_CASE("ga cost path agrees with explain", "[ga]") {
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg);
  auto const terms = load_md_terms(SEQUANT_UNIT_TESTS_SOURCE_DIR
                                   "/ga_pao_ccsd_df.md");
  REQUIRE(terms.size() == 86);

  SECTION("f-only universe") {
    auto const& ref = ga_reference::fonly;
    auto kt = ga::build_key_table(universe(terms, 0), python_extent);
    ga::Fitness F(kt, ga::CostModel::python_parity());
    // the prototype's control schedule is included as a base because it is the
    // one genome known to drive extraction hard (see the coverage CHECKs)
    check_cost_path_agrees(kt, F,
                           {ga::seed_genome(kt),
                            slice_genome(ref.g_ctrl, ref.h_ctrl),
                            slice_genome(ref.g_opt, ref.h0)},
                           0x5eed'0f'01ull);
  }

  SECTION("f-only universe, Exact producer resolution") {
    // the Exact branch of `resolve` is the one `operator()` now enters with
    // pick_out == nullptr, so it gets its own pass
    auto const& ref = ga_reference::fonly;
    auto kt = ga::build_key_table(universe(terms, 0), python_extent);
    ga::Fitness F(kt, ga::CostModel::python_parity(),
                  ga::ProducerResolution::Exact);
    check_cost_path_agrees(
        kt, F,
        {ga::seed_genome(kt), slice_genome(ref.g_ctrl, ref.h_ctrl)},
        0x5eed'0f'02ull);
  }

  SECTION("g-only universe") {
    auto const& ref = ga_reference::gonly;
    auto kt = ga::build_key_table(universe(terms, 1), python_extent);
    ga::Fitness F(kt, ga::CostModel::python_parity());
    check_cost_path_agrees(kt, F,
                           {ga::seed_genome(kt),
                            slice_genome(ref.winner_g, ref.winner_h),
                            slice_genome(ref.g_opt, ref.h0)},
                           0x5eed'09'01ull);
  }

  SECTION("81-term R1+R2 universe at DZ extents") {
    // The universe the GATE-B fingerprint and the mpqc4 CSV-CCSD job are taken
    // on: R2's wide terms are where the deep Fx chains and the CSV composite
    // faces live, so it is the only one that covers them.
    //
    // The base genome is deliberately NOT `seed_genome`: its per-term DP costs
    // ~21 s here and this tag runs before every commit. All-zeros is a valid
    // deterministic genome and the kicks reach the same code paths -- what is
    // under test is the L2 recursion, not the quality of the trees. The
    // coverage CHECKs at the bottom of `check_cost_path_agrees` are what keep
    // that claim honest.
    auto kt = ga::build_key_table(universe(terms, 3), dz_extent);
    REQUIRE(kt.terms.size() == 81);
    ga::Fitness F(kt);
    // The RNG seed is chosen, not arbitrary: multi-tensor extraction is rare
    // in this universe (most seeds give `fx_wide == 0` over 201 genomes), and
    // the coverage CHECK is what stops the seed from silently going stale.
    check_cost_path_agrees(kt, F, {zero_genome(F.layout())}, 0x5eed'd2'08ull);
  }
}

// The bit tables the L2 pass runs on since T-A6, and above all the ORDER one of
// them carries.
//
// `Fitness::find_beta` enumerates candidate face bijections per kind with
// `std::next_permutation` and returns the FIRST one satisfying (Sigma1-3);
// which bijection that is decides how emission renames the second operand of
// every extraction. It used to enumerate over `std::sort`ed `Index` objects and
// now enumerates over index bits ordered by `TermTable::index_rank`. Bit order
// is first-encounter order while walking the term's tensors and agrees with
// `Index::operator<` only by accident, so getting `index_rank` wrong would keep
// every flop count -- and quite possibly the GATE-B fingerprint -- identical
// while rewriting the emitted schedule. Hence this pins the order directly,
// against an independently computed `std::sort`, rather than trusting the
// emission tests to notice.
//
// (`cmake-build-release` sets SEQUANT_ASSERT_BEHAVIOR=IGNORE, so the same check
// inside `build_term` does NOT run under the [ga] tag. This does.)
TEST_CASE("ga index bit tables", "[ga]") {
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg);
  auto const terms = load_md_terms(SEQUANT_UNIT_TESTS_SOURCE_DIR
                                   "/ga_pao_ccsd_df.md");
  REQUIRE(terms.size() == 86);

  // counters, not per-item REQUIREs: the sweeps below run over every subset of
  // every term
  std::size_t bad_rank = 0, bad_proto = 0, bad_kind = 0, bad_face = 0,
              bad_group = 0, bad_canon = 0;
  std::size_t composites = 0, multi_kind_faces = 0, ranks_off_bit = 0,
              canon_slots = 0, canon_off_bit = 0;

  auto check = [&](ga::KeyTable const& kt) {
    for (auto const& T : kt.terms) {
      const std::size_t ni = T.index_list.size();
      REQUIRE(T.index_rank.size() == ni);
      REQUIRE(T.proto_mask.size() == ni);
      REQUIRE(T.kind_mask.size() == kt.kind_ids.size());

      // (1) index_rank IS the position in `std::sort`-by-`Index` order --
      //     recomputed here from the real `Index` objects, independently of
      //     how build_term derived it
      container::svector<Index> sorted(T.index_list.begin(),
                                       T.index_list.end());
      std::sort(sorted.begin(), sorted.end());
      for (std::size_t b = 0; b < ni; ++b) {
        if (!(sorted[T.index_rank[b]] == T.index_list[b])) ++bad_rank;
        if (T.index_rank[b] != b) ++ranks_off_bit;  // coverage, see below
      }

      // (2) proto_mask / kind_mask reproduce the Index-level tables
      for (std::size_t b = 0; b < ni; ++b) {
        ga::IndexMask pm = 0;
        for (auto const& p : T.index_list[b].proto_indices())
          pm |= ga::IndexMask{1} << T.bit_of(p);
        if (T.proto_mask[b] != pm) ++bad_proto;
        if (pm) ++composites;
        if (T.kind[b] != kt.kind_of(T.index_list[b]) ||
            !((T.kind_mask[static_cast<std::size_t>(T.kind[b])] >> b) & 1))
          ++bad_kind;
      }
      for (std::size_t k = 0; k < T.kind_mask.size(); ++k)
        for (ga::IndexMask m = T.kind_mask[k]; m; m &= m - 1)
          if (T.kind[std::countr_zero(m)] != static_cast<int>(k)) ++bad_kind;

      // per-node effective index sets, at `Index` level -- the input to the
      // independent recomputation of F(S) below
      const std::size_t nv = T.n();
      container::svector<ga::FaceSet> eff(nv);
      for (std::size_t v = 0; v < nv; ++v)
        for (auto const& ix : T.tensors[v]->as<Tensor>().const_indices()) {
          for (auto const& p : ix.proto_indices()) eff[v].emplace(p);
          eff[v].emplace(ix);
        }

      // (3) per subset: the masks ARE the face. Since T-A7 `face_mask` is the
      //     table's ONLY record of F(S) -- `container::vector<FaceSet> face`
      //     is gone -- so F(S) is recomputed here straight from the
      //     definition (eff-closure members that are external or also carried
      //     outside S) and compared against what `face_set` reconstructs.
      //     Then: the by-kind grouping find_beta builds out of the mask is
      //     exactly the one the old `container::map<int, svector<Index>>` +
      //     `std::sort` produced.
      for (std::size_t S = 1; S < T.canon_face_bits.size(); ++S) {
        ga::FaceSet ref;  // F(S), from the definition
        for (std::size_t v = 0; v < nv; ++v) {
          if (!(S & (std::size_t{1} << v))) continue;
          for (auto const& ix : eff[v]) {
            bool keep = T.ext.find(ix) != T.ext.end();
            for (std::size_t w = 0; !keep && w < nv; ++w)
              if (!(S & (std::size_t{1} << w)) && eff[w].find(ix) != eff[w].end())
                keep = true;
            if (keep) ref.emplace(ix);
          }
        }
        auto const fs = T.face_set(static_cast<ga::NodeMask>(S));
        if (fs != ref) ++bad_face;
        const ga::IndexMask fm = T.face_mask(static_cast<ga::NodeMask>(S));
        if (static_cast<std::size_t>(std::popcount(fm)) != fs.size())
          ++bad_face;
        for (auto const& ix : fs)
          if (!((fm >> T.bit_of(ix)) & 1)) ++bad_face;

        container::map<int, container::svector<Index>> by;  // the old grouping
        for (auto const& ix : fs) by[kt.kind_of(ix)].push_back(ix);
        for (auto& [k, v] : by) std::sort(v.begin(), v.end());
        if (by.size() > 1) ++multi_kind_faces;
        std::size_t ki = 0;
        for (std::size_t k = 0; k < T.kind_mask.size(); ++k) {  // the new one
          container::svector<std::uint8_t> bits;
          for (ga::IndexMask m = fm & T.kind_mask[k]; m; m &= m - 1)
            bits.push_back(static_cast<std::uint8_t>(std::countr_zero(m)));
          if (bits.empty()) continue;
          std::sort(bits.begin(), bits.end(),
                    [&](std::uint8_t a, std::uint8_t b) {
                      return T.index_rank[a] < T.index_rank[b];
                    });
          // kinds must come out ascending and in the same blocks...
          if (ki >= by.size() ||
              (by.begin() + ki)->first != static_cast<int>(k))
            ++bad_group;
          else {
            auto const& v = (by.begin() + ki)->second;
            if (bits.size() != v.size()) {
              ++bad_group;
            } else {  // ...and each block in the same ORDER
              for (std::size_t i = 0; i < bits.size(); ++i)
                if (!(T.index_list[bits[i]] == v[i])) ++bad_group;
            }
          }
          ++ki;
        }
        if (ki != by.size()) ++bad_group;

        // (4) canon_face_bits is the canonicalizer's slot list, IN ORDER.
        //     Slot position is the array axis (emit.cpp's `named_leaf` gives
        //     every slot Nonsymm/aux so nothing may reorder them), and since
        //     T-A7 build_term stores bit ids rather than the `Index` objects
        //     themselves -- so re-run canonicalize_slots here, independently,
        //     and compare the reconstructed `Index` list ELEMENT FOR ELEMENT.
        //     A permutation of this list changes every emitted array layout
        //     while leaving all costs, and hence the GATE-B fingerprint,
        //     untouched.
        container::vector<ExprPtr> ts;
        for (std::size_t v = 0; v < nv; ++v)
          if (S & (std::size_t{1} << v))
            ts.emplace_back(T.tensors[v]->as<Tensor>().clone());
        auto tn = TensorNetwork{ts};
        auto meta = tn.canonicalize_slots(
            TensorCanonicalizer::cardinal_tensor_labels(), &ref);
        auto const cf = T.canon_face_indices(static_cast<ga::NodeMask>(S));
        if (cf.size() != meta.named_indices_canonical.size()) {
          ++bad_canon;
        } else {
          for (std::size_t i = 0; i < cf.size(); ++i)
            if (!(cf[i] == *meta.named_indices_canonical[i])) ++bad_canon;
          canon_slots += cf.size();
          // coverage: the canonical order is a REORDERING of the face bits,
          // so "store the face bits ascending" would not reproduce it
          auto const& cb = T.canon_face_bits[S];
          if (!std::is_sorted(cb.begin(), cb.end())) ++canon_off_bit;
        }
      }
    }
  };

  check(ga::build_key_table(universe(terms, 0), python_extent));
  check(ga::build_key_table(universe(terms, 1), python_extent));

  CHECK(bad_rank == 0);
  CHECK(bad_proto == 0);
  CHECK(bad_kind == 0);
  CHECK(bad_face == 0);
  CHECK(bad_group == 0);
  CHECK(bad_canon == 0);
  // Coverage, so the sweep cannot go vacuous: composite (CSV/PAO) slots with
  // proto-indices are present, faces with more than one kind are present, and
  // -- the point of the whole test case -- `index_rank` really does disagree
  // with bit order, so (1) and (3) are not tautologies.
  CHECK(composites > 0);
  CHECK(multi_kind_faces > 0);
  CHECK(ranks_off_bit > 0);
  // ...and the slot-order sweep (4) saw real slots, in orders that are not
  // simply the ascending bit order, so it is not a tautology either.
  CHECK(canon_slots > 0);
  CHECK(canon_off_bit > 0);
}

// T-C1. build_key_table canonicalizes subsets in parallel but hands out key ids
// serially, in ascending-subset order, so the table it produces must not depend
// on the thread count at all. The GATE-B fingerprint is the end-to-end check on
// that; this one localizes a failure to the table itself -- and to WHICH of
// key / kind / canonical face order moved -- rather than to "the winner
// changed". Ids are dense and first-encounter-ordered, so a thread-count
// dependence here is not a cosmetic renumbering: it reorders `resolve`'s
// ascending-key walk, hence the floating-point summation order.
TEST_CASE("ga key table does not depend on the thread count", "[ga]") {
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg);
  auto const terms = load_md_terms(SEQUANT_UNIT_TESTS_SOURCE_DIR
                                   "/ga_pao_ccsd_df.md");

  // put back whatever the rest of the suite runs at, on the failure path too
  struct RestoreThreads {
    int n;
    ~RestoreThreads() { set_num_threads(n); }
  } const restore{num_threads()};

  std::size_t subsets_seen = 0, slots_seen = 0, parallel_terms = 0;

  auto compare = [&](int mode) {
    auto const targets = universe(terms, mode);
    set_num_threads(1);
    auto const a = ga::build_key_table(targets, python_extent);
    set_num_threads(8);
    auto const b = ga::build_key_table(targets, python_extent);

    REQUIRE(a.n_keys == b.n_keys);
    REQUIRE(a.terms.size() == b.terms.size());
    // kind ids are handed out in first-encounter order too (register_kind), so
    // the whole map must match key for key AND value for value
    REQUIRE(a.kind_ids.size() == b.kind_ids.size());
    CHECK(std::equal(a.kind_ids.begin(), a.kind_ids.end(), b.kind_ids.begin()));

    std::size_t bad_key = 0, bad_canon = 0, bad_kind = 0;
    for (std::size_t t = 0; t < a.terms.size(); ++t) {
      auto const& A = a.terms[t];
      auto const& B = b.terms[t];
      REQUIRE(A.key.size() == B.key.size());
      REQUIRE(A.canon_face_bits.size() == A.key.size());
      REQUIRE(B.canon_face_bits.size() == B.key.size());
      REQUIRE(A.kind.size() == B.kind.size());
      for (std::size_t i = 0; i < A.kind.size(); ++i)
        if (A.kind[i] != B.kind[i]) ++bad_kind;
      for (std::size_t S = 0; S < A.key.size(); ++S) {
        if (A.key[S] != B.key[S]) ++bad_key;
        if (!std::equal(A.canon_face_bits[S].begin(),
                        A.canon_face_bits[S].end(),
                        B.canon_face_bits[S].begin(),
                        B.canon_face_bits[S].end()))
          ++bad_canon;
        slots_seen += A.canon_face_bits[S].size();
      }
      subsets_seen += A.key.size() - 1;
      // coverage: this term really did take the threaded path (the sweep falls
      // back to a serial loop below canon_min_parallel_subsets = 32 subsets)
      if (A.key.size() - 1 >= 32) ++parallel_terms;
    }
    CHECK(bad_key == 0);
    CHECK(bad_canon == 0);
    CHECK(bad_kind == 0);
  };

  compare(0);
  compare(1);

  // ...and the sweep is not vacuous
  CHECK(subsets_seen > 0);
  CHECK(slots_seen > 0);
  CHECK(parallel_terms > 0);
}

// T-C4 evaluates one block's whole NNI neighbourhood concurrently and then
// picks the winner in a separate serial scan. The failure mode that a
// thread-count comparison alone can MISS is the interesting one: the rule is
// `if (c < best)`, strict, so among equal minimal costs the lowest move index
// wins, and a version that reduced inside the parallel region would agree with
// the serial answer on every run where no tie happened to occur -- and then
// disagree, nondeterministically, on the run that matters.
//
// So the oracle here is not "another thread count", it is the sequential loop
// itself, transcribed from the pre-T-C4 body of `hill_climb` (ga.cpp at
// 822e01be) and left un-refactored on purpose: splice candidate, evaluate,
// keep on strict `<`, otherwise roll back. It also counts how often a block's
// minimum was attained by more than one candidate AND was an improvement, i.e.
// how many times the first-wins rule actually decided something -- asserted
// nonzero below, so the tie coverage cannot silently lapse.
double reference_hill_climb(ga::Fitness const& F, ga::Genome& genome,
                            std::size_t max_sweeps, std::size_t& improvements,
                            std::size_t& tie_decisions) {
  auto const& kt = F.table();
  auto const& lay = F.layout();
  double best = F(genome);
  auto try_block = [&](container::svector<int>& code,
                       std::pair<int, int> slice, int n) {
    if (n < 3) return false;
    bool improved = false;
    ga::TreeCode block(code.begin() + slice.first, code.begin() + slice.second);
    const ga::Laminar fam = ga::decode_tree(block, n);
    const double entry = best;
    container::svector<double> costs;
    for (auto const& mv : ga::nni_moves(fam)) {
      auto cand = ga::encode_tree(mv, n);
      std::copy(cand.begin(), cand.end(), code.begin() + slice.first);
      const double c = F(genome);
      costs.push_back(c);
      if (c < best) {
        best = c;
        block = cand;
        improved = true;
      } else {
        std::copy(block.begin(), block.end(), code.begin() + slice.first);
      }
    }
    if (improved) ++improvements;
    // did a tie at the minimum have to be broken, on an improving block?
    if (!costs.empty()) {
      const double mn = *std::min_element(costs.begin(), costs.end());
      if (mn < entry &&
          std::count(costs.begin(), costs.end(), mn) > 1)
        ++tie_decisions;
    }
    return improved;
  };
  for (std::size_t sweep = 0; sweep < max_sweeps; ++sweep) {
    bool improved = false;
    for (std::size_t d = 0; d < kt.terms.size(); ++d)
      improved |= try_block(genome.g, lay.g_slice[d],
                            static_cast<int>(kt.terms[d].n()));
    for (std::size_t t = 0; t < kt.targets.size(); ++t)
      improved |= try_block(genome.h, lay.h_slice[t],
                            static_cast<int>(kt.targets[t].terms.size()));
    if (!improved) break;
  }
  return best;
}

TEST_CASE("ga hill_climb is the sequential argmin, first-wins", "[ga]") {
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg);
  auto const terms = load_md_terms(SEQUANT_UNIT_TESTS_SOURCE_DIR
                                   "/ga_pao_ccsd_df.md");

  // put back whatever the rest of the suite runs at, on the failure path too
  struct RestoreThreads {
    int n;
    ~RestoreThreads() { set_num_threads(n); }
  } const restore{num_threads()};

  // three universes rather than one: a tie at an improving block's minimum is
  // not common (1-3 per hill climb here), so the coverage is pooled
  std::size_t improvements = 0, tie_decisions = 0;
  for (int mode : {0, 1, 3}) {  // f-only, g-only, the whole CCSD residual set
    INFO("universe mode = " << mode);
    auto const targets = universe(terms, mode);
    auto const kt = ga::build_key_table(targets, python_extent);
    ga::Fitness F(kt);
    auto const seed = ga::seed_genome(kt);
    const double seed_cost = F(seed);
    const std::size_t sweeps = 40;  // the default; it converges well before

    ga::Genome ref = seed;
    const double ref_cost =
        reference_hill_climb(F, ref, sweeps, improvements, tie_decisions);

    for (int nthreads : {1, 2, 8}) {
      INFO("threads = " << nthreads);
      set_num_threads(nthreads);
      ga::Genome g = seed;
      const double c = ga::hill_climb(F, g, sweeps);
      // exact equality: this is one genome's cost from one evaluation, never a
      // reduction over threads
      CHECK(c == ref_cost);
      CHECK(g.g == ref.g);
      CHECK(g.h == ref.h);
      CHECK(F(g) == ref_cost);
    }

    // ...and none of the above is vacuous
    CHECK(ref_cost < seed_cost);
    CHECK(ref.g != seed.g);
  }
  CHECK(improvements > 0);
  CHECK(tie_decisions > 0);
}

// T-C3 made `ga_once` breed serially and evaluate in parallel. Two things can
// break that, and only one of them is a race: (a) the parallel phase is not
// actually independent per kid, or (b) phase 1 adds, drops or reorders a draw
// on the single mt19937_64 the whole search shares -- which changes every
// decision after it. Both surface the same way, as the search landing
// somewhere else, so this pins the OUTCOME (cost AND genome) across thread
// counts.
//
// The GATE-B fingerprint on the 81-term universe is the end-to-end oracle for
// the same property; this is the fast one, and it localizes a failure to
// `run_ga` rather than to setup.
TEST_CASE("ga run_ga does not depend on the thread count", "[ga]") {
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg);
  auto const terms = load_md_terms(SEQUANT_UNIT_TESTS_SOURCE_DIR
                                   "/ga_pao_ccsd_df.md");

  // put back whatever the rest of the suite runs at, on the failure path too
  struct RestoreThreads {
    int n;
    ~RestoreThreads() { set_num_threads(n); }
  } const restore{num_threads()};

  auto const targets = universe(terms, 1);  // g-only
  auto const kt = ga::build_key_table(targets, python_extent);
  ga::Fitness F(kt);
  auto const seed = ga::seed_genome(kt);
  const double seed_cost = F(seed);

  ga::GAOptions opts;
  opts.population = 16;
  opts.generations = 5;
  opts.restarts = 1;
  opts.hill_climb_sweeps = 2;

  auto run_at = [&](int nthreads) {
    set_num_threads(nthreads);
    return ga::run_ga(F, seed, opts);
  };

  auto const [c1, g1] = run_at(1);
  auto const [c2, g2] = run_at(2);
  auto const [c8, g8] = run_at(8);

  // the search must have gone somewhere, or every thread count agrees on the
  // seed for reasons that have nothing to do with T-C3
  CHECK(c1 < seed_cost);
  CHECK(g1.g != seed.g);

  // exact equality, not a tolerance: the winner is one genome's cost, taken
  // from one evaluation, not a reduction over the threads
  CHECK(c2 == c1);
  CHECK(c8 == c1);
  CHECK(g2.g == g1.g);
  CHECK(g2.h == g1.h);
  CHECK(g8.g == g1.g);
  CHECK(g8.h == g1.h);

  // ...and the parallel phase had something to spread over 8 threads
  CHECK(opts.reproduction * opts.population >= 8);
}

TEST_CASE("ga emission expands to the original terms", "[ga]") {
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg);

  auto const terms = load_md_terms(SEQUANT_UNIT_TESTS_SOURCE_DIR
                                   "/ga_pao_ccsd_df.md");

  auto require_expansion = [](container::svector<ga::TargetInput> const&
                                  targets,
                              ga::Genome const& genome) {
    auto kt = ga::build_key_table(targets, python_extent);
    ga::Fitness F(kt, ga::CostModel::python_parity(),
                  ga::ProducerResolution::Greedy);
    auto sch = F.explain(genome);
    auto emitted = ga::emit(F, sch);
    REQUIRE(emitted.size() == targets.size());
    for (std::size_t t = 0; t < targets.size(); ++t) {
      Sum original;
      for (std::size_t s = 0; s < targets[t].summands.size(); ++s) {
        // emission names externals after the target's first term
        REQUIRE(targets[t].ext[s] == targets[t].ext.front());
        original.append(targets[t].summands[s]->clone());
      }
      REQUIRE_THAT(emitted[t], EquivalentTo(ex<Sum>(std::move(original))));
      // the emitted trees binarize: SeQuant's evaluation machinery can
      // consume the unrolled DAG directly
      auto const& ext = targets[t].ext.front();
      ResultExpr res =
          ext.empty()
              ? ResultExpr(Variable(L"Z"), emitted[t])
              : ResultExpr(Tensor(L"Z", bra(container::svector<Index>(
                                             ext.begin(), ext.end())),
                                  ket{}),
                           emitted[t]);
      auto node = binarize(res);
      CHECK(!node.leaf());
    }
  };

  SECTION("f-only, control schedule") {
    auto const& ref = ga_reference::fonly;
    require_expansion(universe(terms, 0),
                      slice_genome(ref.g_ctrl, ref.h_ctrl));
  }
  SECTION("f-only, seed") {
    auto const& ref = ga_reference::fonly;
    require_expansion(universe(terms, 0),
                      slice_genome(ref.g_opt, ref.h0));
  }
  SECTION("g-only, prototype winner") {
    auto const& ref = ga_reference::gonly;
    require_expansion(universe(terms, 1),
                      slice_genome(ref.winner_g, ref.winner_h));
  }
}

// Replay weighting (CostModel::volatile_weight) rescales the flops of merges
// whose cluster contains an amplitude leaf. It is a reweighting of the SAME
// flop objective, not a second metric -- but because Fitness::resolve picks
// each key's producer by weighted cost, a different weight can select a
// different producer and therefore emit a DIFFERENT expression. What no weight
// may change is the mathematics: every emission, at every weight, has to expand
// back to the original sum. This also pins the two properties that make the
// knob safe to ship: it is inert without a volatile-leaf predicate (so every
// pre-existing number is untouched), and the cost is monotone in the weight.
//
// Monotone, NOT affine: each candidate producer's cost is nondecreasing in the
// weight, so the per-key minimum is too, but resolve() re-picks producers at
// every weight, making the total a minimum over affine pieces (concave
// piecewise-linear) rather than a single line.
TEST_CASE("ga replay weighting preserves the mathematics", "[ga]") {
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg);

  auto const terms = load_md_terms(SEQUANT_UNIT_TESTS_SOURCE_DIR
                                   "/ga_pao_ccsd_df.md");
  // g-only: carries amplitudes, and its terms have amplitude-FREE subclusters
  // (the g*C half-transforms), so the weight genuinely discriminates rather
  // than rescaling everything uniformly.
  auto const targets = universe(terms, 1);
  auto const& ref = ga_reference::gonly;
  auto const genome = slice_genome(ref.winner_g, ref.winner_h);

  std::function<bool(Tensor const&)> const amplitude = [](Tensor const& t) {
    return t.label() == L"t";
  };
  std::function<bool(Tensor const&)> const no_predicate{};

  // Evaluate the fixed genome under one weight; optionally hand back the
  // emitted (fully inlined) expressions. The KeyTable must outlive the
  // Fitness that borrows it, so everything stays inside the lambda.
  auto schedule_cost = [&](double weight,
                           std::function<bool(Tensor const&)> const& pred,
                           container::svector<ExprPtr>* emitted = nullptr) {
    auto kt = ga::build_key_table(targets, python_extent, pred);
    auto cm = ga::CostModel::python_parity();
    cm.volatile_weight = weight;
    ga::Fitness F(kt, cm, ga::ProducerResolution::Greedy);
    auto sch = F.explain(genome);
    if (emitted) *emitted = ga::emit(F, sch);
    return sch.total;
  };

  SECTION("inert without a volatile-leaf predicate") {
    // No predicate => volatile_mask is 0 => nothing is ever charged twice.
    // Exact equality is the claim: weighting by 1.0 is an exact no-op in IEEE.
    const double blind = schedule_cost(1.0, no_predicate);
    CHECK(schedule_cost(20.0, no_predicate) == blind);
    CHECK(schedule_cost(1.0, amplitude) == blind);
  }

  SECTION("monotone in the replay count, and actually engaged") {
    const double c1 = schedule_cost(1.0, amplitude);
    const double c2 = schedule_cost(2.0, amplitude);
    const double c20 = schedule_cost(20.0, amplitude);
    CHECK(c2 >= c1);
    CHECK(c20 >= c2);
    // there IS amplitude-dependent work here, so the knob must bite
    CHECK(c20 > c1);
  }

  SECTION("every weight emits the same mathematics") {
    for (double w : {1.0, 2.0, 20.0, 100.0}) {
      container::svector<ExprPtr> emitted;
      schedule_cost(w, amplitude, &emitted);
      REQUIRE(emitted.size() == targets.size());
      for (std::size_t t = 0; t < targets.size(); ++t) {
        Sum original;
        for (std::size_t s = 0; s < targets[t].summands.size(); ++s) {
          // emission names externals after the target's first term
          REQUIRE(targets[t].ext[s] == targets[t].ext.front());
          original.append(targets[t].summands[s]->clone());
        }
        INFO("volatile_weight = " << w << ", target " << t);
        REQUIRE_THAT(emitted[t], EquivalentTo(ex<Sum>(std::move(original))));
      }
    }
  }
}

// The named-emission runtime contract: every shared array is a named tensor
// whose definition is emitted once, in dependency order; each use leaf's
// slots pair position-for-position with the definition head's slots (same
// space, same proto rank); and -- what the mpqc yielder's fast path relies
// on -- the canonical leaf order (EvalExpr construction) is the same
// positional pattern for the definition head and every use tensor, so the
// stored definition value serves each use leaf without a permute.
TEST_CASE("ga named emission runtime contract", "[ga]") {
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg);

  auto const terms = load_md_terms(SEQUANT_UNIT_TESTS_SOURCE_DIR
                                   "/ga_pao_ccsd_df.md");

  auto verify = [](container::svector<ga::TargetInput> const& targets,
                   ga::Genome const& genome) {
    auto kt = ga::build_key_table(targets, python_extent);
    ga::Fitness F(kt, ga::CostModel::python_parity(),
                  ga::ProducerResolution::Greedy);
    auto sch = F.explain(genome);
    auto em = ga::emit_named(F, sch);
    REQUIRE(em.targets.size() == targets.size());

    container::map<std::wstring, std::size_t> def_pos;
    for (std::size_t i = 0; i < em.definitions.size(); ++i) {
      auto const& d = em.definitions[i];
      REQUIRE(ga::is_named_intermediate(d.label()));
      REQUIRE(def_pos.emplace(d.label(), i).second);  // unique names
    }

    container::map<std::wstring, std::size_t> n_uses;
    std::function<void(ExprPtr const&, std::size_t)> scan =
        [&](ExprPtr const& e, std::size_t host_pos) {
          if (e->is<Product>()) {
            for (auto const& f : e->as<Product>().factors()) scan(f, host_pos);
            return;
          }
          if (e->is<Sum>()) {
            for (auto const& s : e->as<Sum>().summands()) scan(s, host_pos);
            return;
          }
          if (!e->is<Tensor>()) return;
          auto const& t = e->as<Tensor>();
          if (!ga::is_named_intermediate(t.label())) return;
          ++n_uses[std::wstring(t.label())];
          auto it = def_pos.find(std::wstring(t.label()));
          REQUIRE(it != def_pos.end());
          // dependency order: a definition body references earlier ones only
          REQUIRE(it->second < host_pos);
          auto const& def = em.definitions[it->second];
          auto const& ds = def.aux();
          auto const& us = t.aux();
          REQUIRE(t.bra_rank() == 0);
          REQUIRE(t.ket_rank() == 0);
          REQUIRE(us.size() == ds.size());
          container::map<Index, Index> u2d;
          for (std::size_t k = 0; k < ds.size(); ++k) {
            CHECK(us[k].space() == ds[k].space());
            REQUIRE(us[k].proto_indices().size() ==
                    ds[k].proto_indices().size());
            u2d.emplace(us[k], ds[k]);
            for (std::size_t j = 0; j < us[k].proto_indices().size(); ++j)
              u2d.try_emplace(us[k].proto_indices()[j],
                              ds[k].proto_indices()[j]);
          }
          // canonical-leaf-order pattern equality: translate the use leaf's
          // canonical order into definition labels and compare with the
          // definition head's own canonical annotation
          EvalExpr ue{t};
          EvalExpr de{def.result_as_tensor()};
          Index::index_vector mapped;
          for (auto const& ix : ue.canon_indices()) {
            auto found = u2d.find(ix);
            REQUIRE(found != u2d.end());
            mapped.push_back(found->second);
          }
          CHECK(mapped == de.canon_indices());
        };
    for (auto const& tgt : em.targets) scan(tgt, em.definitions.size());
    for (std::size_t i = 0; i < em.definitions.size(); ++i)
      scan(em.definitions[i].expression(), i);

    // every definition is genuinely shared (>= 2 use sites)
    for (auto const& d : em.definitions) {
      auto it = n_uses.find(d.label());
      REQUIRE(it != n_uses.end());
      CHECK(it->second >= 2);
    }
    std::printf("  named emission: %zu definitions over %zu targets\n",
                em.definitions.size(), em.targets.size());
  };

  SECTION("f-only, control schedule") {
    auto const& ref = ga_reference::fonly;
    verify(universe(terms, 0), slice_genome(ref.g_ctrl, ref.h_ctrl));
  }
  SECTION("g-only, prototype winner") {
    auto const& ref = ga_reference::gonly;
    verify(universe(terms, 1), slice_genome(ref.winner_g, ref.winner_h));
  }
}

// The runtime gate: binarize the emitted forest and account flops the way
// evaluate() + CacheManager do. The unique contraction nodes must reproduce
// the schedule's shared L1 flops exactly, the sum layer charged per
// occurrence must reproduce L2, and the CacheManager must register exactly
// the multi-use subtrees.
TEST_CASE("ga emitted forest reproduces the schedule under the runtime cache",
          "[ga]") {
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg);

  auto const terms = load_md_terms(SEQUANT_UNIT_TESTS_SOURCE_DIR
                                   "/ga_pao_ccsd_df.md");

  auto verify = [](container::svector<ga::TargetInput> const& targets,
                   ga::Genome const& genome) {
    auto kt = ga::build_key_table(targets, python_extent);
    ga::Fitness F(kt, ga::CostModel::native(), ga::ProducerResolution::Greedy);
    auto sch = F.explain(genome);
    auto emitted = ga::emit(F, sch);

    container::svector<EvalNode<EvalExpr>> nodes;
    for (std::size_t t = 0; t < emitted.size(); ++t) {
      auto const& ext = targets[t].ext.front();
      ResultExpr res =
          ext.empty()
              ? ResultExpr(Variable(L"Z"), emitted[t])
              : ResultExpr(Tensor(L"Z", bra(container::svector<Index>(
                                             ext.begin(), ext.end())),
                                  ket{}),
                           emitted[t]);
      nodes.push_back(binarize(res));
    }

    RuntimeCoster rc{python_extent};
    for (auto const& n : nodes) rc.walk_cached(n);
    for (auto const& n : nodes) rc.walk_l2(n);

    // runtime unique-array accounting == the schedule's shared L1 flops:
    // producer substitution makes subtree identity coincide with array
    // identity, key for key
    CHECK(rc.l1_unique == sch.l1);
    // the sum layer, charged per occurrence, is the schedule's L2
    CHECK(rc.l2_every == sch.l2);
    // what evaluate() + cache would spend; below the schedule total exactly
    // by the structurally duplicate sums the cache shares anyway
    CHECK(rc.cache_total == sch.l1 + rc.l2_unique);
    CHECK(rc.cache_total <= sch.total);

    // the CacheManager registers exactly the multi-use subtrees, each with
    // max_life equal to its cache-semantics use count
    auto cm = cache_manager(nodes);
    std::size_t entries = 0;
    cm.for_each_key([&](RuntimeCoster::Node const& k) {
      ++entries;
      auto it = rc.uses.find(&k);
      REQUIRE(it != rc.uses.end());
      CHECK(cm.max_life(k) == static_cast<int>(it->second));
    });
    std::size_t multi = 0;
    for (auto const& [n, c] : rc.uses) multi += c >= 2;
    CHECK(entries == multi);

    std::printf("  schedule %.6e = l1 %.6e + l2 %.6e; cache-aware %.6e "
                "(dup-sum savings %.6e); %zu picked arrays, %zu cache "
                "entries\n",
                sch.total, sch.l1, sch.l2, rc.cache_total,
                sch.total - rc.cache_total, sch.pick.size(), entries);
  };

  SECTION("f-only, control schedule") {
    auto const& ref = ga_reference::fonly;
    verify(universe(terms, 0), slice_genome(ref.g_ctrl, ref.h_ctrl));
  }
  SECTION("g-only, prototype winner") {
    auto const& ref = ga_reference::gonly;
    verify(universe(terms, 1), slice_genome(ref.winner_g, ref.winner_h));
  }
}
