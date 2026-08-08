#include <SeQuant/core/optimize/ga/fitness.hpp>

#include <SeQuant/core/bliss.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/algorithm/find.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <map>
#include <mutex>
#include <numeric>
#include <unordered_map>

namespace sequant::opt::ga {

// `lex_less` (the child/fibre ordering) now lives in genome.hpp, next to the
// codec, because the decode memo in fitness.hpp needs it too.

GenomeLayout GenomeLayout::of(KeyTable const& kt) {
  GenomeLayout out;
  int p = 0;
  for (auto const& T : kt.terms) {
    out.g_slice.emplace_back(p, p + static_cast<int>(T.n()) - 1);
    for (int k = 2; k <= static_cast<int>(T.n()); ++k)
      out.g_ranges.push_back(2 * k - 3);
    p = out.g_slice.back().second;
  }
  p = 0;
  for (auto const& tgt : kt.targets) {
    out.h_slice.emplace_back(p, p + static_cast<int>(tgt.terms.size()) - 1);
    for (int k = 2; k <= static_cast<int>(tgt.terms.size()); ++k)
      out.h_ranges.push_back(2 * k - 3);
    p = out.h_slice.back().second;
  }
  return out;
}

NodeMask const* ForestState::Term::children_of(NodeMask S) const {
  for (auto const& e : ch)
    if (e[0] == S) return &e[1];
  return nullptr;
}

// The beta cache key, flattened into a word buffer (T-A6). It used to hold four
// `Cluster`s plus two `svector<pair<Index, AmbientTag>>` -- ~350 B per entry at
// `sizeof(Index) == 168`, compared tuple-wise by `Index::operator<` down a
// `std::map`, on EVERY call rather than every miss. Everything in it is now
// small integers:
//
//   [0..3]  x1.S, v1.S, x2.S, v2.S
//   [4..5]  the four term ids
//   [6]     the two ambient terms
//   [7]     the two ambient sizes
//   [8..]   a1's entries, then a2's, each `(bit << 56) | tag` (tags leave the
//           top byte free, see AmbientTag)
//
// The two ambient terms are redundant -- a1.d is always x1.d and a2.d always
// x2.d, since an ambient's keys are the face of the value it belongs to -- but
// including them costs one word and makes the key's fidelity a property of the
// buffer rather than of that argument.
//
// Fidelity vs. the old key: within a term, bit <-> `Index` is a bijection and
// the term is in the buffer, so no two distinct old keys collide here. (The
// converse is not quite true -- an `Index`-valued tag from two different terms
// used to compare equal and now does not -- which can only cost a recomputation
// of a pure function, never change an answer.)
struct Fitness::BetaKey {
  /// 24 words of inline storage: the fixed part is 8 and the two ambients are
  /// face-sized, so the buffer stays on the stack for the whole search on
  /// C4H10/DZ. (Measured: 32 is not faster, so nothing is spilling.)
  container::svector<std::uint64_t, 24> w;
  /// `mix(w)`, computed once by `seal()`. Precomputed rather than derived in
  /// the hasher because T-C2 needs the hash BEFORE the map is reached (it
  /// picks the shard), and hashing the buffer twice per lookup would be a
  /// measurable tax on the hottest cache in the evaluation.
  std::uint64_t h = 0;
  /// The old `BetaKeyHash` body, verbatim -- same hash VALUE, so keys bucket
  /// exactly as before.
  void seal() {
    std::uint64_t x2 = 0x9e3779b97f4a7c15ull;
    for (std::uint64_t x : w) {
      x2 ^= x;
      x2 *= 0x100000001b3ull;
      x2 ^= x2 >> 29;
    }
    x2 ^= x2 >> 32;
    x2 *= 0xc4ceb9fe1a85ec53ull;
    x2 ^= x2 >> 32;
    h = x2;
  }
  friend bool operator==(BetaKey const& l, BetaKey const& r) {
    return l.w == r.w;
  }
};

struct Fitness::BetaKeyHash {
  std::size_t operator()(BetaKey const& k) const {
    return static_cast<std::size_t>(k.h);
  }
};

// --- T-C2: the memoization caches are shared across evaluating threads ------
//
// These three (four, counting the emission-only `Index` form) maps are the
// ONLY shared-mutable state on the evaluation path -- `EvalScratch` is
// per-thread by construction (see its docs) and the `KeyTable` is read-only
// after `build_key_table`. So making them concurrent is the whole of what
// Group C needs before it may run evaluations in parallel.
//
// Each is a fixed 64-way shard array, one `std::mutex` per shard, and the
// protocol at every call site is the same three steps:
//
//   1. lock the key's shard, look up, unlock;
//   2. on a miss compute the value OUTSIDE any lock;
//   3. lock the shard again and `try_emplace` -- FIRST INSERTION WINS, the
//      loser drops its copy and uses the winner's.
//
// Two threads racing on the same key therefore both compute it, and that is
// deliberately benign rather than merely tolerated: every one of these values
// is a pure function of its key (see the per-cache notes below), so the two
// computations agree bit-for-bit and it cannot matter which one is stored.
// Holding the lock across the computation would instead serialize `face_perms`
// (bliss) and `find_beta` (`std::next_permutation` over the face), i.e. the
// exact work Group C exists to spread out.
//
// REFERENCE STABILITY IS PART OF THE CONTRACT, not an implementation detail:
// `correspondences`/`correspondences_bits`/`face_perms` return references into
// their map and `find_beta` returns a pointer into the beta cache (T-A6),
// which callers hold across further calls. Both `std::map` and
// `std::unordered_map` are node-based -- a rehash relinks nodes, it does not
// move them -- and NOTHING here ever erases, so a reference handed out under
// the shard lock stays valid for the life of the `Fitness`, whatever any other
// thread does to any shard afterwards. That is also why the shard count is
// fixed at construction: the shard array itself is never resized.
//
// Sharding cannot change WHICH value is memoized. It changes only where a
// (key, value) pair is stored and, on a race, which of two identical values is
// stored -- never how a value is computed from its key.
namespace {

template <typename Map>
class ShardedCache {
 public:
  using key_type = typename Map::key_type;
  using mapped_type = typename Map::mapped_type;

  /// 64 shards: enough that 8-12 evaluating threads collide rarely, small
  /// enough that the array (64 mutexes + 64 empty maps, one cache line each)
  /// is noise next to a `Fitness`.
  static constexpr std::size_t n_shards = 64;

  /// The stored value for \p k, or null. \p h is any hash of \p k; it must be
  /// the same one `put` is given for that key.
  mapped_type const* find(key_type const& k, std::uint64_t h) {
    Shard& s = shard(h);
    std::lock_guard<std::mutex> lock(s.m);
    auto it = s.map.find(k);
    return it == s.map.end() ? nullptr : &it->second;
  }

  /// Stores \p v under \p k if the key is absent, and returns the stored
  /// value -- which is the one that got there first.
  mapped_type const& put(key_type&& k, mapped_type&& v, std::uint64_t h) {
    Shard& s = shard(h);
    std::lock_guard<std::mutex> lock(s.m);
    return s.map.try_emplace(std::move(k), std::move(v)).first->second;
  }

 private:
  /// Cache-line aligned so two shards' mutexes never share a line.
  struct alignas(64) Shard {
    std::mutex m;
    Map map;
  };
  /// Fibonacci hashing on the TOP bits: the per-cache hashes below feed their
  /// low bits to the shard's own map (`unordered_map`'s bucket index, `map`'s
  /// comparisons are unaffected), so the shard must not consume the same bits.
  Shard& shard(std::uint64_t h) {
    return shards_[(h * 0x9e3779b97f4a7c15ull) >> 58];
  }
  std::array<Shard, n_shards> shards_;
};

/// Shard selector for the `Cluster`-keyed caches (`corr`, `corr_ix`, `auts`).
/// Only used to pick a shard -- the maps inside a shard are still ordered by
/// `Cluster`'s `<=>`, so nothing about their contents depends on this.
inline std::uint64_t cluster_hash(Cluster c) {
  std::uint64_t h = c.S * 0x9e3779b97f4a7c15ull;
  h ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(c.d)) *
       0xff51afd7ed558ccdull;
  h ^= h >> 29;
  return h * 0xc4ceb9fe1a85ec53ull;
}
inline std::uint64_t cluster_pair_hash(std::pair<Cluster, Cluster> const& k) {
  return cluster_hash(k.first) ^ (cluster_hash(k.second) * 3);
}

}  // namespace

struct Fitness::Caches {
  /// Pure memoization of `find_beta`, so the container is free to change; what
  /// is NOT free to change is the enumeration order inside `find_beta`, which
  /// decides which beta is memoized. That order is a function of the key
  /// alone (`TermTable::index_rank` groups the two faces, `next_permutation`
  /// walks them, `correspondences_bits(v1, v2)` -- itself keyed and pure --
  /// supplies Sigma3), so which thread computes an entry, and whether it is
  /// computed once or twice, is unobservable.
  ///
  /// `std::unordered_map` is node-based, so `find_beta` may hand out a pointer
  /// into it and its callers may hold that pointer across further calls; the
  /// sharded wrapper preserves that (see above).
  ShardedCache<
      std::unordered_map<BetaKey, std::optional<BitPairs>, BetaKeyHash>>
      beta;
  /// Also pure in its key: the canonical face orders come from the immutable
  /// key table and the permutations from `face_perms(c1)`.
  ShardedCache<
      std::map<std::pair<Cluster, Cluster>, container::svector<BitPairs>>>
      corr;
  /// The `Index` form of `corr`, built on demand for emission.
  ShardedCache<
      std::map<std::pair<Cluster, Cluster>,
               container::svector<container::svector<std::pair<Index, Index>>>>>
      corr_ix;
  /// bliss automorphisms of one cut, closed into a group -- deterministic for
  /// a given graph, and the graph is rebuilt from the immutable key table, so
  /// pure in `c` like the rest. This is the cache that must NOT go
  /// thread-local: it is the expensive one, and its whole value is that the
  /// bliss run amortizes ACROSS threads.
  ShardedCache<std::map<Cluster, container::svector<container::svector<int>>>>
      auts;
};

Fitness::Fitness(KeyTable const& kt, CostModel cost,
                 ProducerResolution resolution)
    : kt_(&kt),
      cost_(cost),
      resolution_(resolution),
      layout_(GenomeLayout::of(kt)),
      caches_(std::make_shared<Caches>()),
      scratch_(kt) {
  // Per-term leaf data: pure functions of the KeyTable, so hoisted out of the
  // ~164k-evaluation inner loop (T-A5). Body is verbatim the old
  // `leaf_summand`; only the moment it runs changed.
  leaf_ambient_.reserve(kt.terms.size());
  leaf_val_.reserve(kt.terms.size());
  container::svector<int> rank;
  for (std::size_t d = 0; d < kt.terms.size(); ++d) {
    TermTable const& T = kt.terms[d];
    // eta: external index -> (kind, rank-within-kind) slot of the target. The
    // rank is handed out in ASCENDING `Index` order of the externals, which in
    // bit form is ascending `index_rank` (T-A6).
    AmbientMap ambient;
    ambient.d = static_cast<int>(d);
    container::svector<std::uint8_t> ext;
    for (auto const& ix : T.ext)
      ext.push_back(static_cast<std::uint8_t>(T.bit_of(ix)));
    std::sort(ext.begin(), ext.end(), [&](std::uint8_t a, std::uint8_t b) {
      return T.index_rank[a] < T.index_rank[b];
    });
    rank.assign(kt.kind_ids.size(), 0);
    for (std::uint8_t b : ext) {
      const int k = T.kind[b];
      ambient.push(b, eta_tag(k * 4096 + rank[static_cast<std::size_t>(k)]++));
    }
    ambient.seal();
    leaf_ambient_.push_back(std::move(ambient));
    leaf_val_.push_back(std::make_shared<Val const>(
        Val{Val::Cl, static_cast<int>(d), T.full(), 0, {}, {}, {}}));
  }
}

ForestState Fitness::decode_forest(Genome const& genome,
                                   EvalScratch& sc) const {
  ForestState st;
  for (std::size_t d = 0; d < kt_->terms.size(); ++d) {
    // the slice is [lo, lo + n - 1); the memo derives the length from n
    const int lo = layout_.g_slice[d].first;
    ForestState::Term term;
    // memo hit or fresh decode_tree + build_children -- identical either way
    sc.trees.decode(genome.g.data() + lo, static_cast<int>(kt_->terms[d].n()),
                    term.fam, term.ch);
    st.terms.push_back(std::move(term));
  }
  // Producer fibres: internal clusters grouped by key, in (term d ascending,
  // Term::ch order) within a key. Two passes over the SAME cluster sequence --
  // count then place -- reproduce that push order exactly (see Fibres).
  sc.fib.begin();
  for (std::size_t d = 0; d < st.terms.size(); ++d)
    for (auto const& e : st.terms[d].ch) sc.fib.count(kt_->terms[d].key[e[0]]);
  sc.fib.seal();
  for (std::size_t d = 0; d < st.terms.size(); ++d)
    for (auto const& e : st.terms[d].ch)
      sc.fib.place(kt_->terms[d].key[e[0]], {static_cast<int>(d), e[0]});
  sc.fib.close();
  return st;
}

Summand Fitness::leaf_summand(int d) const {
  Summand s;
  s.coeff = kt_->terms[d].scalar;
  s.val = leaf_val_[d];         // shared, immutable (see Fitness::leaf_val_)
  s.ambient = leaf_ambient_[d];  // precomputed in the ctor
  return s;
}

// The ways to read a value as X * V (proto_csv_opt.options): a cluster splits
// through either of its children; an already-factored value re-presents its
// stored split with the combined residual as the X part.
namespace {
struct Option {
  std::size_t key;
  Cluster x, v;
  Summand r;  // residual value (coeff/ambient filled by build_node)
};
/// The cost path's `Option`: the same (key, x, v) decision, but the residual
/// is one index into the per-evaluation `CostVal` arena.
struct CostOption {
  std::size_t key;
  Cluster x, v;
  int r;
};
}  // namespace

Summand Fitness::build_node(ForestState const& st, Summand s1,
                            Summand s2) const {
  auto options = [&](Summand const& s) {
    container::svector<Option> out;
    Val const& v = *s.val;
    if (v.kind == Val::Cl && std::popcount(v.S) >= 2) {
      auto const* cc = st.terms[v.d].children_of(v.S);
      for (int i = 0; i < 2; ++i) {
        const NodeMask keep = cc[i], other = cc[1 - i];
        out.push_back({kt_->terms[v.d].key[keep],
                       {v.d, other},
                       {v.d, keep},
                       Summand{1,
                               std::make_shared<Val>(Val{Val::Cl, v.d, other,
                                                         0, {}, {}, {}}),
                               {}}});
      }
    } else if (v.kind == Val::Fx) {
      out.push_back({kt_->terms[v.d].key[v.V],
                     {v.d, v.S},
                     {v.d, v.V},
                     Summand{1, v.inner, {}}});
    }
    return out;
  };

  for (auto const& o1 : options(s1))
    for (auto const& o2 : options(s2)) {
      if (o1.key != o2.key) continue;
      auto const* beta =
          find_beta(o1.x, o1.v, s1.ambient, o2.x, o2.v, s2.ambient);
      if (!beta) continue;
      // chain deeper: the residuals' ambient is X1's own face, identified
      // through beta on the second side
      Summand r1 = o1.r, r2 = o2.r;
      r1.coeff = s1.coeff;
      r2.coeff = s2.coeff;
      r1.ambient.d = o1.x.d;
      r2.ambient.d = o2.x.d;
      for (auto const& [x, y] : *beta) {
        r1.ambient.push(x, face_tag(o1.x.d, x));
        r2.ambient.push(y, face_tag(o1.x.d, x));
      }
      r1.ambient.seal();
      r2.ambient.seal();
      Summand inner = build_node(st, std::move(r1), std::move(r2));
      Summand out;
      out.coeff = 1;
      out.val = std::make_shared<Val>(
          Val{Val::Fx, o1.x.d, o1.x.S, o1.v.S, inner.val, {}, {}});
      out.ambient = s1.ambient;
      return out;
    }

  // no extraction possible: a plain addition on the face of s1
  Val const& v1 = *s1.val;
  const NodeMask face = v1.kind == Val::Fx ? (v1.S | v1.V) : v1.S;
  Summand out;
  out.coeff = 1;
  out.ambient = s1.ambient;
  out.val = std::make_shared<Val>(
      Val{Val::Sum, v1.d, face, 0, {}, std::move(s1), std::move(s2)});
  return out;
}

// The cost path's twin of `build_node` above (T-A5). READ THEM SIDE BY SIDE:
// they are the same algorithm, and the whole safety argument of having two
// evaluation paths rests on their staying that way.
//
//   * `options` enumerates the same splits in the same order, so the (o1, o2)
//     pairs are visited in the same order;
//   * `find_beta` is therefore called with the same arguments in the same
//     sequence -- which is what matters, because the beta cache turns the
//     DECISION SEQUENCE into the shared state of the two paths;
//   * the resulting value tree is node-for-node the same shape, so
//     `cost_of_cost_val` accumulates exactly what `cost_of_value` does, in the
//     same floating-point order.
//
// What is dropped: the `make_shared<Val>` per node, the `Summand::coeff`
// (`Complex<cpp_rational>`), the `AmbientMap` copies (borrowed here -- every
// map is written once and then only read), the retained `beta` and the
// retained Sum operands. None of them is read by the cost.
Fitness::CostSummand Fitness::build_node_cost(ForestState const& st,
                                              EvalScratch& sc, CostSummand s1,
                                              CostSummand s2) const {
  auto options = [&](CostSummand const& s,
                     container::svector<CostOption>& out) {
    out.clear();
    // by value: `new_val` below can reallocate the arena
    CostVal const v = sc.vals[static_cast<std::size_t>(s.val)];
    if (v.kind == Val::Cl && std::popcount(v.S) >= 2) {
      auto const* cc = st.terms[v.d].children_of(v.S);
      for (int i = 0; i < 2; ++i) {
        const NodeMask keep = cc[i], other = cc[1 - i];
        out.push_back({kt_->terms[v.d].key[keep],
                       {v.d, other},
                       {v.d, keep},
                       sc.new_val(Val::Cl, v.d, other, 0)});
      }
    } else if (v.kind == Val::Fx) {
      out.push_back(
          {kt_->terms[v.d].key[v.V], {v.d, v.S}, {v.d, v.V}, v.a});
    }
  };

  container::svector<CostOption> o1s, o2s;
  // `build_node` re-evaluates `options(s2)` once per outer iteration; hoisting
  // it here is behavior-neutral because `options` is a pure function of
  // (st, s.val) up to which arena slot the fresh Cl residual lands in, and at
  // most one (o1, o2) pair is ever consumed.
  options(s1, o1s);
  options(s2, o2s);
  for (auto const& o1 : o1s)
    for (auto const& o2 : o2s) {
      if (o1.key != o2.key) continue;
      auto const* beta = find_beta(o1.x, o1.v, *s1.amb, o2.x, o2.v, *s2.amb);
      if (!beta) continue;
      // chain deeper: the residuals' ambient is X1's own face, identified
      // through beta on the second side
      AmbientMap& a1 = sc.new_ambient();
      AmbientMap& a2 = sc.new_ambient();
      a1.d = o1.x.d;
      a2.d = o2.x.d;
      for (auto const& [x, y] : *beta) {
        a1.push(x, face_tag(o1.x.d, x));
        a2.push(y, face_tag(o1.x.d, x));
      }
      a1.seal();
      a2.seal();
      const CostSummand inner = build_node_cost(st, sc, CostSummand{o1.r, &a1},
                                                CostSummand{o2.r, &a2});
      return CostSummand{
          sc.new_val(Val::Fx, o1.x.d, o1.x.S, o1.v.S, inner.val), s1.amb};
    }

  // no extraction possible: a plain addition on the face of s1
  CostVal const v1 = sc.vals[static_cast<std::size_t>(s1.val)];
  const NodeMask face = v1.kind == Val::Fx ? (v1.S | v1.V) : v1.S;
  return CostSummand{sc.new_val(Val::Sum, v1.d, face, 0, s1.val, s2.val),
                     s1.amb};
}

// (Sigma1)(Sigma2)(Sigma3) of the prototype's Definition 5.4d: a bijection
// beta : F(x1) -> F(x2), kind-respecting, mapping proto-index structure onto
// proto-index structure (S1), agreeing with the ambient identification on the
// axes that survive into the parent face (S2), and reading the shared factor
// through a valid axis correspondence of its cut (S3).
//
// T-A6: every question below is a bit question about the two terms, and is now
// asked that way. F(S) is `T.face_mask(S)` (since T-A7 that is the only form
// the key table stores it in), `kt_->kind_of(ix)` is
// `T.kind[bit]`, `FP1.find(x) != end` is a bit test, and the proto relation is
// a comparison of two `proto_mask`s. What did NOT change is the enumeration
// order, which is the one thing here that is observable in the emitted
// schedule: see the `index_rank` note below.
BitPairs const* Fitness::find_beta(Cluster x1, Cluster v1, AmbientMap const& a1,
                                   Cluster x2, Cluster v2,
                                   AmbientMap const& a2) const {
  BetaKey key;
  auto& w = key.w;
  w.push_back(x1.S);
  w.push_back(v1.S);
  w.push_back(x2.S);
  w.push_back(v2.S);
  auto pack2 = [](int lo32_a, int lo32_b) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(lo32_a))
            << 32) |
           static_cast<std::uint32_t>(lo32_b);
  };
  w.push_back(pack2(x1.d, v1.d));
  w.push_back(pack2(x2.d, v2.d));
  w.push_back(pack2(a1.d, a2.d));
  w.push_back(pack2(static_cast<int>(a1.e.size()),
                    static_cast<int>(a2.e.size())));
  for (auto const& [b, t] : a1.e)
    w.push_back((static_cast<std::uint64_t>(b) << 56) | t);
  for (auto const& [b, t] : a2.e)
    w.push_back((static_cast<std::uint64_t>(b) << 56) | t);

  key.seal();

  // T-C2: look the memo up, and on a miss compute the beta with no lock held
  // and publish it. `out` below is the LOCAL answer; the pointer this function
  // returns always points into the cache, which is where the first writer's
  // answer lives -- and every writer's answer is the same one, since the whole
  // body below is a pure function of the key.
  const std::uint64_t kh = key.h;
  if (auto const* hit = caches_->beta.find(key, kh))
    return *hit ? &**hit : nullptr;
  std::optional<BitPairs> out;
  // Publishes `out` under the shard lock and hands back a pointer into the
  // cache node -- never into `out`, which dies with this call. Note that the
  // "no beta exists" answer is memoized too (as `nullopt`), exactly as the
  // serial version's default-constructed slot did.
  auto publish = [&]() -> BitPairs const* {
    auto const& stored = caches_->beta.put(std::move(key), std::move(out), kh);
    return stored ? &*stored : nullptr;
  };

  TermTable const& T1 = kt_->terms[x1.d];
  TermTable const& T2 = kt_->terms[x2.d];
  const IndexMask F1 = T1.face_mask(x1.S);
  const IndexMask F2 = T2.face_mask(x2.S);

  // Group both faces by kind, ASCENDING `index_rank` within a kind and
  // ascending kind between them. This is the load-bearing line of the whole
  // task: it is what `std::sort` on the `Index` objects and a
  // `container::map<int, ...>` gave, hence the order `std::next_permutation`
  // walks the candidate bijections, hence WHICH valid beta is found first,
  // hence how emission renames the second operand of every extraction.
  // Ordering by bit instead would leave every flop count identical and change
  // the emitted schedule.
  container::svector<int> kinds;
  container::svector<container::svector<std::uint8_t>> by1, by2;
  auto bits_by_rank = [](TermTable const& T, IndexMask m,
                         container::svector<std::uint8_t>& out_bits) {
    out_bits.clear();
    for (; m; m &= m - 1)
      out_bits.push_back(static_cast<std::uint8_t>(std::countr_zero(m)));
    std::sort(out_bits.begin(), out_bits.end(),
              [&](std::uint8_t a, std::uint8_t b) {
                return T.index_rank[a] < T.index_rank[b];
              });
  };
  for (std::size_t k = 0; k < T1.kind_mask.size(); ++k) {
    const IndexMask m1 = F1 & T1.kind_mask[k];
    const IndexMask m2 = F2 & T2.kind_mask[k];
    // covers the old `F1set.size() != F2set.size()`, `by1.size() !=
    // by2.size()` and per-kind size checks in one pass: the kinds partition
    // both faces, so equal per-kind counts is exactly equal kind sets with
    // equal blocks (and hence equal totals).
    if (std::popcount(m1) != std::popcount(m2)) return publish();
    if (!m1) continue;
    kinds.push_back(static_cast<int>(k));
    by1.emplace_back();
    by2.emplace_back();
    bits_by_rank(T1, m1, by1.back());
    bits_by_rank(T2, m2, by2.back());
  }

  const IndexMask FP1 = T1.face_mask(x1.S | v1.S);
  const IndexMask FP2 = T2.face_mask(x2.S | v2.S);
  const IndexMask FV1 = T1.face_mask(v1.S);
  const IndexMask FV2 = T2.face_mask(v2.S);
  auto const& sigmas = correspondences_bits(v1, v2);

  // enumerate kind-respecting bijections: per-kind permutations of by2
  container::svector<container::svector<int>> perms(kinds.size());
  // beta as a direct map over term-1 bits; `absent` marks a bit outside F1
  static constexpr std::uint8_t absent = 0xff;
  std::uint8_t bmap[64];
  // lambda-with-self, not `std::function`: this is called once per (x1, v1,
  // x2, v2, a1, a2) cache miss and the type erasure was pure overhead.
  auto enumerate = [&](auto&& self, std::size_t ki) -> bool {
    if (ki == kinds.size()) {
      std::memset(bmap, absent, sizeof bmap);
      for (std::size_t i = 0; i < kinds.size(); ++i)
        for (std::size_t j = 0; j < by1[i].size(); ++j)
          bmap[by1[i][j]] = by2[i][static_cast<std::size_t>(perms[i][j])];
      // (Sigma1) proto-index relation. Mapping x's protos through beta and
      // comparing the resulting MASK against `proto_mask[beta(x)]` is the old
      // "sort both proto lists and compare": beta is injective and an index's
      // protos are distinct, so the two sets have the sizes the two vectors
      // had.
      for (IndexMask m = F1; m; m &= m - 1) {
        const int x = std::countr_zero(m);
        const IndexMask px = T1.proto_mask[static_cast<std::size_t>(x)];
        if (!px) continue;  // == !x.has_proto_indices()
        IndexMask mapped = 0;
        for (IndexMask q = px; q; q &= q - 1) {
          const std::uint8_t t = bmap[std::countr_zero(q)];
          if (t == absent) return false;  // proto not in face
          mapped |= IndexMask{1} << t;
        }
        if (mapped != T2.proto_mask[bmap[x]]) return false;
      }
      // (Sigma2) surviving axes agree with the ambient identification
      for (IndexMask m = F1; m; m &= m - 1) {
        const int x = std::countr_zero(m);
        const std::uint8_t y = bmap[x];
        const bool sp1 = (FP1 >> x) & 1;
        const bool sp2 = (FP2 >> y) & 1;
        if (sp1 != sp2) return false;
        if (sp1) {  // == `tag_of(a1, x) != tag_of(a2, y)` over optionals
          AmbientTag const* t1 = a1.find(static_cast<std::uint8_t>(x));
          AmbientTag const* t2 = a2.find(y);
          if ((t1 == nullptr) != (t2 == nullptr)) return false;
          if (t1 && *t1 != *t2) return false;
        }
      }
      // (Sigma3) existential over the shared factor's axis correspondences
      for (auto const& sig : sigmas) {
        bool ok = true;
        for (IndexMask m = F1; m; m &= m - 1) {
          const int x = std::countr_zero(m);
          const std::uint8_t y = bmap[x];
          const bool c1 = (FV1 >> x) & 1;
          const bool c2 = (FV2 >> y) & 1;
          if (c1 != c2) {
            ok = false;
            break;
          }
          if (c1) {
            auto it = std::find_if(sig.begin(), sig.end(), [&](auto const& p) {
              return p.first == x;
            });
            if (it == sig.end() || it->second != y) {
              ok = false;
              break;
            }
          }
        }
        if (ok) return true;
      }
      return false;
    }
    auto& perm = perms[ki];
    perm.resize(by1[ki].size());
    std::iota(perm.begin(), perm.end(), 0);
    do {
      if (self(self, ki + 1)) return true;
    } while (std::next_permutation(perm.begin(), perm.end()));
    return false;
  };

  // Sigma1-3 are checked inside `enumerate`, which stops at the first valid
  // bijection; reconstruct it from the perms it left behind.
  if (enumerate(enumerate, 0)) {
    BitPairs beta;
    for (std::size_t i = 0; i < kinds.size(); ++i)
      for (std::size_t j = 0; j < by1[i].size(); ++j)
        beta.emplace_back(by1[i][j],
                          by2[i][static_cast<std::size_t>(perms[i][j])]);
    out = std::move(beta);
  }
  return publish();
}

void Fitness::cost_of_value(ValPtr const& val, double& l2,
                            container::svector<Cluster>& demanded) const {
  Val const& v = *val;
  switch (v.kind) {
    case Val::Cl:
      if (std::popcount(v.S) >= 2) demanded.push_back({v.d, v.S});
      break;
    case Val::Fx:
      cost_of_value(v.inner, l2, demanded);
      if (std::popcount(v.V) >= 2) demanded.push_back({v.d, v.V});
      // L2 structure is never named: emission keeps Fx products and Sums inside
      // the TARGET expression (only keyed L1 clusters become IGA<n>
      // definitions), and every target tree is re-evaluated from scratch each
      // replay. So L2 work is replayed regardless of whether it happens to be
      // amplitude-free.
      l2 += cost_.merge(kt_->terms[v.d], v.S, v.V, kt_->volatility_aware);
      break;
    case Val::Sum:
      cost_of_value(v.s1.val, l2, demanded);
      cost_of_value(v.s2.val, l2, demanded);
      l2 += cost_.addition(kt_->terms[v.d], v.S, kt_->volatility_aware);
      break;
  }
}

// Twin of `cost_of_value` over the `CostVal` arena -- same switch, same order
// of the two recursive calls, same `l2 +=` sequence, so the floating-point
// summation order is the same. Edit the two together.
void Fitness::cost_of_cost_val(EvalScratch const& sc, int val, double& l2,
                               container::svector<Cluster>& demanded) const {
  CostVal const& v = sc.vals[static_cast<std::size_t>(val)];
  switch (v.kind) {
    case Val::Cl:
      if (std::popcount(v.S) >= 2) demanded.push_back({v.d, v.S});
      break;
    case Val::Fx:
      cost_of_cost_val(sc, v.a, l2, demanded);
      if (std::popcount(v.V) >= 2) demanded.push_back({v.d, v.V});
      l2 += cost_.merge(kt_->terms[v.d], v.S, v.V, kt_->volatility_aware);
      break;
    case Val::Sum:
      cost_of_cost_val(sc, v.a, l2, demanded);
      cost_of_cost_val(sc, v.b, l2, demanded);
      l2 += cost_.addition(kt_->terms[v.d], v.S, kt_->volatility_aware);
      break;
  }
}

// The working sets below are dense arrays over `[0, kt_->n_keys)` living in
// the scratch (T-A4); they used to be `container::map`/`container::set`, i.e.
// one O(k) memmove per insert over the 240-580 keys a genome touches. Every
// walk order that the ordered containers used to supply for free is
// reproduced explicitly, because those orders ARE the floating-point
// summation order of this function:
//   * the fibre loop takes `sc.fib.keys()`, sorted ascending == flat_map order;
//   * the two stack walks are seeded from `seeds` (ascending), pop from the
//     back and are gated by first insertion -- unchanged, only the gate's
//     storage differs, and `KeyStamps::insert` returns exactly what
//     `flat_set::insert(k).second` returned;
//   * `pick_out` iterates the pass-2 keys ascending, so they are collected in
//     visit order and sorted before the copy.
double Fitness::resolve(ForestState const& st, EvalScratch& sc,
                        container::svector<Cluster> const& demanded,
                        container::map<std::size_t, Cluster>* pick_out) const {
  auto& seeds = sc.seeds;
  seeds.clear();
  for (auto const& c : demanded)
    if (std::popcount(c.S) >= 2) seeds.push_back(kt_->terms[c.d].key[c.S]);
  std::sort(seeds.begin(), seeds.end());
  seeds.erase(std::unique(seeds.begin(), seeds.end()), seeds.end());

  auto merge_of = [&](Cluster c) {
    auto const* cc = st.terms[c.d].children_of(c.S);
    return cost_.merge(kt_->terms[c.d], cc[0], cc[1]);
  };
  // An amplitude-free array is only paid ONCE if something actually holds it
  // between replays, and emission (emit_named) names -- hence stores -- exactly
  // the keys used at more than one site; a single-use key is inlined into its
  // consumer and rebuilt whenever that consumer is. So the replay count of a
  // key is `1` only when it is BOTH amplitude-free AND shared. Weighting every
  // amplitude-free merge by 1 instead promises reuse the runtime never
  // provides, and the search answers by manufacturing single-use amplitude-free
  // work that is then inlined and recomputed every iteration (measured: R2's
  // inlined amplitude-free flops went 7.1e8 -> 7.2e11 at volatile_weight 20).
  // This mirrors the eval cache's own rule, `count >= min_repeats ||
  // persistent`: reuse needs a holder, not just amplitude independence.
  auto& uses = sc.uses;
  auto& uses_set = sc.uses_set;
  auto& walked = sc.walked;
  auto& seen = sc.seen;
  auto& seen_list = sc.seen_list;
  auto& pick = sc.pick;
  auto& pick_set = sc.pick_set;
  auto& stack = sc.stack;
  // `uses[k]` is only defined where `uses_set` says so -- the flat_map's
  // default-constructed 0 becomes an explicit first-touch store.
  auto bump = [&](std::size_t k) {
    if (uses_set.insert(k))
      uses[k] = 1;
    else
      ++uses[k];
  };
  auto producer_of = [&](std::size_t k) {
    assert(pick_set.contains(k));  // every walked key has a fibre, hence a pick
    return pick[k];
  };
  // Leaves the pass-2 keys in `seen_list`, in visit order (the caller sorts).
  auto dag_cost = [&] {
    // pass 1: how many sites consume each key -- L2 demands (with
    // multiplicity, so `demanded` rather than the deduplicated `seeds`) plus
    // one per referencing parent slot.
    uses_set.clear();
    for (auto const& c : demanded)
      if (std::popcount(c.S) >= 2) bump(kt_->terms[c.d].key[c.S]);
    {
      stack.assign(seeds.begin(), seeds.end());
      walked.clear();
      while (!stack.empty()) {
        const std::size_t k = stack.back();
        stack.pop_back();
        if (!walked.insert(k)) continue;
        Cluster const c = producer_of(k);
        auto const* cc = st.terms[c.d].children_of(c.S);
        for (int i = 0; i < 2; ++i)
          if (std::popcount(cc[i]) >= 2) {
            const std::size_t ck = kt_->terms[c.d].key[cc[i]];
            bump(ck);
            stack.push_back(ck);
          }
      }
    }
    // pass 2: cost, each key charged once per replay it actually incurs
    stack.assign(seeds.begin(), seeds.end());
    seen.clear();
    seen_list.clear();
    double cost = 0;
    while (!stack.empty()) {
      const std::size_t k = stack.back();
      stack.pop_back();
      if (!seen.insert(k)) continue;
      seen_list.push_back(static_cast<std::uint32_t>(k));
      Cluster const c = producer_of(k);
      auto const* cc = st.terms[c.d].children_of(c.S);
      auto const& T = kt_->terms[c.d];
      const bool shared = (uses_set.contains(k) ? uses[k] : 0) >= 2;
      // Without volatility tracking there is nothing to replay-weight: the
      // objective is the historical single-shot DAG count, every distinct key
      // charged exactly once.
      const bool replayed = kt_->volatility_aware &&
                            (CostModel::is_volatile(T, c.S) || !shared);
      cost += cost_.merge(T, cc[0], cc[1], replayed);
      for (int i = 0; i < 2; ++i)
        if (std::popcount(cc[i]) >= 2) stack.push_back(T.key[cc[i]]);
    }
    return cost;
  };

  // `sc.fib.keys()` is ascending, so this is the old `for (auto const& [k, fl]
  // : st.fib)`; the inner scan still re-scores `fl.front()` first, so the
  // first-wins tie-break over the fibre's push order is unchanged.
  auto const& fkeys = sc.fib.keys();
  pick_set.clear();
  for (std::uint32_t k : fkeys) {
    Cluster const* fl = sc.fib.data(k);
    const std::size_t nf = sc.fib.size(k);
    Cluster best = fl[0];
    double bc = merge_of(best);
    for (std::size_t i = 0; i < nf; ++i) {
      const double cc = merge_of(fl[i]);
      if (cc < bc) {
        bc = cc;
        best = fl[i];
      }
    }
    pick_set.insert(k);
    pick[k] = best;
  }

  auto copy_picks = [&] {
    std::sort(seen_list.begin(), seen_list.end());  // == flat_set order
    for (std::uint32_t k : seen_list) (*pick_out)[k] = pick[k];
  };

  if (resolution_ == ProducerResolution::Exact) {
    // enumerate one-producer-per-key assignments over multi-candidate fibres
    container::svector<std::size_t> multi;
    double space = 1;
    for (std::uint32_t k : fkeys)
      if (sc.fib.size(k) > 1) {
        multi.push_back(k);
        space *= static_cast<double>(sc.fib.size(k));
      }
    if (space <= 1e6) {
      double best_cost = std::numeric_limits<double>::max();
      // The candidate assignment is applied in place: every `multi` key is
      // overwritten on every iteration, and no other key is ever touched, so
      // `pick` is exactly the old `cand = pick; cand[multi[i]] = ...`. The
      // winner is remembered as its odometer reading and replayed at the end.
      container::svector<std::size_t> idx(multi.size(), 0),
          best_idx(multi.size(), 0);
      auto apply = [&](container::svector<std::size_t> const& ix) {
        for (std::size_t i = 0; i < multi.size(); ++i)
          pick[multi[i]] = sc.fib.data(multi[i])[ix[i]];
      };
      while (true) {
        apply(idx);
        const double c = dag_cost();
        if (c < best_cost) {
          best_cost = c;
          best_idx = idx;
        }
        std::size_t i = 0;
        for (; i < multi.size(); ++i) {
          if (++idx[i] < sc.fib.size(multi[i])) break;
          idx[i] = 0;
        }
        if (i == multi.size()) break;
      }
      if (pick_out) {
        apply(best_idx);
        dag_cost();
        copy_picks();
      }
      return best_cost;
    }
  }

  const double c = dag_cost();
  if (pick_out) copy_picks();
  return c;
}

double Fitness::operator()(Genome const& genome) const {
  return (*this)(genome, scratch_);
}

// The cost-only path (T-A5). Step for step `explain` below, minus everything
// only emission reads: no `Schedule`, no `Summand`/`Val`/`shared_ptr`, no
// coefficient, no retained roots, and `resolve` is asked for the number
// without the `pick` map (whose flat_map copy alone is ~273 O(k) inserts).
double Fitness::operator()(Genome const& genome, EvalScratch& sc) const {
  ForestState const forest = decode_forest(genome, sc);
  sc.begin_l2();
  sc.l2_roots.clear();
  for (std::size_t t = 0; t < kt_->targets.size(); ++t) {
    auto const& tgt = kt_->targets[t];
    const int lo = layout_.h_slice[t].first;
    sc.trees.decode(genome.h.data() + lo, static_cast<int>(tgt.terms.size()),
                    sc.l2_fam);
    auto rec = [&](auto&& self, NodeMask S) -> CostSummand {
      if (std::popcount(S) == 1) {
        const int d = static_cast<int>(tgt.terms[std::countr_zero(S)]);
        return CostSummand{sc.new_val(Val::Cl, d, kt_->terms[d].full(), 0),
                           &leaf_ambient_[d]};
      }
      auto [c1, c2] = tree_children(sc.l2_fam, S);
      if (lex_less(c2, c1)) std::swap(c1, c2);
      // sequenced explicitly, and identically in `explain`: the two paths
      // must present the same `find_beta` sequence, and function-argument
      // evaluation order is unspecified in C++.
      const CostSummand a = self(self, c1);
      const CostSummand b = self(self, c2);
      return build_node_cost(forest, sc, a, b);
    };
    sc.l2_roots.push_back(
        rec(rec, tgt.terms.empty() ? 0 : (NodeMask{1} << tgt.terms.size()) - 1)
            .val);
  }
  double l2 = 0;
  auto& demanded = sc.l2_demanded;
  demanded.clear();
  for (int root : sc.l2_roots) cost_of_cost_val(sc, root, l2, demanded);
  const double l1 = resolve(forest, sc, demanded, nullptr);
  return l1 + l2;
}

Schedule Fitness::explain(Genome const& genome) const {
  return explain(genome, scratch_);
}

Schedule Fitness::explain(Genome const& genome, EvalScratch& sc) const {
  Schedule sch;
  sch.forest = decode_forest(genome, sc);
  Laminar fam;
  for (std::size_t t = 0; t < kt_->targets.size(); ++t) {
    auto const& tgt = kt_->targets[t];
    const int lo = layout_.h_slice[t].first;
    sc.trees.decode(genome.h.data() + lo,
                    static_cast<int>(tgt.terms.size()), fam);
    auto rec = [&](auto&& self, NodeMask S) -> Summand {
      if (std::popcount(S) == 1)
        return leaf_summand(
            static_cast<int>(tgt.terms[std::countr_zero(S)]));
      auto [c1, c2] = tree_children(fam, S);
      if (lex_less(c2, c1)) std::swap(c1, c2);
      // see the identical comment in operator(): both paths sequence the two
      // child recursions left to right, explicitly.
      Summand a = self(self, c1);
      Summand b = self(self, c2);
      return build_node(sch.forest, std::move(a), std::move(b));
    };
    sch.roots.push_back(rec(
        rec, tgt.terms.empty() ? 0 : (NodeMask{1} << tgt.terms.size()) - 1));
  }
  container::svector<Cluster> demanded;
  for (auto const& root : sch.roots) cost_of_value(root.val, sch.l2, demanded);
  sch.l1 = resolve(sch.forest, sc, demanded, &sch.pick);
  sch.total = sch.l1 + sch.l2;
  return sch;
}

// Automorphism group of the cut of cluster c, as permutations of the
// positions of its canonical face order; the identity is always present.
// Mirrors the prototype's tie-order enumeration (Universe.axes), capped at
// the same 48 elements.
container::svector<container::svector<int>> const& Fitness::face_perms(
    Cluster c) const {
  const std::uint64_t h = cluster_hash(c);
  if (auto const* hit = caches_->auts.find(c, h)) return *hit;
  // T-C2: built into a local and published at the end, so the bliss run below
  // holds no lock. See the ShardedCache comment for why racing runs are
  // benign.
  container::svector<container::svector<int>> out;
  auto publish = [&]() -> container::svector<container::svector<int>> const& {
    return caches_->auts.put(Cluster{c}, std::move(out), h);
  };
  TermTable const& T = kt_->terms[c.d];
  // the only two places left that want real `Index` objects out of the key
  // table; both are per-cluster and cached, never per evaluation (T-A7)
  auto const face = T.canon_face_indices(c.S);
  const std::size_t nf = face.size();
  container::svector<int> id(nf);
  std::iota(id.begin(), id.end(), 0);
  out.push_back(id);
  if (nf < 2) return publish();

  // rebuild the cut's colored graph and collect automorphism generators
  container::vector<ExprPtr> ts;
  for (std::size_t v = 0; v < T.n(); ++v)
    if (c.S & (NodeMask{1} << v))
      ts.emplace_back(T.tensors[v]->as<Tensor>().clone());
  auto tn = TensorNetwork{ts};
  const FaceSet fs = T.face_set(c.S);  // must outlive create_graph
  auto graph = tn.create_graph({.named_indices = &fs,
                                .distinct_named_indices = false,
                                .make_idx_to_vertex = true});
  container::svector<std::size_t> face_vertex(nf);
  container::map<std::size_t, int> vertex_pos;
  for (std::size_t i = 0; i < nf; ++i) {
    face_vertex[i] = graph.idx_to_vertex.at(face[i]);
    vertex_pos.emplace(face_vertex[i], static_cast<int>(i));
  }
  container::svector<container::svector<int>> gens;
  auto save_aut = [&](unsigned /* n */, const unsigned* aut) {
    container::svector<int> p(nf);
    for (std::size_t i = 0; i < nf; ++i) {
      auto it = vertex_pos.find(aut[face_vertex[i]]);
      if (it == vertex_pos.end()) return;  // moves face out of face: skip
      p[i] = it->second;
    }
    gens.push_back(std::move(p));
  };
  bliss::Stats stats;
  graph.bliss_graph->find_automorphisms(
      stats, &bliss::aut_hook<decltype(save_aut)>, &save_aut);

  // close the generators into a group over face positions (cap: 48)
  for (std::size_t i = 0; i < out.size() && out.size() < 48; ++i)
    for (auto const& g : gens) {
      container::svector<int> comp(nf);
      for (std::size_t p = 0; p < nf; ++p) comp[p] = g[out[i][p]];
      if (ranges::find(out, comp) == ranges::end(out)) {
        out.push_back(std::move(comp));
        if (out.size() >= 48) break;
      }
    }
  return publish();
}

container::svector<BitPairs> const& Fitness::correspondences_bits(
    Cluster c1, Cluster c2) const {
  auto key = std::pair{c1, c2};
  const std::uint64_t h = cluster_pair_hash(key);
  if (auto const* hit = caches_->corr.find(key, h)) return *hit;
  // T-C2: local + publish, so no lock is held across `face_perms` (which takes
  // the `auts` shard lock of its own -- these two caches are never locked at
  // the same time, so no lock ordering exists to get wrong).
  container::svector<BitPairs> out;
  auto publish = [&]() -> container::svector<BitPairs> const& {
    return caches_->corr.put(std::move(key), std::move(out), h);
  };
  // position-wise alignment of the canonical face orders, composed with the
  // automorphisms of the first cut: sigma_pi(x) = sigma0(pi(x)).
  TermTable const& T1 = kt_->terms[c1.d];
  TermTable const& T2 = kt_->terms[c2.d];
  // The canonical faces ARE bit lists since T-A7, in slot order -- the
  // interning T-A6 did here per cluster pair now happens once, at build time.
  auto const& b1 = T1.canon_face_bits[c1.S];
  auto const& b2 = T2.canon_face_bits[c2.S];
  if (b1.size() != b2.size()) return publish();
  for (auto const& pi : face_perms(c1)) {
    BitPairs sigma;
    for (std::size_t i = 0; i < b1.size(); ++i)
      sigma.emplace_back(b1[static_cast<std::size_t>(pi[i])], b2[i]);
    // == the old `std::sort` on the `Index` pairs: within a term the `Index`
    // order IS the `index_rank` order, and a sigma's first components are
    // distinct (they are a permutation of the canonical face), so the first
    // component already decides every comparison.
    std::sort(sigma.begin(), sigma.end(), [&](auto const& p, auto const& q) {
      if (T1.index_rank[p.first] != T1.index_rank[q.first])
        return T1.index_rank[p.first] < T1.index_rank[q.first];
      return T2.index_rank[p.second] < T2.index_rank[q.second];
    });
    if (ranges::find(out, sigma) == ranges::end(out))
      out.push_back(std::move(sigma));
  }
  return publish();
}

container::svector<container::svector<std::pair<Index, Index>>> const&
Fitness::correspondences(Cluster c1, Cluster c2) const {
  auto key = std::pair{c1, c2};
  const std::uint64_t h = cluster_pair_hash(key);
  if (auto const* hit = caches_->corr_ix.find(key, h)) return *hit;
  auto const& bits = correspondences_bits(c1, c2);
  auto const& il1 = kt_->terms[c1.d].index_list;
  auto const& il2 = kt_->terms[c2.d].index_list;
  container::svector<container::svector<std::pair<Index, Index>>> out;
  for (auto const& sig : bits) {
    container::svector<std::pair<Index, Index>> sigma;
    for (auto const& [a, b] : sig) sigma.emplace_back(il1[a], il2[b]);
    out.push_back(std::move(sigma));
  }
  return caches_->corr_ix.put(std::move(key), std::move(out), h);
}

}  // namespace sequant::opt::ga
