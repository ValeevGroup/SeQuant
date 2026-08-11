#ifndef SEQUANT_CORE_OPTIMIZE_GA_FITNESS_HPP
#define SEQUANT_CORE_OPTIMIZE_GA_FITNESS_HPP

// Genome evaluation: decode the L1 forests and L2 summation trees, run the
// always-extract distributive pass (factor a shared multiplicand out of two
// summands whenever the arrays match and a (Sigma1)(Sigma2)(Sigma3) index
// bijection exists), then resolve one producer per demanded array key and sum
// the DAG's merge flops, each key counted once. Ported from proto_csv_opt.

#include <SeQuant/core/optimize/ga/cost.hpp>
#include <SeQuant/core/optimize/ga/genome.hpp>
#include <SeQuant/core/optimize/ga/key_table.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace sequant::opt::ga {

struct Genome {
  container::svector<int> g;  ///< per-term L1 codes, concatenated
  container::svector<int> h;  ///< per-target L2 codes, concatenated
};

/// Gene offsets and ranges induced by a KeyTable.
struct GenomeLayout {
  container::svector<std::pair<int, int>> g_slice;  ///< per term [lo, hi)
  container::svector<std::pair<int, int>> h_slice;  ///< per target [lo, hi)
  container::svector<int> g_ranges;  ///< gene k of a block has range 2k-3
  container::svector<int> h_ranges;

  static GenomeLayout of(KeyTable const& kt);
};

/// A cluster of one term: (term id, subset of its tensors).
struct Cluster {
  int d;
  NodeMask S;
  friend auto operator<=>(Cluster const&, Cluster const&) = default;
};

/// Ambient axis tag: which target axis (or enclosing-face index) an open index
/// of an L2 value corresponds to; two summands may be added or share a factored
/// multiplicand only if their faces agree tag-by-tag. Two forms, told apart by
/// `ambient_eta_flag`: an eta slot `ambient_eta_flag | (kind*4096 + rank)`
/// naming a TARGET axis (leaf ambients only), or a face index
/// `(term << 8) | bit`. Tags are only compared for equality; the top byte stays
/// free so a (bit, tag) pair packs into one word in `Fitness::BetaKey`.
using AmbientTag = std::uint64_t;
inline constexpr AmbientTag ambient_eta_flag = AmbientTag{1} << 55;
/// \p slot is the `kind*4096 + rank` eta slot of a target axis.
inline constexpr AmbientTag eta_tag(int slot) {
  return ambient_eta_flag | static_cast<AmbientTag>(slot);
}
/// The index at \p bit of term \p d.
inline constexpr AmbientTag face_tag(int d, unsigned bit) {
  return (static_cast<AmbientTag>(static_cast<unsigned>(d)) << 8) | bit;
}

/// The ambient identification of one L2 value's face: for each index BIT of
/// term `d`, the axis that index corresponds to. Entries are kept sorted
/// ascending by bit, i.e. in a CANONICAL order -- that is what lets
/// `Fitness::BetaKey` compare two ambients as a flat word buffer.
struct AmbientMap {
  int d = 0;  ///< the term whose index bits are the keys
  container::svector<std::pair<std::uint8_t, AmbientTag>> e;

  void clear() {
    d = 0;
    e.clear();
  }
  /// Appends an entry; the batch must be closed with `seal()`.
  void push(std::uint8_t bit, AmbientTag tag) { e.emplace_back(bit, tag); }
  /// Restores the canonical ascending-bit order after a batch of `push`es.
  void seal() { std::sort(e.begin(), e.end()); }
  /// The tag of \p bit, or null if this ambient does not cover it.
  AmbientTag const* find(std::uint8_t bit) const {
    for (auto const& kv : e)
      if (kv.first == bit) return &kv.second;
    return nullptr;
  }
};

/// A face bijection as bit pairs: `find_beta`'s beta (bits of term x1.d ->
/// bits of term x2.d) and `correspondences`' sigma alike.
using BitPairs = container::svector<std::pair<std::uint8_t, std::uint8_t>>;

/// L2 value: 'cl' = an L1 array; 'fx' = (inner) * A_{d,V}, a factored sum
/// whose parent face is F(d, X|V); 'sum' = a plain elementwise addition.
struct Val;
using ValPtr = std::shared_ptr<Val const>;

struct Summand {
  Constant::scalar_type coeff = 1;
  ValPtr val;
  AmbientMap ambient;
};

struct Val {
  enum Kind { Cl, Fx, Sum };
  Kind kind;
  int d;        ///< owning term of the face
  NodeMask S;   ///< Cl: the cluster; Fx: X (residual); Sum: the face cluster
  NodeMask V;   ///< Fx only: the extracted shared cluster
  ValPtr inner; ///< Fx only: combined residual value
  Summand s1, s2;  ///< Sum only (ambients retained for emission)
};

/// The cost path's L2 value: `Val` stripped to what the L2 cost recursion and
/// `find_beta` read. `Fitness::operator()` builds the SAME value tree
/// `explain` builds, node for node and in the same order, out of these.
struct CostVal {
  Val::Kind kind;
  int d;       ///< owning term of the face
  NodeMask S;  ///< Cl: the cluster; Fx: X (residual); Sum: the face cluster
  NodeMask V;  ///< Fx only: the extracted shared cluster
  int a;       ///< Fx: inner; Sum: s1.val -- index into EvalScratch::vals
  int b;       ///< Sum: s2.val -- index into EvalScratch::vals
};

/// The decoded per-genome state: one laminar family per term and the internal
/// clusters' children. The key fibres are per-evaluation scratch and live in
/// `EvalScratch::fib`.
struct ForestState {
  struct Term {
    Laminar fam;
    /// per internal cluster: (S, c1, c2) with c1 < c2 lexicographically
    ChildTable ch;
    /// The stored (c1, c2) pair of internal cluster \p S, in STORED lex order
    /// (`tree_children` returns the same pair canon-larger first; the two
    /// orders are not interchangeable -- each caller's order feeds emission).
    /// Throws `std::logic_error` if \p S is not an internal cluster.
    NodeMask const* children_of(NodeMask S) const;
  };
  container::svector<Term> terms;
};

/// Memo for `decode_tree` keyed on the FULL code slice plus its leaf count.
/// The key must be the whole slice: a hash-only or truncated key would let a
/// collision silently substitute a different tree and shift the objective, so
/// the 64-bit hash is only a probe accelerator and every candidate slot is
/// confirmed gene by gene. Storage is three flat arenas addressed by offsets;
/// the children table is filled lazily.
class TreeMemo {
 public:
  /// Entries kept before the memo is cleared wholesale; measured, not guessed
  /// (reuse is temporal, so a bounded window costs no hit rate and keeps the
  /// memo off the peak-RSS budget). Override with `GAOptions::memo_capacity`.
  static constexpr std::size_t default_capacity = 1u << 15;

  TreeMemo() : TreeMemo(default_capacity) {}
  explicit TreeMemo(std::size_t max_entries)
      : cap_(max_entries ? max_entries : 1), tab_(1024) {}

  /// fam := decode_tree(code[0, n-1), n).
  void decode(int const* code, int n, Laminar& fam) {
    lookup(code, n, fam, nullptr);
  }
  /// fam := decode_tree(code[0, n-1), n), ch := build_children(fam).
  void decode(int const* code, int n, Laminar& fam, ChildTable& ch) {
    lookup(code, n, fam, &ch);
  }

  /// Record a family whose code the caller already holds. REQUIRES
  /// decode_tree(code, n) == fam, i.e. `code` came from `encode_tree(fam, n)`;
  /// a bogus pair would poison every later decode. No-op if already known.
  void seed(int const* code, int n, Laminar const& fam) {
    if (n < 2) return;  // see lookup(): nothing to remember
    const std::uint64_t h = hash_of(code, n);
    const std::size_t i = probe(h, code, n);
    if (tab_[i].n != 0) return;
    store(i, h, code, n, fam, nullptr);
  }

  std::size_t hits() const { return hits_; }
  std::size_t misses() const { return misses_; }
  std::size_t clears() const { return clears_; }

 private:
  static constexpr std::uint32_t no_children = ~std::uint32_t{0};
  /// Arena offsets are 32-bit; keep them below the `no_children` sentinel.
  static constexpr std::size_t arena_max = (std::size_t{1} << 32) - 1024;
  struct Slot {
    std::uint64_t h = 0;
    std::uint32_t key = 0;  ///< offset into keys_ (n-1 genes)
    std::uint32_t fam = 0;  ///< offset into fams_ (2n-1 masks)
    std::uint32_t ch = no_children;  ///< offset into chs_ (n-1 entries)
    std::int32_t n = 0;              ///< 0 = empty slot
  };

  static std::size_t fam_len(int n) { return std::size_t(2 * n - 1); }
  static std::size_t ch_len(int n) { return std::size_t(n - 1); }
  static std::size_t key_len(int n) { return std::size_t(n - 1); }

  /// FNV-1a-ish over (n, genes); equality still compares everything.
  static std::uint64_t hash_of(int const* code, int n) {
    std::uint64_t h = 0x9e3779b97f4a7c15ull ^
                      (static_cast<std::uint64_t>(n) * 0xff51afd7ed558ccdull);
    for (std::size_t i = 0, e = key_len(n); i < e; ++i) {
      h ^= static_cast<std::uint32_t>(code[i]);
      h *= 0x100000001b3ull;
      h ^= h >> 29;
    }
    h ^= h >> 32;
    h *= 0xc4ceb9fe1a85ec53ull;
    h ^= h >> 32;
    return h;
  }

  /// Index of the matching slot, or of the empty slot where it would go.
  std::size_t probe(std::uint64_t h, int const* code, int n) const {
    const std::size_t mask = tab_.size() - 1;
    std::size_t i = h & mask;
    while (true) {
      Slot const& s = tab_[i];
      if (s.n == 0) return i;
      if (s.h == h && s.n == n &&
          std::equal(code, code + key_len(n), keys_.data() + s.key))
        return i;
      i = (i + 1) & mask;
    }
  }

  void lookup(int const* code, int n, Laminar& fam, ChildTable* ch) {
    if (n < 2) {  // 0 or 1 leaves: decode_tree's initial family, no internals.
      fam.assign(std::size_t{1}, NodeMask{1});  // (and the length formulas
      if (ch) ch->clear();                      //  below would underflow)
      return;
    }
    const std::uint64_t h = hash_of(code, n);
    const std::size_t i = probe(h, code, n);
    if (tab_[i].n != 0) {
      ++hits_;
      Slot& s = tab_[i];
      NodeMask const* f = fams_.data() + s.fam;
      fam.assign(f, f + fam_len(n));
      if (ch) {
        if (s.ch == no_children) {  // seeded (or L2-only) entry: fill it in
          build_children(fam, *ch);
          if (chs_.size() + ch_len(n) <= arena_max) {
            s.ch = static_cast<std::uint32_t>(chs_.size());
            chs_.insert(chs_.end(), ch->begin(), ch->end());
          }
        } else {
          auto const* c = chs_.data() + s.ch;
          ch->assign(c, c + ch_len(n));
        }
      }
      return;
    }
    ++misses_;
    TreeCode c(code, code + key_len(n));
    fam = decode_tree(c, n);
    if (ch) build_children(fam, *ch);
    store(i, h, code, n, fam, ch);
  }

  void store(std::size_t i, std::uint64_t h, int const* code, int n,
             Laminar const& fam, ChildTable const* ch) {
    if (n_entries_ >= cap_ || keys_.size() + key_len(n) > arena_max ||
        fams_.size() + fam_len(n) > arena_max ||
        chs_.size() + ch_len(n) > arena_max) {
      clear();
      i = probe(h, code, n);
    } else if ((n_entries_ + 1) * 2 > tab_.size()) {
      rehash(tab_.size() * 2);
      i = probe(h, code, n);
    }
    Slot& s = tab_[i];
    s.h = h;
    s.n = n;
    s.key = static_cast<std::uint32_t>(keys_.size());
    keys_.insert(keys_.end(), code, code + key_len(n));
    s.fam = static_cast<std::uint32_t>(fams_.size());
    fams_.insert(fams_.end(), fam.begin(), fam.end());
    if (ch) {
      s.ch = static_cast<std::uint32_t>(chs_.size());
      chs_.insert(chs_.end(), ch->begin(), ch->end());
    } else {
      s.ch = no_children;
    }
    ++n_entries_;
  }

  void rehash(std::size_t n_slots) {
    std::vector<Slot> next(n_slots);
    const std::size_t mask = n_slots - 1;
    for (Slot const& s : tab_) {
      if (s.n == 0) continue;
      std::size_t i = s.h & mask;
      while (next[i].n != 0) i = (i + 1) & mask;
      next[i] = s;
    }
    tab_.swap(next);
  }

  /// Wholesale eviction (reuse is temporal, so no LRU index is worth its
  /// cost). Arena capacity is retained so the refill does not re-grow.
  void clear() {
    std::fill(tab_.begin(), tab_.end(), Slot{});
    keys_.clear();
    fams_.clear();
    chs_.clear();
    n_entries_ = 0;
    ++clears_;
  }

  std::size_t cap_;
  std::vector<Slot> tab_;  ///< power-of-two, linear probing, load factor <= .5
  std::vector<int> keys_;
  std::vector<NodeMask> fams_;
  std::vector<std::array<NodeMask, 3>> chs_;
  std::size_t n_entries_ = 0;
  std::size_t hits_ = 0, misses_ = 0, clears_ = 0;
};

/// A set of key ids with O(1) membership AND O(1) clear: dense stamps over
/// `[0, n_keys)` against a running epoch. Used for `resolve`'s first-insertion
/// gates; it is UNORDERED, so wherever `resolve` also consumes an order the
/// caller collects the ids and sorts them explicitly.
class KeyStamps {
 public:
  void resize(std::size_t n) {
    s_.assign(n, 0);
    cur_ = 1;  // stamps are 0, so nothing is a member yet
  }
  /// Empties the set. The epoch must never alias a stale stamp, hence the
  /// rewind once every 2^32 clears.
  void clear() {
    if (++cur_ == 0) {
      std::fill(s_.begin(), s_.end(), 0);
      cur_ = 1;
    }
  }
  bool contains(std::size_t k) const { return s_[k] == cur_; }
  /// True iff `k` was not already present.
  bool insert(std::size_t k) {
    if (s_[k] == cur_) return false;
    s_[k] = cur_;
    return true;
  }

 private:
  std::vector<std::uint32_t> s_;
  std::uint32_t cur_ = 1;
};

/// The key fibres of one decoded forest: every internal cluster of every term,
/// grouped by the key id of the array it produces, filled in two passes over
/// the same cluster sequence (count, then place).
///
/// **Walk order is load bearing: it IS the floating-point summation order of
/// `resolve`.** `keys()` is ascending, and within a key clusters keep their
/// PUSH order (term `d` ascending, then `ForestState::Term::ch` order).
class Fibres {
 public:
  void resize(std::size_t n) {
    size_.assign(n, 0);
    off_.assign(n, 0);
  }

  /// Pass 0: drop the previous evaluation's fibres (touched keys only).
  void begin() {
    for (std::uint32_t k : keys_) size_[k] = 0;
    keys_.clear();
    cl_.clear();
  }
  /// Pass 1: one call per (cluster, key), in push order.
  void count(std::size_t k) {
    if (size_[k]++ == 0) keys_.push_back(static_cast<std::uint32_t>(k));
  }
  /// Between the passes: fix the ascending key order and the arena layout.
  void seal() {
    std::sort(keys_.begin(), keys_.end());
    std::uint32_t run = 0;
    for (std::uint32_t k : keys_) {
      off_[k] = run;
      run += size_[k];
    }
    cl_.resize(run);
  }
  /// Pass 2: one call per (cluster, key), in the SAME order as `count`.
  void place(std::size_t k, Cluster c) { cl_[off_[k]++] = c; }
  /// After pass 2: rewind the cursors so `off_[k]` is the fibre's start again.
  void close() {
    for (std::uint32_t k : keys_) off_[k] -= size_[k];
  }

  /// Touched key ids, ascending.
  std::vector<std::uint32_t> const& keys() const { return keys_; }
  std::size_t size(std::size_t k) const { return size_[k]; }
  /// The fibre of `k`, in push order.
  Cluster const* data(std::size_t k) const { return cl_.data() + off_[k]; }

 private:
  std::vector<std::uint32_t> size_;  ///< per key: fibre length (0 = untouched)
  std::vector<std::uint32_t> off_;   ///< per key: arena start (cursor in pass 2)
  std::vector<std::uint32_t> keys_;  ///< touched keys, ascending after seal()
  std::vector<Cluster> cl_;          ///< fibres, concatenated in key order
};

/// Per-evaluation workspace: the decoded-tree memo and `resolve`'s dense
/// per-key working sets, which is why the constructor takes the table.
///
/// **It is pure scratch**: nothing in it can change the number an evaluation
/// produces, so N workspaces give N bit-identical answers -- which is what
/// makes the parallel evaluation in ga.cpp correctness-neutral by
/// construction. `Fitness` owns one `mutable` instance for its serial API.
struct EvalScratch {
  explicit EvalScratch(KeyTable const& kt,
                       std::size_t memo_capacity = TreeMemo::default_capacity)
      : trees(memo_capacity), n_keys(kt.n_keys) {
    // Unconditional (asserts are compiled out of the shipping build): this is
    // the only guard on the uint32 truncation Fibres/KeyStamps depend on.
    if (n_keys > std::size_t{0xffffffffu})
      throw std::length_error(
          "ga: key count exceeds the uint32 key-id representation");
    fib.resize(n_keys);
    uses.assign(n_keys, 0);
    uses_set.resize(n_keys);
    walked.resize(n_keys);
    seen.resize(n_keys);
    pick.assign(n_keys, Cluster{0, 0});
    pick_set.resize(n_keys);
  }
  EvalScratch(EvalScratch&&) = default;
  EvalScratch& operator=(EvalScratch&&) = default;
  EvalScratch(EvalScratch const&) = default;
  EvalScratch& operator=(EvalScratch const&) = default;

  TreeMemo trees;
  std::size_t n_keys = 0;

  /// `resolve`'s working sets, dense over `[0, n_keys)` and emptied by bumping
  /// an epoch. Payload arrays (`uses`, `pick`) are only meaningful where the
  /// paired `KeyStamps` says so.
  Fibres fib;                       ///< produced by decode_forest
  std::vector<std::int32_t> uses;   ///< consumer count per key (dag_cost pass 1)
  KeyStamps uses_set;
  KeyStamps walked;                 ///< dag_cost pass 1 first-visit gate
  KeyStamps seen;                   ///< dag_cost pass 2 first-visit gate
  std::vector<std::uint32_t> seen_list;  ///< pass-2 keys, in visit order
  std::vector<Cluster> pick;             ///< chosen producer per key
  KeyStamps pick_set;
  std::vector<std::size_t> seeds;  ///< resolve: demanded keys, sorted unique
  std::vector<std::size_t> stack;  ///< dag_cost: the LIFO walk stack

  /// The cost path's L2 workspace: rewound, never freed, per evaluation.
  std::vector<CostVal> vals;  ///< arena: the whole L2 value forest of one eval
  /// Ambient maps for the extracted residuals. A `std::deque` because a
  /// `CostSummand` holds a pointer into it that must survive the pool growing.
  std::deque<AmbientMap> amb_pool;
  std::size_t amb_used = 0;              ///< high-water mark within `amb_pool`
  std::vector<int> l2_roots;             ///< one `vals` index per target
  container::svector<Cluster> l2_demanded;  ///< clusters demanded by L2
  Laminar l2_fam;                        ///< the decoded L2 tree of one target

  /// Rewinds the L2 arena and the ambient pool. Once per EVALUATION, not per
  /// target: a target's root still points into the pool while the next is built.
  void begin_l2() {
    vals.clear();
    amb_used = 0;
  }
  /// The next free ambient map, emptied.
  AmbientMap& new_ambient() {
    if (amb_used == amb_pool.size()) amb_pool.emplace_back();
    AmbientMap& m = amb_pool[amb_used++];
    m.clear();
    return m;
  }
  int new_val(Val::Kind kind, int d, NodeMask S, NodeMask V, int a = -1,
              int b = -1) {
    vals.push_back(CostVal{kind, d, S, V, a, b});
    return static_cast<int>(vals.size()) - 1;
  }
};

enum class ProducerResolution { Greedy, Exact };

/// Result of one evaluation, sufficient to emit the schedule.
struct Schedule {
  double total = 0, l1 = 0, l2 = 0;
  container::svector<Summand> roots;  ///< one per target
  container::map<std::size_t, Cluster> pick;  ///< key -> producing cluster
  /// key -> how many sites `resolve` charged it to: the L2 demands (with
  /// multiplicity) plus one per parent slot in the resolved DAG. `>= 2` is the
  /// `shared` bit the objective feeds `CostModel::runtime_amortized`, and it
  /// must equal emission's render-site count (`ga::emission_use_counts`) key
  /// for key -- the [ga] "persistent naming matches the objective" test pins
  /// that, because the matched pair is only matched if the two counts agree.
  /// Filled by `explain` only (the cost path never materializes it).
  container::map<std::size_t, int> uses;
  ForestState forest;
};

class Fitness {
 public:
  Fitness(KeyTable const& kt, CostModel cost = CostModel::native(),
          ProducerResolution resolution = ProducerResolution::Greedy,
          std::size_t memo_capacity = TreeMemo::default_capacity);

  /// The cost of \p genome. Bit-for-bit `explain(genome).total` without ever
  /// materializing the `Schedule`. The two paths stay in lockstep by
  /// construction (`build_node_cost` mirrors `build_node` node for node) and by
  /// the exact-equality tripwire in the `[ga]` "cost path agrees with explain"
  /// test, which must be extended, never relaxed.
  double operator()(Genome const& genome) const;
  Schedule explain(Genome const& genome) const;

  /// Same evaluations against a caller-owned workspace. Identical results --
  /// the scratch is pure scratch -- so a thread that owns one needs no
  /// synchronization for anything reachable from here except the `Caches`.
  double operator()(Genome const& genome, EvalScratch& scratch) const;
  Schedule explain(Genome const& genome, EvalScratch& scratch) const;

  /// A fresh workspace sized for this table; one per concurrent evaluation.
  EvalScratch make_scratch() const {
    return EvalScratch(*kt_, memo_capacity_);
  }
  /// The workspace behind the serial public API above.
  EvalScratch& scratch() const { return scratch_; }

  GenomeLayout const& layout() const { return layout_; }
  KeyTable const& table() const { return *kt_; }
  CostModel const& cost() const { return cost_; }

  /// All valid axis correspondences F(c1) -> F(c2) between two key-equal
  /// clusters (the prototype's sigma set): position-wise alignment of the
  /// canonical face orders composed with the automorphisms of the cut. The
  /// `Index` form of `correspondences_bits`.
  container::svector<container::svector<std::pair<Index, Index>>> const&
  correspondences(Cluster c1, Cluster c2) const;

  /// The eta ambient of term \p d's leaf value (external index bit -> target
  /// axis slot). Emission needs the same map to rename a target's externals.
  AmbientMap const& leaf_ambient(int d) const {
    return leaf_ambient_[static_cast<std::size_t>(d)];
  }

 private:
  /// The cost path's `Summand`: a value (index into `EvalScratch::vals`) plus
  /// the ambient `find_beta` needs. The map is BORROWED -- a per-term leaf map
  /// or a pooled residual map, each written once and thereafter only read.
  struct CostSummand {
    int val;
    AmbientMap const* amb;
  };

  ForestState decode_forest(Genome const& genome, EvalScratch& scratch) const;
  Summand leaf_summand(int d) const;
  Summand build_node(ForestState const& st, Summand s1, Summand s2) const;
  /// The cost-path twin of `build_node`. Same options in the same order, same
  /// `find_beta` arguments in the same sequence, same resulting value tree --
  /// only the representation differs. Edit the two together.
  CostSummand build_node_cost(ForestState const& st, EvalScratch& scratch,
                              CostSummand s1, CostSummand s2) const;
  /// The FIRST (Sigma1)(Sigma2)(Sigma3)-valid bijection F(x1) -> F(x2) in the
  /// enumeration order pinned by `TermTable::index_rank`; null if there is
  /// none. Memoized on (x1, v1, x2, v2, a1, a2) in a node-based cache, so the
  /// returned pointer outlives any number of later insertions.
  BitPairs const* find_beta(Cluster x1, Cluster v1, AmbientMap const& a1,
                            Cluster x2, Cluster v2, AmbientMap const& a2) const;
  /// `correspondences` in the bit form the evaluation path uses.
  container::svector<BitPairs> const& correspondences_bits(Cluster c1,
                                                           Cluster c2) const;
  container::svector<container::svector<int>> const& face_perms(
      Cluster c) const;
  /// Accumulates \p val's L2 cost and the clusters it demands. All L2 work is
  /// charged at the replay weight: only keyed L1 clusters become named
  /// definitions, so Fx/Sum structure is rebuilt on every replay.
  void cost_of_value(ValPtr const& val, double& l2,
                     container::svector<Cluster>& demanded) const;
  /// The cost-path twin of `cost_of_value`, over the `CostVal` arena. Must run
  /// as a SECOND pass over the finished forest, as `explain` does: accumulating
  /// during construction would reorder the sum and over-count.
  void cost_of_cost_val(EvalScratch const& scratch, int val, double& l2,
                        container::svector<Cluster>& demanded) const;
  /// \param pick_out if non-null, receives the producer of every walked key;
  ///        \param uses_out (only read when \p pick_out is) receives the
  ///        matching per-key charge counts (\ref Schedule::uses).
  double resolve(ForestState const& st, EvalScratch& scratch,
                 container::svector<Cluster> const& demanded,
                 container::map<std::size_t, Cluster>* pick_out,
                 container::map<std::size_t, int>* uses_out = nullptr) const;

  KeyTable const* kt_;
  CostModel cost_;
  ProducerResolution resolution_;
  GenomeLayout layout_;
  /// Per-term leaf data, precomputed once: pure functions of the KeyTable. The
  /// `Val` is SHARED by every use of term `d`'s leaf -- sound because `ValPtr`
  /// is `shared_ptr<Val const>` and nothing keys off a `Val`'s address.
  container::svector<AmbientMap> leaf_ambient_;
  container::svector<ValPtr> leaf_val_;
  // caches (semantically const; keyed on data that outlives them)
  struct BetaKey;
  struct BetaKeyHash;
  struct Caches;
  std::shared_ptr<Caches> caches_;
  /// Memo capacity handed to every workspace this object makes.
  std::size_t memo_capacity_ = TreeMemo::default_capacity;
  /// Serial workspace. Declared last: its constructor reads `*kt_`.
  mutable EvalScratch scratch_;
};

}  // namespace sequant::opt::ga

#endif  // SEQUANT_CORE_OPTIMIZE_GA_FITNESS_HPP
