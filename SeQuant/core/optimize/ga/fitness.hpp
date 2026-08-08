#ifndef SEQUANT_CORE_OPTIMIZE_GA_FITNESS_HPP
#define SEQUANT_CORE_OPTIMIZE_GA_FITNESS_HPP

// Genome evaluation: decode the L1 forests and L2 summation trees, run the
// always-extract distributive pass (factor a shared multiplicand out of two
// summands whenever the arrays match and an index bijection satisfying the
// prototype's (Sigma1)(Sigma2)(Sigma3) conditions exists), then resolve one
// producer per demanded array key and sum the DAG's merge flops, each key
// counted once. Faithful port of proto_csv_flat.Flat / proto_csv_opt.

#include <SeQuant/core/optimize/ga/cost.hpp>
#include <SeQuant/core/optimize/ga/genome.hpp>
#include <SeQuant/core/optimize/ga/key_table.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
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

/// Ambient axis tag: which target axis (or which enclosing-face index) an
/// open index of an L2 value corresponds to. Two summands may be added or
/// share a factored multiplicand only if their faces agree tag-by-tag.
///
/// One uint64 rather than a `std::variant<int, Index>` (T-A6). Two forms, told
/// apart by `ambient_eta_flag`:
///   * eta slot -- `ambient_eta_flag | (kind*4096 + rank)`: which axis of the
///     TARGET this index is (leaf ambients only);
///   * face index -- `(term << 8) | bit`: which index of the enclosing face it
///     was identified with by the parent extraction's beta.
/// Tags are only ever compared for equality, and the top byte is left free so
/// a (bit, tag) pair packs into one word in `Fitness::BetaKey`.
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
/// term `d`, the axis that index corresponds to.
///
/// Was `container::map<Index, AmbientTag>` (T-A6). Both the keys and the
/// `Index`-form tags were per-term index bits all along, and every consumer
/// either looks one up, compares two maps for equality, or -- at emission only
/// -- needs the real `Index` back, which is `kt.terms[d].index_list[bit]`.
///
/// Entries are kept sorted ascending by bit, i.e. in a CANONICAL order: that
/// is what lets `Fitness::BetaKey` compare two ambients as a flat word buffer.
/// (The old map was sorted by `Index` instead; the order is not otherwise
/// observable -- every consumer either searches by key or searches for a tag,
/// and tags are unique within a map.)
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
  // (T-A6 dropped a `beta` member here: the face bijection was stored on every
  // Fx node and read by nobody -- emission renames the second residual through
  // the ambient tags the beta INDUCES, which live on the residuals' own
  // `Summand::ambient`, never through the beta itself.)
};

/// The cost path's L2 value (T-A5): `Val` stripped to exactly what the L2 cost
/// recursion and `find_beta` read.
///
/// `Fitness::operator()` builds the SAME value tree `explain` builds, node for
/// node and in the same order, but out of these: no `shared_ptr` (children are
/// indices into the per-evaluation arena `EvalScratch::vals`), no
/// `Constant::scalar_type` coefficient (a `Complex<cpp_rational>` -- 3 % of the
/// search profile just to construct, and neither `cost_of_value` nor
/// `find_beta` ever reads it), no retained `beta`, and no `Summand` operands on
/// a Sum. What is left is a 24-byte POD.
struct CostVal {
  Val::Kind kind;
  int d;       ///< owning term of the face
  NodeMask S;  ///< Cl: the cluster; Fx: X (residual); Sum: the face cluster
  NodeMask V;  ///< Fx only: the extracted shared cluster
  int a;       ///< Fx: inner; Sum: s1.val -- index into EvalScratch::vals
  int b;       ///< Sum: s2.val -- index into EvalScratch::vals
};

// `ChildTable` and `build_children` moved to genome.hpp (the NNI slot walk is
// driven by the table now); the definitions are unchanged and everything here
// consumes them through that header.

/// The decoded per-genome state: one laminar family per term and the internal
/// clusters' children. The key fibres that used to live here are now
/// `EvalScratch::fib` (T-A4): they are pure per-evaluation scratch, consumed
/// only by `resolve`, and as a `container::map` they were the single most
/// expensive structure in the eval.
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
///
/// Why it pays: `cross` hands every block of a child genome intact from a
/// parent, `nni_kick` perturbs only k blocks and `hill_climb` changes exactly
/// one -- so of the 83 blocks an evaluation decodes (81 L1 terms + 2 L2
/// targets on the C4H10/DZ job) all but a handful were decoded before.
///
/// Why the key is the whole slice: `decode_tree(code, n)` is injective on
/// codes, so two distinct slices are two distinct trees. A hash-only or
/// truncated key would let one collision silently substitute a different tree
/// and shift the objective by an amount no test would attribute to the memo.
/// The 64-bit hash is only a probe accelerator; every candidate slot is
/// confirmed by comparing `n` and all n-1 genes.
///
/// Storage is three flat arenas (keys / families / children tables) addressed
/// by offsets, because every length is a function of `n`: the code has n-1
/// genes, the family 2n-1 masks, the children table n-1 entries. That keeps
/// one entry to a single contiguous run instead of three small_vectors, makes
/// copy-out a memcpy, and makes the cap enforceable by a single counter.
/// The children table is filled lazily: L2 blocks and the `nni_moves` decodes
/// in ga.cpp never ask for it, and a seeded entry (below) does not have it.
class TreeMemo {
 public:
  /// Entries kept before the memo is cleared wholesale. Reuse is temporal (a
  /// generation re-decodes its own blocks), so a bounded window costs almost
  /// no hit rate, and the cap is what keeps the memo off the peak-RSS budget
  /// -- which matters twice over once Group C gives every thread a scratch.
  /// Measured on C4H10/DZ (10.79 M lookups; bench peak RSS 1.636 GB without
  /// the memo), fingerprint identical at every setting:
  ///
  ///     cap        run_ga   peak RSS   hit rate   memo bytes   entries
  ///     unbounded  36.92 s   1.731 GB   100.00%      71.5 MB    64899
  ///     1<<15      36.95 s   1.681 GB    99.99%      35.8 MB    (2 clears)
  ///     8192       37.55 s   1.644 GB    99.94%       9.5 MB   (14 clears)
  ///
  /// so 1<<15 buys the whole speedup for 40 % less memory than unbounded.
  /// Override with SEQUANT_GA_MEMO_CAP (entries).
  static constexpr std::size_t default_capacity = 1u << 15;

  TreeMemo() : TreeMemo(env_capacity()) {}
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

  /// Record a family whose code the caller already holds -- `nni_kick` and
  /// `hill_climb` call `encode_tree(mv, n)` with `mv` in hand, so the entry is
  /// free. Requires decode_tree(code, n) == fam, i.e. `code` came from
  /// `encode_tree(fam, n)` (pinned exhaustively by the [ga] codec test).
  /// A no-op if the code is already known.
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
  /// Arena offsets are 32-bit (the Slot is then 24 B); keep every reachable
  /// offset strictly below the `no_children` sentinel.
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

  static std::size_t env_capacity() {
    if (char const* s = std::getenv("SEQUANT_GA_MEMO_CAP")) {
      const long v = std::atol(s);
      if (v > 0) return static_cast<std::size_t>(v);
    }
    return default_capacity;
  }

  /// FNV-1a-ish over (n, genes); `n` is mixed in so trees of different leaf
  /// counts do not share probe chains. Equality still compares everything.
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

  /// Wholesale eviction: an LRU would need a second index and the reuse here
  /// is temporal anyway (a generation's blocks are re-decoded within that
  /// generation), so dropping the window is both cheaper and simpler. Arena
  /// capacity is retained so the refill does not re-grow.
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
/// `[0, n_keys)` compared against a running epoch, so clearing is one
/// increment instead of touching 43k entries.
///
/// It replaces the `container::set<std::size_t>` working sets of `resolve`,
/// whose only use is the first-insertion gate. `insert` returns exactly what
/// `flat_set::insert(k).second` returned, so the gated walks keep their shape;
/// what is gone is the ordered storage (an O(k) memmove per insert) that the
/// gate never needed. Where an ORDER was also consumed -- `seen` in
/// `pick_out` -- the caller collects the ids and sorts them (see `resolve`).
class KeyStamps {
 public:
  /// Sizes the set for `n` key ids and leaves it empty.
  void resize(std::size_t n) {
    s_.assign(n, 0);
    cur_ = 1;  // stamps are 0, so nothing is a member yet
  }
  /// Empties the set. O(1) except once every 2^32 clears, when the stamps are
  /// rewound (the epoch must never alias a stale stamp).
  void clear() {
    if (++cur_ == 0) {
      std::fill(s_.begin(), s_.end(), 0);
      cur_ = 1;
    }
  }
  bool contains(std::size_t k) const { return s_[k] == cur_; }
  /// == `flat_set::insert(k).second`: true iff `k` was not already present.
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
/// grouped by the key id of the array it produces. Replaces
/// `container::map<std::size_t, container::svector<Cluster>>`, which took one
/// O(k) memmove per insert over the ~240-580 keys a genome touches.
///
/// Filled in two passes over the same cluster sequence -- count, then place --
/// so the clusters of one key land contiguously in a flat arena and `keys()`
/// is a short list of only the keys actually touched. Both dense arrays are
/// `n_keys` uint32s; the arena holds one `Cluster` per internal cluster
/// (~570 on C4H10/DZ), and neither is reallocated after the first evaluation.
///
/// **Walk order is load bearing** -- it is the floating-point summation order
/// of `resolve` -- and is reproduced exactly:
///   * `keys()` is sorted ascending, which is what iterating the `flat_map`
///     gave; `resolve` walks it to build `pick` and to enumerate `multi`.
///   * Within a key, clusters keep their PUSH order (term `d` ascending, then
///     `ForestState::Term::ch` order), because `place()` revisits the clusters
///     in exactly the order `count()` saw them and appends each at its key's
///     running cursor.
class Fibres {
 public:
  void resize(std::size_t n) {
    size_.assign(n, 0);
    off_.assign(n, 0);
  }

  /// Pass 0: drop the previous evaluation's fibres. Only the keys that were
  /// touched are reset -- that is what keeps this O(fibres), not O(n_keys).
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

  /// Touched key ids, ascending (== the old flat_map's iteration order).
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

/// Per-evaluation workspace.
///
/// Everything an evaluation needs beyond the immutable `KeyTable` lives here:
/// the decoded-tree memo (T-A3) and the dense per-key working sets of
/// `resolve` (T-A4: `fib`, `uses`, `walked`, `seen`, `pick`), which is why the
/// constructor takes the table -- those arrays are sized `n_keys`.
///
/// Two properties make this the seam the threading work (Group C) builds on:
///   * it is *pure scratch*. Nothing in it can change the number an evaluation
///     produces: the memo caches `decode_tree`, a pure function of (code, n),
///     and T-A4's arrays are cleared-by-epoch per call. So N scratches give N
///     bit-identical answers, and a per-thread scratch is correctness-neutral
///     by construction rather than by review.
///   * it is threaded by reference from the public entry points downwards, so
///     going thread-local is a call-site change (`make_scratch()` per worker),
///     never a signature change.
/// `Fitness` owns one `mutable` instance for its serial public API; that is
/// the ONLY shared-mutable state added here, and it disappears the moment a
/// caller passes its own.
struct EvalScratch {
  explicit EvalScratch(KeyTable const& kt) : n_keys(kt.n_keys) {
    // Key lists are held as uint32 (43k keys on C4H10/DZ; the subset tables
    // themselves would be astronomically larger long before this binds). Runs
    // once per scratch, so it is unconditional rather than an `assert`: it is
    // the only guard on the truncation `Fibres`/`KeyStamps` depend on, and the
    // shipping build compiles asserts out.
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

  /// `resolve`'s working sets, dense over `[0, n_keys)`. Each is logically
  /// empty at the start of the call that uses it and is emptied by bumping an
  /// epoch, never by touching n_keys entries. Payload arrays (`uses`, `pick`)
  /// are only meaningful where the paired `KeyStamps` says so.
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

  /// The cost path's L2 workspace (T-A5). All of it is rewound, never freed,
  /// per evaluation, so the steady state allocates nothing at all for L2.
  std::vector<CostVal> vals;  ///< arena: the whole L2 value forest of one eval
  /// Ambient maps the cost recursion builds for the extracted residuals. A
  /// `std::deque` because a `CostSummand` holds a pointer into it that must
  /// survive the pool growing underneath a deeper recursion level -- and
  /// because, unlike a vector of `unique_ptr`, it leaves `EvalScratch`
  /// copyable. Entries are reused across evaluations (`flat_map::clear()`
  /// keeps its buffer), which is what takes the residual ambients off the
  /// allocator entirely.
  std::deque<AmbientMap> amb_pool;
  std::size_t amb_used = 0;              ///< high-water mark within `amb_pool`
  std::vector<int> l2_roots;             ///< one `vals` index per target
  container::svector<Cluster> l2_demanded;  ///< clusters demanded by L2
  Laminar l2_fam;                        ///< the decoded L2 tree of one target

  /// Rewinds the L2 arena and the ambient pool. Called once per evaluation --
  /// NOT per target: a target's root `CostSummand` still points into the pool
  /// while the next target is built.
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
  /// Appends a value to the arena and returns its index.
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
  ForestState forest;
};

class Fitness {
 public:
  Fitness(KeyTable const& kt, CostModel cost = CostModel::native(),
          ProducerResolution resolution = ProducerResolution::Greedy);

  /// The cost of \p genome. Bit-for-bit `explain(genome).total`, but it never
  /// materializes the `Schedule`: the L2 pass runs on `CostVal`s in the
  /// scratch arena and `resolve` is asked for the number only (T-A5). The two
  /// paths are kept in lockstep by construction -- `build_node_cost` mirrors
  /// `build_node` node for node, so the `find_beta` call sequence and the L2
  /// recursion order are identical -- and by the exact-equality tripwire in
  /// the `[ga]` "cost path agrees with explain" test case, which must be
  /// extended, never relaxed, by any future change to either path.
  double operator()(Genome const& genome) const;
  Schedule explain(Genome const& genome) const;

  /// Same evaluations against a caller-owned workspace. Identical results --
  /// the scratch is pure scratch (see EvalScratch) -- so a thread that owns
  /// one needs no synchronization for anything reachable from here except the
  /// `Caches` below (T-C2).
  double operator()(Genome const& genome, EvalScratch& scratch) const;
  Schedule explain(Genome const& genome, EvalScratch& scratch) const;

  /// A fresh workspace sized for this table; Group C gives each worker one.
  EvalScratch make_scratch() const { return EvalScratch(*kt_); }
  /// The workspace behind the serial public API above (also used by the search
  /// in ga.cpp to decode blocks and to seed the memo).
  EvalScratch& scratch() const { return scratch_; }

  GenomeLayout const& layout() const { return layout_; }
  KeyTable const& table() const { return *kt_; }
  CostModel const& cost() const { return cost_; }

  /// All valid axis correspondences F(c1) -> F(c2) between two key-equal
  /// clusters (the prototype's sigma set): position-wise alignment of the
  /// canonical face orders, composed with the automorphisms of the cut.
  ///
  /// The evaluation path uses the bit form (`correspondences_bits`); this is
  /// the `Index` form emission and the sigma-soundness tests want, materialized
  /// from it on demand and cached separately.
  container::svector<container::svector<std::pair<Index, Index>>> const&
  correspondences(Cluster c1, Cluster c2) const;

  /// The eta ambient of term \p d's leaf value (external index bit -> target
  /// axis slot). Emission needs the same map to rename a target's externals.
  AmbientMap const& leaf_ambient(int d) const {
    return leaf_ambient_[static_cast<std::size_t>(d)];
  }

 private:
  /// The cost path's `Summand`: a value (index into `EvalScratch::vals`) plus
  /// the ambient identification `find_beta` needs. The map is BORROWED -- it
  /// is either a precomputed per-term leaf map or one of the scratch's pooled
  /// residual maps -- because every one of them is written once and then only
  /// read, so the copy `Summand` makes buys nothing but `Index` copies.
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
  /// The first (Sigma1)(Sigma2)(Sigma3)-valid bijection F(x1) -> F(x2) in the
  /// enumeration order pinned by `TermTable::index_rank`, as bit pairs owned by
  /// the beta cache; null if there is none. Memoized on (x1, v1, x2, v2, a1,
  /// a2) -- the cache is node-based, so the returned pointer outlives any
  /// number of later insertions.
  BitPairs const* find_beta(Cluster x1, Cluster v1, AmbientMap const& a1,
                            Cluster x2, Cluster v2, AmbientMap const& a2) const;
  /// `correspondences` in the bit form the evaluation path uses.
  container::svector<BitPairs> const& correspondences_bits(Cluster c1,
                                                           Cluster c2) const;
  container::svector<container::svector<int>> const& face_perms(
      Cluster c) const;
  /// Accumulates \p val's L2 cost and the clusters it demands. All L2 work is
  /// charged at the replay weight: emission leaves Fx/Sum structure inside the
  /// target expression (only keyed L1 clusters become named definitions), and
  /// targets are rebuilt every replay.
  void cost_of_value(ValPtr const& val, double& l2,
                     container::svector<Cluster>& demanded) const;
  /// The cost-path twin of `cost_of_value`, over the `CostVal` arena. Run as a
  /// SECOND pass over the finished value forest, exactly as `explain` runs
  /// `cost_of_value` after building all roots -- accumulating during
  /// construction would both reorder the floating-point sum and over-count,
  /// since a successful extraction discards its operands' `Cl` values.
  void cost_of_cost_val(EvalScratch const& scratch, int val, double& l2,
                        container::svector<Cluster>& demanded) const;
  double resolve(ForestState const& st, EvalScratch& scratch,
                 container::svector<Cluster> const& demanded,
                 container::map<std::size_t, Cluster>* pick_out) const;

  KeyTable const* kt_;
  CostModel cost_;
  ProducerResolution resolution_;
  GenomeLayout layout_;
  /// Per-term leaf data, precomputed once (T-A5). `leaf_summand` rebuilt both
  /// of these for every target leaf of every evaluation -- a `std::sort` of the
  /// externals, a rank map, one `kind_of` per external, one `make_shared<Val>`
  /// -- although they are pure functions of the KeyTable. The `Val` is now
  /// SHARED by every use of term `d`'s leaf: `ValPtr` is `shared_ptr<Val
  /// const>` so it cannot be mutated, and neither the evaluation nor emission
  /// keys off a `Val`'s address.
  container::svector<AmbientMap> leaf_ambient_;
  container::svector<ValPtr> leaf_val_;
  // caches (semantically const; keyed on data that outlives them)
  struct BetaKey;
  struct BetaKeyHash;
  struct Caches;
  std::shared_ptr<Caches> caches_;
  /// Serial workspace. Declared last: its constructor reads `*kt_`.
  mutable EvalScratch scratch_;
};

}  // namespace sequant::opt::ga

#endif  // SEQUANT_CORE_OPTIMIZE_GA_FITNESS_HPP
