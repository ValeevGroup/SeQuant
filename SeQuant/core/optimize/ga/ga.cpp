#include <SeQuant/core/optimize/ga/ga.hpp>

#include <SeQuant/core/optimize/ga/cost.hpp>
#include <SeQuant/core/runtime.hpp>

#include <range/v3/view/iota.hpp>

#include <bit>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <limits>
#include <mutex>
#include <random>
#include <vector>

namespace sequant::opt::ga {

namespace {

/// Per-subset DP scratch, reused across terms: `assign` refills in place, so
/// the buffers only grow when a term is wider than every term before it.
struct SeedScratch {
  container::vector<double> cost;
  container::vector<NodeMask> split;
  container::svector<NodeMask> stack;
};

/// The exact single-term contraction-order DP, run directly on the KeyTable's
/// per-subset index masks (T-B1).
///
/// This replaces a per-term call to
/// `opt::detail::single_term_opt<DenseFLOPs>(TensorNetwork{T.tensors}, ...)`,
/// which rebuilt an `Index` set per candidate merge and cost 20.6 s on the
/// 81-term C4H10/cc-pVDZ universe; `merge_volume` reads the masks
/// `build_key_table` already computed, and the whole universe now seeds in
/// ~0.1 s. It is the SAME recurrence, and the reproduction is exact -- see the
/// four points below, each pinned to the code it mirrors:
///
///  1. **cost.** `AdditiveModel::relax` (cost_model.hpp:124-160) charges
///     `w * flops_counter(F(a), F(b), F(a|b)) + fp + cost[a] + cost[b]`.
///     Seeding passes default `CostParams`, so `volatile_mask` is 0 (w == 1)
///     and `footprint_weight` is 0 (fp == 0) -- the seed is deliberately
///     volatility-BLIND even when the Fitness cost is volatility-aware, and
///     that is preserved here. `merge_volume` + the `v == 1. -> 0.` scalar rule
///     is `flops_counter` (single_term_detail.hpp:112-129) exactly.
///  2. **subset order.** Ascending `S`, matching `solve_single_term`'s `n` loop
///     (cost_model.hpp:23-42).
///  3. **bipartition order.** `bits::bipartitions(S)` (algorithm.hpp:179-215)
///     is the first half of the ordered bipartitions, i.e. the pairs
///     `(a, S^a)` with `a` running over the submasks of `S` in ASCENDING
///     numeric order while `S^a` descends -- so the half is exactly `a < S^a`,
///     and the driver drops the `a == 0` pair. `(a - S) & S` is the ascending
///     submask successor; `a > S^a` ends the half.
///  4. **tie-break.** `relax` accepts on `new_cost <= acc.ops`
///     (cost_model.hpp:161), so the LAST equal-cost bipartition wins.
///
/// `AdditiveModel::reconstruct`'s child ordering (cost_model.hpp:192-213) has
/// no analogue here and needs none: the laminar family -- and hence
/// `encode_tree`'s output -- is determined by the SET of chosen splits.
TreeCode seed_tree(TermTable const& T, SeedScratch& scr) {
  const int n = static_cast<int>(T.n());
  if (n < 3) return TreeCode(n - 1, 0);  // 1 or 2 leaves: the only tree
  const NodeMask full = T.full();
  const std::size_t N = std::size_t{1} << n;
  scr.cost.assign(N, 0.);  // singletons (and the unused empty set) cost 0
  scr.split.assign(N, 0);
  for (NodeMask S = 1; S < N; ++S) {
    if (std::popcount(S) < 2) continue;
    double best = std::numeric_limits<double>::max();
    NodeMask best_a = 0;
    for (NodeMask a = S & (~S + 1); a; a = (a - S) & S) {
      const NodeMask b = S ^ a;
      if (a > b) break;
      const double v = merge_volume(T, a, b);
      const double c = scr.cost[a] + scr.cost[b] + (v == 1. ? 0. : v);
      if (c <= best) {  // last equal-cost bipartition wins
        best = c;
        best_a = a;
      }
    }
    scr.cost[S] = best;
    scr.split[S] = best_a;
  }
  // the laminar family of the optimal tree: every singleton, plus every
  // cluster reachable from `full` through the stored splits
  Laminar fam;
  fam.reserve(2 * n - 1);
  for (int v = 0; v < n; ++v) fam.push_back(NodeMask{1} << v);
  scr.stack.assign(1, full);
  while (!scr.stack.empty()) {
    const NodeMask S = scr.stack.back();
    scr.stack.pop_back();
    fam.push_back(S);
    const NodeMask a = scr.split[S], b = S ^ a;
    if (std::popcount(a) >= 2) scr.stack.push_back(a);
    if (std::popcount(b) >= 2) scr.stack.push_back(b);
  }
  canon_sort(fam);
  return encode_tree(std::move(fam), n);
}

struct Rng {
  std::mt19937_64 eng;
  explicit Rng(std::uint64_t seed) : eng(seed) {}
  double uniform() { return std::uniform_real_distribution<>{}(eng); }
  std::size_t below(std::size_t n) {
    return std::uniform_int_distribution<std::size_t>{0, n - 1}(eng);
  }
  int exponential(double mean) {
    return static_cast<int>(std::exponential_distribution<>{1.0 / mean}(eng));
  }
};

// One planned NNI kick: which layer, which block, which move index. The move
// FAMILY is deliberately absent -- since `nni_move_count` became a function of
// the leaf count alone, the rng draw a kick consumes no longer needs the
// kicked block decoded, so a kick can be DRAWN serially and APPLIED later, on
// another thread, against whatever code the block holds by then (kicks to one
// block compose in recorded order; see ga_once).
struct KickPlan {
  std::uint32_t blk;  ///< term id (l1) or target id (!l1)
  std::uint32_t mv;   ///< move index in [0, nni_move_count(n))
  bool l1;
};

// Apply one planned kick to `x`, using (and warming) the caller's memo.
// This is the old `nni_kick` body minus every rng draw: decode the block
// (memo, with its cached ChildTable), build ONLY the drawn neighbour (T-A2c),
// encode it, seed the memo with the fresh code, splice it in.
void apply_kick(Fitness const& F, Genome& x, KickPlan const& kk,
                EvalScratch& sc, Laminar& fam, ChildTable& ch) {
  auto const& kt = F.table();
  auto const& lay = F.layout();
  auto& code = kk.l1 ? x.g : x.h;
  const auto slice = kk.l1 ? lay.g_slice[kk.blk] : lay.h_slice[kk.blk];
  const int n = kk.l1
                    ? static_cast<int>(kt.terms[kk.blk].n())
                    : static_cast<int>(kt.targets[kk.blk].terms.size());
  sc.trees.decode(code.data() + slice.first, n, fam, ch);
  const Laminar mv = nni_move_at(fam, ch, kk.mv);
  auto next = encode_tree(mv, n);
  // free warm entry: `next` is exactly the code the next decode of this
  // block will present, and `mv` is exactly what it would decode to.
  sc.trees.seed(next.data(), n, mv);
  std::copy(next.begin(), next.end(), code.begin() + slice.first);
}

// whole-block crossover with the coins already drawn: each term's L1 tree and
// each target's L2 tree comes intact from one parent. `coin` has one entry
// per block, g blocks first -- the order the old `cross` drew its
// `rng.uniform() < .5` per block.
Genome cross_apply(Fitness const& F, Genome const& a, Genome const& b,
                   std::uint8_t const* coin) {
  auto const& lay = F.layout();
  Genome c = a;
  for (auto const& [lo, hi] : lay.g_slice)
    if (*coin++)
      std::copy(b.g.begin() + lo, b.g.begin() + hi, c.g.begin() + lo);
  for (auto const& [lo, hi] : lay.h_slice)
    if (*coin++)
      std::copy(b.h.begin() + lo, b.h.begin() + hi, c.h.begin() + lo);
  return c;
}

using Scored = std::pair<double, Genome>;

/// The evaluation workspaces the parallel phase of `ga_once` hands out, one
/// per concurrent evaluation (T-C3).
///
/// `sequant::for_each` exposes no thread id, so the obvious `thread_local`
/// cannot be indexed by one -- and, less obviously, it would also be the wrong
/// LIFETIME: `for_each` creates a fresh set of `std::thread`s on every call
/// and joins them, so a `thread_local` scratch is destroyed and rebuilt once
/// per generation. That is ~2100 constructions of a ~1.9 MB dense workspace
/// over a default search, and it throws away the decoded-tree memo (T-A3)
/// every generation -- precisely the structure whose whole value is that it
/// accumulates. Both variants were built and measured on C4H10/DZ at 8
/// threads, fingerprint identical either way:
///
///     phase 2 (evaluation)   thread_local 0.90 s   this pool 0.73 s
///     ga_once total          thread_local 3.43 s   this pool 3.25 s
///     bench peak RSS         thread_local 0.564    this pool 0.562-0.579 GB
///
/// i.e. the pool buys 19 % of the phase it governs, for a memory difference
/// that is inside the run-to-run spread of the peak (the pool holds
/// num_threads() workspaces per restart: 1.89 MB dense + ~7-10 MB of memo
/// each, freed between restarts). The memo cap never binds at either setting
/// -- 8.5k entries against the 32768 default, 0 clears -- so there is nothing
/// to buy back by lowering SEQUANT_GA_MEMO_CAP here.
///
/// A lease is taken per EVALUATION rather than per thread, so no thread
/// identity is needed at all: two uncontended mutex operations against a
/// ~46 us evaluation. WHICH scratch a given evaluation gets is therefore
/// nondeterministic -- and irrelevant, because `EvalScratch` is pure scratch
/// (see its docs): nothing in it can change the number an evaluation
/// produces, so N workspaces give N bit-identical answers.
///
/// `std::deque` because a lease holds a reference into it while another
/// evaluation may cause it to grow.
class ScratchPool {
 public:
  explicit ScratchPool(Fitness const& F) : F_(&F) {}
  EvalScratch& take() {
    std::lock_guard<std::mutex> lock(m_);
    if (free_.empty()) {
      slots_.push_back(F_->make_scratch());
      return slots_.back();
    }
    EvalScratch& s = *free_.back();
    free_.pop_back();
    return s;
  }
  void give(EvalScratch& s) {
    std::lock_guard<std::mutex> lock(m_);
    free_.push_back(&s);
  }

 private:
  Fitness const* F_;
  std::mutex m_;
  std::deque<EvalScratch> slots_;
  std::vector<EvalScratch*> free_;
};

/// RAII lease on one of the pool's workspaces.
class Lease {
 public:
  explicit Lease(ScratchPool& p) : p_(p), s_(p.take()) {}
  ~Lease() { p_.give(s_); }
  Lease(Lease const&) = delete;
  Lease& operator=(Lease const&) = delete;
  EvalScratch& operator*() const { return s_; }

 private:
  ScratchPool& p_;
  EvalScratch& s_;
};

/// Optional `hill_climb` / `ga_once` wall-clock split, printed by `run_ga`
/// when SEQUANT_GA_PHASE_STATS is set. `run_ga` is the sum of the two and only
/// `ga_once` is parallel (T-C3), so the split is what any judgement about the
/// speedup has to be made against.
struct PhaseStats {
  double hill_climb = 0, ga_once = 0;
  /// ga_once's own three phases: serial breeding, parallel evaluation, serial
  /// sort/select. `breed` is the Amdahl floor of T-C3 and the number T-C4 and
  /// anything after it has to be sized against.
  double breed = 0, eval = 0, select = 0;
  static bool enabled() {
    static const bool on = std::getenv("SEQUANT_GA_PHASE_STATS") != nullptr;
    return on;
  }
};
double now_s() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

// --- T-C3: one generation is plan / apply+evaluate / assemble ---------------
//
// The whole search consumes ONE std::mt19937_64, and every decision after a
// draw depends on every draw before it, so the only restructuring that is
// allowed here is one that leaves the draw SEQUENCE untouched. The fused loop
// this evolved from did, per child: two `rng.below` (parents), `cross` (one
// `rng.uniform` per L1 and per L2 block), one `rng.exponential`, then per
// kick: one `rng.uniform` (layer), one `rng.below` (block), and -- iff the
// block has n >= 3 leaves -- one `rng.below(nni_move_count)` for the move.
//
// Since `nni_move_count` became a closed form of the leaf count, EVERY value
// in that sequence is now a function of earlier draws and the immutable
// layout alone -- nothing needs a genome decoded, or even materialized. So
// phase 1 makes exactly those calls in exactly that order and merely RECORDS
// the outcomes (parent indices, crossover coins, kick plans); building the
// children moved into the parallel phase, where each task materializes its
// own kid -- cross_apply, then its kicks in recorded order on its leased
// workspace -- and scores it. The kid a task builds is byte-identical to the
// one the serial loop built: cross_apply consumes the same coins, kicks to
// one block compose in the same order, and every function involved is a pure
// function of (genome bytes, plan) -- the memo is a cache of the pure
// decode_tree/build_children pair, so WHOSE memo answers is unobservable.
//
// The n < 3 guard mirrors the old kick_block exactly: such a kick consumed
// its layer/block draws but neither drew a move nor changed anything, so it
// is drawn here and simply not recorded. (n >= 3 implies the move count
// 2(n-2) >= 2 > 0, so the old `n_moves == 0` early-out cannot fire there.)
//
// Phase 3 rebuilds `all` in the ORIGINAL child order -- kids ascending, then
// the incumbent population -- because `std::sort` is not stable: the pre-sort
// order decides which of two equal-cost genomes survives. `std::shuffle` is
// untouched and consumes the same number of draws, since that depends only on
// the element count.
Scored ga_once(Fitness const& F, Genome const& g0, GAOptions const& opts,
               std::uint64_t seed, PhaseStats* stats) {
  Rng rng(seed);
  auto const& kt = F.table();
  auto const& lay = F.layout();
  // member pointer, not a reference: `stats` may be null and the argument
  // would then be a dereference of it at the call site
  auto tick = [&](double PhaseStats::*into, double t0) {
    if (stats) stats->*into += now_s() - t0;
  };
  double t0 = stats ? now_s() : 0.;
  ScratchPool pool(F);

  // the recorded plans of one batch of children (reused across generations)
  const std::size_t n_blocks = lay.g_slice.size() + lay.h_slice.size();
  std::vector<std::pair<std::uint32_t, std::uint32_t>> parents;  // per child
  std::vector<std::uint8_t> coins;  // n_blocks per child, g blocks first
  std::vector<KickPlan> kicks;      // flat arena
  std::vector<std::pair<std::uint32_t, std::uint32_t>> kick_span;  // per child
  auto plan_reset = [&] {
    parents.clear();
    coins.clear();
    kicks.clear();
    kick_span.clear();
  };
  // draw one child's kicks; records only the ones that do something
  auto plan_kicks = [&](int k) {
    const auto lo = static_cast<std::uint32_t>(kicks.size());
    for (int i = 0; i < k; ++i) {
      if (rng.uniform() < opts.l1_move_prob) {
        const auto d = rng.below(kt.terms.size());
        const int n = static_cast<int>(kt.terms[d].n());
        if (n < 3) continue;
        kicks.push_back({static_cast<std::uint32_t>(d),
                         static_cast<std::uint32_t>(
                             rng.below(nni_move_count(n))),
                         true});
      } else {
        const auto t = rng.below(kt.targets.size());
        const int n = static_cast<int>(kt.targets[t].terms.size());
        if (n < 3) continue;
        kicks.push_back({static_cast<std::uint32_t>(t),
                         static_cast<std::uint32_t>(
                             rng.below(nni_move_count(n))),
                         false});
      }
    }
    kick_span.emplace_back(lo, static_cast<std::uint32_t>(kicks.size()));
  };

  container::svector<Genome> kids;
  container::svector<double> scores;
  container::svector<Scored> pop;
  // phase 2: materialize child i from its plan and score it, in parallel.
  // Shared with sibling tasks: the `Fitness` (sharded Caches, read-only
  // KeyTable), the read-only plan arrays and `pop`/`base`, and the pre-sized
  // output slots `kids[i]`/`scores[i]` each task writes and no one else reads.
  auto materialize_and_evaluate = [&](Genome const* base) {
    const std::size_t nk = kick_span.size();
    kids.resize(nk);
    scores.assign(nk, 0.);
    if (nk == 0) return;
    auto ids = ranges::views::iota(std::size_t{0}, nk);
    sequant::for_each(ids, [&](std::size_t i) {
      Lease scratch(pool);
      Genome kid = base ? *base
                        : cross_apply(F, pop[parents[i].first].second,
                                      pop[parents[i].second].second,
                                      coins.data() + i * n_blocks);
      Laminar fam;  // per-task temps, reused across this child's kicks
      ChildTable ch;
      for (auto j = kick_span[i].first; j < kick_span[i].second; ++j)
        apply_kick(F, kid, kicks[j], *scratch, fam, ch);
      scores[i] = F(kid, *scratch);
      kids[i] = std::move(kid);
    });
  };

  pop.reserve(opts.population);
  pop.emplace_back(F(g0), g0);
  for (std::size_t i = 1; i < opts.population; ++i)
    plan_kicks(1 + static_cast<int>((i - 1) % 12));
  materialize_and_evaluate(&g0);
  for (std::size_t i = 0; i < kids.size(); ++i)
    pop.emplace_back(scores[i], std::move(kids[i]));
  auto by_cost = [](Scored const& a, Scored const& b) {
    return a.first < b.first;
  };
  std::sort(pop.begin(), pop.end(), by_cost);
  Scored best = pop.front();
  const std::size_t n_elite = std::max<std::size_t>(
      1, static_cast<std::size_t>((1 - opts.mixing) * opts.population));
  for (std::size_t gen = 0; gen < opts.generations; ++gen) {
    const double frac =
        1.0 - static_cast<double>(gen) /
                  std::max<std::size_t>(1, opts.generations - 1);
    // ---- phase 1 (serial): draw, consuming `rng` in exactly the old order
    if (stats) t0 = now_s();
    plan_reset();
    for (std::size_t i = 0; i < opts.reproduction * opts.population; ++i) {
      // sequenced explicitly: function-argument evaluation order is
      // unspecified, and parent a's draw must precede parent b's
      const auto pa = static_cast<std::uint32_t>(rng.below(pop.size()));
      const auto pb = static_cast<std::uint32_t>(rng.below(pop.size()));
      parents.emplace_back(pa, pb);
      for (std::size_t bcoin = 0; bcoin < n_blocks; ++bcoin)
        coins.push_back(rng.uniform() < .5 ? 1 : 0);
      const int k = 1 + rng.exponential(std::max(0.3, opts.kick0 * frac));
      plan_kicks(k);
    }
    tick(&PhaseStats::breed, t0);

    // ---- phase 2 (parallel): materialize and score them
    if (stats) t0 = now_s();
    materialize_and_evaluate(nullptr);
    tick(&PhaseStats::eval, t0);

    // ---- phase 3 (serial): assemble in the original order, then select
    if (stats) t0 = now_s();
    container::svector<Scored> all;
    all.reserve(opts.reproduction * opts.population + pop.size());
    for (std::size_t i = 0; i < kids.size(); ++i)
      all.emplace_back(scores[i], std::move(kids[i]));
    for (auto& s : pop) all.push_back(std::move(s));
    std::sort(all.begin(), all.end(), by_cost);
    // elites survive; the rest of the next population is sampled uniformly
    // from the remainder to retain genetic variability (EP&H)
    container::svector<Scored> next(all.begin(), all.begin() + n_elite);
    std::shuffle(all.begin() + n_elite, all.end(), rng.eng);
    for (std::size_t i = 0; next.size() < opts.population; ++i)
      next.push_back(std::move(all[n_elite + i]));
    std::sort(next.begin(), next.end(), by_cost);
    pop = std::move(next);
    if (pop.front().first < best.first) best = pop.front();
    tick(&PhaseStats::select, t0);
  }
  return best;
}

}  // namespace

Genome seed_genome(KeyTable const& kt) {
  Genome out;
  SeedScratch scr;
  for (auto const& T : kt.terms) {
    const TreeCode code = seed_tree(T, scr);
    out.g.insert(out.g.end(), code.begin(), code.end());
  }
  for (auto const& tgt : kt.targets)
    out.h.insert(out.h.end(), tgt.terms.size() - 1, 0);
  return out;
}

// --- T-C4: one block's NNI neighbourhood is evaluated concurrently ----------
//
// Why the parallel argmin is EXACT, i.e. reproduces the sequential loop this
// replaces bit for bit.
//
//  1. *The candidates are mutually independent.* `nni_moves` is computed ONCE,
//     from the block as it stands when `try_block` is entered, so every
//     candidate is a family in that one list. The old loop spliced candidate i
//     into the slice, called `F(genome)`, and -- when it did not improve --
//     restored the slice before candidate i+1. It never touched any OTHER
//     block. So the i-th evaluation was always
//     `F(entry genome with this slice replaced by cand_i)`, whatever happened
//     for j < i: the incumbent `block` the loop restores to is itself one of
//     the cand_j, so the value that comes back for cand_i does not depend on
//     it. That is exactly what each task computes here, against the shared,
//     unmutated `genome`.
//  2. *The reduction is an argmin.* `block` was reassigned only under
//     `c < best`, and `best` only ever decreases, so the loop selected the
//     minimum over the whole move list -- and left `code` holding it, because
//     an improving candidate is the one splice that is not rolled back.
//  3. *Ties are first-wins.* The comparison is strict `<`, so among equal
//     minimal costs the LOWEST move index is kept. The serial scan below runs
//     `i` ascending with the identical `if (c < best)` rule and therefore
//     makes the identical choice -- this, not the evaluation, is where a
//     parallel version could silently change the answer, which is why the
//     scan is a separate serial pass over a pre-sized array rather than a
//     reduction inside the parallel region. Task completion order is
//     irrelevant: no task reads `best`.
//  4. *`best` still threads through the blocks.* It is the running best of the
//     whole hill climb, updated in the serial scan, and blocks and sweeps stay
//     sequential -- accepting block d's move changes the genome block d+1 is
//     evaluated against.
//
// The evaluations run on leased pool workspaces (T-C3), never on
// `F.scratch()`. The one behavioural detail that has to move with them is the
// memo seeding: the old loop handed `encode_tree`'s output straight to
// `F.scratch().trees` so the very next `F(genome)` decode of the block was a
// hit. Each task now seeds the workspace IT is about to evaluate on, which is
// the same trade in the same place; the accepted block is additionally seeded
// into `F.scratch()`, whose only remaining job here is the one `decode` per
// block per sweep. Measured on C4H10/DZ: the serial workspace's memo goes from
// 2.80 M lookups at 99.98 % to 0.71 M at 99.90 % (the 2.08 M that left are
// exactly hill_climb's evaluations), and the workspaces they landed on run at
// 99.85-100.00 % with ~110 misses each. No scratch loses its memo, and the
// total miss count over all of them (~1.5 k at 8 threads against 692 before)
// is 0.05 % of the lookups either way.
//
// No serial threshold is applied, because the batch sizes never get small
// enough to need one. `sequant::for_each` costs ~65 us per call at 8 threads
// on this machine (measured: 79 us at batch 1, 64 us at batch 106, i.e. it is
// spawn/join latency, not per-task work) against a ~36 us evaluation, so a
// block breaks even at 2.1 candidates and every block that reaches this code
// has at least 2 (n >= 3, and an NNI move list is empty or has >= 2 entries).
// The L1 blocks -- mean 14 candidates -- are where the win is thin, not
// negative: ~4x predicted, and the measured whole-`hill_climb` speedup is
// dominated by them since they are 88 % of the evaluations.
double hill_climb(Fitness const& F, Genome& genome, std::size_t max_sweeps) {
  auto const& kt = F.table();
  auto const& lay = F.layout();
  double best = F(genome);
  auto& memo = F.scratch().trees;
  ScratchPool pool(F);
  Laminar fam;
  ChildTable ch;
  // reused across blocks and sweeps: the neighbour list, move i's code (at
  // cands[i*(n-1) .. (i+1)*(n-1))) and move i's cost
  container::svector<Laminar> moves;
  std::vector<int> cands;
  std::vector<double> costs;
  // `layer` rather than a reference to the code: a task needs to splice into
  // the same layer of its OWN genome copy.
  auto try_block = [&](container::svector<int> Genome::*layer,
                       std::pair<int, int> slice, int n) {
    if (n < 3) return false;
    auto& code = genome.*layer;
    const int lo = slice.first;
    const int w = slice.second - slice.first;  // == n - 1 genes per block
    memo.decode(code.data() + lo, n, fam, ch);
    moves = nni_moves(fam, ch);
    if (moves.empty()) return false;

    // ---- parallel: every candidate against the same entry genome. `genome`
    // is only read here; the two outputs are pre-sized and every task writes
    // its own disjoint slot and reads no one else's.
    cands.resize(moves.size() * static_cast<std::size_t>(w));
    costs.assign(moves.size(), 0.);
    auto ids = ranges::views::iota(std::size_t{0}, moves.size());
    sequant::for_each(ids, [&](std::size_t i) {
      Lease scratch(pool);
      const TreeCode cand = encode_tree(moves[i], n);
      // free warm entry, into the workspace that is about to decode it
      (*scratch).trees.seed(cand.data(), n, moves[i]);
      Genome x = genome;
      std::copy(cand.begin(), cand.end(), (x.*layer).begin() + lo);
      std::copy(cand.begin(), cand.end(),
                cands.begin() + static_cast<std::ptrdiff_t>(i) * w);
      costs[i] = F(x, *scratch);
    });

    // ---- serial: today's rule, in ascending move order (first wins)
    std::size_t win = moves.size();
    for (std::size_t i = 0; i < moves.size(); ++i)
      if (costs[i] < best) {
        best = costs[i];
        win = i;
      }
    if (win == moves.size()) return false;
    int const* won = cands.data() + win * static_cast<std::size_t>(w);
    std::copy(won, won + w, code.begin() + lo);
    memo.seed(won, n, moves[win]);
    return true;
  };
  for (std::size_t sweep = 0; sweep < max_sweeps; ++sweep) {
    bool improved = false;
    for (std::size_t d = 0; d < kt.terms.size(); ++d)
      improved |= try_block(&Genome::g, lay.g_slice[d],
                            static_cast<int>(kt.terms[d].n()));
    for (std::size_t t = 0; t < kt.targets.size(); ++t)
      improved |= try_block(&Genome::h, lay.h_slice[t],
                            static_cast<int>(kt.targets[t].terms.size()));
    if (!improved) break;
  }
  return best;
}

std::pair<double, Genome> run_ga(Fitness const& F, Genome seed,
                                 GAOptions const& opts) {
  // The restarts stay sequential: they nest badly with T-C3's per-generation
  // parallel phase, and each one's rng stream is a function of the one before
  // it only through `opts.seed + 977 * s`, so there is nothing to gain here
  // that the generation loop does not already give.
  PhaseStats stats;
  const bool timed = PhaseStats::enabled();
  double t0 = timed ? now_s() : 0.;
  hill_climb(F, seed, opts.hill_climb_sweeps);
  if (timed) stats.hill_climb += now_s() - t0;
  Scored best{F(seed), seed};
  for (std::size_t s = 0; s < opts.restarts; ++s) {
    if (timed) t0 = now_s();
    auto b = ga_once(F, seed, opts, opts.seed + 977 * s,
                     timed ? &stats : nullptr);
    if (timed) stats.ga_once += now_s() - t0;
    if (b.first < best.first) best = std::move(b);
  }
  if (timed) t0 = now_s();
  hill_climb(F, best.second, opts.hill_climb_sweeps);
  if (timed) stats.hill_climb += now_s() - t0;
  best.first = F(best.second);
  if (timed)
    std::fprintf(stderr,
                 "[ga phases] threads=%d hill_climb=%.2f s ga_once=%.2f s"
                 " (breed=%.2f eval=%.2f select=%.2f)\n",
                 num_threads(), stats.hill_climb, stats.ga_once, stats.breed,
                 stats.eval, stats.select);
  return best;
}

}  // namespace sequant::opt::ga
