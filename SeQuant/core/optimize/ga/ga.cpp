#include <SeQuant/core/optimize/ga/ga.hpp>

#include <SeQuant/core/optimize/ga/cost.hpp>
#include <SeQuant/core/runtime.hpp>

#include <range/v3/view/iota.hpp>

#include <bit>
#include <chrono>
#include <cstdio>
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

/// The exact single-term contraction-order DP, on the KeyTable's per-subset
/// index masks. Same recurrence as `opt::detail::single_term_opt<DenseFLOPs>`,
/// and four things must stay exact because they decide which tree comes out:
/// the cost (`merge_volume` + the `v == 1. -> 0.` scalar rule is
/// `flops_counter`, under default `CostParams`, so the seed is deliberately
/// volatility-BLIND even when the Fitness cost is not); the ascending subset
/// order; the bipartition order (the `a < S^a` half, `a` ascending, `(a-S)&S`
/// being the submask successor); and the `<=` TIE-BREAK, so the LAST
/// equal-cost bipartition wins. Child ordering needs no analogue: the laminar
/// family is determined by the SET of chosen splits.
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
// FAMILY is deliberately absent, so a kick can be drawn serially and applied
// later on another thread against whatever code the block holds by then;
// kicks to one block compose in recorded order (see ga_once).
struct KickPlan {
  std::uint32_t blk;  ///< term id (l1) or target id (!l1)
  std::uint32_t mv;   ///< move index in [0, nni_move_count(n))
  bool l1;
};

// Apply one planned kick to `x`, using (and warming) the caller's memo:
// decode the block, build ONLY the drawn neighbour, encode it, seed the memo,
// splice it in. Consumes no rng draw.
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
// each target's L2 tree comes intact from one parent. `coin` has one entry per
// block, g blocks first -- the order the coins were drawn in.
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

/// The evaluation workspaces the parallel phase hands out, one per concurrent
/// evaluation. A pool rather than `thread_local` because `sequant::for_each`
/// rebuilds its threads per call, which would destroy the decoded-tree memo.
/// A lease is taken per EVALUATION, so no thread identity is needed; WHICH
/// workspace an evaluation gets is nondeterministic and irrelevant, since
/// `EvalScratch` is pure scratch. `std::deque` because a lease holds a
/// reference into it while another evaluation may make it grow.
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
/// when `GAOptions::report_phases` is set.
struct PhaseStats {
  double hill_climb = 0, ga_once = 0;
  /// ga_once's own three phases: serial breeding (the Amdahl floor), parallel
  /// evaluation, serial sort/select.
  double breed = 0, eval = 0, select = 0;
};
double now_s() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

// One generation is plan / apply+evaluate / assemble, and it is EXACT: the
// whole search consumes one std::mt19937_64, so the DRAW SEQUENCE is the
// invariant. Phase 1 (serial) makes exactly the calls the fused serial loop
// made, in that order -- per child two `rng.below` (parents), one
// `rng.uniform` per block (coins), one `rng.exponential`, then per kick one
// `rng.uniform` (layer), one `rng.below` (block) and, iff n >= 3, one
// `rng.below(nni_move_count)` -- and only RECORDS the outcomes; an n < 3 kick
// still consumes its draws and is simply not recorded. Phase 2 (parallel)
// materializes each kid from its plan, which is byte-identical to building it
// serially because everything involved is a pure function of (genome bytes,
// plan), so WHOSE memo answers is unobservable. Phase 3 rebuilds `all` in the
// ORIGINAL child order because `std::sort` is not stable: the pre-sort order
// decides which of two equal-cost genomes survives.
Scored ga_once(Fitness const& F, Genome const& g0, GAOptions const& opts,
               std::uint64_t seed, PhaseStats* stats) {
  Rng rng(seed);
  auto const& kt = F.table();
  auto const& lay = F.layout();
  // member pointer, not a reference: `stats` may be null
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
  // Everything shared is read-only except the pre-sized output slots
  // `kids[i]`/`scores[i]`, which each task writes and no one else reads.
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

// One block's NNI neighbourhood is evaluated concurrently, and the parallel
// argmin is EXACT -- it reproduces the sequential scan bit for bit. The
// candidates are mutually independent (`nni_moves` is computed once from the
// block as it stands on entry, and each evaluation is the entry genome with
// this slice replaced by one candidate). The reduction is an argmin whose TIES
// ARE FIRST-WINS on strict `<`, which is why the scan is a separate SERIAL
// pass over a pre-sized array and not a reduction inside the parallel region:
// task completion order must not decide between equal-cost moves. Blocks and
// sweeps stay sequential, since accepting block d's move changes the genome
// block d+1 is evaluated against.
double hill_climb(Fitness const& F, Genome& genome, std::size_t max_sweeps) {
  auto const& kt = F.table();
  auto const& lay = F.layout();
  double best = F(genome);
  auto& memo = F.scratch().trees;
  ScratchPool pool(F);
  Laminar fam;
  ChildTable ch;
  // reused across blocks and sweeps: neighbour list, move i's code (at
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

    // ---- parallel: every candidate against the same entry genome, which is
    // only read here; each task writes its own disjoint output slots.
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

    // ---- serial argmin in ascending move order: ties go to the first
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
                                 GAOptions const& opts, SearchTrace* trace) {
  // Restarts stay sequential: they nest badly with the per-generation parallel
  // phase, and each one's rng stream is `opts.seed + 977 * s`.
  PhaseStats stats;
  const bool timed = opts.report_phases;
  // The seed's own cost is the per-term-equivalent baseline, and it is only
  // observable HERE: hill_climb rewrites `seed` in place, so after the next
  // line the genome that entered is gone. Costs one fitness evaluation, hence
  // only on request.
  if (trace) trace->seed_cost = F(seed);
  double t0 = timed ? now_s() : 0.;
  hill_climb(F, seed, opts.hill_climb_sweeps);
  if (timed) stats.hill_climb += now_s() - t0;
  Scored best{F(seed), seed};
  if (trace) trace->hill_climbed_cost = best.first;
  for (std::size_t s = 0; s < opts.restarts; ++s) {
    if (timed) t0 = now_s();
    auto b = ga_once(F, seed, opts, opts.seed + 977 * s,
                     timed ? &stats : nullptr);
    if (timed) stats.ga_once += now_s() - t0;
    if (trace) trace->restart_costs.push_back(b.first);
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
