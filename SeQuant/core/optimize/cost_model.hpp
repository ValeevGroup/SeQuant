#ifndef SEQUANT_CORE_OPTIMIZE_COST_MODEL_HPP
#define SEQUANT_CORE_OPTIMIZE_COST_MODEL_HPP

#include <SeQuant/core/optimize/single_term_detail.hpp>  // helpers + EvalSequence + OptRes

#include <range/v3/view/concat.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <utility>

namespace sequant::opt::detail {

/// \brief Fills the per-subset State table bottom-up via the model's hooks.
///
/// The driver owns only what every single-term DP shares: the subset lattice
/// (subsets in increasing order) and the bipartition enumeration. The model
/// supplies the per-objective recurrence through leaf/init/relax/finalize.
template <class Model, typename TIdxs>
container::vector<typename Model::State> solve_single_term(
    Model const& m, TensorNetwork const& network, TIdxs const& tidxs,
    typename Model::Context& ctx) {
  auto const nt = network.tensors().size();
  container::vector<typename Model::State> st(size_t{1} << nt);
  auto const connected =
      outer_product_connectivity(network, tidxs, m.prune_outer_products);
  for (size_t n = 1; n < st.size(); ++n) {
    if (std::popcount(n) == 1) {
      st[n] = m.leaf(ctx, n);
      continue;
    }
    if (!connected[n]) continue;  // never form a disconnected subset
    typename Model::State acc = m.init(ctx, n);
    for (auto&& [lp, rp] : bits::bipartitions(n))
      if (lp != 0 && rp != 0 && connected[lp] && connected[rp])
        m.relax(ctx, n, lp, rp, st[lp], st[rp], acc);
    st[n] = std::move(acc);
    m.finalize(ctx, n, st);
  }
  return st;
}

/// \brief Generic single-term optimization: build context, solve, reconstruct.
///
/// \tparam Model A type satisfying the CostModel concept (a built-in objective
///         such as \ref AdditiveModel, or a user-defined model).
/// \param m The model instance carrying the objective's parameters.
/// \param network The TensorNetwork containing the tensors to be contracted.
/// \param tidxs The set of indices that should remain open in the final result.
/// \return The optimal EvalSequence under the model's objective.
template <class Model, typename TIdxs>
EvalSequence run_single_term_opt(Model const& m, TensorNetwork const& network,
                                 TIdxs const& tidxs) {
  auto const nt = network.tensors().size();
  if (nt == 1) return EvalSequence{0};
  if (nt == 2) return EvalSequence{0, 1, -1};
  typename Model::Context ctx = m.build_context(network, tidxs);
  auto st = solve_single_term(m, network, tidxs, ctx);
  return m.reconstruct(ctx, st);
}

/// \brief Companion to \ref run_single_term_opt that also reports, for each
/// contraction (\c -1) node of the returned \c EvalSequence in emission order,
/// the sliced-set of batchable \c Index values realized at that node. Requires
/// \p Model to additionally expose \c reconstruct_axes (currently only \ref
/// PeakBatchedModel does); see its doc comment for the precise per-node
/// convention (RPN / post-order, left-first, matching the shared \c build
/// recursion used by \ref reconstruct).
///
/// \return The optimal EvalSequence, paired with one \c container::svector
///         of sliced \c Index per \c -1 token of that sequence, in the same
///         left-first post-order the sequence itself was emitted in. For the
///         nt==1 shortcut (no contractions) the axes vector is empty; for the
///         nt==2 shortcut (single contraction, no DP context is built) the
///         axes vector holds one empty entry (no batching info available).
template <class Model, typename TIdxs>
std::pair<EvalSequence,
          container::vector<container::svector<std::pair<Index, AxisKind>>>>
run_single_term_opt_axes(Model const& m, TensorNetwork const& network,
                         TIdxs const& tidxs) {
  auto const nt = network.tensors().size();
  if (nt == 1) return {EvalSequence{0}, {}};
  if (nt == 2)
    return {EvalSequence{0, 1, -1},
            container::vector<container::svector<std::pair<Index, AxisKind>>>{
                container::svector<std::pair<Index, AxisKind>>{}}};
  typename Model::Context ctx = m.build_context(network, tidxs);
  auto st = solve_single_term(m, network, tidxs, ctx);
  return m.reconstruct_axes(ctx, st);
}

/// \brief Additive single-term cost model (FLOPs or operand storage size).
///
/// Implements the additive single-term DP, factored into the CostModel hooks
/// driven by \ref run_single_term_opt. The per-contraction cost is supplied by
/// \p CostFn (e.g. \ref flops_counter or \ref memsize_counter); a
/// per-intermediate footprint penalty, volatile replay weighting, and subnet
/// common-subexpression elimination (CSE) are all handled by the hooks.
///
/// \tparam CostFn A callable <tt>(lhs, rhs, result) -> double</tt>.
/// \tparam FootprintFn A callable <tt>(result) -> double</tt>.
template <typename CostFn, typename FootprintFn>
struct AdditiveModel {
  CostFn cost_fn;
  FootprintFn footprint_fn;
  size_t volatile_mask;
  double volatile_weight;
  double footprint_weight;
  bool subnet_cse;
  /// Prune disconnected (outer-product) subsets from the DP (see
  /// OptimizeOptions::prune_outer_products). Default true.
  bool prune_outer_products = true;

  /// Per-subset DP cell: the relevant OptRes fields the additive DP mutates.
  struct State {
    double ops = std::numeric_limits<double>::max();
    size_t lp = 0;
    size_t rp = 0;
    container::vector<size_t> subnets;
  };

  /// Precomputed tables AND mutable DP scratch, built once by build_context.
  struct Context {
    /// Per-subset open indices (and scratch) from init_results.
    container::vector<OptRes> results;
    /// Canonical-subnet ids per subset (CSE only).
    container::vector<size_t> meta_ids;
    /// Optimal cost per canonical subnet id, populated during the DP (CSE
    /// only).
    container::vector<double> unique_meta_costs;
  };

  template <typename TIdxs>
  Context build_context(TensorNetwork const& network,
                        TIdxs const& tidxs) const {
    Context ctx;
    ctx.results.resize(size_t{1} << network.tensors().size());
    // Outer-product pruning: skip building index sets for disconnected subsets
    // the DP will never form (solve_single_term also skips them).
    auto const connected =
        outer_product_connectivity(network, tidxs, prune_outer_products);
    init_results(network, tidxs, ctx.results, &connected);
    if (subnet_cse) {
      auto md = build_subnet_metadata(network, ctx.results, connected);
      ctx.meta_ids = std::move(md.meta_ids);
      ctx.unique_meta_costs = std::move(md.unique_meta_costs);
    }
    return ctx;
  }

  State leaf(Context const& /*ctx*/, size_t /*n*/) const {
    // ops 0; the singleton sequence is implicit in reconstruct.
    return State{0.0, 0, 0, {}};
  }

  State init(Context const& /*ctx*/, size_t /*n*/) const {
    return State{std::numeric_limits<double>::max(), 0, 0, {}};
  }

  void relax(Context& ctx, size_t n, size_t lp, size_t rp, State const& lp_st,
             State const& rp_st, State& acc) const {
    // A subset is volatile iff it contains any volatile leaf; the contraction
    // that forms it is then re-executed on every replay. volatile_mask == 0
    // (no predicate / DenseSize / volatile_weight<=1) makes w == 1 everywhere.
    double const w = (volatile_mask & n) ? volatile_weight : 1.0;

    // Per-intermediate memory-footprint penalty (storage of THIS result),
    // added once and NOT scaled by the replay weight w (peak footprint is a
    // one-time materialization cost). Zero when footprint_weight == 0.
    double const fp =
        footprint_weight != 0.0
            ? footprint_weight * footprint_fn(ctx.results[n].indices)
            : 0.0;

    double new_cost = 0;
    container::vector<size_t> combined_subnets;
    if (subnet_cse) {
      // subnets is always kept sorted; set_union requires sorted inputs and
      // produces sorted output -- this invariant is maintained throughout.
      std::set_union(lp_st.subnets.begin(), lp_st.subnets.end(),
                     rp_st.subnets.begin(), rp_st.subnets.end(),
                     std::back_inserter(combined_subnets));
      new_cost = w * cost_fn(ctx.results[lp].indices,  //
                             ctx.results[rp].indices,  //
                             ctx.results[n].indices)   //
                 + fp;
      for (auto id : combined_subnets) {
        new_cost += ctx.unique_meta_costs[id];
      }
    } else {
      new_cost = w * cost_fn(ctx.results[lp].indices,  //
                             ctx.results[rp].indices,  //
                             ctx.results[n].indices)   //
                 + fp + lp_st.ops + rp_st.ops;
    }

    if (new_cost <= acc.ops) {
      acc.ops = new_cost;
      acc.lp = lp;
      acc.rp = rp;
      if (subnet_cse) {
        acc.subnets = std::move(combined_subnets);
      }
    }
  }

  void finalize(Context& ctx, size_t n, container::vector<State>& st) const {
    if (!subnet_cse) return;
    auto mid = ctx.meta_ids[n];
    // Canonically equivalent subnetworks share the same topology and index
    // sizes, so their cost is identical. Overwriting with a later bitmask's
    // cost is intentional and benign.
    // Recompute w exactly as relax does: a subset is volatile iff it contains
    // any volatile leaf, so the stored cost must use the same scaling.
    double const w = (volatile_mask & n) ? volatile_weight : 1.0;
    ctx.unique_meta_costs[mid] =
        w * cost_fn(ctx.results[st[n].lp].indices,
                    ctx.results[st[n].rp].indices, ctx.results[n].indices) +
        (footprint_weight != 0.0
             ? footprint_weight * footprint_fn(ctx.results[n].indices)
             : 0.0);
    auto it = std::lower_bound(st[n].subnets.begin(), st[n].subnets.end(), mid);
    if (it == st[n].subnets.end() || *it != mid) {
      st[n].subnets.insert(it, mid);
    }
  }

  EvalSequence reconstruct(Context const& /*ctx*/,
                           container::vector<State> const& st) const {
    using ranges::views::concat;
    // Rebuild per-subset sequences bottom-up: singletons emit their tensor
    // index, internal nodes emit the
    // children ordered by lseq[0] < rseq[0] followed by -1. Subsets are visited
    // in increasing order, so st[lp]/st[rp] (lp,rp < n) are already built.
    container::vector<EvalSequence> seq(st.size());
    for (size_t n = 0; n < st.size(); ++n) {
      auto const pc = std::popcount(n);
      if (pc == 1) {
        seq[n].emplace_back(std::countr_zero(n));
      } else if (pc >= 2) {
        // Pruned (outer-product) subsets are never relaxed, so lp/rp stay at
        // their default 0/0 sentinel (bits::bipartitions guarantees a real
        // split always assigns both nonzero together); skip them here too --
        // they are never referenced as a child by a subset actually on the
        // reconstructed path, since solve_single_term only relaxes into a
        // parent through children that are themselves connected.
        if (st[n].lp == 0 && st[n].rp == 0) continue;
        auto const& lseq = seq[st[n].lp];
        auto const& rseq = seq[st[n].rp];
        seq[n] = (lseq[0] < rseq[0] ? concat(lseq, rseq) : concat(rseq, lseq)) |
                 ranges::to<EvalSequence>;
        seq[n].push_back(-1);
      }
    }
    return seq.back();
  }
};

/// \brief Insert a (peak, flops) trade-off into a Pareto frontier with
/// domination pruning: skip the new point if an existing one is no worse in
/// both objectives; otherwise drop every existing point the new one dominates
/// and append it. \tparam FP any struct with \c peak and \c flops members.
template <typename FP>
void pareto_insert(container::vector<FP>& f, FP p) {
  for (auto const& e : f)
    if (e.peak <= p.peak && e.flops <= p.flops) return;  // dominated -> skip
  f.erase(std::remove_if(f.begin(), f.end(),
                         [&](FP const& e) {
                           return p.peak <= e.peak && p.flops <= e.flops;
                         }),
          f.end());
  f.push_back(p);
}

/// \brief Index of the lexicographic (peak, then flops) optimum on a frontier.
template <typename FP>
int pareto_best(container::vector<FP> const& f) {
  int best = 0;
  for (int i = 1; i < static_cast<int>(f.size()); ++i)
    if (f[i].peak < f[best].peak ||
        (f[i].peak == f[best].peak && f[i].flops < f[best].flops))
      best = i;
  return best;
}

/// \brief Per-contraction roofline secondary cost (tie-break wall-time proxy).
///
/// Returns \c max(flops, beta * Q), with data movement
/// \c Q = max(traffic, kappa * flops / sqrt(M / c0)) combining compulsory
/// single-pass traffic with the finite-cache (Hong-Kung) re-read bound. With
/// \c beta (machine_balance) <= 0 this is exactly \c flops (pure-flop
/// tie-break, no behavior change). \c traffic is the operand+result footprint
/// (elements), \c M is fast_mem_elems, \c c0 is block_tiles, \c kappa is
/// block_prefactor. See doc/dev/specs/2026-06-23-roofline-tiebreak-cost.md.
inline double roofline_op_cost(double flops, double traffic,
                               double machine_balance, double fast_mem_elems,
                               double block_tiles,
                               double block_prefactor) noexcept {
  if (machine_balance <= 0.0) return flops;
  double Q = traffic;
  if (fast_mem_elems > 0.0 && block_tiles > 0.0)
    Q = std::max(
        Q, block_prefactor * flops / std::sqrt(fast_mem_elems / block_tiles));
  return std::max(flops, machine_balance * Q);
}

/// \brief Peak-memory single-term cost model (DensePeakSize objective).
///
/// Implements the all-co-resident pebble-game DP, factored into the CostModel
/// hooks driven by \ref run_single_term_opt. The recurrence minimizes peak
/// memory: while one child evaluates, the other child's leaf inputs sit
/// resident. A Pareto frontier of (peak, flops) per subset lets the
/// lexicographic (peak, then flops) optimum be reached. No CSE.
///
/// \tparam IdxToSz A callable mapping an Index to its extent.
template <typename IdxToSz>
struct PeakModel {
  IdxToSz idxsz;
  /// Optional k-aware inner (CSV/PNO composite) extent; see footprint_counter.
  std::function<double(Index const&, std::size_t)> inner_pow = {};
  /// Predicate marking a leaf tensor as volatile (amplitude-dependent). Used
  /// ONLY to weight the secondary flop tie-break: a volatile contraction is
  /// replayed every iteration, so its flops are scaled by \c volatile_weight.
  std::function<bool(Tensor const&)> is_volatile_leaf = {};
  /// Replay weight applied to volatile contractions in the flop tie-break.
  double volatile_weight = 1.0;
  /// Roofline parameters for the secondary (tie-break) cost; see
  /// \ref RooflineParams and \ref roofline_op_cost. machine_balance == 0
  /// (default) => pure-flop tie-break (no behavior change).
  double machine_balance = 0.0;
  double fast_mem_elems = 0.0;
  double block_tiles = 3.0;
  double block_prefactor = 1.0;
  /// Relative peak tolerance for the final (root) selection: among frontier
  /// points whose peak is within (1 + peak_flops_tolerance) of the minimum
  /// peak, pick the one with the fewest flops. 0 (default) = strict peak-min
  /// (exact-tie flop tie-break only). A small positive value trades a bounded
  /// peak increase for a potentially large flop reduction (e.g. forming a
  /// persistent 4-PNO integral instead of recomputing a ladder).
  double peak_flops_tolerance = 0.0;
  /// Perf-first / peak-second selection: when true, `reconstruct` selects the
  /// root-frontier point by (flops, then peak) instead of (peak, then flops),
  /// bypassing `peak_flops_tolerance`. Default false = peak-first (unchanged).
  bool perf_first = false;
  /// Prune disconnected (outer-product) subsets from the DP (see
  /// OptimizeOptions::prune_outer_products). Default true.
  bool prune_outer_products = true;

  /// One non-dominated (peak, flops) trade-off for a subset, with the
  /// bipartition / order / child-frontier-indices needed to reconstruct it.
  /// \c lp_idx / \c rp_idx select which frontier point of each child was used.
  struct FrontPoint {
    double peak = std::numeric_limits<double>::max();
    double flops = std::numeric_limits<double>::max();
    size_t lp = 0;
    size_t rp = 0;
    bool lp_first = true;
    int lp_idx = -1;
    int rp_idx = -1;
  };

  /// Per-subset DP cell: the PARETO FRONTIER of non-dominated (peak, flops)
  /// trade-offs for building the subtree rooted at this subset. A pure
  /// peak-min DP is degenerate when the peak is set by an unavoidable
  /// leaf/intermediate (e.g. a DF integral) that dominates many factorizations:
  /// a single peak-min cell with a *local* flop tie-break does not give the
  /// global flop-min among peak-optimal schedules (the max-recurrence lacks
  /// optimal substructure for the secondary objective). Carrying the frontier
  /// lets a parent combine a child's slightly-higher-peak/lower-flops point
  /// when that peak is hidden under the parent's peak-determining term, so the
  /// final lexicographic (peak, then flops) optimum is reachable -- e.g. it can
  /// form a persistent 4-PNO integral (cheap flops) instead of recomputing the
  /// whole ladder, when both are peak-equal.
  using State = container::vector<FrontPoint>;

  /// Precomputed tables: S[n] = footprint of subset n's result tensor,
  /// L[n] = sum of leaf (singleton) sizes in subset n, idx[n] = subset n's
  /// open (result) indices (for the per-contraction flop tie-break), and
  /// flops_of(lhs, rhs, result) the flop count of one binary contraction.
  struct Context {
    container::vector<double> S;
    container::vector<double> L;
    container::vector<IndexSet> idx;
    std::function<double(IndexSet const&, IndexSet const&, IndexSet const&)>
        flops_of;
    /// Bitmask of volatile leaf tensors (for the flop tie-break weight).
    std::size_t volatile_mask = 0;
  };

  template <typename TIdxs>
  Context build_context(TensorNetwork const& network,
                        TIdxs const& tidxs) const {
    // CSE is not supported for DensePeakSize.
    Context ctx;
    auto const nt = network.tensors().size();
    auto const sz = size_t{1} << nt;
    container::vector<OptRes> results(sz);
    // Outer-product pruning: skip building tables for disconnected subsets the
    // DP will never form (solve_single_term also skips them).
    auto const connected =
        outer_product_connectivity(network, tidxs, prune_outer_products);
    init_results(network, tidxs, results, &connected);
    auto fp = footprint_counter(idxsz, inner_pow);
    ctx.S.assign(sz, 0.0);
    ctx.idx.resize(sz);
    for (size_t n = 0; n < sz; ++n) {
      ctx.S[n] = (n == 0 || !connected[n]) ? 0.0 : fp(results[n].indices);
      ctx.idx[n] = std::move(results[n].indices);
    }
    ctx.L.assign(sz, 0.0);
    for (size_t n = 0; n < sz; ++n)
      for (size_t b = 0; b < nt; ++b)
        if (n & (size_t{1} << b)) ctx.L[n] += ctx.S[size_t{1} << b];
    ctx.flops_of = [c = flops_counter(idxsz, inner_pow)](
                       IndexSet const& lhs, IndexSet const& rhs,
                       IndexSet const& result) { return c(lhs, rhs, result); };
    ctx.volatile_mask = leaf_volatile_mask(network, is_volatile_leaf);
    return ctx;
  }

  State leaf(Context const& ctx, size_t n) const {
    return State{
        FrontPoint{ctx.S[n], 0.0, 0, 0, true, -1, -1}};  // load: 0 flops
  }

  State init(Context const& /*ctx*/, size_t /*n*/) const {
    return State{};  // empty frontier; relax fills it
  }

  void relax(Context& ctx, size_t n, size_t lp, size_t rp, State const& lp_st,
             State const& rp_st, State& acc) const {
    double const both = ctx.S[lp] + ctx.S[rp] + ctx.S[lp | rp];
    // A volatile contraction (subset n contains a volatile leaf) is replayed
    // every iteration, so its flops are scaled by volatile_weight -- matching
    // the DenseFLOPs convention. The contraction flops are order-independent.
    double const w = (ctx.volatile_mask & n) ? volatile_weight : 1.0;
    // Secondary (tie-break) cost: roofline wall-time proxy per replay, charged
    // w times for volatile (replayed) contractions. traffic = operand+result
    // footprint. With machine_balance==0 this reduces to flops (no change).
    double const cflops =
        w * roofline_op_cost(ctx.flops_of(ctx.idx[lp], ctx.idx[rp], ctx.idx[n]),
                             ctx.S[lp] + ctx.S[rp] + ctx.S[n], machine_balance,
                             fast_mem_elems, block_tiles, block_prefactor);
    // Cross every (peak,flops) trade-off of the two children.
    for (int li = 0; li < static_cast<int>(lp_st.size()); ++li)
      for (int ri = 0; ri < static_cast<int>(rp_st.size()); ++ri) {
        double const pL = lp_st[li].peak, pR = rp_st[ri].peak;
        double const lp_first_cand =
            std::max({ctx.L[rp] + pL, ctx.S[lp] + pR, both});
        double const rp_first_cand =
            std::max({ctx.L[lp] + pR, ctx.S[rp] + pL, both});
        pareto_insert(
            acc, FrontPoint{std::min(lp_first_cand, rp_first_cand),
                            lp_st[li].flops + rp_st[ri].flops + cflops, lp, rp,
                            lp_first_cand <= rp_first_cand, li, ri});
      }
  }

  void finalize(Context& /*ctx*/, size_t /*n*/,
                container::vector<State>& /*st*/) const {}

  EvalSequence reconstruct(Context const& /*ctx*/,
                           container::vector<State> const& st) const {
    size_t const full = st.size() - 1;
    if (perf_first) {
      // Perf-first / peak-second (non-batched): min flops, ties by lower peak,
      // bypassing the peak_flops_tolerance epsilon band (a peak-first knob).
      auto const& rootf = st[full];
      int pbest = 0;
      for (int i = 1; i < static_cast<int>(rootf.size()); ++i)
        if (rootf[i].flops < rootf[pbest].flops ||
            (rootf[i].flops == rootf[pbest].flops &&
             rootf[i].peak < rootf[pbest].peak))
          pbest = i;
      // Reuse the existing back-pointer walk with the chosen root index.
      std::function<EvalSequence(size_t, int)> pbuild =
          [&](size_t n, int idx) -> EvalSequence {
        if (std::popcount(n) == 1)
          return EvalSequence{static_cast<int>(std::countr_zero(n))};
        FrontPoint const& fp = st[n][idx];
        size_t const fs = fp.lp_first ? fp.lp : fp.rp;
        int const fi = fp.lp_first ? fp.lp_idx : fp.rp_idx;
        size_t const ss = fp.lp_first ? fp.rp : fp.lp;
        int const si = fp.lp_first ? fp.rp_idx : fp.lp_idx;
        EvalSequence s = pbuild(fs, fi);
        EvalSequence b = pbuild(ss, si);
        s.insert(s.end(), b.begin(), b.end());
        s.push_back(-1);
        return s;
      };
      return pbuild(full, pbest);
    }
    // ε-tolerant selection: among frontier points within
    // (1 + peak_flops_tolerance) of the minimum peak, take the fewest flops
    // (ties broken by lower peak). tolerance == 0 recovers strict peak-min.
    auto const& root = st[full];
    double minpeak = std::numeric_limits<double>::max();
    for (auto const& fp : root) minpeak = std::min(minpeak, fp.peak);
    double const thresh = minpeak * (1.0 + peak_flops_tolerance);
    int best = -1;
    for (int i = 0; i < static_cast<int>(root.size()); ++i)
      if (root[i].peak <= thresh &&
          (best < 0 || root[i].flops < root[best].flops ||
           (root[i].flops == root[best].flops &&
            root[i].peak < root[best].peak)))
        best = i;
    // Follow back-pointers (which child + which child frontier point).
    std::function<EvalSequence(size_t, int)> build =
        [&](size_t n, int idx) -> EvalSequence {
      if (std::popcount(n) == 1)
        return EvalSequence{static_cast<int>(std::countr_zero(n))};
      FrontPoint const& fp = st[n][idx];
      size_t const fs = fp.lp_first ? fp.lp : fp.rp;
      int const fi = fp.lp_first ? fp.lp_idx : fp.rp_idx;
      size_t const ss = fp.lp_first ? fp.rp : fp.lp;
      int const si = fp.lp_first ? fp.rp_idx : fp.lp_idx;
      EvalSequence s = build(fs, fi);
      EvalSequence b = build(ss, si);
      s.insert(s.end(), b.begin(), b.end());
      s.push_back(-1);
      return s;
    };
    return build(full, best);
  }
};

/// \brief Multi-mode batched peak-memory single-term cost model
/// (DensePeakSizeBatched objective).
///
/// Implements the per-batchable-index all-co-resident pebble-game DP, factored
/// into the CostModel hooks driven by \ref run_single_term_opt. Follows the
/// batch-aware cost model design (section 6.2, model A): each DP cell is
/// indexed by both a subset \c n and a sliced-set context \c B over the
/// batchable indices, so a model State is the \c [B]-vector of per-context
/// \ref BatchedRes. No CSE; persistence-gated batching.
///
/// \tparam IdxToSz A callable mapping an Index to its extent.
template <typename IdxToSz>
struct PeakBatchedModel {
  IdxToSz idxsz;
  std::function<bool(Index const&)> is_batchable;
  std::function<std::size_t(Index const&)> batch;
  std::function<bool(Tensor const&)> is_volatile_leaf;
  /// Optional k-aware inner (CSV/PNO composite) extent; see footprint_counter.
  std::function<double(Index const&, std::size_t)> inner_pow = {};
  /// Replay weight applied to volatile contractions in the flop tie-break.
  double volatile_weight = 1.0;
  /// Roofline parameters for the secondary (tie-break) cost; see
  /// \ref RooflineParams and \ref roofline_op_cost. machine_balance == 0
  /// (default) => pure-flop tie-break (no behavior change). Uses full
  /// (unsliced) operand+result footprints (total per-replay traffic; slicing
  /// reduces peak, not total work).
  double machine_balance = 0.0;
  double fast_mem_elems = 0.0;
  double block_tiles = 3.0;
  double block_prefactor = 1.0;
  /// If true, batch only persistent subtrees (decline any subset containing a
  /// volatile leaf). Default false = batch across the board. See
  /// BatchPolicy::persistent_only.
  bool batch_persistent_only = false;
  /// Unused by \ref reconstruct (superseded by the threshold-gated selection
  /// below, driven by \ref peak_threshold / \ref numeric_size); retained for
  /// source compatibility. See PeakModel::peak_flops_tolerance, which is
  /// still consulted by the (unbatched) DensePeakSize model.
  double peak_flops_tolerance = 0.0;
  /// In-flight batch-contribution footprint multiplier; see
  /// BatchPolicy::accumulation_factor. Charged only on nodes that contract a
  /// batchable index (Ap != 0), into the all-co-resident peak term, to price
  /// the accumulator + contribution co-residency of K += contribution.
  double accumulation_factor = 0.0;
  /// Peak-memory budget (BYTES) for threshold-gated selection; see
  /// BatchPolicy::peak_threshold. +infinity (default) => min-flops (no
  /// batching).
  double peak_threshold = std::numeric_limits<double>::infinity();
  /// Bytes per stored element, to compare the model's element-count peak to
  /// peak_threshold (bytes). Default 8 (double / TensorD).
  double numeric_size = 8.0;
  /// Perf-first / peak-second selection: when true, `select_root` selects the
  /// root-frontier point by (flops, then peak) and does NOT consult
  /// `peak_threshold` as a feasibility gate (it can no longer force a
  /// FLOPS-catastrophic factorization for its sliceability). Default false =
  /// peak-first threshold-gated selection (unchanged).
  bool perf_first = false;

  /// If true (the default), charge the batch RECOMPUTATION cost on the flops/
  /// exec axis. The batched evaluator re-executes each contraction per tile of
  /// the ancestor batch axes its result does NOT carry (across-batch work is
  /// recomputed; within-batch sharing is cached -- see eval.hpp "replays the
  /// build of every compatible persistent final"). A node at ancestor-sliced-
  /// set B is charged nbatch(b) for each b in B not open in the node, so a
  /// schedule that slices many axes it must recompute across pays for it. The
  /// alternative (false) assumes WORK PARITY (batching is free on flops), which
  /// under-costs heavily-sliced families and does not reflect the true cost of
  /// batching; kept only as an escape hatch for comparison.
  bool charge_batch_recompute = true;
  /// Prune disconnected (outer-product) subsets from the DP (see
  /// OptimizeOptions::prune_outer_products). Default true.
  bool prune_outer_products = true;

  /// One non-dominated (peak, flops) trade-off for a (subset, sliced-set \c B)
  /// cell. \c aprime is the sliced-set chosen at this node; the children are
  /// read at context \c C = B | aprime, at frontier indices \c lp_idx /
  /// \c rp_idx. See \ref PeakModel::FrontPoint for why a frontier (not a single
  /// peak-min cell) is needed.
  struct BFrontPoint {
    double peak = std::numeric_limits<double>::max();
    double flops = std::numeric_limits<double>::max();
    size_t lp = 0;
    size_t rp = 0;
    bool lp_first = true;
    std::size_t aprime = 0;
    int lp_idx = -1;
    int rp_idx = -1;
  };

  /// Per-subset DP cell: a \c [B]-vector (size \c nB = 2^m) of Pareto
  /// frontiers, one non-dominated (peak, flops) set per sliced-set context B.
  using State = container::vector<container::vector<BFrontPoint>>;

  /// Precomputed tables and per-(subset, sliced-set) lookup parameters built
  /// once by build_context.
  struct Context {
    /// Ordered, deduplicated batchable indices (bit \c k maps to \c aux[k]).
    container::vector<Index> aux;
    /// Number of batchable indices (= aux.size()).
    std::size_t m = 0;
    /// Number of sliced-sets (= 2^m).
    std::size_t nB = 1;
    /// Number of tensors in the network.
    std::size_t nt = 0;
    /// tables[B][n] = footprint of subset n under sliced-set B.
    container::vector<container::vector<double>> tables;
    /// open_aux[n] = bitmask of batchable indices open in subset n.
    container::vector<std::size_t> open_aux;
    /// Bitmask of volatile leaf tensors.
    std::size_t volatile_mask = 0;
    /// idx[n] = subset n's open (result) indices, for the flop tie-break.
    container::vector<IndexSet> idx;
    /// flops_of(lhs, rhs, result) = flop count of one binary contraction.
    /// Retained as the reference for the fast_flops parity test; the relax
    /// hot loop uses fast_flops (see below).
    std::function<double(IndexSet const&, IndexSet const&, IndexSet const&)>
        flops_of;
    /// nbatch[k] = number of batch tiles of aux[k] = ceil(extent / target),
    /// clamped to >= 1. Used to charge batch recomputation (see
    /// charge_batch_recompute): a node inside an ancestor batch loop over
    /// aux[k] that does not carry aux[k] is re-executed nbatch[k] times.
    container::vector<double> nbatch;

    // --- fast per-subset flop precompute (relax tie-break hot path) ---
    // Per subset, sorted (FullLabelCompare-ordered) atom IDs: outer atoms are
    // the plain indices + the proto-indices of composites; inner atoms are the
    // composites themselves -- exactly tot_indices()'s outer/inner split. IDs
    // are assigned in FullLabelCompare order so the union merge multiplies in
    // the SAME order as inner_aware_volume, giving a bit-identical result.
    // fast_flops(lp, rp) == flops_of(idx[lp], idx[rp], idx[lp|rp]) by
    // construction (gated by the [fast-flops-parity] test). use_fast_flops is
    // false only in that test, to exercise the flops_of reference path.
    bool use_fast_flops = true;
    container::vector<container::svector<std::uint32_t>> f_outer;
    container::vector<container::svector<std::uint32_t>> f_inner;
    container::vector<double> fo_ext;         // outer atom -> extent
    container::vector<double> fi_ext;         // inner atom -> extent (fallback)
    container::vector<Index> fi_index;        // inner atom -> composite Index
    container::vector<std::uint32_t> fi_grp;  // inner atom -> proto-group id
    bool f_inner_engaged = false;
    std::function<double(Index const&, std::size_t)> f_inner_pow;
    mutable container::vector<std::uint32_t> f_grpcount;   // scratch (ngrp)
    mutable container::svector<std::uint32_t> f_union_in;  // scratch

    /// Flop count of the binary contraction (subset \c lp) x (subset \c rp),
    /// equivalent to flops_of(idx[lp], idx[rp], idx[n]) with n = lp|rp, but
    /// reading precomputed atom-ID lists instead of rebuilding index sets per
    /// call. (idx[n] is a subset of idx[lp] U idx[rp], so it adds nothing.)
    double fast_flops(std::size_t lp, std::size_t rp) const {
      double mem = 1.0;
      {
        auto const& A = f_outer[lp];
        auto const& B = f_outer[rp];
        std::size_t i = 0, j = 0;
        while (i < A.size() && j < B.size()) {
          if (A[i] < B[j])
            mem *= fo_ext[A[i++]];
          else if (B[j] < A[i])
            mem *= fo_ext[B[j++]];
          else {
            mem *= fo_ext[A[i]];
            ++i;
            ++j;
          }
        }
        while (i < A.size()) mem *= fo_ext[A[i++]];
        while (j < B.size()) mem *= fo_ext[B[j++]];
      }
      auto const& AI = f_inner[lp];
      auto const& BI = f_inner[rp];
      if (!AI.empty() || !BI.empty()) {
        f_union_in.clear();
        std::size_t i = 0, j = 0;
        while (i < AI.size() && j < BI.size()) {
          std::uint32_t id;
          if (AI[i] < BI[j])
            id = AI[i++];
          else if (BI[j] < AI[i])
            id = BI[j++];
          else {
            id = AI[i];
            ++i;
            ++j;
          }
          f_union_in.push_back(id);
        }
        while (i < AI.size()) f_union_in.push_back(AI[i++]);
        while (j < BI.size()) f_union_in.push_back(BI[j++]);
        if (f_inner_engaged) {
          for (auto id : f_union_in) ++f_grpcount[fi_grp[id]];
          for (auto id : f_union_in)
            mem *= f_inner_pow(fi_index[id], f_grpcount[fi_grp[id]]);
          for (auto id : f_union_in) f_grpcount[fi_grp[id]] = 0;
        } else {
          for (auto id : f_union_in) mem *= fi_ext[id];
        }
      }
      return mem == 1.0 ? 0.0 : mem;
    }

    /// Context-restricted size of subset s under sliced-set ctx (the table is
    /// indexed by the part of ctx actually open in s; mirrors the oracle).
    double sz(std::size_t s, std::size_t ctx) const {
      return tables[ctx & open_aux[s]][s];
    }
    /// Per-context leaf-sum of subset s (sum of singleton sizes under ctx).
    double Lof(std::size_t s, std::size_t ctx) const {
      double r = 0.0;
      for (std::size_t b = 0; b < nt; ++b)
        if (s & (std::size_t{1} << b)) r += sz(std::size_t{1} << b, ctx);
      return r;
    }
  };

  template <typename TIdxs>
  Context build_context(TensorNetwork const& network,
                        TIdxs const& tidxs) const {
    // CSE is not supported for DensePeakSizeBatched.
    Context ctx;
    ctx.nt = network.tensors().size();
    ctx.aux = batchable_index_list(network, is_batchable);
    ctx.m = ctx.aux.size();
    // Per-axis batch-tile count for the recompute charge: ceil(extent/target),
    // >= 1. batch() returns the target tile size; a 0/absent target => 1 tile
    // (no batching of that axis => no recompute).
    ctx.nbatch.assign(ctx.m, 1.0);
    for (std::size_t k = 0; k < ctx.m; ++k) {
      double const ext = static_cast<double>(idxsz(ctx.aux[k]));
      double const tgt = static_cast<double>(batch(ctx.aux[k]));
      ctx.nbatch[k] = (tgt > 0.0) ? std::max(1.0, std::ceil(ext / tgt)) : 1.0;
    }
    // accumulation_factor is charged per accumulation node (Ap != 0) and is
    // valid for any number of batchable indices: with nested accumulation the
    // per-node charges co-exist at the peak. Validated by the identity
    // peak_cost_batched == reconstructed_batched_peak (test [batched-accum]).
    ctx.nB = std::size_t{1} << ctx.m;
    // Outer-product pruning: skip building tables for disconnected subsets the
    // DP will never form (solve_single_term also skips them). connected[n]==1
    // for singletons/empty and for connected subsets; the (~2x connected)
    // needed-mask is derived internally where a complement lookup requires it.
    auto const connected =
        outer_product_connectivity(network, tidxs, prune_outer_products);
    ctx.tables = sliced_footprints(network, tidxs, idxsz, is_batchable, batch,
                                   ctx.aux, inner_pow, &connected);
    // open_aux is NOT pruned: is_spectator_axis scans open_aux over the FULL
    // subset lattice (including disconnected subsets) to verify an axis is
    // never contracted, so every subset's open-axis bitmask must be real.
    ctx.open_aux = subset_open_aux(network, tidxs, ctx.aux);
    ctx.volatile_mask = leaf_volatile_mask(network, is_volatile_leaf);
    // Per-subset open indices + a flop counter, for the lexicographic
    // (peak, then flops) tie-break (mirrors PeakModel). The flop tie-break uses
    // the unbatched contraction flops; total work summed over batches matches
    // it (work parity), so it consistently orders equal-peak schedules.
    auto const sz = std::size_t{1} << ctx.nt;
    container::vector<OptRes> results(sz);
    init_results(network, tidxs, results, &connected);
    ctx.idx.resize(sz);
    for (std::size_t n = 0; n < sz; ++n)
      ctx.idx[n] = std::move(results[n].indices);
    ctx.flops_of = [c = flops_counter(idxsz, inner_pow)](
                       IndexSet const& lhs, IndexSet const& rhs,
                       IndexSet const& result) { return c(lhs, rhs, result); };

    // Fast-flops atom tables (see Context::fast_flops). Atom IDs are assigned
    // in FullLabelCompare order so the union-merge multiply order matches
    // inner_aware_volume exactly (bit-identical flops). Disconnected subsets
    // have empty ctx.idx (pruned) and contribute no atoms; fast_flops is only
    // ever called for connected bipartitions.
    {
      std::set<Index, Index::FullLabelCompare> outer_atoms, inner_atoms;
      for (std::size_t n = 0; n < sz; ++n)
        for (auto const& ix : ctx.idx[n]) {
          if (ix.has_proto_indices()) {
            inner_atoms.insert(ix);
            for (auto const& p : ix.proto_indices()) outer_atoms.insert(p);
          } else {
            outer_atoms.insert(ix);
          }
        }
      std::map<Index, std::uint32_t, Index::FullLabelCompare> oid, iid;
      ctx.fo_ext.reserve(outer_atoms.size());
      for (auto const& a : outer_atoms) {
        oid.emplace(a, static_cast<std::uint32_t>(ctx.fo_ext.size()));
        ctx.fo_ext.push_back(static_cast<double>(idxsz(a)));
      }
      std::map<std::wstring, std::uint32_t> gid;  // proto-signature -> group id
      ctx.fi_ext.reserve(inner_atoms.size());
      for (auto const& c : inner_atoms) {
        iid.emplace(c, static_cast<std::uint32_t>(ctx.fi_index.size()));
        ctx.fi_index.push_back(c);
        ctx.fi_ext.push_back(static_cast<double>(idxsz(c)));
        std::wstring sig;
        for (auto const& p : c.proto_indices()) {
          sig += p.full_label();
          sig += L',';
        }
        auto it = gid.find(sig);
        ctx.fi_grp.push_back(
            it != gid.end()
                ? it->second
                : gid.emplace(sig, static_cast<std::uint32_t>(gid.size()))
                      .first->second);
      }
      ctx.f_grpcount.assign(gid.size(), 0);
      ctx.f_outer.resize(sz);
      ctx.f_inner.resize(sz);
      for (std::size_t n = 0; n < sz; ++n) {
        auto& oo = ctx.f_outer[n];
        auto& ii = ctx.f_inner[n];
        for (auto const& ix : ctx.idx[n]) {
          if (ix.has_proto_indices()) {
            ii.push_back(iid.at(ix));
            for (auto const& p : ix.proto_indices()) oo.push_back(oid.at(p));
          } else {
            oo.push_back(oid.at(ix));
          }
        }
        std::sort(oo.begin(), oo.end());
        oo.erase(std::unique(oo.begin(), oo.end()), oo.end());
        std::sort(ii.begin(), ii.end());
        ii.erase(std::unique(ii.begin(), ii.end()), ii.end());
      }
      ctx.f_inner_engaged = option_engaged(inner_pow);
      ctx.f_inner_pow = inner_pow;
    }
    return ctx;
  }

  State leaf(Context const& ctx, size_t n) const {
    State s(ctx.nB);
    for (std::size_t B = 0; B < ctx.nB; ++B)
      s[B].push_back(BFrontPoint{ctx.sz(n, B), 0.0, 0, 0, true, 0, -1, -1});
    return s;
  }

  State init(Context const& ctx, size_t /*n*/) const {
    return State(ctx.nB);  // nB empty frontiers; relax fills them
  }

  void relax(Context& ctx, size_t n, size_t lp, size_t rp, State const& lp_st,
             State const& rp_st, State& acc) const {
    // Secondary (tie-break) cost: roofline wall-time proxy per replay, charged
    // volatile_weight times for volatile (replayed) contractions. Uses the full
    // (unsliced) operand+result footprint as the per-replay traffic; slicing
    // reduces peak (primary axis), not total work. machine_balance==0 => flops.
    double const w = (ctx.volatile_mask & n) ? volatile_weight : 1.0;
    double const cflops =
        w * roofline_op_cost(
                ctx.use_fast_flops
                    ? ctx.fast_flops(lp, rp)
                    : ctx.flops_of(ctx.idx[lp], ctx.idx[rp], ctx.idx[n]),
                ctx.sz(lp, 0) + ctx.sz(rp, 0) + ctx.sz(n, 0), machine_balance,
                fast_mem_elems, block_tiles, block_prefactor);
    for (std::size_t B = 0; B < ctx.nB; ++B) {
      // Batch recomputation charge: this node sits inside the ancestor batch
      // loops over the axes in B. For each b in B whose axis this node's result
      // does NOT carry (b not open in n), the node is re-executed nbatch(b)
      // times (across-batch recompute); axes it carries are partitioned (x1).
      // Default off => rf==1 => historical work-parity cost.
      double rf = 1.0;
      if (charge_batch_recompute) {
        std::size_t const esc =
            B & ~ctx.open_aux[n];  // B axes n does not carry
        for (std::size_t k = 0; k < ctx.m; ++k)
          if (esc & (std::size_t{1} << k)) rf *= ctx.nbatch[k];
      }
      double const cflops_B = cflops * rf;
      // Batchable indices contracted at THIS node: open at children but not at
      // the parent. By default batching is applied ACROSS THE BOARD: slicing
      // the batch axis shrinks any intermediate carrying it regardless of
      // volatility (footprint objective) while leaving flops unchanged, so the
      // persistence gate would only ever raise the modelled peak. Set
      // batch_persistent_only to restore the persistent-only gate (decline to
      // slice subsets that contain a volatile leaf).
      std::size_t const Acand =
          (batch_persistent_only && (ctx.volatile_mask & n))
              ? std::size_t{0}
              : ((ctx.open_aux[lp] | ctx.open_aux[rp]) & ~ctx.open_aux[n]);
      // Enumerate every subset A' of Acand (including the empty set).
      std::size_t Ap = Acand;
      while (true) {
        std::size_t const C = B | Ap;
        double const szlp = ctx.sz(lp, C), szrp = ctx.sz(rp, C),
                     szn = ctx.sz(n, B);
        // A node that contracts a batchable index (Ap != 0) is accumulated over
        // the aux batches (K += contribution); the in-flight contribution (same
        // index set as the result, size szn) co-resides with the accumulator.
        // Charge it once, on the all-co-resident moment only -- the pre-result
        // staged terms (Lrp+pl, szlp+prr) exclude it since szn is not yet
        // built.
        double const contrib = (Ap != 0) ? accumulation_factor * szn : 0.0;
        double const both = szlp + szrp + szn + contrib;
        double const Lrp = ctx.Lof(rp, C), Llp = ctx.Lof(lp, C);
        // Cross every (peak,flops) trade-off of the two children at context C.
        for (int li = 0; li < static_cast<int>(lp_st[C].size()); ++li)
          for (int ri = 0; ri < static_cast<int>(rp_st[C].size()); ++ri) {
            double const pl = lp_st[C][li].peak, prr = rp_st[C][ri].peak;
            double const lpf = std::max({Lrp + pl, szlp + prr, both});
            double const rpf = std::max({Llp + prr, szrp + pl, both});
            pareto_insert(acc[B], BFrontPoint{std::min(lpf, rpf),
                                              lp_st[C][li].flops +
                                                  rp_st[C][ri].flops + cflops_B,
                                              lp, rp, lpf <= rpf, Ap, li, ri});
          }
        if (Ap == 0) break;
        Ap = (Ap - 1) & Acand;
      }
    }
  }

  void finalize(Context& /*ctx*/, size_t /*n*/,
                container::vector<State>& /*st*/) const {}

  /// \brief True iff batchable axis bit \p k is a genuine external spectator of
  /// this network: it is OPEN on the root result AND is contracted at NO node.
  ///
  /// A spectator is carried unchanged from the leaves up to the root: whenever
  /// any tensor in a subset carries the axis, the axis stays OPEN in that
  /// subset, so it never appears in any node's contracted-at-node set
  /// (\c Acand = (open_aux[lp]|open_aux[rp]) & ~open_aux[n]) and slicing it is
  /// purely a footprint change with identical work. The forest-batching seed
  /// path asserts this before honoring a seed axis: seeding an axis that is
  /// actually contracted somewhere would mis-size the nodes that contract it.
  bool is_spectator_axis(Context const& ctx, std::size_t k) const {
    std::size_t const root = (std::size_t{1} << ctx.nt) - 1;
    // Must be carried on the root result (a genuine external index).
    if (!((ctx.open_aux[root] >> k) & 1u)) return false;
    // Leaves (single-tensor subsets) that carry axis k.
    std::size_t leafmask = 0;
    for (std::size_t b = 0; b < ctx.nt; ++b)
      if ((ctx.open_aux[std::size_t{1} << b] >> k) & 1u)
        leafmask |= (std::size_t{1} << b);
    // Whenever the axis is AVAILABLE in a subset (some carrying leaf is in it)
    // it must be OPEN in that subset -- else it is contracted at the node that
    // forms that subset, i.e. not a pure spectator.
    for (std::size_t n = 1; n <= root; ++n)
      if ((leafmask & n) && !((ctx.open_aux[n] >> k) & 1u)) return false;
    return true;
  }

  /// Threshold-gated root-frontier selection shared by \ref reconstruct and
  /// \ref reconstruct_axes: among points whose peak (bytes) fits
  /// peak_threshold, pick fewest flops (ties by lower peak). If none fit, pick
  /// min peak (best effort). peak_threshold == +inf => all feasible => min
  /// flops => the non-batched schedule. Returns the chosen index into
  /// \c st[root][root_B].
  ///
  /// \param root_B The batch context read at the ROOT. The default 0 (the
  ///        empty sliced-set) is the historical behavior, kept byte-identical.
  ///        The forest-batching seed path passes a non-zero mask to SEED a
  ///        spectator external axis into the root frontier, so the whole tree
  ///        is sized (and its peak reported) with that axis sliced --
  ///        \c st[root][root_B] is already computed by the DP's B-loop; only
  ///        which B the root selection reads changes.
  int select_root(Context const& ctx, container::vector<State> const& st,
                  std::size_t root_B = 0) const {
    std::size_t const root = (std::size_t{1} << ctx.nt) - 1;
    auto const& rootf = st[root][root_B];
    auto peak_bytes = [this](double peak_elems) {
      return peak_elems * numeric_size;
    };
    if (perf_first) {
      // Perf-first / peak-second: min flops, ties by lower peak. The frontier
      // keeps one min-peak point per distinct flops value (pareto_insert prunes
      // equal-flops higher-peak points), so this both picks the cheapest
      // factorization and takes its fully-sliced (min-peak) realization.
      // peak_threshold is deliberately NOT consulted here.
      int pbest = 0;
      for (int i = 1; i < static_cast<int>(rootf.size()); ++i)
        if (rootf[i].flops < rootf[pbest].flops ||
            (rootf[i].flops == rootf[pbest].flops &&
             rootf[i].peak < rootf[pbest].peak))
          pbest = i;
      return pbest;
    }
    int best = -1;
    bool any_feasible = false;
    for (int i = 0; i < static_cast<int>(rootf.size()); ++i)
      if (peak_bytes(rootf[i].peak) <= peak_threshold) {
        any_feasible = true;
        if (best < 0 || rootf[i].flops < rootf[best].flops ||
            (rootf[i].flops == rootf[best].flops &&
             rootf[i].peak < rootf[best].peak))
          best = i;
      }
    if (!any_feasible) {
      // Infeasible: no schedule fits the budget. Fall back to min peak.
      double minpeak = std::numeric_limits<double>::max();
      for (int i = 0; i < static_cast<int>(rootf.size()); ++i)
        if (rootf[i].peak < minpeak) {
          minpeak = rootf[i].peak;
          best = i;
        }
    }
    // DIAGNOSTIC (gated): the DP's OWN chosen frontier cost (.flops is the
    // roofline exec cost, volatile-weighted -- the true objective), plus the
    // global-min .flops over the WHOLE frontier (feasible or not) so a
    // non-min-feasible selection is visible. Sum the printed chosen_flops
    // across terms to test threshold monotonicity of the DP's real objective.
    if (best >= 0 && std::getenv("SEQUANT_SELROOT_DEBUG")) {
      double gmin = std::numeric_limits<double>::max();
      for (auto const& p : rootf) gmin = std::min(gmin, p.flops);
      std::cerr << "[selroot] chosen_flops=" << rootf[best].flops
                << " chosen_peak_gb=" << (peak_bytes(rootf[best].peak) / 1e9)
                << " feasible=" << (any_feasible ? 1 : 0)
                << " nfront=" << rootf.size() << " global_min_flops=" << gmin
                << "\n";
    }
    return best;
  }

  EvalSequence reconstruct(Context const& ctx,
                           container::vector<State> const& st) const {
    std::size_t const root = (std::size_t{1} << ctx.nt) - 1;
    int const best = select_root(ctx, st);
    // Recursive back-pointer walk: at (n, B, idx) read the chosen front point,
    // form child context C = B | aprime, recurse in lp_first order.
    std::function<EvalSequence(std::size_t, std::size_t, int)> build =
        [&](std::size_t n, std::size_t B, int idx) -> EvalSequence {
      if (std::popcount(n) == 1)
        return EvalSequence{static_cast<int>(std::countr_zero(n))};
      BFrontPoint const& r = st[n][B][idx];
      std::size_t const C = B | r.aprime;
      std::size_t const fs = r.lp_first ? r.lp : r.rp;
      int const fi = r.lp_first ? r.lp_idx : r.rp_idx;
      std::size_t const ss = r.lp_first ? r.rp : r.lp;
      int const si = r.lp_first ? r.rp_idx : r.lp_idx;
      EvalSequence s = build(fs, C, fi);
      EvalSequence b = build(ss, C, si);
      s.insert(s.end(), b.begin(), b.end());
      s.push_back(-1);
      return s;
    };
    return build(root, 0, best);
  }

  /// Companion to \ref reconstruct that additionally reports, for each \c -1
  /// (contraction) entry emitted in the returned \c EvalSequence in emission
  /// order, the vector of \c Index sliced at that node (\c ctx.aux[bit] for
  /// each set bit of that node's \c aprime). Leaf entries contribute nothing.
  /// Does not change \ref reconstruct's own output; the two walks are kept in
  /// lock-step so the RPN order and the per-node axes line up.
  std::pair<EvalSequence,
            container::vector<container::svector<std::pair<Index, AxisKind>>>>
  reconstruct_axes(Context const& ctx,
                   container::vector<State> const& st) const {
    std::size_t const root = (std::size_t{1} << ctx.nt) - 1;
    int const best = select_root(ctx, st);  // shared helper (see below)
    container::vector<container::svector<std::pair<Index, AxisKind>>> node_axes;
    std::function<EvalSequence(std::size_t, std::size_t, int)> build =
        [&](std::size_t n, std::size_t B, int idx) -> EvalSequence {
      if (std::popcount(n) == 1)
        return EvalSequence{static_cast<int>(std::countr_zero(n))};
      BFrontPoint const& r = st[n][B][idx];
      std::size_t const C = B | r.aprime;
      std::size_t const fs = r.lp_first ? r.lp : r.rp;
      int const fi = r.lp_first ? r.lp_idx : r.rp_idx;
      std::size_t const ss = r.lp_first ? r.rp : r.lp;
      int const si = r.lp_first ? r.rp_idx : r.lp_idx;
      EvalSequence s = build(fs, C, fi);
      EvalSequence b = build(ss, C, si);
      s.insert(s.end(), b.begin(), b.end());
      container::svector<std::pair<Index, AxisKind>> axes;
      for (std::size_t k = 0; k < ctx.m; ++k)
        if (r.aprime & (std::size_t{1} << k))
          axes.push_back({ctx.aux[k], AxisKind::Contracted});
      node_axes.push_back(std::move(axes));  // one entry per -1, in RPN order
      s.push_back(-1);
      return s;
    };
    auto seq = build(root, 0, best);
    return {std::move(seq), std::move(node_axes)};
  }
};

/// \brief Achieved minimum peak memory (the DensePeakSize objective value) for
/// the whole network under its optimal order. Builds \ref PeakModel, runs the
/// generic driver's \ref solve_single_term, and returns the root subset's peak.
/// Used by tests to compare against the brute-force oracle.
template <typename TIdxs, typename IdxToSz>
double peak_cost(TensorNetwork const& network, TIdxs const& tidxs,
                 IdxToSz&& idxsz) {
  PeakModel<std::decay_t<IdxToSz>> model{std::forward<IdxToSz>(idxsz)};
  auto ctx = model.build_context(network, tidxs);
  auto st = solve_single_term(model, network, tidxs, ctx);
  // root subset == full set == last element; its frontier's smallest peak is
  // the achieved minimum peak memory.
  double mn = std::numeric_limits<double>::max();
  for (auto const& fp : st.back()) mn = std::min(mn, fp.peak);
  return mn;
}

/// \brief Achieved minimum batched peak memory: the peak[root][B=0] objective
/// of the multi-mode batched DP. Builds \ref PeakBatchedModel, runs
/// \ref solve_single_term, and returns the root subset's B=0 peak.
template <typename TIdxs, typename IdxToSz>
double peak_cost_batched(
    TensorNetwork const& network, TIdxs const& tidxs, IdxToSz&& idxsz,
    std::function<bool(Index const&)> const& is_batchable,
    std::function<std::size_t(Index const&)> const& batch_target_size,
    std::function<bool(Tensor const&)> const& is_volatile_leaf,
    double accumulation_factor = 0.0) {
  PeakBatchedModel<std::decay_t<IdxToSz>> model{
      std::forward<IdxToSz>(idxsz),
      is_batchable,
      batch_target_size,
      is_volatile_leaf,
      /* inner_pow */ {},
      /* volatile_weight */ 1.0,
      /* machine_balance */ 0.0,
      /* fast_mem_elems */ 0.0,
      /* block_tiles */ 3.0,
      /* block_prefactor */ 1.0,
      /* batch_persistent_only */ false,
      /* peak_flops_tolerance */ 0.0,
      accumulation_factor};
  auto ctx = model.build_context(network, tidxs);
  auto st = solve_single_term(model, network, tidxs, ctx);
  // root subset's B=0 frontier; its smallest peak is the achieved minimum.
  double mn = std::numeric_limits<double>::max();
  for (auto const& fp : st.back()[0]) mn = std::min(mn, fp.peak);
  return mn;
}

/// \brief Result of seeding a spectator external axis into the batched DP's
/// ROOT frontier (forest batching, P1). The seed axis is contracted at NO node
/// (a pure spectator carried on the root and every subtree holding a composite
/// leg that bears it), so slicing it does NOT change total work (flops) but
/// shrinks every intermediate that carries it -- in particular the perf-first
/// PPL W giant -- by ~occ_block/occ. Carries both the unseeded (B=0) and seeded
/// (B={k_seed}) root selections so a caller can report the drop, plus the
/// chosen axis + block size for Task 4 to thread through the public optimize
/// API.
struct SeededBatchedResult {
  double unseeded_peak_bytes = 0;
  double seeded_peak_bytes = 0;
  double unseeded_flops = 0;
  double seeded_flops = 0;
  std::optional<Index> seeded_axis;
  std::size_t occ_block = 0;
  bool spectator_ok = false;
};

/// \brief Runs the batched peak DP once and reports the ROOT-frontier selection
/// at B=0 (unseeded) and at B={k_seed} (the seed axis sliced into the root
/// batch context), under the model's own (perf-first or peak-first) selector.
///
/// The spectator external axis (e.g. the residual's own external occupied
/// index, carried only as a composite protoindex) is admitted to \c ctx.aux by
/// \ref batchable_index_list's proto pass regardless of the base
/// \c is_batchable predicate; but \ref sliced_footprints only shrinks an axis
/// for which \c is_batchable is true. This routine therefore extends the
/// predicate (and gives the seed axis a block-sized batch target \p occ_block)
/// so the DP tables actually slice the seed axis when its bit is set. Every
/// other axis is untouched, and B=0 never sets the seed bit (the spectator is
/// contracted nowhere, so it never enters any node's \c Acand and hence never
/// enters any child context at B=0), so the unseeded (B=0) selection is
/// byte-identical to \p base_model's.
///
/// \p seed_axis MUST be a genuine external spectator (asserted via
/// \ref PeakBatchedModel::is_spectator_axis): carried on the root and
/// contracted at no node. Seeding a contracted axis would mis-size the nodes
/// that contract it, so this asserts rather than silently mis-report.
template <typename TIdxs, typename IdxToSz>
SeededBatchedResult seeded_root_peak_batched(
    PeakBatchedModel<IdxToSz> base_model, TensorNetwork const& network,
    TIdxs const& tidxs, Index const& seed_axis, std::size_t occ_block) {
  SeededBatchedResult out;
  out.occ_block = occ_block;
  std::wstring const seed_label(seed_axis.full_label());

  auto base_batchable = base_model.is_batchable;
  auto base_batch = base_model.batch;
  base_model.is_batchable = [base_batchable, seed_label](Index const& ix) {
    return (base_batchable && base_batchable(ix)) ||
           std::wstring(ix.full_label()) == seed_label;
  };
  base_model.batch = [base_batch, seed_label, occ_block](Index const& ix) {
    if (std::wstring(ix.full_label()) == seed_label) return occ_block;
    return base_batch ? base_batch(ix) : std::size_t{0};
  };

  auto ctx = base_model.build_context(network, tidxs);
  // Locate the seed axis's bit in ctx.aux.
  std::size_t k_seed = ctx.m;  // sentinel == not found
  for (std::size_t k = 0; k < ctx.m; ++k)
    if (std::wstring(ctx.aux[k].full_label()) == seed_label) {
      k_seed = k;
      break;
    }
  SEQUANT_ASSERT(k_seed < ctx.m &&
                 "seed axis not recognized as a batchable axis");
  // SAFEGUARD: only honor a genuine external spectator.
  out.spectator_ok = base_model.is_spectator_axis(ctx, k_seed);
  SEQUANT_ASSERT(out.spectator_ok &&
                 "seed axis is not a pure external spectator "
                 "(open on the root AND contracted at no node)");
  out.seeded_axis = ctx.aux[k_seed];

  auto st = solve_single_term(base_model, network, tidxs, ctx);
  std::size_t const root = (std::size_t{1} << ctx.nt) - 1;
  std::size_t const seed = std::size_t{1} << k_seed;

  int const best0 = base_model.select_root(ctx, st, 0);
  int const bestS = base_model.select_root(ctx, st, seed);
  SEQUANT_ASSERT(best0 >= 0 && bestS >= 0);
  auto const& p0 = st[root][0][best0];
  auto const& pS = st[root][seed][bestS];
  out.unseeded_peak_bytes = p0.peak * base_model.numeric_size;
  out.seeded_peak_bytes = pS.peak * base_model.numeric_size;
  out.unseeded_flops = p0.flops;
  out.seeded_flops = pS.flops;
  return out;
}

/// \brief Independent memory-simulation recomputation of the chosen batched
/// reconstruction's model-A peak. Builds \ref PeakBatchedModel, runs
/// \ref solve_single_term for the back-pointer table, and recomputes the
/// subtree peak by direct memory simulation (NOT the DP's max/+ formula). Must
/// EQUAL \ref peak_cost_batched; a mismatch signals a DP/reconstruction bug.
template <typename TIdxs, typename IdxToSz>
double reconstructed_batched_peak(
    TensorNetwork const& network, TIdxs const& tidxs, IdxToSz&& idxsz,
    std::function<bool(Index const&)> const& is_batchable,
    std::function<std::size_t(Index const&)> const& batch_target_size,
    std::function<bool(Tensor const&)> const& is_volatile_leaf,
    double accumulation_factor = 0.0) {
  PeakBatchedModel<std::decay_t<IdxToSz>> model{
      std::forward<IdxToSz>(idxsz),
      is_batchable,
      batch_target_size,
      is_volatile_leaf,
      /* inner_pow */ {},
      /* volatile_weight */ 1.0,
      /* machine_balance */ 0.0,
      /* fast_mem_elems */ 0.0,
      /* block_tiles */ 3.0,
      /* block_prefactor */ 1.0,
      /* batch_persistent_only */ false,
      /* peak_flops_tolerance */ 0.0,
      accumulation_factor};
  auto ctx = model.build_context(network, tidxs);
  auto st = solve_single_term(model, network, tidxs, ctx);
  auto const nt = network.tensors().size();

  // Simulate the peak of evaluating subtree n at ancestor context B by walking
  // the chosen back-pointers. A leaf is resident at its own size. For an
  // internal node, with child context C = B | aprime, evaluate the lp_first
  // child fully (peak: its own simulated peak), then hold its result (sized at
  // C) while evaluating the second child (whose inputs co-reside at Lof), then
  // both results co-reside while the parent result (sized at B) is formed.
  // Re-derives the chosen reconstruction's peak by following the back-pointer
  // tree (contexts/orders chosen by the DP) and recomputing each child's peak
  // via recursion, rather than reading the DP's minimized st[*].peak table.
  // The per-node combination (stage_first/stage_second/stage_form) uses the
  // same staged-peak formula as the DP's lpf. What this validates independently
  // is the back-pointer walk itself (which children, order, context). The
  // Task-2/Task-3 batched oracle is the independent guard on the staged-peak
  // algebra.
  auto sim = [&](auto&& self, std::size_t n, std::size_t B, int idx) -> double {
    if (std::popcount(n) == 1) return ctx.sz(n, B);
    auto const& r = st[n][B][idx];
    std::size_t const C = B | r.aprime;
    std::size_t const f = r.lp_first ? r.lp : r.rp;  // evaluated first
    int const fi = r.lp_first ? r.lp_idx : r.rp_idx;
    std::size_t const s = r.lp_first ? r.rp : r.lp;  // evaluated second
    int const si = r.lp_first ? r.rp_idx : r.lp_idx;
    double const peak_f = self(self, f, C, fi);
    double const peak_s = self(self, s, C, si);
    // While the first child evaluates, the second child's leaf inputs sit
    // resident (Lof(s, C)). While the second child evaluates, the first
    // child's result sits resident (sz(f, C)). When both results exist, the
    // parent result (sz(n, B)) is materialized alongside them.
    double const stage_first = ctx.Lof(s, C) + peak_f;
    double const stage_second = ctx.sz(f, C) + peak_s;
    // Mirror the DP's relax() contribution charge: a node that contracts a
    // batchable index (aprime != 0) is accumulated over the aux batches, so
    // the in-flight contribution (sized like the result, ctx.sz(n, B))
    // co-resides with the accumulator at the all-co-resident moment.
    double const contrib =
        (r.aprime != 0) ? model.accumulation_factor * ctx.sz(n, B) : 0.0;
    double const stage_form =
        ctx.sz(f, C) + ctx.sz(s, C) + ctx.sz(n, B) + contrib;
    return std::max({stage_first, stage_second, stage_form});
  };

  std::size_t const root = (std::size_t{1} << nt) - 1;
  return sim(sim, root, 0, pareto_best(st[root][0]));
}

/// \brief Compile-time concept for a single-term-DP cost model.
///
/// A type \c M satisfies \c CostModel if it provides two associated types
/// (\c State and \c Context) and the six methods that \ref
/// run_single_term_opt calls.  The four built-in models -- \ref AdditiveModel
/// (FLOPs and Size variants), \ref PeakModel, and \ref PeakBatchedModel --
/// all satisfy this concept.  Users may also implement \c CostModel directly
/// and pass the instance to \ref run_single_term_opt to obtain an
/// \ref EvalSequence under any custom objective.
///
/// Requirements:
/// - \c M::State  -- the per-subset DP cell (the driver never inspects it).
/// - \c M::Context -- precomputed tables and mutable scratch.
/// - \c build_context(net, tidxs) \c const -> \c Context
/// - \c leaf(ctx, n) \c const -> \c State
/// - \c init(ctx, n) \c const -> \c State
/// - \c relax(ctx, n, lp, rp, lp_st, rp_st, acc) \c const -- updates \c acc
/// - \c finalize(ctx, n, states) \c const -- per-subset post-processing hook
/// - \c reconstruct(ctx, states) \c const -> \c EvalSequence
template <class M>
concept CostModel =
    requires {
      typename M::State;
      typename M::Context;
    } && requires(M const& m, TensorNetwork const& net,
                  container::svector<Index> const& tidxs,
                  typename M::Context& ctx, typename M::Context const& cctx,
                  size_t n, typename M::State const& cst, typename M::State& st,
                  container::vector<typename M::State>& sts,
                  container::vector<typename M::State> const& csts) {
      { m.build_context(net, tidxs) } -> std::same_as<typename M::Context>;
      { m.leaf(cctx, n) } -> std::same_as<typename M::State>;
      { m.init(cctx, n) } -> std::same_as<typename M::State>;
      m.relax(ctx, n, n, n, cst, cst, st);
      m.finalize(ctx, n, sts);
      { m.reconstruct(cctx, csts) } -> std::same_as<EvalSequence>;
    };

}  // namespace sequant::opt::detail

#endif  // SEQUANT_CORE_OPTIMIZE_COST_MODEL_HPP
