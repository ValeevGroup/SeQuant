#ifndef SEQUANT_CORE_OPTIMIZE_COST_MODEL_HPP
#define SEQUANT_CORE_OPTIMIZE_COST_MODEL_HPP

#include <SeQuant/core/eval/node_batch_annotation.hpp>
#include <SeQuant/core/optimize/single_term_detail.hpp>  // helpers + EvalSequence + OptRes

#include <range/v3/view/concat.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <string>
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
/// \p Model to additionally expose \c reconstruct_batched_modes (currently only
/// \ref PeakBatchedModel does); see its doc comment for the precise per-node
/// convention (RPN / post-order, left-first, matching the shared \c build
/// recursion used by \ref reconstruct).
///
/// \return The optimal EvalSequence, paired with one \c container::svector
///         of sliced \c Index per \c -1 token of that sequence, in the same
///         left-first post-order the sequence itself was emitted in. For the
///         nt==1 shortcut (no contractions) the modes vector is empty; for the
///         nt==2 shortcut (single contraction, no DP context is built) the
///         modes vector holds one empty entry (no batching info available).
template <class Model, typename TIdxs>
std::pair<EvalSequence, container::vector<NodeBatchAnnotation>>
run_single_term_opt_axes(Model const& m, TensorNetwork const& network,
                         TIdxs const& tidxs) {
  auto const nt = network.tensors().size();
  if (nt == 1) return {EvalSequence{0}, {}};
  if (nt == 2)
    return {EvalSequence{0, 1, -1},
            container::vector<NodeBatchAnnotation>{NodeBatchAnnotation{}}};
  typename Model::Context ctx = m.build_context(network, tidxs);
  auto st = solve_single_term(m, network, tidxs, ctx);
  return m.reconstruct_batched_modes(ctx, st);
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

/// \brief Slice-count-aware Pareto insert for the perf-first batched frontier.
///
/// With \p use_nsl == false this is byte-identical to \ref pareto_insert (plain
/// (peak, flops) domination; \c FP::nsl is not consulted). With
/// \p use_nsl == true (perf-first + a finite peak_threshold) the cumulative
/// sliced-mode count \c FP::nsl becomes a third Pareto objective, so a point
/// dominates only when it is no worse in peak, flops, and slice count. Because
/// contracted slicing is flops-neutral, the unsliced realization (higher peak,
/// \c nsl == 0) and a sliced one (lower peak, \c nsl > 0) are then Pareto-
/// incomparable: both survive on the frontier, all the way to the root, so an
/// ancestor that needs the lower-peak (sliced) subtree still has it and the
/// root selection still has the unsliced one to pick when it fits the budget.
/// \ref select_root then chooses the least-sliced feasible schedule -- no free
/// slicing below the ceiling. (A per-flops peak/nsl trade-off keeps at most one
/// point per distinct slice count, so the frontier stays bounded: slicing more
/// modes monotonically lowers peak.) The peak-first path passes false and is
/// unchanged.
template <typename FP>
void pareto_insert_ceiling(container::vector<FP>& f, FP p, bool use_nsl) {
  auto const dominates = [use_nsl](FP const& e, FP const& q) {
    if (!use_nsl) return e.peak <= q.peak && e.flops <= q.flops;
    return e.peak <= q.peak && e.flops <= q.flops && e.nsl <= q.nsl;
  };
  for (auto const& e : f)
    if (dominates(e, p)) return;  // p dominated -> skip
  f.erase(std::remove_if(f.begin(), f.end(),
                         [&](FP const& e) { return dominates(p, e); }),
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

/// \brief Peak-memory single-term cost model (DenseSpaceTime objective).
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
  /// k-aware inner (CSV/PNO composite) extent; see footprint_counter. Required
  /// whenever the network has composite indices (empty => inner_aware_volume
  /// throws); pass an explicit no-op only for composite-free networks. No
  /// default: omitting it silently mis-sized composites (4-PAO-integral bug).
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
    // CSE is not supported for DenseSpaceTime.
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
/// (DenseSpaceTimeBatched objective).
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
  std::function<std::size_t(Index const&)> batch;
  std::function<bool(Tensor const&)> is_volatile_leaf;
  /// k-aware inner (CSV/PNO composite) extent; see footprint_counter. Required
  /// whenever the network has composite indices (empty => inner_aware_volume
  /// throws); pass an explicit no-op only for composite-free networks. No
  /// default: omitting it silently mis-sized composites (4-PAO-integral bug).
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
  /// still consulted by the (unbatched) DenseSpaceTime model.
  double peak_flops_tolerance = 0.0;
  /// In-flight batch-contribution footprint multiplier; see
  /// BatchPolicy::accumulation_factor. Charged only on nodes that contract a
  /// batchable index (Ap != 0), into the all-co-resident peak term, to price
  /// the accumulator + contribution co-residency of K += contribution.
  double accumulation_factor = 0.0;
  /// Peak-memory budget (bytes) for threshold-gated selection; see
  /// BatchPolicy::peak_threshold. +infinity (default) => min-flops (no
  /// batching).
  double peak_threshold = std::numeric_limits<double>::infinity();
  /// Bytes per stored element, to compare the model's element-count peak to
  /// peak_threshold (bytes). Default 8 (double / TensorD).
  double numeric_size = 8.0;
  /// Perf-first / peak-second selection: when true, `select_root` selects the
  /// root-frontier point by (flops, then peak) and does not consult
  /// `peak_threshold` as a feasibility gate, so the threshold cannot force a
  /// flops-catastrophic factorization for its sliceability. Default false =
  /// peak-first, threshold-gated selection.
  bool perf_first = false;

  /// If true (the default), charge the batch recomputation cost on the flops/
  /// exec mode. The batched evaluator re-executes each contraction per tile of
  /// the ancestor batch modes its result does not carry (across-batch work is
  /// recomputed; within-batch sharing is cached -- see eval.hpp "replays the
  /// build of every compatible persistent final"). A node at ancestor-sliced-
  /// set B is charged nbatches(b) for each b in B not open in the node, so a
  /// schedule that slices many modes it must recompute across pays for it. The
  /// alternative (false) assumes work parity (batching is free on flops), which
  /// under-costs heavily-sliced families and does not reflect the true cost of
  /// batching; kept only as an escape hatch for comparison.
  ///
  /// The charge is order-aware: DP cells are ordered nests (see \ref
  /// Context::ordered), so a node is billed only for the enclosing loops outer
  /// to its innermost-carried placement (\ref Context::escaped_outer). Loops
  /// inner to that placement, and every loop when the node carries none of the
  /// enclosing modes, hoist above for free (rf == 1). An order-blind set charge
  /// would systematically over-bill the hoistable case.
  bool charge_batch_recompute = true;
  /// If true (the default), a value bound to a batch loop is charged as often
  /// as the loop runs. A loop runs once per production of the value that
  /// closes it (the node contracting the mode in batches), so a non-volatile
  /// subtree sliced by a loop that a volatile node closes is rebuilt every
  /// replay -- it cannot persist across evaluations (the runtime's persistence
  /// rule, explicit-cells design section 12) -- and takes \ref volatile_weight
  /// exactly as an amplitude-dependent value does. A non-volatile nest closed
  /// by a non-volatile node (a persistent intermediate summed over the DF
  /// index from two persistent factors) runs once and is charged once, as
  /// before. Without this the model assumed every amplitude-independent
  /// subtree is computed once, ever, regardless of slicing, and priced slicing
  /// a persistent operand of an amplitude contraction as free -- part of the
  /// water-20 occupied-batching cost gap (3x measured vs 1x modelled). Applied
  /// at the closing node, where its volatility is known: for each child
  /// subtree that is non-volatile and carries a mode sliced here, the child's
  /// whole cost is scaled (an over-charge only for descendants of that child
  /// which do not carry the mode and hoist above the loop). An external loop
  /// the DP opens here (\ref BFrontPoint::eopen) closes here too -- the full
  /// result is assembled outside it, per production of this node -- so it is
  /// charged the same way: the scale keys on the whole sliced set S (the
  /// contracted \c aprime or the external opens), not on \c aprime alone.
  bool charge_bound_persistence = true;
  /// Spaces batchable in the contracted role -- i.e. a mode of such a space is
  /// batchable where it is summed at some node. Building block; companion to
  /// \ref is_batchable_external_index. Declared here (outside the
  /// positionally-initialized prefix) and adjacent to its external companion so
  /// the two roles read together; set it by member assignment. Defaults to
  /// decline every index => no mode is batchable in the contracted role.
  std::function<bool(Index const&)> is_batchable_contracted_index =
      [](Index const&) { return false; };
  /// Spaces batchable in the external role -- i.e. a mode of such a space is
  /// batchable when it is open on the term root (a spectator carried to the
  /// result), not when it is contracted. Companion to \ref
  /// is_batchable_contracted_index, which admits spaces batchable in the
  /// contracted role. Keeping the two roles as separate caller-supplied space
  /// sets is what lets this layer stay domain-generic: the caller decides which
  /// spaces are batchable in which role, and \ref build_context drops every
  /// mode whose role's predicate rejects it (so a space batchable only as
  /// external never bloats the 2^m search with its contracted occurrences).
  /// Defaults to decline every index; a caller that wants external batching
  /// sets it explicitly (there is no fallback to the contracted-role
  /// predicate). Declared here (outside the positionally-initialized prefix) so
  /// existing aggregate initializations are unaffected; set it by member
  /// assignment.
  std::function<bool(Index const&)> is_batchable_external_index =
      [](Index const&) { return false; };

  /// Derived "batchable in any role": true iff the index is batchable in the
  /// contracted or the external role. This is not a settable field; the DP's
  /// role filters consume the individual building blocks, never this union.
  /// The building blocks default-decline, so both are always callable here.
  bool is_batchable(Index const& ix) const {
    return is_batchable_contracted_index(ix) || is_batchable_external_index(ix);
  }
  /// Prune disconnected (outer-product) subsets from the DP (see
  /// OptimizeOptions::prune_outer_products). Default true.
  bool prune_outer_products = true;
  /// Term-level gate for \ref reconstruct_batched_modes emitting \c
  /// BatchModeType::External entries (genuine external modes; see \ref
  /// is_external_mode), threaded from \ref CostParams::batch_spectator_indices
  /// / BatchPolicy::batch_spectator_indices. Default false so every other
  /// PeakBatchedModel construction (peak_cost_batched, compute_external_batch_
  /// axis's own model, existing tests) is unaffected and emits no External
  /// entries -- byte-identical to before this member existed. With it on, an
  /// external batch loop is a per-node DP choice (\ref BFrontPoint::eopen):
  /// external bits become nestable cell modes and \ref
  /// reconstruct_batched_modes stamps \c External exactly where the DP opened
  /// or inherited them. There is no post-DP placement pass.
  bool batch_spectator_indices = false;

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
    /// Cumulative count of batchable modes sliced anywhere in this realization
    /// (this node's \c aprime popcount plus both children's \c nsl). Used only
    /// by the perf-first ceiling's threshold-aware frontier domination
    /// (\ref pareto_insert_ceiling) to break peak-below-budget ties toward the
    /// least-sliced realization, so the unsliced schedule survives whenever it
    /// fits the budget. Zero for leaves and for the unbatched / peak-first
    /// paths, where it is never consulted.
    std::size_t nsl = 0;
    /// External modes whose batch loop this node opens: carried by the node,
    /// absent from its enclosing set, nested outer of \c aprime. The node is
    /// produced per batch and scattered into its full result outside the loop;
    /// the children see the mode in their enclosing set, non-carriers under it
    /// pay the recompute charge, and the opener is the loop's closer for
    /// \ref charge_bound_persistence. Children are read at
    /// descend(descend(B, eopen), aprime).
    std::size_t eopen = 0;
  };

  /// Per-subset DP cell: a \c [B]-vector (size \c nB = 2^m) of Pareto
  /// frontiers, one non-dominated (peak, flops) set per sliced-set context B.
  using State = container::vector<container::vector<BFrontPoint>>;

  /// Precomputed tables and per-(subset, sliced-set) lookup parameters built
  /// once by build_context.
  struct Context {
    /// Ordered, deduplicated batchable indices (bit \c k maps to \c
    /// batchable_modes[k]).
    container::vector<Index> batchable_modes;
    /// Number of batchable indices (= batchable_modes.size()).
    std::size_t m = 0;
    /// Number of sliced-sets (= 2^m).
    std::size_t nB = 1;
    /// Number of tensors in the network.
    std::size_t nt = 0;
    /// tables[B][n] = footprint of subset n under sliced-set B.
    container::vector<container::vector<double>> tables;
    /// open_modes[n] = bitmask of batchable indices open in subset n.
    container::vector<std::size_t> open_modes;
    /// Bitmask of volatile leaf tensors.
    std::size_t volatile_mask = 0;
    /// idx[n] = subset n's open (result) indices, for the flop tie-break.
    container::vector<IndexSet> idx;
    /// flops_of(lhs, rhs, result) = flop count of one binary contraction.
    /// Retained as the reference for the fast_flops parity test; the relax
    /// hot loop uses fast_flops (see below).
    std::function<double(IndexSet const&, IndexSet const&, IndexSet const&)>
        flops_of;
    /// nbatches[k] = number of batch tiles of batchable_modes[k] = ceil(extent
    /// / target), clamped to >= 1. Used to charge batch recomputation (see
    /// charge_batch_recompute): a node inside an ancestor batch loop over
    /// batchable_modes[k] that does not carry batchable_modes[k] is re-executed
    /// nbatches[k] times.
    container::vector<double> nbatches;

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

    /// Footprint of subset s under an explicit sliced-set union mask (the
    /// order-independent bitmask, as returned by \ref cell_union). Slicing is a
    /// pure footprint change, so the size depends only on which modes are in
    /// the union, never on the cell's nesting order -- this is the primitive
    /// that
    /// \ref sz and the external-placement re-price (\ref subtree_peak) share,
    /// so an external mode can be injected into the union directly without a
    /// precomputed ordered cell (external bits are excluded from build_cells).
    double sz_u(std::size_t U, std::size_t s) const {
      return tables[U & open_modes[s]][s];
    }
    /// Per-context leaf-sum of subset s under an explicit union mask.
    double Lof_u(std::size_t U, std::size_t s) const {
      double r = 0.0;
      for (std::size_t b = 0; b < nt; ++b)
        if (s & (std::size_t{1} << b)) r += sz_u(U, std::size_t{1} << b);
      return r;
    }

    /// Context-restricted size of subset s under sliced-set ctx (the table is
    /// indexed by the part of ctx actually open in s; mirrors the oracle).
    double sz(std::size_t s, std::size_t id) const {
      return sz_u(cell_union(id), s);
    }
    /// Per-context leaf-sum of subset s (sum of singleton sizes under ctx).
    double Lof(std::size_t s, std::size_t id) const {
      return Lof_u(cell_union(id), s);
    }

    // --- ordered-cell layer ------------------------------------------------
    // A DP cell is identified by `id`. Normally (`ordered` true) `id` indexes
    // an ordered sequence of batched modes (outer to inner): cell_union
    // recovers the bitmask for the (order-independent) footprint tables,
    // descend appends a contracted-here set as inner modes, and escaped_outer
    // charges only enclosing modes outer to a node's innermost-carried
    // placement. `ordered` is cleared only by build_cells' enumeration-blowup
    // guard, after which `id` is the sliced-set bitmask B (identity;
    // nCells == nB) and every helper here reduces to the plain bitmask ops --
    // a bounded-nest schedule is still correct, only optimality is lost.
    bool ordered = true;
    std::size_t cap = 3;     // max sequence length when ordered
    std::size_t nCells = 1;  // number of DP cells (== nB when !ordered)
    // Bitmask of external batchable modes (bit k set iff is_external_mode(k)):
    // open on the root, contracted at no node (built from open_modes, which
    // build_context assigns before build_cells()). Unless the term batches
    // spectators, build_cells skips these bits when enumerating sequences, so a
    // cell's union/sequence never contains an external mode -- nesting order is
    // meaningless for a mode that is never a contracted-here set at any node,
    // and admitting it would blow up the enumeration for nothing. With
    // spectator batching on they are nestable: the DP opens an external loop at
    // a node of its choosing (BFrontPoint::eopen).
    std::size_t external_mask = 0;
    container::vector<std::size_t> cell_union_;  // id -> union bitmask
    container::vector<container::svector<std::uint8_t>>
        cell_seq_;  // id -> ordered mode-bit indices
    container::vector<std::size_t>
        cell_descend_;  // id*m + k -> child id / SIZE_MAX

    std::size_t cell_union(std::size_t id) const {
      return ordered ? cell_union_[id] : id;
    }
    // Append the modes of `Ap` (ascending bit order = canonical co-contracted
    // order) as inner positions. SIZE_MAX if that exceeds cap or repeats a
    // mode.
    std::size_t descend(std::size_t id, std::size_t Ap) const {
      if (!ordered) return id | Ap;
      for (std::size_t k = 0; k < m; ++k)
        if (Ap & (std::size_t{1} << k)) {
          id = cell_descend_[id * m + k];
          if (id == std::numeric_limits<std::size_t>::max()) return id;
        }
      return id;
    }
    // The children's context of a frontier point: its external opens
    // (outer) then its contracted-here modes (inner).
    std::size_t descend_pt(std::size_t id, std::size_t eopen,
                           std::size_t Ap) const {
      if (eopen) {
        id = descend(id, eopen);
        if (id == std::numeric_limits<std::size_t>::max()) return id;
      }
      return descend(id, Ap);
    }
    // Enclosing modes charged as recompute for a node carrying `carried`:
    // those outer to the node's innermost-carried placement it does not carry.
    // Carr == 0 (no carried enclosing mode) hoists above every loop => none.
    std::size_t escaped_outer(std::size_t id, std::size_t carried) const {
      if (!ordered) return id & ~carried;
      auto const& seq = cell_seq_[id];
      int placement = -1;
      for (int p = static_cast<int>(seq.size()) - 1; p >= 0; --p)
        if ((carried >> seq[p]) & 1u) {
          placement = p;
          break;
        }
      if (placement < 0) return 0;
      std::size_t esc = 0;
      for (int p = 0; p < placement; ++p)
        if (!((carried >> seq[p]) & 1u)) esc |= (std::size_t{1} << seq[p]);
      return esc;
    }

    // Fill the cell tables: enumerate every ordered sequence of batched modes
    // up to `cap` length (id 0 = the empty sequence = the term root);
    // cell_union_/cell_seq_/cell_descend_ index them. `ext_nestable` admits the
    // external bits into the enumeration (spectator batching; see
    // external_mask). If the enumeration would blow up, `ordered` is cleared,
    // nCells == nB and the tables stay empty (helpers use the bitmask
    // directly).
    void build_cells(bool ext_nestable) {
      if (!ordered) {
        nCells = nB;
        return;
      }
      // Enumeration blowup guard: estimate Sum_{k<=cap} P(m,k) and fall back to
      // set-keyed if it would be huge (a bounded-nest schedule is still
      // correct; only optimality is lost). With cap==3 this is ~m^3/6, so it
      // only trips for pathologically many batchable indices.
      {
        std::size_t const limc = std::min(m, cap);
        std::size_t est = 1, term = 1;
        for (std::size_t k = 1; k <= limc; ++k) {
          term *= (m - k + 1);
          est += term;
        }
        if (est > 100000) {
          ordered = false;
          nCells = nB;
          return;
        }
      }
      cell_seq_.assign(1, container::svector<std::uint8_t>{});
      cell_union_.assign(1, std::size_t{0});
      std::map<container::svector<std::uint8_t>, std::size_t> seq_id;
      seq_id.emplace(cell_seq_[0], std::size_t{0});
      std::size_t const lim = std::min(m, cap);
      for (std::size_t id = 0; id < cell_seq_.size(); ++id) {
        if (cell_seq_[id].size() >= lim) continue;
        for (std::uint8_t k = 0; k < m; ++k) {
          if (!ext_nestable && ((external_mask >> k) & 1u))
            continue;  // external: not nestable unless the DP opens them
          if (cell_union_[id] & (std::size_t{1} << k)) continue;
          auto ns = cell_seq_[id];
          ns.push_back(k);
          if (seq_id.find(ns) == seq_id.end()) {
            seq_id.emplace(ns, cell_seq_.size());
            cell_union_.push_back(cell_union_[id] | (std::size_t{1} << k));
            cell_seq_.push_back(std::move(ns));
          }
        }
      }
      nCells = cell_seq_.size();
      cell_descend_.assign(nCells * m, std::numeric_limits<std::size_t>::max());
      for (std::size_t id = 0; id < nCells; ++id) {
        if (cell_seq_[id].size() >= lim) continue;
        for (std::uint8_t k = 0; k < m; ++k) {
          if (!ext_nestable && ((external_mask >> k) & 1u)) continue;
          if (cell_union_[id] & (std::size_t{1} << k)) continue;
          auto ns = cell_seq_[id];
          ns.push_back(k);
          cell_descend_[id * m + k] = seq_id.at(ns);
        }
      }
    }
  };

  template <typename TIdxs>
  Context build_context(TensorNetwork const& network,
                        TIdxs const& tidxs) const {
    // CSE is not supported for DenseSpaceTimeBatched.
    Context ctx;
    ctx.nt = network.tensors().size();
    // Candidates from both batchability roles (contracted / external); the role
    // filter below keeps each mode only if its actual role admits it.
    ctx.batchable_modes = batchable_mode_list(
        network, is_batchable_contracted_index, is_batchable_external_index);
    // open_modes must be assigned before build_cells() so the ordered-cell
    // enumeration can identify and exclude external modes (see
    // Context::external_mask below). Not pruned: is_external_mode scans it
    // over the full subset lattice (including disconnected subsets) to verify
    // a mode is never contracted, so every subset's open-mode bitmask must be
    // real. Computed here first because the role filter needs the root open
    // set, then recomputed if the filter shrinks the mode list.
    ctx.open_modes = subset_open_aux(network, tidxs, ctx.batchable_modes);
    {
      // Role filter: a mode open on the root occurs in the external role and is
      // batchable only if is_batchable_external_index admits it; otherwise it
      // occurs contracted and is batchable only if
      // is_batchable_contracted_index admits it. Dropping the rest keeps 2^m
      // free of modes that can never be batched in the role they actually occur
      // in (e.g. the contracted members of a space that is only
      // external-batchable). Each role consults only its own building block --
      // there is no fallback from the external role to the contracted one, so a
      // space batchable only in the contracted role never admits its external
      // occurrences. Both building blocks default-decline, hence are always
      // callable here.
      std::size_t const root = (std::size_t{1} << ctx.nt) - 1;
      container::vector<Index> kept;
      kept.reserve(ctx.batchable_modes.size());
      for (std::size_t k = 0; k < ctx.batchable_modes.size(); ++k) {
        Index const& ix = ctx.batchable_modes[k];
        bool const ext = (ctx.open_modes[root] >> k) & 1u;
        bool const keep = ext ? is_batchable_external_index(ix)
                              : is_batchable_contracted_index(ix);
        if (keep) kept.push_back(ix);
      }
      if (kept.size() != ctx.batchable_modes.size()) {
        ctx.batchable_modes = std::move(kept);
        ctx.open_modes = subset_open_aux(network, tidxs, ctx.batchable_modes);
      }
    }
    ctx.m = ctx.batchable_modes.size();
    // Per-mode batch-tile count for the recompute charge: ceil(extent/target),
    // >= 1. batch() returns the target tile size; a 0/absent target => 1 tile
    // (no batching of that mode => no recompute).
    ctx.nbatches.assign(ctx.m, 1.0);
    for (std::size_t k = 0; k < ctx.m; ++k) {
      double const ext = static_cast<double>(idxsz(ctx.batchable_modes[k]));
      double const tgt = static_cast<double>(batch(ctx.batchable_modes[k]));
      ctx.nbatches[k] = (tgt > 0.0) ? std::max(1.0, std::ceil(ext / tgt)) : 1.0;
    }
    // accumulation_factor is charged per accumulation node (Ap != 0) and is
    // valid for any number of batchable indices: with nested accumulation the
    // per-node charges co-exist at the peak. Validated by the identity
    // peak_cost_batched == reconstructed_batched_peak (test [batched-accum]).
    ctx.nB = std::size_t{1} << ctx.m;
    // A mode open on the root and contracted at no node (is_external_mode) is
    // carried unchanged from leaves to root: it never appears in any node's
    // contracted-here set. Without spectator batching the ordered enumeration
    // excludes it (cells hold contracted modes only); with it, the DP may open
    // its loop at a node of its choosing (see ext_nestable below).
    for (std::size_t k = 0; k < ctx.m; ++k)
      if (is_external_mode(ctx, k)) ctx.external_mask |= (std::size_t{1} << k);

    // Order-aware cell layer: DP cells index ordered sequences of batched
    // modes. cap == m => full enumeration (never worse than the set-keyed DP);
    // build_cells falls back to set-keyed if m is large.
    // Cap the ordered nest depth (sequence length), not m: the number of modes
    // that co-nest in a single term (m_B) is small (<=3 for C60 contracted),
    // while m (all batchable indices) can be large -- cap==m gives Sum P(m,k)
    // cells (13700 at m=7). Depth 3 => Sum_{k<=3} P(m,k) (260 at m=7),
    // polynomial in m, so it engages on large-m terms too. A term needing a
    // deeper nest loses only optimality (a bounded-nest schedule is still
    // correct), never correctness. build_cells still guards a hard cell-count
    // blowup.
    ctx.cap = std::min<std::size_t>(ctx.m, 3);
    // Spectator batching: external modes become nestable cell modes (their
    // loop is opened at a node of the DP's choosing); allow them room in the
    // nest cap. Decided here, after `cap` is set and before build_cells.
    bool const ext_nestable = batch_spectator_indices && ctx.external_mask != 0;
    if (ext_nestable)
      ctx.cap = std::min<std::size_t>(
          ctx.m, ctx.cap + std::popcount(ctx.external_mask));
    ctx.build_cells(ext_nestable);
    // Outer-product pruning: skip building tables for disconnected subsets the
    // DP will never form (solve_single_term also skips them). connected[n]==1
    // for singletons/empty and for connected subsets; the (~2x connected)
    // needed-mask is derived internally where a complement lookup requires it.
    auto const connected =
        outer_product_connectivity(network, tidxs, prune_outer_products);
    ctx.tables =
        sliced_footprints(network, tidxs, idxsz, is_batchable_contracted_index,
                          batch, ctx.batchable_modes, inner_pow, &connected);
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
    State s(ctx.nCells);
    for (std::size_t B = 0; B < ctx.nCells; ++B)
      s[B].push_back(BFrontPoint{ctx.sz(n, B), 0.0, 0, 0, true, 0, -1, -1});
    return s;
  }

  State init(Context const& ctx, size_t /*n*/) const {
    return State(ctx.nCells);  // empty frontiers; relax fills them
  }

  void relax(Context& ctx, size_t n, size_t lp, size_t rp, State const& lp_st,
             State const& rp_st, State& acc) const {
    // Secondary (tie-break) cost: roofline wall-time proxy per replay, charged
    // volatile_weight times for volatile (replayed) contractions. Uses the full
    // (unsliced) operand+result footprint as the per-replay traffic; slicing
    // reduces peak (primary mode), not total work. machine_balance==0 => flops.
    bool const is_volatile = (ctx.volatile_mask & n) != 0;
    double const w = is_volatile ? volatile_weight : 1.0;
    double const cflops =
        w * roofline_op_cost(
                ctx.use_fast_flops
                    ? ctx.fast_flops(lp, rp)
                    : ctx.flops_of(ctx.idx[lp], ctx.idx[rp], ctx.idx[n]),
                ctx.sz(lp, 0) + ctx.sz(rp, 0) + ctx.sz(n, 0), machine_balance,
                fast_mem_elems, block_tiles, block_prefactor);
    // Perf-first ceiling gate: under a finite budget, peak below the budget is
    // free, so flops-neutral contracted slicing must not be applied merely to
    // lower a sub-budget peak. Enabling nsl as a third Pareto objective keeps
    // the unsliced realization incomparable-to (hence co-resident with) the
    // sliced ones on the frontier; select_root then declines to slice below the
    // ceiling. Off for peak-first, or for a +inf budget where relax already
    // forces contracted_here == 0 (nothing to slice): byte-identical.
    bool const ceiling_nsl = perf_first && std::isfinite(peak_threshold);
    for (std::size_t B = 0; B < ctx.nCells; ++B) {
      // Batch recomputation charge: this node sits inside the ancestor batch
      // loops over the modes in B. For each b in B whose mode this node's
      // result does not carry (b not open in n), the node is re-executed
      // nbatches(b) times (across-batch recompute); modes it carries are
      // partitioned (x1). Default off => rf==1 => work-parity cost.
      double rf = 1.0;
      if (charge_batch_recompute) {
        // Recompute charge over the enclosing modes this node is re-executed
        // across. !ordered: B & ~carried (every escaped mode -- the set
        // charge). ordered: only escaped modes outer to the node's
        // innermost- carried placement (escaped modes inner to it hoist above
        // for free; Carr == 0 hoists above the whole nest => none, subsuming
        // A3a). escaped_outer collapses to the set charge when !ordered,
        // byte-identical.
        std::size_t const esc = ctx.escaped_outer(B, ctx.open_modes[n]);
        for (std::size_t k = 0; k < ctx.m; ++k)
          if (esc & (std::size_t{1} << k)) rf *= ctx.nbatches[k];
      }
      double const cflops_B = cflops * rf;
      // Batchable indices contracted at THIS node: open at children but not at
      // the parent. By default batching is applied ACROSS THE BOARD: slicing
      // the batch mode shrinks any intermediate carrying it regardless of
      // volatility (footprint objective) while leaving flops unchanged, so the
      // persistence gate would only ever raise the modelled peak. Set
      // batch_persistent_only to restore the persistent-only gate (decline to
      // slice subsets that contain a volatile leaf).
      // Revert-to-no-batching gate: an infinite peak_threshold is an unlimited
      // budget -- no term can ever be over budget, so nothing should batch.
      // Because contracted-index slicing is flops-neutral, the min-flops
      // frontier would otherwise always keep its fully-sliced (min-peak)
      // realization and slice unconditionally (a free peak reduction nobody
      // asked for). Forcing contracted_here == 0 here makes the batched cost
      // model produce the same schedule as the unbatched model when the budget
      // is unlimited: one model, cleanly reverting, rather than two models that
      // might disagree on the factorization. (No batchable axes -> ctx.m == 0
      // -> open_modes are all 0 -> contracted_here is already 0, so that revert
      // case needs no special handling.) A finite budget still enumerates
      // slicing; whether a term actually needs it is the select_root ceiling's
      // job.
      std::size_t const contracted_here =
          (!std::isfinite(peak_threshold) ||
           (batch_persistent_only && (ctx.volatile_mask & n)))
              ? std::size_t{0}
              : ((ctx.open_modes[lp] | ctx.open_modes[rp]) &
                 ~ctx.open_modes[n]);
      // External loops this node may open: carried here, not already sliced by
      // an enclosing loop. Under an unlimited budget nothing batches (the same
      // revert as contracted_here above). External bits are nestable cell modes
      // only under spectator batching, so this is 0 without it.
      std::size_t const ext_here =
          (batch_spectator_indices && std::isfinite(peak_threshold))
              ? (ctx.external_mask & ctx.open_modes[n] & ~ctx.cell_union(B))
              : std::size_t{0};
      // Enumerate every subset S of (contracted_here | ext_here), including
      // the empty set: Ap = its contracted part, E = its external part.
      std::size_t const choices = contracted_here | ext_here;
      std::size_t S = choices;
      while (true) {
        std::size_t const Ap = S & ~ctx.external_mask;
        std::size_t const E = S & ctx.external_mask;
        std::size_t const C = ctx.descend_pt(B, E, Ap);
        // ordered: skip an over-cap / repeat descent (SIZE_MAX). !ordered:
        // descend == B|Ap and never returns SIZE_MAX, so this is
        // byte-identical.
        if (C != std::numeric_limits<std::size_t>::max()) {
          double const szlp = ctx.sz(lp, C), szrp = ctx.sz(rp, C),
                       szn = ctx.sz(n, B);
          // A node that contracts a batchable index (Ap != 0) is accumulated
          // over the batches of that mode (K += contribution); the in-flight
          // contribution (same index set as the result, size szn) co-resides
          // with the accumulator. Charge it once, on the all-co-resident moment
          // only
          // -- the pre-result staged terms (Lrp+pl, szlp+prr) exclude it since
          // szn is not yet built.
          // An external loop opened here (E != 0) produces the node per
          // batch (size sliced on E) and scatters it into the full result
          // (szn, pre-sized, held across the loop); with a contracted mode
          // also sliced here the per-batch block accumulates over its
          // batches (accumulation_factor).
          double const contrib =
              (E != 0) ? ((Ap != 0) ? accumulation_factor : 1.0) *
                             ctx.sz_u(ctx.cell_union(B) | E, n)
                       : ((Ap != 0) ? accumulation_factor * szn : 0.0);
          double const both = szlp + szrp + szn + contrib;
          double const Lrp = ctx.Lof(rp, C), Llp = ctx.Lof(lp, C);
          // Resident-scan: a node that batches (Ap !=
          // 0) allocates its accumulator (szn) up front and holds it across the
          // batch loop, so it co-resides with the children as they evaluate per
          // batch -- the pre-result staged terms (which exclude szn as "not yet
          // built") must include it. Composed recursively, each ancestor
          // batching node's szn stacks onto every descendant's peak: the
          // Sum-of-enclosing-residents scan. S == 0 (unbatched, built once)
          // keeps szn out.
          double const res = (S != 0) ? szn : 0.0;
          // Cross every (peak,flops) trade-off of the two children at context
          // C.
          // Bound-persistence charge (charge_bound_persistence): this node
          // closes the loops over Ap. If it is volatile, every batch of those
          // loops re-runs each replay, so a non-volatile child sliced by them
          // (it carries a mode of Ap) is rebuilt each replay -- scale that
          // child's cost by volatile_weight (a volatile child already is).
          // (An external loop opened here is closed here too: the full result
          // is assembled outside it, per production of this node.)
          double const lp_scale =
              (charge_bound_persistence && S != 0 && is_volatile &&
               !(ctx.volatile_mask & lp) && (ctx.open_modes[lp] & S) != 0)
                  ? volatile_weight
                  : 1.0;
          double const rp_scale =
              (charge_bound_persistence && S != 0 && is_volatile &&
               !(ctx.volatile_mask & rp) && (ctx.open_modes[rp] & S) != 0)
                  ? volatile_weight
                  : 1.0;
          for (int li = 0; li < static_cast<int>(lp_st[C].size()); ++li)
            for (int ri = 0; ri < static_cast<int>(rp_st[C].size()); ++ri) {
              double const pl = lp_st[C][li].peak, prr = rp_st[C][ri].peak;
              double const lpf =
                  std::max({Lrp + pl + res, szlp + prr + res, both});
              double const rpf =
                  std::max({Llp + prr + res, szrp + pl + res, both});
              // Cumulative sliced-mode count of this realization: the modes
              // sliced at this node (popcount Ap) plus both children's. Fed to
              // the perf-first ceiling's threshold-aware domination so the
              // unsliced realization survives when it fits the budget.
              std::size_t const nsl =
                  lp_st[C][li].nsl + rp_st[C][ri].nsl +
                  static_cast<std::size_t>(std::popcount(S));
              pareto_insert_ceiling(
                  acc[B],
                  BFrontPoint{std::min(lpf, rpf),
                              lp_scale * lp_st[C][li].flops +
                                  rp_scale * rp_st[C][ri].flops + cflops_B,
                              lp, rp, lpf <= rpf, Ap, li, ri, nsl, E},
                  ceiling_nsl);
            }
        }  // C != SIZE_MAX
        if (S == 0) break;
        S = (S - 1) & choices;
      }
    }
  }

  void finalize(Context& /*ctx*/, size_t /*n*/,
                container::vector<State>& /*st*/) const {}

  /// \brief True iff batchable mode bit \p k is a genuine external mode of
  /// this network: it is open on the root result and is contracted at no node.
  ///
  /// An external mode is carried unchanged from the leaves up to the root:
  /// whenever any tensor in a subset carries the mode, the mode stays open in
  /// that subset, so it never appears in any node's contracted-at-node set
  /// (\c contracted_here = (open_modes[lp]|open_modes[rp]) & ~open_modes[n])
  /// and slicing it is purely a footprint change with identical work. This is
  /// what \c Context::external_mask is built from: such a mode is a nestable
  /// cell mode only under \ref batch_spectator_indices, and only then may the
  /// DP open its loop at a node (\ref BFrontPoint::eopen).
  bool is_external_mode(Context const& ctx, std::size_t k) const {
    std::size_t const root = (std::size_t{1} << ctx.nt) - 1;
    // Must be carried on the root result (a genuine external index).
    if (!((ctx.open_modes[root] >> k) & 1u)) return false;
    // Leaves (single-tensor subsets) that carry mode k.
    std::size_t leafmask = 0;
    for (std::size_t b = 0; b < ctx.nt; ++b)
      if ((ctx.open_modes[std::size_t{1} << b] >> k) & 1u)
        leafmask |= (std::size_t{1} << b);
    // Whenever the mode is available in a subset (some carrying leaf is in it)
    // it must be open in that subset -- else it is contracted at the node that
    // forms that subset, i.e. not a pure external mode.
    for (std::size_t n = 1; n <= root; ++n)
      if ((leafmask & n) && !((ctx.open_modes[n] >> k) & 1u)) return false;
    return true;
  }

  /// Threshold-gated root-frontier selection shared by \ref reconstruct and
  /// \ref reconstruct_batched_modes: among points whose peak (bytes) fits
  /// peak_threshold, pick fewest flops (ties by lower peak). If none fit, pick
  /// min peak (best effort). peak_threshold == +inf => all feasible => min
  /// flops => the non-batched schedule. Returns the chosen index into
  /// \c st[root][0] -- the term root is always read at cell 0 (the empty
  /// sequence): every loop the schedule uses, contracted or external, is
  /// opened inside the tree by the DP, so the root has no enclosing nest.
  int select_root(Context const& ctx,
                  container::vector<State> const& st) const {
    std::size_t const root = (std::size_t{1} << ctx.nt) - 1;
    auto const& rootf = st[root][0];
    auto peak_bytes = [this](double peak_elems) {
      return peak_elems * numeric_size;
    };
    if (perf_first) {
      // Perf-first / peak-second: min flops, ties by lower peak. The frontier
      // keeps one min-peak point per distinct flops value (pareto_insert prunes
      // equal-flops higher-peak points), so this both picks the cheapest
      // factorization and takes its fully-sliced (min-peak) realization.
      //
      // peak_threshold acts as a ceiling, not a min-peak objective: among the
      // points whose byte peak fits the budget, take the fewest flops (ties by
      // lower peak). A term is therefore batched only when its cheapest
      // (unbatched) schedule would exceed the budget -- flops are never traded
      // for peak below the ceiling, so this cannot force the flops-catastrophic
      // factorization that peak-first (the !perf_first branch) can. If no point
      // fits (the budget is below even the min-peak schedule) keep the
      // perf-first character and fall back to global min flops (best effort,
      // accepting the overage) rather than peak-first's min-peak fallback.
      // peak_threshold == +inf (the default when no budget is set) makes every
      // point feasible, reducing this to the pure min-flops selection --
      // byte-identical to before this ceiling existed.
      // Among the feasible frontier, min flops; ties broken toward the
      // least-sliced realization (fewer nsl), then lower peak. The nsl tiebreak
      // is what realizes the ceiling's intent: when the unsliced and a sliced
      // schedule both fit the budget at equal flops (kept distinct by
      // pareto_insert_ceiling), pick the unsliced one -- no free slicing below
      // the ceiling. With +inf / peak-first the frontier carries a single point
      // per flops, so nsl never discriminates: byte-identical.
      auto better = [&](int i, int j) {
        return rootf[i].flops < rootf[j].flops ||
               (rootf[i].flops == rootf[j].flops &&
                (rootf[i].nsl < rootf[j].nsl ||
                 (rootf[i].nsl == rootf[j].nsl &&
                  rootf[i].peak < rootf[j].peak)));
      };
      int pbest = -1;
      for (int i = 0; i < static_cast<int>(rootf.size()); ++i)
        if (peak_bytes(rootf[i].peak) <= peak_threshold &&
            (pbest < 0 || better(i, pbest)))
          pbest = i;
      bool const fit = pbest >= 0;  // a schedule met the ceiling
      if (!fit) {
        // Nothing fits the budget: best-effort min flops, ties by min peak (the
        // most-sliced realization). Not min nsl -- the nsl "don't slice below
        // the ceiling" tiebreak in `better` applies only among feasible points;
        // for an over-budget term it must be inverted, or the least-sliced
        // (max-peak) schedule is chosen and blows the budget by the most (a
        // giant unbatched integral -> OOM, e.g. the water-20 2.6 TB composite).
        // This restores the pre-nsl fallback: min flops, then min peak. (A
        // min-peak fallback was tried and rejected: for the infeasible terms it
        // takes the most-sliced, flops-catastrophic factorizations -- 8x the
        // FLOPs on water-20 at 100 GB.)
        auto const better_peak = [&](int i, int j) {
          return rootf[i].flops < rootf[j].flops ||
                 (rootf[i].flops == rootf[j].flops &&
                  rootf[i].peak < rootf[j].peak);
        };
        for (int i = 0; i < static_cast<int>(rootf.size()); ++i)
          if (pbest < 0 || better_peak(i, pbest)) pbest = i;
      }
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
      std::size_t const C = ctx.descend_pt(B, r.eopen, r.aprime);
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

  /// Modeled peak (in elements) of the chosen back-pointer subtree rooted at
  /// subset \p n, read at DP schedule cell \p Bsched and sized under the
  /// explicit union mask \p Usize. This is the reusable, node-level re-price of
  /// the peak reconstruction in \ref reconstructed_batched_peak: it follows the
  /// same back-pointer walk (children/order/aprime chosen by the DP) and the
  /// same staged-peak algebra (\c stage_first / \c stage_second / \c stage_form
  /// with the resident-scan \c res and accumulation \c contrib terms), but
  /// sizes every subset from \p Usize instead of from the schedule cell's
  /// union. With
  /// \p Usize == \c cell_union(Bsched) the two agree exactly.
  double subtree_peak(Context const& ctx, container::vector<State> const& st,
                      std::size_t n, std::size_t Bsched, std::size_t Usize,
                      int idx) const {
    if (std::popcount(n) == 1) return ctx.sz_u(Usize, n);
    auto const& r = st[n][Bsched][idx];
    std::size_t const C = ctx.descend_pt(Bsched, r.eopen, r.aprime);
    std::size_t const Uc = Usize | r.aprime;
    std::size_t const f = r.lp_first ? r.lp : r.rp;
    int const fi = r.lp_first ? r.lp_idx : r.rp_idx;
    std::size_t const s = r.lp_first ? r.rp : r.lp;
    int const si = r.lp_first ? r.rp_idx : r.lp_idx;
    double const peak_f = subtree_peak(ctx, st, f, C, Uc, fi);
    double const peak_s = subtree_peak(ctx, st, s, C, Uc, si);
    double const res = (r.aprime != 0) ? ctx.sz_u(Usize, n) : 0.0;
    double const stage_first = ctx.Lof_u(Uc, s) + peak_f + res;
    double const stage_second = ctx.sz_u(Uc, f) + peak_s + res;
    double const contrib =
        (r.aprime != 0) ? accumulation_factor * ctx.sz_u(Usize, n) : 0.0;
    double const stage_form =
        ctx.sz_u(Uc, f) + ctx.sz_u(Uc, s) + ctx.sz_u(Usize, n) + contrib;
    return std::max({stage_first, stage_second, stage_form});
  }

  /// Companion to \ref reconstruct that additionally reports, for each \c -1
  /// (contraction) entry emitted in the returned \c EvalSequence in emission
  /// order, the vector of \c Index sliced at that node (\c
  /// ctx.batchable_modes[bit] for each set bit of that node's \c aprime). Leaf
  /// entries contribute nothing. Does not change \ref reconstruct's own output;
  /// the two walks are kept in lock-step so the RPN order and the per-node
  /// modes line up.
  ///
  /// \p out_root_peak_bytes (when non-null) receives the term's reported root
  /// peak in bytes: the selected root frontier point's footprint, which already
  /// reflects every loop -- contracted and external alike -- the DP chose to
  /// open, since external opens are priced inside the DP (\ref
  /// BFrontPoint::eopen) rather than stamped on afterwards.
  std::pair<EvalSequence, container::vector<NodeBatchAnnotation>>
  reconstruct_batched_modes(Context const& ctx,
                            container::vector<State> const& st,
                            double* out_root_peak_bytes = nullptr) const {
    std::size_t const root = (std::size_t{1} << ctx.nt) - 1;
    int const best = select_root(ctx, st);  // shared helper (see below)
    // The selected root frontier point's own footprint is the reported peak:
    // every external loop the schedule uses was opened by the DP inside that
    // cost (BFrontPoint::eopen), so there is nothing to re-size afterwards.
    if (out_root_peak_bytes)
      *out_root_peak_bytes = st[root][0][best].peak * numeric_size;
    // Shared child-extraction for the back-pointer walk below: given a node
    // subset `n`, its enclosing cell `B`, and the chosen frontier index `idx`,
    // fetch the frontier point, descend to the children's cell `C`, and return
    // the two child subsets/indices in the canonical `lp_first` order.
    struct ChildFrontier {
      std::size_t f;
      int fi;
      std::size_t s;
      int si;
      std::size_t C;
    };
    auto child_frontier = [&](std::size_t n, std::size_t B,
                              int idx) -> ChildFrontier {
      auto const& r = st[n][B][idx];
      std::size_t const C = ctx.descend_pt(B, r.eopen, r.aprime);
      return ChildFrontier{
          r.lp_first ? r.lp : r.rp, r.lp_first ? r.lp_idx : r.rp_idx,
          r.lp_first ? r.rp : r.lp, r.lp_first ? r.rp_idx : r.lp_idx, C};
    };
    container::vector<NodeBatchAnnotation> node_axes;
    std::function<EvalSequence(std::size_t, std::size_t, int)> build =
        [&](std::size_t n, std::size_t B, int idx) -> EvalSequence {
      if (std::popcount(n) == 1)
        return EvalSequence{static_cast<int>(std::countr_zero(n))};
      BFrontPoint const& r = st[n][B][idx];
      auto const [fs, fi, ss, si, C] = child_frontier(n, B, idx);
      EvalSequence s = build(fs, C, fi);
      EvalSequence b = build(ss, C, si);
      s.insert(s.end(), b.begin(), b.end());
      container::svector<std::pair<Index, BatchModeType>> modes;
      // External modes first (outer), then Contracted: `modes` becomes
      // `ann.axes`, the node's own realized-loop order == the runtime pick
      // order (eval.hpp), so emitting External before Contracted realizes a
      // co-carried external mode outer of the contracted loop, preventing the
      // scatter from widening the external mode to full extent per contracted
      // block. Emit follows selection: stamp the external loops the DP itself
      // opened -- every external mode of this node's enclosing cell (or opened
      // right here) that the node carries -- and nothing else. Without
      // spectator batching no external bit ever enters a cell or an `eopen`,
      // so this is 0 and the push-order is a no-op.
      std::size_t const emit_mask_n =
          batch_spectator_indices ? ((ctx.cell_union(B) | r.eopen) &
                                     ctx.external_mask & ctx.open_modes[n])
                                  : std::size_t{0};
      for (std::size_t k = 0; k < ctx.m; ++k)
        if (emit_mask_n & (std::size_t{1} << k))
          modes.push_back({ctx.batchable_modes[k], BatchModeType::External});
      for (std::size_t k = 0; k < ctx.m; ++k)
        if (r.aprime & (std::size_t{1} << k))
          modes.push_back({ctx.batchable_modes[k], BatchModeType::Contracted});
      NodeBatchAnnotation ann;
      ann.axes = std::move(modes);
      // Loop-open emission: the subset of axes for which this node
      // introduces the physical batch loop (vs a node that only carries the
      // sliced mode). peak_profile builds its enclosing-loop nest (ectx) from
      // this, so one physical loop counts once instead of
      // once-per-carrying-node.
      //  - External: the DP's own open site, this frontier point's `eopen`.
      //  - Contracted: a contracted batch loop contracts (and so opens) at its
      //    unique node = this node's aprime, mirroring the Contracted axes
      //    above.
      {
        std::size_t const opened_ext_mask =
            batch_spectator_indices ? r.eopen : std::size_t{0};
        container::svector<std::pair<Index, BatchModeType>> opened;
        for (std::size_t k = 0; k < ctx.m; ++k)
          if (opened_ext_mask & (std::size_t{1} << k))
            opened.push_back({ctx.batchable_modes[k], BatchModeType::External});
        for (std::size_t k = 0; k < ctx.m; ++k)
          if (r.aprime & (std::size_t{1} << k))
            opened.push_back(
                {ctx.batchable_modes[k], BatchModeType::Contracted});
        ann.opened_here = std::move(opened);
      }
      // Order-aware placement/lifetime bridge: emit this node's order-aware
      // gate and effective use count for a later runtime hoist pass. `B` is the
      // node's ordered enclosing cell, so all the needed pieces (open_modes,
      // escaped_outer, nbatches) are in hand at this frame -- no second
      // structure or pass. The node's runtime residency (its
      // enclosing-contracted and external modes alike) is not emitted here: it
      // is the all-batched-modes cross-occurrence meet stamped onto
      // EvalExpr::sliced_modes by stamp_lifetime_masks.
      {
        // Order-aware gate for per-level placement: this node participates in
        // hoist placement even when its residency is empty (a whole-nest
        // invariant -> root scope).
        ann.order_aware = true;
        // effective_count = rf = prod of nbatches[k] over the escaped-outer set
        // (enclosing loops the node does not carry). With no CSE on this path
        // the back-pointer object is a strict tree (single consumer), so this
        // per-node rf is the effective use count.
        std::size_t const esc = ctx.escaped_outer(B, ctx.open_modes[n]);
        double rf = 1.0;
        for (std::size_t k = 0; k < ctx.m; ++k)
          if (esc & (std::size_t{1} << k)) rf *= ctx.nbatches[k];
        ann.effective_count = static_cast<std::size_t>(std::llround(rf));
      }
      node_axes.push_back(std::move(ann));  // one entry per -1, in RPN order
      s.push_back(-1);
      return s;
    };
    auto seq = build(root, 0, best);
    return {std::move(seq), std::move(node_axes)};
  }
};

/// \brief Achieved minimum peak memory (the DenseSpaceTime objective value) for
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
  model.is_batchable_contracted_index = is_batchable;
  auto ctx = model.build_context(network, tidxs);
  auto st = solve_single_term(model, network, tidxs, ctx);
  // root subset's B=0 frontier; its smallest peak is the achieved minimum.
  double mn = std::numeric_limits<double>::max();
  for (auto const& fp : st.back()[0]) mn = std::min(mn, fp.peak);
  return mn;
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
  model.is_batchable_contracted_index = is_batchable;
  auto ctx = model.build_context(network, tidxs);
  auto st = solve_single_term(model, network, tidxs, ctx);
  auto const nt = network.tensors().size();

  // Simulate the peak of evaluating the chosen root subtree by walking the
  // chosen back-pointers (which children, order, context chosen by the DP) and
  // recomputing each node's peak via the staged-peak algebra
  // (stage_first/stage_second/stage_form with the resident-scan res term and
  // the accumulation contrib term), rather than reading the DP's minimized
  // st[*].peak table. This is exactly PeakBatchedModel::subtree_peak sized at
  // the schedule cell's own union (no external injection): sizing under
  // cell_union(B) is what makes the two agree exactly. What this
  // validates independently is the back-pointer walk itself; the batched oracle
  // is the independent guard on the staged-peak algebra.
  std::size_t const root = (std::size_t{1} << nt) - 1;
  int const best = pareto_best(st[root][0]);
  return model.subtree_peak(ctx, st, root, /*Bsched=*/0,
                            /*Usize=*/ctx.cell_union(0), best);
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
