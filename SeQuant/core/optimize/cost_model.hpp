#ifndef SEQUANT_CORE_OPTIMIZE_COST_MODEL_HPP
#define SEQUANT_CORE_OPTIMIZE_COST_MODEL_HPP

#include <SeQuant/core/optimize/single_term_detail.hpp>  // helpers + EvalSequence + OptRes

#include <range/v3/view/concat.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <concepts>
#include <functional>
#include <limits>
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
  (void)tidxs;
  auto const nt = network.tensors().size();
  container::vector<typename Model::State> st(size_t{1} << nt);
  for (size_t n = 1; n < st.size(); ++n) {
    if (std::popcount(n) == 1) {
      st[n] = m.leaf(ctx, n);
      continue;
    }
    typename Model::State acc = m.init(ctx, n);
    for (auto&& [lp, rp] : bits::bipartitions(n))
      if (lp != 0 && rp != 0) m.relax(ctx, n, lp, rp, st[lp], st[rp], acc);
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
    init_results(network, tidxs, ctx.results);
    if (subnet_cse) {
      auto md = build_subnet_metadata(network, ctx.results);
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
    init_results(network, tidxs, results);
    auto fp = footprint_counter(idxsz, inner_pow);
    ctx.S.assign(sz, 0.0);
    ctx.idx.resize(sz);
    for (size_t n = 0; n < sz; ++n) {
      ctx.S[n] = (n == 0) ? 0.0 : fp(results[n].indices);
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
  /// Relative peak tolerance for the final (root) selection; see
  /// PeakModel::peak_flops_tolerance. 0 (default) = strict peak-min.
  double peak_flops_tolerance = 0.0;
  /// In-flight batch-contribution footprint multiplier; see
  /// BatchPolicy::accumulation_factor. Charged only on nodes that contract a
  /// batchable index (Ap != 0), into the all-co-resident peak term, to price
  /// the accumulator + contribution co-residency of K += contribution.
  double accumulation_factor = 0.0;

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
    std::function<double(IndexSet const&, IndexSet const&, IndexSet const&)>
        flops_of;

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
    // The accumulation_factor charge is per accumulation node (charged on each
    // node that contracts a batchable index). Its semantics are only
    // well-defined for a single batch axis; with multiple batchable indices the
    // per-node, once-per-node charge would conflate independent accumulations.
    SEQUANT_ASSERT(
        (accumulation_factor == 0.0 || ctx.m <= 1) &&
        "DensePeakSizeBatched: accumulation_factor != 0 requires at most one "
        "batchable index");
    ctx.nB = std::size_t{1} << ctx.m;
    ctx.tables = sliced_footprints(network, tidxs, idxsz, is_batchable, batch,
                                   ctx.aux, inner_pow);
    ctx.open_aux = subset_open_aux(network, tidxs, ctx.aux);
    ctx.volatile_mask = leaf_volatile_mask(network, is_volatile_leaf);
    // Per-subset open indices + a flop counter, for the lexicographic
    // (peak, then flops) tie-break (mirrors PeakModel). The flop tie-break uses
    // the unbatched contraction flops; total work summed over batches matches
    // it (work parity), so it consistently orders equal-peak schedules.
    auto const sz = std::size_t{1} << ctx.nt;
    container::vector<OptRes> results(sz);
    init_results(network, tidxs, results);
    ctx.idx.resize(sz);
    for (std::size_t n = 0; n < sz; ++n)
      ctx.idx[n] = std::move(results[n].indices);
    ctx.flops_of = [c = flops_counter(idxsz, inner_pow)](
                       IndexSet const& lhs, IndexSet const& rhs,
                       IndexSet const& result) { return c(lhs, rhs, result); };
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
        w * roofline_op_cost(ctx.flops_of(ctx.idx[lp], ctx.idx[rp], ctx.idx[n]),
                             ctx.sz(lp, 0) + ctx.sz(rp, 0) + ctx.sz(n, 0),
                             machine_balance, fast_mem_elems, block_tiles,
                             block_prefactor);
    for (std::size_t B = 0; B < ctx.nB; ++B) {
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
                                                  rp_st[C][ri].flops + cflops,
                                              lp, rp, lpf <= rpf, Ap, li, ri});
          }
        if (Ap == 0) break;
        Ap = (Ap - 1) & Acand;
      }
    }
  }

  void finalize(Context& /*ctx*/, size_t /*n*/,
                container::vector<State>& /*st*/) const {}

  EvalSequence reconstruct(Context const& ctx,
                           container::vector<State> const& st) const {
    std::size_t const root = (std::size_t{1} << ctx.nt) - 1;
    // ε-tolerant selection on the root's B=0 frontier: among points within
    // (1 + peak_flops_tolerance) of the minimum peak, fewest flops (ties broken
    // by lower peak). tolerance == 0 recovers strict peak-min.
    auto const& rootf = st[root][0];
    double minpeak = std::numeric_limits<double>::max();
    for (auto const& fp : rootf) minpeak = std::min(minpeak, fp.peak);
    double const thresh = minpeak * (1.0 + peak_flops_tolerance);
    int best = -1;
    for (int i = 0; i < static_cast<int>(rootf.size()); ++i)
      if (rootf[i].peak <= thresh &&
          (best < 0 || rootf[i].flops < rootf[best].flops ||
           (rootf[i].flops == rootf[best].flops &&
            rootf[i].peak < rootf[best].peak)))
        best = i;
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
    std::function<bool(Tensor const&)> const& is_volatile_leaf) {
  PeakBatchedModel<std::decay_t<IdxToSz>> model{std::forward<IdxToSz>(idxsz),
                                                is_batchable, batch_target_size,
                                                is_volatile_leaf};
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
    std::function<bool(Tensor const&)> const& is_volatile_leaf) {
  PeakBatchedModel<std::decay_t<IdxToSz>> model{std::forward<IdxToSz>(idxsz),
                                                is_batchable, batch_target_size,
                                                is_volatile_leaf};
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
    double const stage_form = ctx.sz(f, C) + ctx.sz(s, C) + ctx.sz(n, B);
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
