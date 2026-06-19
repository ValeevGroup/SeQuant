#ifndef SEQUANT_CORE_OPTIMIZE_COST_MODEL_HPP
#define SEQUANT_CORE_OPTIMIZE_COST_MODEL_HPP

#include <SeQuant/core/optimize/single_term.hpp>  // helpers + EvalSequence + OptRes

#include <range/v3/view/concat.hpp>

#include <algorithm>
#include <bit>
#include <concepts>
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
/// Reproduces \ref single_term_opt_impl exactly, factored into the
/// CostModel hooks driven by \ref run_single_term_opt. The per-contraction
/// cost is supplied by \p CostFn (e.g. \ref flops_counter or
/// \ref memsize_counter); a per-intermediate footprint penalty, volatile
/// replay weighting, and subnet common-subexpression elimination (CSE) are
/// all handled identically to the original implementation.
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
    // Rebuild per-subset sequences bottom-up, exactly as single_term_opt_impl's
    // final pass: singletons emit their tensor index, internal nodes emit the
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

/// \brief Peak-memory single-term cost model (DensePeakSize objective).
///
/// Reproduces \ref single_term_opt_peak_impl / \ref peak_dp exactly, factored
/// into the CostModel hooks driven by \ref run_single_term_opt. The recurrence
/// is the all-co-resident pebble game: while one child evaluates, the other
/// child's leaf inputs sit resident. No CSE; no volatile weighting.
///
/// \tparam IdxToSz A callable mapping an Index to its extent.
template <typename IdxToSz>
struct PeakModel {
  IdxToSz idxsz;

  /// Per-subset DP cell: minimum peak to build the subtree rooted at this
  /// subset, plus the winning bipartition and evaluation order.
  struct State {
    double peak = std::numeric_limits<double>::max();
    size_t lp = 0;
    size_t rp = 0;
    bool lp_first = true;
  };

  /// Precomputed tables: S[n] = footprint of subset n's result tensor,
  /// L[n] = sum of leaf (singleton) sizes in subset n.
  struct Context {
    container::vector<double> S;
    container::vector<double> L;
  };

  template <typename TIdxs>
  Context build_context(TensorNetwork const& network,
                        TIdxs const& tidxs) const {
    // CSE is not supported for DensePeakSize.
    Context ctx;
    ctx.S = subset_footprints(network, tidxs, idxsz);
    auto const sz = ctx.S.size();
    auto const nt = network.tensors().size();
    ctx.L.assign(sz, 0.0);
    for (size_t n = 0; n < sz; ++n)
      for (size_t b = 0; b < nt; ++b)
        if (n & (size_t{1} << b)) ctx.L[n] += ctx.S[size_t{1} << b];
    return ctx;
  }

  State leaf(Context const& ctx, size_t n) const {
    return State{ctx.S[n], 0, 0, true};
  }

  State init(Context const& /*ctx*/, size_t /*n*/) const {
    return State{std::numeric_limits<double>::max(), 0, 0, true};
  }

  void relax(Context& ctx, size_t /*n*/, size_t lp, size_t rp,
             State const& lp_st, State const& rp_st, State& acc) const {
    double const both = ctx.S[lp] + ctx.S[rp] + ctx.S[lp | rp];
    double const lp_first_cand =
        std::max({ctx.L[rp] + lp_st.peak, ctx.S[lp] + rp_st.peak, both});
    double const rp_first_cand =
        std::max({ctx.L[lp] + rp_st.peak, ctx.S[rp] + lp_st.peak, both});
    double const cand = std::min(lp_first_cand, rp_first_cand);
    if (cand < acc.peak) {
      acc.peak = cand;
      acc.lp = lp;
      acc.rp = rp;
      acc.lp_first = (lp_first_cand <= rp_first_cand);
    }
  }

  void finalize(Context& /*ctx*/, size_t /*n*/,
                container::vector<State>& /*st*/) const {}

  EvalSequence reconstruct(Context const& /*ctx*/,
                           container::vector<State> const& st) const {
    using ranges::views::concat;
    // Bottom-up emit: a subset's sequence is (first-child seq)(second-child
    // seq) -1, where first/second are chosen by lp_first (the lower-peak child
    // is emitted last, i.e. evaluated second / built first in the sequence).
    container::vector<EvalSequence> seq(st.size());
    for (size_t n = 0; n < st.size(); ++n) {
      if (std::popcount(n) == 1) {
        seq[n] = EvalSequence{static_cast<int>(std::countr_zero(n))};
      } else if (std::popcount(n) >= 2) {
        auto const& a = st[n].lp_first ? seq[st[n].lp] : seq[st[n].rp];
        auto const& b = st[n].lp_first ? seq[st[n].rp] : seq[st[n].lp];
        seq[n] = concat(a, b) | ranges::to<EvalSequence>;
        seq[n].push_back(-1);
      }
    }
    return seq.back();
  }
};

/// \brief Multi-mode batched peak-memory single-term cost model
/// (DensePeakSizeBatched objective).
///
/// Reproduces \ref single_term_opt_peak_batched_impl / \ref peak_dp_batched
/// exactly, factored into the CostModel hooks driven by \ref
/// run_single_term_opt. The recurrence is the per-batchable-index all-co-
/// resident pebble game of the batch-aware cost model design (section 6.2,
/// model A): each DP cell is indexed by both a subset \c n and a sliced-set
/// context \c B over the batchable indices, so a model State is the \c [B]
/// -vector of per-context \ref BatchedRes. No CSE; persistence-gated batching.
///
/// \tparam IdxToSz A callable mapping an Index to its extent.
template <typename IdxToSz>
struct PeakBatchedModel {
  IdxToSz idxsz;
  std::function<bool(Index const&)> is_batchable;
  std::function<std::size_t(Index const&)> batch;
  std::function<bool(Tensor const&)> is_volatile_leaf;

  /// Per-subset DP cell: the full \c [B]-vector (size \c nB = 2^m) of per-
  /// sliced-set \ref BatchedRes back-pointers.
  using State = container::vector<BatchedRes>;

  /// Precomputed tables and per-(subset, sliced-set) lookup parameters built
  /// once by build_context, mirroring \ref peak_dp_batched's locals.
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
    ctx.nB = std::size_t{1} << ctx.m;
    ctx.tables =
        sliced_footprints(network, tidxs, idxsz, is_batchable, batch, ctx.aux);
    ctx.open_aux = subset_open_aux(network, tidxs, ctx.aux);
    ctx.volatile_mask = leaf_volatile_mask(network, is_volatile_leaf);
    return ctx;
  }

  State leaf(Context const& ctx, size_t n) const {
    State s(ctx.nB);
    for (std::size_t B = 0; B < ctx.nB; ++B) s[B].peak = ctx.sz(n, B);
    return s;
  }

  State init(Context const& ctx, size_t /*n*/) const {
    return State(ctx.nB);  // BatchedRes defaults to peak == max
  }

  void relax(Context& ctx, size_t n, size_t lp, size_t rp, State const& lp_st,
             State const& rp_st, State& acc) const {
    for (std::size_t B = 0; B < ctx.nB; ++B) {
      auto& curr = acc[B];
      // Batchable indices contracted at THIS node: open at children but not at
      // parent, gated by persistence.
      std::size_t const Acand =
          (ctx.volatile_mask & n)
              ? std::size_t{0}
              : ((ctx.open_aux[lp] | ctx.open_aux[rp]) & ~ctx.open_aux[n]);
      // Enumerate every subset A' of Acand (including the empty set).
      std::size_t Ap = Acand;
      while (true) {
        std::size_t const C = B | Ap;
        double const pl = lp_st[C].peak;
        double const prr = rp_st[C].peak;
        double const both = ctx.sz(lp, C) + ctx.sz(rp, C) + ctx.sz(n, B);
        double const lpf =
            std::max({ctx.Lof(rp, C) + pl, ctx.sz(lp, C) + prr, both});
        double const rpf =
            std::max({ctx.Lof(lp, C) + prr, ctx.sz(rp, C) + pl, both});
        double const cand = std::min(lpf, rpf);
        if (cand < curr.peak) {
          curr.peak = cand;
          curr.lp = lp;
          curr.rp = rp;
          curr.lp_first = (lpf <= rpf);
          curr.aprime = Ap;
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
    using ranges::views::concat;
    // Recursive back-pointer walk from the root subset at the empty ancestor
    // context (root, B=0). At each node (n, B) read the winning aprime, form
    // the child context C = B | aprime, and recurse into the children at C in
    // lp_first order (the lower-peak child is emitted last). A singleton emits
    // its tensor index; an internal node emits concat(first, second), -1. Each
    // st[child] is a [B]-vector, so children are read at index C.
    auto build = [&](auto&& self, std::size_t n,
                     std::size_t B) -> EvalSequence {
      if (std::popcount(n) == 1)
        return EvalSequence{static_cast<int>(std::countr_zero(n))};
      auto const& r = st[n][B];
      std::size_t const C = B | r.aprime;
      auto first = self(self, r.lp_first ? r.lp : r.rp, C);
      auto second = self(self, r.lp_first ? r.rp : r.lp, C);
      auto seq = concat(first, second) | ranges::to<EvalSequence>;
      seq.push_back(-1);
      return seq;
    };

    std::size_t const root = (std::size_t{1} << ctx.nt) - 1;
    return build(build, root, 0);
  }
};

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
