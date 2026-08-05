#ifndef SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_COST_PROFILE_HPP
#define SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_COST_PROFILE_HPP

#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/lifetime_mask.hpp>
#include <SeQuant/core/eval/placement_router.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/optimize/single_term_detail.hpp>

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace sequant::eval::dryrun {

/// Configuration for a faithful (gated) dry-run cache: the same footprint gate
/// and cross-occurrence batch-variant veto the real batched eval loop applies,
/// so a batch-variant giant (a `mu~`/`K`-carrying DF intermediate) is NOT
/// cached whole but recomputed sliced under each consumer's batch trigger.
///
/// The element types mirror the gated \c sequant::cache_manager overload
/// (\c cache_manager.hpp): \c is_volatile is invoked on every \c TreeNode
/// (deduced as \c EvalNodeDryRun for the dry-run backend).
///
/// This struct lives here (rather than only in the test) because Task 4's
/// \c cost_profile() entry point consumes it.
struct CacheConfig {
  /// Footprint gate (bytes): a node whose result footprint exceeds this is not
  /// cached. 0 (default) disables the gate.
  double max_footprint = 0.;
  /// Minimum non-persistent repeats to cache an internal node (CSE rule).
  std::size_t min_repeats = 2;
  /// `bool(EvalNodeDryRun const&)`: true if the node is intrinsically volatile
  /// (typically the amplitude leaves). Empty => nothing is volatile.
  std::function<bool(EvalNodeDryRun const&)> is_volatile;
};

/// Builds a gated dry-run cache from an eval-node range, a \p cfg, and a
/// \p regime that supplies the moment-aware node-size model used for the
/// footprint gate.
///
/// The footprint functor sizes a node's result (its \c canon_indices()) with
/// the SAME moment-aware counter the DryRun \c Result uses
/// (\c memsize_counter over \c regime.idx_to_extent()/inner_pow_fn()), scaled
/// to bytes, so the gate compares like-for-like against \c cfg.max_footprint.
///
/// Unlike the SIMPLE \c cache_manager(nodes) factory the ad-hoc dry-run test
/// sites use, this routes through the GATED overload so free-batchable-axis
/// giants are vetoed (matching the real run). The returned cache is used across
/// the WHOLE forest without a per-summand reset (matching a real solve's
/// whole-iteration cache scope; \c cost_profile relies on the lifetime mask to
/// release each value after its last cross-term use, so cross-summand-shared
/// values are reused, not rebuilt).
///
/// \param nodes the evaluation forest (a range of \c EvalNodeDryRun).
/// \param cfg   footprint/repeat/volatility/batchability configuration.
/// \param regime the size regime supplying extents and CSV moment tables.
/// \return a \c CacheManager over \c EvalNodeDryRun.
template <typename NodeRange>
auto build_dryrun_cache(NodeRange const& nodes, CacheConfig const& cfg,
                        SizeRegime const& regime) {
  auto memsize = sequant::opt::detail::memsize_counter(regime.idx_to_extent(),
                                                       regime.inner_pow_fn());

  // Footprint (bytes) of a node's RESULT: canon_indices() fed to the
  // moment-aware counter (as the counter's `result` slot; the empty lhs/rhs
  // contribute nothing) times 8 bytes/element. Same arithmetic as the DryRun
  // Result::size_in_bytes(), so the gate is faithful.
  auto footprint_of =
      [memsize = std::move(memsize)](EvalNodeDryRun const& n) -> double {
    std::vector<Index> const result(n->canon_indices().begin(),
                                    n->canon_indices().end());
    return memsize(std::vector<Index>{}, std::vector<Index>{}, result) * 8.0;
  };

  // Default the volatility predicate so the gated factory never invokes an
  // empty std::function (nothing volatile leaves that gate inert, matching the
  // factory's own default).
  std::function<bool(EvalNodeDryRun const&)> is_volatile =
      cfg.is_volatile ? cfg.is_volatile
                      : std::function<bool(EvalNodeDryRun const&)>(
                            [](EvalNodeDryRun const&) { return false; });

  // Note: the gated sequant::cache_manager() overload below stamps the
  // cross-occurrence lifetime mask on `nodes` itself before its DAG walk /
  // veto (cache_manager.hpp), so this call site does not need to do so.
  return sequant::cache_manager(nodes, std::move(is_volatile), cfg.min_repeats,
                                std::move(footprint_of), cfg.max_footprint);
}

/// One recomputed value's avoidable-recompute breakdown (see
/// \c CostProfile::avoidable_nodes). \c label is the value's signature (a
/// `result;lhs;rhs` full-label string built in \c DryRunOps::prod). \c flops is
/// the avoidable recompute in FLOPs -- `total_flops - full_flops`, the
/// arithmetic the batched replay repeated beyond building this value ONCE at
/// full extent (the batching-free / unlimited-memory ideal). \c count is the
/// equivalent number of extra full rebuilds (`total_flops/full_flops - 1`): 0
/// when the builds tile the value once (disjoint slices), ~N-1 when the value
/// is rebuilt full N times (an un-hoisted invariant). Only values with a
/// positive avoidable FLOP count are recorded.
struct AvoidableNode {
  std::string label;
  double count = 0;
  double flops = 0;
};

/// Roll a replay's per-value build tally (a \c CostSink populated by
/// \c CostModel::tally_node during a \c Trace::On replay) into the avoidable-
/// recompute breakdown: per distinct value (label), avoidable FLOPs =
/// max(0, total_flops - full_flops), sorted by avoidable FLOPs descending,
/// keeping only values with a positive amount. Shared by \c cost_profile()
/// (whole-forest rollup) and any caller that attaches its OWN sink to a replay
/// it drives (e.g. the schedule-dump test, which emits these per-value numbers
/// so the visualizer consumes cost_profile's avoidable instead of recomputing
/// it). Single-threaded caller expected: the replay has finished and the sink
/// is quiescent, so no lock is taken.
inline std::vector<AvoidableNode> avoidable_nodes_from_sink(
    CostSink const& sink) {
  // Per value (keyed by label): avoidable FLOPs = total_flops - full_flops --
  // the arithmetic the batched replay repeated beyond building the value ONCE
  // at full extent (the batching-free ideal). FLOPs is additive across slices,
  // so disjoint per-block slices that tile the value sum to exactly full_flops
  // (0 avoidable), while a value rebuilt full per block sums to N*full
  // ((N-1)*full avoidable). The per-value max(0, .) keeps every entry -- and
  // hence the sum
  // -- in [0, total_flops], so avoidable_time() stays in [0, 1] by
  // construction.
  std::vector<AvoidableNode> out;
  for (auto const& [label, nc] : sink.per_node) {
    if (nc.full_flops <= 0.0) continue;
    double const avoidable = nc.total_flops - nc.full_flops;
    if (avoidable <= 0.0) continue;
    out.push_back({label,
                   nc.total_flops / nc.full_flops - 1.0 /*extra full rebuilds*/,
                   avoidable});
  }
  std::sort(out.begin(), out.end(),
            [](AvoidableNode const& a, AvoidableNode const& b) {
              return a.flops > b.flops;
            });

  // DIAGNOSTIC (SEQUANT_AVOIDABLE_DEBUG): per-value FLOP accounting, worst
  // first. total_flops == full_flops => builds tiled the value once (0
  // avoidable); total_flops == N*full_flops => value rebuilt full N times.
  if (std::getenv("SEQUANT_AVOIDABLE_DEBUG")) {
    double sum_total = 0, sum_avoid = 0;
    for (auto const& [label, nc] : sink.per_node) {
      sum_total += nc.total_flops;
      if (nc.full_flops > 0.0 && nc.total_flops > nc.full_flops)
        sum_avoid += nc.total_flops - nc.full_flops;
    }
    std::fprintf(stderr,
                 "[avoidable-debug] dryrun_flops=%.6g avoidable_flops=%.6g "
                 "frac=%.4f n_values=%zu\n",
                 sum_total, sum_avoid,
                 sum_total > 0 ? sum_avoid / sum_total : 0,
                 sink.per_node.size());
    std::size_t shown = 0;
    for (auto const& an : out) {
      if (shown++ >= 20) break;
      auto const& nc = sink.per_node.at(an.label);
      std::fprintf(stderr,
                   "  builds=%zu total_flops=%.4g full_flops=%.4g "
                   "avoid_flops=%.4g (=%.3g full-rebuilds) | %s\n",
                   nc.builds, nc.total_flops, nc.full_flops, an.flops, an.count,
                   an.label.c_str());
    }
  }
  return out;
}

/// Summary of the modeled cost of a factorized dry-run eval forest, as produced
/// by \c cost_profile(). All quantities are summed/maxed over every summand
/// tree in the forest.
struct CostProfile {
  /// Predicted peak working-set (bytes): the max over summands of the
  /// batched-scratch high-watermark folded by the Task-3 \c PeakSink and the
  /// outer gated cache's \c working_set_hwmark().
  ///
  /// This ACCOUNTS FOR co-resident residency across the scope chain, rather
  /// than being the max-of-independent-hwmarks lower bound it was before: each
  /// per-op hwmark folded into a cache's \c working_set_hwmark_ (in eval.hpp)
  /// already adds \c CacheManager::chain_residency() of that cache's
  /// scope-chain ancestors at the instant of the op, so a scratch cache's
  /// high-watermark is the max over its life of (scratch residency +
  /// everything alive up its parent chain at that instant) -- the co-resident
  /// sum when a persistent cross-term cache entry is alive at the same instant
  /// as a batched-inner transient. The outer (root) cache has no parent, so
  /// its own hwmark is unaffected (added term is 0); the \c max() fold here
  /// (over \c peak.load() and \c cache.working_set_hwmark()) is then correct in
  /// both regimes -- batched (peak.load() already carries the co-resident
  /// scratch peak) and unbatched (peak.load() == 0, the outer hwmark alone is
  /// the peak).
  ///
  /// Each live buffer is counted exactly once: the per-op operand guards skip
  /// an operand whose buffer any cache on the scope chain already holds
  /// (\c CacheManager::chain_holds(), pointer identity), so a value read full
  /// from an ANCESTOR cache is counted only via \c chain_residency() and not
  /// again as an operand, while a sliced/permuted/phase-shifted read (a
  /// distinct buffer) is correctly added.
  ///
  /// One known deviation remains -- an UNDER-count that keeps this a lower
  /// bound in exactly one place: the external-scatter accumulator (\c dest in
  /// eval.hpp) is a plain local, not a cache entry, so it is invisible to
  /// \c chain_residency() even though it co-resides with every inner block's
  /// working set. Pre-existing, outside this field's scope.
  double peak_bytes = 0;
  /// STATIC per-node DP-MODEL FLOPs: summed unweighted static contraction FLOPs
  /// over all internal nodes, from the order-/batching-blind static walk. NOT
  /// CSE-aware across summands (a cross-term shared intermediate is walked, and
  /// its FLOPs counted, once per occurrence, not once overall), and NOT
  /// replay-aware: each node is priced exactly once, so this never reflects the
  /// per-occ-block recompute the batched replay incurs. Compare against \c
  /// dryrun_flops: the two are ~equal when no batching engages, and \c
  /// dryrun_flops exceeds this by the recompute factor when it does.
  double model_flops = 0;
  /// STATIC per-node DP-MODEL roofline exec cost, summed over all internal
  /// nodes. Same per-occurrence (not CSE-deduplicated), batching-blind caveat
  /// as \c model_flops.
  double model_exec = 0;
  /// Number of internal (contraction) nodes across the forest (static count).
  std::size_t model_n_ops = 0;
  /// REPLAY-tallied (recompute-aware) FLOPs: the sum, over every ACTUAL product
  /// op executed in the \c Trace::On replay, of that op's SLICED-extent flops
  /// (folded via the \c CostSink attached to the shared \c CostModel; see
  /// result.hpp \c DryRunOps::prod). A batched op re-executed once per occ
  /// block is charged once per block at its sliced size, so occ-DEPENDENT work
  /// stays ~work-neutral while occ-INDEPENDENT recompute (persistent
  /// intermediates / leaf re-materializations re-run at full size per block)
  /// inflates -- making this the metric that PREDICTS batched-replay
  /// overcompute (e.g. occ-batching being ~Nx aux-only). Equal to \c
  /// model_flops (up to the product-op-only tally) when no batching engages.
  double dryrun_flops = 0;
  /// REPLAY-tallied roofline exec cost (traffic-dominated, the better wall-time
  /// proxy), summed per product-op execution. Same recompute semantics as \c
  /// dryrun_flops.
  double dryrun_exec = 0;
  /// Number of product-op EXECUTIONS in the replay (counts re-executions per
  /// batch block), so it grows with recompute -- unlike \c model_n_ops.
  std::size_t dryrun_n_ops = 0;

  /// Per-value avoidable-recompute breakdown: one entry per DISTINCT value
  /// whose batched replay repeated arithmetic beyond building it ONCE at full
  /// extent
  /// (`total_flops > full_flops`) -- the recompute a hoist would have avoided.
  /// Sorted by avoidable FLOPs descending. Empty when batching repeats no
  /// arithmetic (every value is built at most once-worth, e.g. disjoint slices
  /// tiling it). See \c DryRunOps::prod (per-build tally) and the post-replay
  /// rollup in \c cost_profile().
  std::vector<AvoidableNode> avoidable_nodes;
  /// Total avoidable recompute in FLOPs (sum of \c avoidable_nodes[i].flops):
  /// arithmetic the batched replay repeated beyond the build-once ideal.
  /// Compare against \c dryrun_flops for the avoidable FRACTION (see \c
  /// avoidable_time()). Zero when batching repeats no arithmetic.
  double avoidable_flops = 0;
  /// Total avoidable recompute expressed as equivalent extra full rebuilds (sum
  /// of \c avoidable_nodes[i].count).
  double avoidable_ops = 0;

  /// Avoidable FRACTION of replay arithmetic: \c avoidable_flops / \c
  /// dryrun_flops (0 when no ops ran), in [0, 1] by construction. The
  /// single-number "how much of the batched replay's arithmetic was repeated
  /// recompute vs. the unlimited-memory ideal" summary. (FLOPs, not roofline
  /// exec: recompute is repeated WORK, and FLOPs -- being linear in extents --
  /// makes disjoint slicing exactly free and keeps this bounded.)
  [[nodiscard]] double avoidable_time() const {
    return dryrun_flops > 0 ? avoidable_flops / dryrun_flops : 0.0;
  }
};

/// Replays a factorized eval forest zero-data through the real eval loop --
/// with a gated cache built from \p cfg (Task 2) and a \c PeakSink threaded
/// through the batched evaluator (Task 3) -- and, alongside, does a static walk
/// of the forest to accumulate FLOPs / roofline exec cost / op count. This is
/// the single reusable entry point both SeQuant tests and MPQC call.
///
/// \par The printing gate
/// \c CacheManager::working_set_hwmark() only accumulates while
/// \c sequant::eval::log::printing() is true (the hwmark update sits on the
/// trace-printing path). This routine therefore FORCES the eval logger's level
/// > 0 around the replay -- discarding the narrow trace to a null sink when no
/// \p trace is requested -- and restores the previous logger state afterward,
/// so \c peak_bytes is non-zero even with no trace stream.
///
/// \par Global state / threading
/// This routine mutates the process-global \c Logger::instance().eval state
/// (\c level and \c stream) for the duration of the replay (restored on every
/// exit path, including exceptions). Because that state is a singleton shared
/// by the whole process, \c cost_profile() MUST be called single-threaded --
/// e.g. as a pre-flight step before, or a post-hoc step after, the real
/// multi-threaded eval -- never concurrently with other code that reads or
/// writes \c Logger::instance().eval (including another concurrent
/// \c cost_profile() call).
///
/// \par FLOPs / exec accounting (model vs dryrun)
/// \c CostProfile::model_flops and \c CostProfile::model_exec are accumulated
/// by a STATIC walk that sums a contribution per BINARIZED internal node of the
/// forest; they are NOT CSE-aware across summands (a shared intermediate that
/// recurs across summand trees, or multiple times within one, is counted once
/// per occurrence) and NOT replay-aware (each node is priced exactly once,
/// blind to order/batching). \c CostProfile::dryrun_flops / \c dryrun_exec /
/// \c dryrun_n_ops are the recompute-aware counterparts, tallied from the
/// \c Trace::On replay below: every ACTUAL product-op execution folds its
/// SLICED-extent cost into a \c CostSink attached to the shared \c CostModel,
/// so an op re-executed once per batch block is counted once per block. When no
/// batching engages the two agree (up to the product-op-only dryrun tally);
/// when it does, \c dryrun_* exceeds \c model_* by the recompute factor -- a
/// TIME (arithmetic) cost that \c peak_bytes, a SPACE (co-resident working
/// set) measure, does not express.
///
/// \param forest per-summand optimized+binarized eval forest (the real IR).
/// \param policy the batch policy driving the replay evaluator; its accept is
///        the derived role union \c policy.is_batchable_index().
/// \param cfg    gated-cache config (footprint gate, volatile, repeats).
/// \param regime the size regime supplying extents and CSV moment tables;
///        the internal \c CostModel and \c DryRunLeafEvaluator are built from
///        it.
/// \param trace  optional per-op trace sink (nullptr = no trace). When
///        non-null, the eval loop's narrow trace is transcoded (UTF-8) into it.
/// \param router optional placement router (see \c placement_router.hpp),
///        attached to the replay cache right after it is built. Every current
///        caller omits this (nullptr, the default), which leaves the router
///        seam in \c evaluate() inert -- byte-identical to before this
///        parameter existed. Phase 2 wires a router only from test call sites.
/// \return the accumulated \c CostProfile.
inline CostProfile cost_profile(
    std::vector<EvalNodeDryRun> const& forest, BatchPolicy const& policy,
    CacheConfig const& cfg, SizeRegime const& regime,
    std::wostream* trace = nullptr,
    PlacementRouter<EvalNodeDryRun> const* router = nullptr,
    sequant::eval::ScheduleSink* schedule_sink = nullptr) {
  CostProfile profile;

  auto cm = std::make_shared<CostModel const>(regime);
  DryRunLeafEvaluator const leaf{cm};

  // ---- static cost walk (independent of the replay) --------------------
  // For every internal node: flops = flops_counter(left, right, result); the
  // roofline exec cost uses the left operand's footprint as the transferred
  // bytes and the arena convention (4096) the [dryrun-costmodel] test fixes.
  auto const flops_of = sequant::opt::detail::flops_counter(
      regime.idx_to_extent(), regime.inner_pow_fn());
  std::function<void(EvalNodeDryRun const&)> walk =
      [&](EvalNodeDryRun const& n) {
        if (n.leaf()) return;
        profile.model_n_ops += 1;
        double const node_flops =
            flops_of(n.left()->canon_indices(), n.right()->canon_indices(),
                     n->canon_indices());
        profile.model_flops += node_flops;
        container::svector<Index> const left(n.left()->canon_indices().begin(),
                                             n.left()->canon_indices().end());
        profile.model_exec +=
            cm->exec_cost(node_flops, cm->memsize(left), 4096);
        walk(n.left());
        walk(n.right());
      };
  for (auto const& root : forest) walk(root);

  // ---- peak replay through the real eval loop --------------------------
  // The cache's batch-variant veto is driven by the cross-occurrence lifetime
  // mask (stamped inside build_dryrun_cache -> cache_manager); the replay
  // EVALUATOR's accept is the derived role union, applied inside
  // make_evaluator(policy) via policy.is_batchable_index().
  auto cache = build_dryrun_cache(forest, cfg, regime);
  // Null (default) => every existing caller's replay is unaffected, since
  // set_placement_router(nullptr) is exactly the cache's own default.
  cache.set_placement_router(router);
  // Route this replay's SCHEDULE_RUN_EVENT records to the caller's sink (if
  // any). Null (default) => no dump, byte-identical. The wet-run sets an
  // equivalent sink on its own eval cache, so the two batched schedules can be
  // captured and diffed for structural equivalence.
  cache.set_schedule_sink(schedule_sink);

  auto& logger = Logger::instance();
  // RAII guard restoring the process-global Logger::eval state on EVERY exit
  // path from this point on -- normal return, early return, or an exception
  // unwinding out of the replay loop below -- not just the two trailing
  // assignments a plain save/restore would rely on. Without this, a throw
  // from anything in the loop OTHER than evaluate<Trace::On> (e.g.
  // std::bad_alloc from make_evaluator/set_custom_evaluator/
  // working_set_hwmark/cache.reset()) would unwind past the local
  // `trace_capture` destructor while `logger.eval.stream` still points at it,
  // leaving a dangling pointer in the process-global singleton with
  // level == 2 still set.
  struct LoggerEvalGuard {
    decltype(logger.eval)& eval;
    std::size_t const prev_level;
    std::ostream* const prev_stream;
    ~LoggerEvalGuard() {
      eval.level = prev_level;
      eval.stream = prev_stream;
    }
  } logger_eval_guard{logger.eval, logger.eval.level, logger.eval.stream};

  // Force printing() on so working_set_hwmark() accumulates. The eval logger
  // stream is narrow; capture into a narrow buffer only when a (wide) trace
  // sink was requested, else discard to a null stream.
  std::ostringstream trace_capture;
  logger.eval.level = 2;
  logger.eval.stream = trace ? &trace_capture : nullptr;

  // Attach a replay cost sink to the shared CostModel so each product op
  // executed in the Trace::On replay below folds its SLICED-extent cost here
  // (DryRunOps::prod -> CostModel::tally_op). Only DryRun Results built from
  // this same `cm` fold in, and only while the sink is attached, so this
  // records exactly the replay recompute for THIS forest. Detached on every
  // exit path by the guard below (the model outlives the results, but leaving a
  // dangling sink pointer set would be a latent hazard if `cm` were reused).
  CostSink costsink;
  cm->set_cost_sink(&costsink);
  struct CostSinkGuard {
    CostModel const& cm;
    ~CostSinkGuard() { cm.set_cost_sink(nullptr); }
  } costsink_guard{*cm};

  std::atomic<double> peak{0.0};
  for (auto const& root : forest) {
    cache.set_custom_evaluator(sequant::make_evaluator(
        policy, leaf, sequant::make_no_scope_guard{}, &peak));
    try {
      (void)sequant::evaluate<Trace::On>(root, leaf, cache);
    } catch (std::exception const&) {
      // A zero-data DryRun sizing throw must not mask the peak read.
    }
    // Fold the outer cached residency BEFORE reset() (which zeroes the
    // hwmark). `peak` folds every batched scratch high-watermark across all
    // summands via std::max, so its running load() is the global scratch peak.
    profile.peak_bytes = std::max(
        {profile.peak_bytes, peak.load(), double(cache.working_set_hwmark())});
    // NO per-term reset. A real solve's cache spans the whole iteration (all
    // summands + equations) and reuses cross-summand values, evicting only by
    // the lifetime mask stamped over the whole forest. Resetting between terms
    // would instead drop non-persistent scratch after each summand, REBUILDING
    // a cross-summand-shared value in every summand -- over-counting recompute
    // (and mis-estimating peak) vs any real run. So the replay keeps the shared
    // cache across the whole forest; the lifetime mask releases each value
    // after its last cross-term use. (This makes the dry-run schedule match the
    // wet run's; see doc/dev/specs/2026-08-05-...schedule-equivalence.)
  }

  // Read the replay-tallied (recompute-aware) totals the sink accumulated over
  // every product op of every summand's Trace::On replay. (costsink_guard
  // detaches the sink from `cm` on function exit.)
  profile.dryrun_flops = costsink.flops.load(std::memory_order_relaxed);
  profile.dryrun_exec = costsink.exec.load(std::memory_order_relaxed);
  profile.dryrun_n_ops = costsink.n_ops.load(std::memory_order_relaxed);

  // Per-value avoidable-recompute rollup (shared with the schedule-dump
  // emitter): for each DISTINCT value the replay built, avoidable FLOPs =
  // total_flops - full_flops, the arithmetic repeated beyond building it ONCE
  // at full extent (see avoidable_nodes_from_sink()).
  profile.avoidable_nodes = avoidable_nodes_from_sink(costsink);
  for (auto const& an : profile.avoidable_nodes) {
    profile.avoidable_flops += an.flops;
    profile.avoidable_ops += an.count;
  }

  // logger_eval_guard's destructor restores logger.eval.{level,stream} at
  // function exit (see above); no manual restore needed here.

  // If a wide trace sink was requested, transcode the captured narrow (UTF-8)
  // eval trace into it (the eval loop writes only to the narrow logger stream;
  // index labels such as mu~/K are multi-byte, so a plain widen would corrupt
  // them -- decode UTF-8 to code points instead).
  if (trace) {
    std::string const s = trace_capture.str();
    std::wstring w;
    w.reserve(s.size());
    for (std::size_t i = 0; i < s.size();) {
      unsigned char const c = static_cast<unsigned char>(s[i]);
      char32_t cp;
      std::size_t len;
      if (c < 0x80) {
        cp = c;
        len = 1;
      } else if ((c >> 5) == 0x6) {
        cp = c & 0x1Fu;
        len = 2;
      } else if ((c >> 4) == 0xE) {
        cp = c & 0x0Fu;
        len = 3;
      } else if ((c >> 3) == 0x1E) {
        cp = c & 0x07u;
        len = 4;
      } else {
        cp = c;  // invalid lead byte: pass through
        len = 1;
      }
      for (std::size_t k = 1; k < len && i + k < s.size(); ++k)
        cp = (cp << 6) | (static_cast<unsigned char>(s[i + k]) & 0x3Fu);
      w.push_back(static_cast<wchar_t>(cp));
      i += len;
    }
    *trace << w;
  }

  return profile;
}

}  // namespace sequant::eval::dryrun

#endif  // SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_COST_PROFILE_HPP
