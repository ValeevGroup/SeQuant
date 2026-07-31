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
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/optimize/single_term_detail.hpp>

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <functional>
#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace sequant::eval::dryrun {

/// Configuration for a faithful (gated) dry-run cache: the same footprint gate
/// and free-batchable-axis veto the real batched eval loop applies, so a
/// free-batchable giant (a `mu~`/`K`-carrying DF intermediate) is NOT cached
/// whole but recomputed sliced under each consumer's batch trigger.
///
/// The element types mirror the gated \c sequant::cache_manager overload
/// (\c cache_manager.hpp): \c is_volatile is invoked on every \c TreeNode
/// (deduced as \c EvalNodeDryRun for the dry-run backend), \c
/// is_batchable_index on every result \c Index.
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
  /// `bool(Index const&)`: true for an index the runtime batched evaluator
  /// slices over (e.g. DF aux `K` / PAO `mu~`). A node whose result carries
  /// such a free index is vetoed from caching. Empty => nothing is batchable.
  ///
  /// ADVISORY when passed to \c cost_profile(): that entry point OVERWRITES
  /// this field with \c policy.is_batchable_contracted_index (the
  /// CONTRACTED-role building block) before building the cache. The cache veto
  /// is contracted-stamp-only, so it is decoupled from the replay evaluator's
  /// accept (the derived role union). Only \c build_dryrun_cache() called
  /// directly honors this field as-is.
  std::function<bool(Index const&)> is_batchable_index;
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
/// giants are vetoed (matching the real run). Call \c CacheManager::reset() on
/// the returned cache between summands to drop per-term non-persistent scratch
/// while keeping persistent (cross-term) entries.
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

  // Default the predicates so the gated factory never invokes an empty
  // std::function (nothing volatile / nothing batchable leaves those gates
  // inert, matching the factory's own defaults).
  std::function<bool(EvalNodeDryRun const&)> is_volatile =
      cfg.is_volatile ? cfg.is_volatile
                      : std::function<bool(EvalNodeDryRun const&)>(
                            [](EvalNodeDryRun const&) { return false; });
  std::function<bool(Index const&)> is_batchable_index =
      cfg.is_batchable_index ? cfg.is_batchable_index
                             : std::function<bool(Index const&)>(
                                   [](Index const&) { return false; });

  // Note: the gated sequant::cache_manager() overload below stamps the
  // cross-occurrence lifetime mask on `nodes` itself before its DAG walk /
  // veto (cache_manager.hpp), so this call site does not need to do so.
  return sequant::cache_manager(nodes, std::move(is_volatile), cfg.min_repeats,
                                std::move(footprint_of), cfg.max_footprint,
                                std::move(is_batchable_index));
}

/// One recomputed node's avoidable-recompute breakdown (see
/// \c CostProfile::avoidable_nodes). \c label is the op signature (a
/// `result;lhs;rhs` full-label string built in \c DryRunOps::prod). \c count is
/// `builds - necessary`: how many of the node's actual replay builds were
/// redundant given its TOUCHED (dependent) modes -- a node rebuilt once per
/// block of a mode it does not touch is avoidable recompute (a missed hoist).
/// \c exec / \c flops are that count times the per-build roofline exec / FLOPs,
/// i.e. the wall-time-proxy and arithmetic cost of the redundancy. Only nodes
/// with \c count > 0 are recorded.
struct AvoidableNode {
  std::string label;
  double count = 0;
  double exec = 0;
  double flops = 0;
};

/// Roll a replay's per-node build tally (a \c CostSink populated by
/// \c CostModel::tally_node during a \c Trace::On replay) into the avoidable-
/// recompute breakdown: per distinct op signature, \c avoidable = builds -
/// necessary (clamped >= 0), weighted by the per-build exec / flops, sorted by
/// exec descending, keeping only nodes with a positive avoidable count. Shared
/// by \c cost_profile() (whole-forest rollup) and any caller that attaches its
/// OWN sink to a replay it drives (e.g. the schedule-dump test, which emits
/// these per-node numbers so the visualizer consumes cost_profile's avoidable
/// instead of recomputing it). Single-threaded caller expected: the replay has
/// finished and the sink is quiescent, so no lock is taken.
inline std::vector<AvoidableNode> avoidable_nodes_from_sink(
    CostSink const& sink) {
  std::vector<AvoidableNode> out;
  for (auto const& [sig, nc] : sink.per_node) {
    double const builds = static_cast<double>(nc.builds);
    double const avoidable =
        builds > nc.necessary ? builds - nc.necessary : 0.0;
    if (avoidable <= 0.0) continue;
    out.push_back({sig, avoidable, avoidable * nc.exec_per_build,
                   avoidable * nc.flops_per_build});
  }
  std::sort(out.begin(), out.end(),
            [](AvoidableNode const& a, AvoidableNode const& b) {
              return a.exec > b.exec;
            });
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
  /// This is computed as \c max(batched-inner scratch high-watermark,
  /// outer cross-term cached residency), NOT their sum. When a persistent
  /// cross-term cache entry co-resides in memory with a batched-inner
  /// transient at the same instant, the two are actually additive, so this
  /// value is a LOWER BOUND on the true peak in that case. It is exact
  /// whenever one of the two terms dominates the other (e.g. the C60 4-PNO
  /// \c W case, where the batched-inner scratch dwarfs any persistent
  /// cross-term residency).
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

  /// Per-node avoidable-recompute breakdown: one entry per DISTINCT op whose
  /// replay build count exceeded the \c necessary count implied by its TOUCHED
  /// (dependent) modes -- a node rebuilt once per block of a mode it does not
  /// touch (a missed hoist). Sorted by \c exec descending. Empty when no
  /// batching engages (every node then builds exactly as its touched-mode
  /// extent demands). See \c DryRunOps::prod (per-build tally) and the
  /// post-replay rollup in \c cost_profile().
  std::vector<AvoidableNode> avoidable_nodes;
  /// Total avoidable roofline exec across all recomputed nodes (sum of
  /// \c avoidable_nodes[i].exec) -- the wall-time-proxy counterpart of the raw
  /// avoidable build COUNT. Compare against \c dryrun_exec for the avoidable
  /// FRACTION (see \c avoidable_time()). Zero when no batching engages.
  double avoidable_exec = 0;
  /// Total avoidable build count (sum of \c avoidable_nodes[i].count): how many
  /// product-op executions across the whole replay were redundant recompute.
  double avoidable_ops = 0;

  /// Avoidable FRACTION of replay exec: \c avoidable_exec / \c dryrun_exec (0
  /// when no ops ran). The single-number "how much of the batched replay's work
  /// was wasted recompute" summary.
  [[nodiscard]] double avoidable_time() const {
    return dryrun_exec > 0 ? avoidable_exec / dryrun_exec : 0.0;
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
/// when it does, \c dryrun_* exceeds \c model_* by the recompute factor --
/// which is what \c peak_bytes (max, not sum) cannot express.
///
/// \param forest per-summand optimized+binarized eval forest (the real IR).
/// \param policy the batch policy driving the replay evaluator; its
///        \c is_batchable_contracted_index is COPIED over
///        \c cfg.is_batchable_index internally (the cache veto is
///        contracted-stamp-only), while the evaluator accept is the derived
///        role union \c policy.is_batchable_index(). The \c cfg field is
///        advisory here.
/// \param cfg    gated-cache config (footprint gate, axis veto, volatile,
///        repeats).
/// \param regime the size regime supplying extents and CSV moment tables;
///        the internal \c CostModel and \c DryRunLeafEvaluator are built from
///        it.
/// \param trace  optional per-op trace sink (nullptr = no trace). When
///        non-null, the eval loop's narrow trace is transcoded (UTF-8) into it.
/// \return the accumulated \c CostProfile.
inline CostProfile cost_profile(std::vector<EvalNodeDryRun> const& forest,
                                BatchPolicy const& policy,
                                CacheConfig const& cfg,
                                SizeRegime const& regime,
                                std::wostream* trace = nullptr) {
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
  // Route the two consumers by intent (they are decoupled):
  //   - the cache VETO takes the CONTRACTED-role building block: only a mode
  //     actually sliced AT a node (BatchModeType::Contracted) vetoes caching,
  //     so the veto must never see the union (an external-only mode is
  //     batch-invariant and stays cacheable).
  //   - the replay EVALUATOR's accept is the derived role union, applied inside
  //     make_evaluator(policy) via policy.is_batchable_index().
  CacheConfig local_cfg = cfg;
  local_cfg.is_batchable_index = policy.is_batchable_contracted_index;
  auto cache = build_dryrun_cache(forest, local_cfg, regime);

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
    cache.reset();  // drop per-term non-persistent scratch; keep persistent
  }

  // Read the replay-tallied (recompute-aware) totals the sink accumulated over
  // every product op of every summand's Trace::On replay. (costsink_guard
  // detaches the sink from `cm` on function exit.)
  profile.dryrun_flops = costsink.flops.load(std::memory_order_relaxed);
  profile.dryrun_exec = costsink.exec.load(std::memory_order_relaxed);
  profile.dryrun_n_ops = costsink.n_ops.load(std::memory_order_relaxed);

  // Per-node avoidable-recompute rollup (shared with the schedule-dump
  // emitter). For each DISTINCT op signature the replay built, `builds` is how
  // many times it actually ran and `necessary` is the product of block counts
  // of its TOUCHED (dependent) modes -- the minimum a perfect-hoisting
  // evaluator would do; the excess is avoidable recompute (a node rebuilt once
  // per block of a mode it does NOT touch). Because the cost model is DENSE
  // (every touched block visited exactly once, no block sparsity), `necessary`
  // is exact, so avoidable = builds - necessary needs no empirical correction.
  profile.avoidable_nodes = avoidable_nodes_from_sink(costsink);
  for (auto const& an : profile.avoidable_nodes) {
    profile.avoidable_exec += an.exec;
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
