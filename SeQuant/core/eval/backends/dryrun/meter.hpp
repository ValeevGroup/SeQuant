#ifndef SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_METER_HPP
#define SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_METER_HPP

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_profile.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/peak_monitor.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/eval/placement_router.hpp>
#include <SeQuant/core/eval/scope_executor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <memory>
#include <ostream>
#include <ranges>
#include <string>
#include <unordered_map>
#include <vector>

namespace sequant::eval::dryrun {

///
/// \brief One value's build-vs-home fidelity witness for a \c MeterReport
/// (see \c assemble_report): how many times a distinct value was built (over
/// its whole recompute tally) versus WHERE it is homed and WHERE it is used,
/// read off the matching \c RichSchedule::ValueCell (looked up by hash).
///
struct HomeFidelity {
  std::string label;       ///< value signature "idx:hash" (idx disambiguates a
                           ///< 64-bit hash collision between distinct nodes;
                           ///< see cost_profile.hpp's sig_full_flops)
  std::size_t hash = 0;    ///< the value's EvalExpr::hash_value()
  std::size_t builds = 0;  ///< total builds across slices (recompute-aware)
  std::string home;        ///< dag-scope of home_modes ("" == root {} -- a
                           ///< whole-nest invariant)
  std::string uses;        ///< dag-scope list of the value's occurrences
};

///
/// \brief Summary of one metered dry-run (or wet) replay: the hierarchy-wide
/// peak (from a wired \c PeakMonitor), the persistent/volatile FLOPs and
/// CostModel exec-time split (rolled up from a \c CacheManager's recompute
/// tally, classified by \c compute_volatility), the total build count, and
/// the per-value build-vs-home fidelity list (see \c HomeFidelity).
///
struct MeterReport {
  double peak_bytes = 0;         ///< PeakMonitor high-water (dense on dry
                                 ///< run / sparse on wet run)
  std::size_t peak_op_hash = 0;  ///< location (op hash) of the peak

  double flops_persistent = 0, flops_volatile = 0;  ///< dense model; dry-only
  double cost_persistent = 0, cost_volatile = 0;  ///< CostModel exec; dry-only

  std::size_t builds_total = 0;

  std::vector<HomeFidelity> home_fidelity;  ///< sorted: builds desc

  BatchScheduler scheduler =
      BatchScheduler::forest_descent;  ///< which executor this report
                                       ///< describes
};

///
/// \brief Bottom-up memoized volatility over an evaluation \p forest: a node
/// is volatile iff \p is_volatile flags it directly, or (for an internal
/// node) either child is volatile -- the SAME rule the gated
/// \c sequant::cache_manager factory applies while building its NV/V
/// frontier (see cache_manager.hpp's DAG walk). Keyed by
/// \c TreeNode::hash_value() (rather than the node identity itself) so a
/// caller can classify a \c CacheManager::recompute_tally() entry -- keyed by
/// the SAME node identity but not necessarily the SAME node object -- by its
/// hash.
///
/// \param forest the evaluation forest (a range of eval nodes).
/// \param is_volatile `bool(TreeNode const&)`: true if the node is
///        intrinsically volatile. Only its value on leaves matters in
///        practice (volatility propagates up), but it is consulted on every
///        node, matching \c cache_manager's gated factory.
/// \return a map from \c hash_value() to whether that value is volatile.
///
template <class Forest, class IsVolatile>
std::unordered_map<std::size_t, bool> compute_volatility(
    Forest const& forest, IsVolatile const& is_volatile) {
  using Node = std::ranges::range_value_t<Forest>;

  std::unordered_map<std::size_t, bool> volatile_of;

  auto visit = [&](auto&& self, Node const& n) -> bool {
    std::size_t const h = n->hash_value();
    if (auto it = volatile_of.find(h); it != volatile_of.end())
      return it->second;
    bool v;
    if (n.leaf()) {
      v = is_volatile(n);
    } else {
      bool const vl = self(self, n.left());
      bool const vr = self(self, n.right());
      v = is_volatile(n) || vl || vr;
    }
    volatile_of.emplace(h, v);
    return v;
  };

  for (auto const& tree : forest) visit(visit, tree);
  return volatile_of;
}

///
/// \brief Assemble a \c MeterReport from a metered replay: a walked
/// \p cache (its \c recompute_tally() populated by \c CacheManager::
/// tally_build over the replay), the hierarchy-wide \p mon (\c PeakMonitor),
/// the \p rich linearized schedule (\c compute_dag_boulevard over the SAME
/// \p forest, supplying each value's home/use dag-scope), and \p is_volatile
/// (fed to \c compute_volatility to classify each distinct value).
///
/// Per distinct value (one \c cache.recompute_tally() entry): \c builds is
/// the sum, over its slices, of each slice's build count; \c node_flops /
/// \c node_exec are the sum, over its slices, of build-count times that
/// slice's actual (flops, exec). The value is classified persistent/volatile
/// by \c compute_volatility's verdict for its hash and folded into the
/// matching \c MeterReport::flops_*/cost_* accumulator. Its \c HomeFidelity
/// entry's \c home/uses are read off the \p rich cell sharing its hash (empty
/// if the value has no matching cell, e.g. a leaf never realized as its own
/// distinct product build).
///
/// \param cache the (root) cache whose \c recompute_tally() was populated by
///        a \c Trace::On metered replay.
/// \param mon the \c PeakMonitor wired onto \p cache's scope chain during the
///        replay.
/// \param rich the linearized schedule (\c compute_dag_boulevard) over the
///        SAME forest the replay walked.
/// \param forest the evaluation forest (fed to \c compute_volatility).
/// \param is_volatile `bool(TreeNode const&)`: intrinsic volatility
///        predicate, as for \c compute_volatility.
/// \param scheduler which executor this report describes (stashed verbatim
///        into \c MeterReport::scheduler).
/// \return the assembled \c MeterReport.
///
template <class Cache, class Forest, class IsVolatile>
MeterReport assemble_report(Cache const& cache, PeakMonitor const& mon,
                            RichSchedule const& rich, Forest const& forest,
                            IsVolatile const& is_volatile,
                            BatchScheduler scheduler) {
  MeterReport report;
  report.scheduler = scheduler;
  report.peak_bytes = static_cast<double>(mon.hwmark_bytes);
  report.peak_op_hash = mon.peak.op_hash;

  auto const volatility = compute_volatility(forest, is_volatile);

  // hash -> ValueCell* lookup, mirroring make_node_meta's map build
  // (scope_executor.hpp): rich.cells is a flat vector, not keyed by hash.
  std::unordered_map<std::size_t, ValueCell const*> cell_by_hash;
  cell_by_hash.reserve(rich.cells.size());
  for (auto const& cell : rich.cells) cell_by_hash.emplace(cell.hash, &cell);

  // dag-scope formatting: comma-joined IndexSpace base_keys, no trailing
  // comma -- the same convention make_node_meta uses (scope_executor.hpp).
  auto const dag_scope = [](auto const& modes) {
    std::string s;
    for (auto const& m : modes) {
      if (!s.empty()) s += ",";
      s += toUtf8(m.space().base_key());
    }
    return s;
  };

  std::size_t idx = 0;
  for (auto const& [node, tally] : cache.recompute_tally()) {
    std::size_t builds = 0;
    double node_flops = 0.0, node_exec = 0.0;
    for (auto const& [sig, rec] : tally.slices) {
      builds += rec.count;
      node_flops += static_cast<double>(rec.count) * rec.flops;
      node_exec += static_cast<double>(rec.count) * rec.exec;
    }
    report.builds_total += builds;

    std::size_t const hash = node->hash_value();
    bool const is_vol = [&] {
      auto it = volatility.find(hash);
      return it != volatility.end() && it->second;
    }();

    if (is_vol) {
      report.flops_volatile += node_flops;
      report.cost_volatile += node_exec;
    } else {
      report.flops_persistent += node_flops;
      report.cost_persistent += node_exec;
    }

    HomeFidelity hf;
    hf.label = std::to_string(idx++) + ":" + std::to_string(hash);
    hf.hash = hash;
    hf.builds = builds;
    if (auto it = cell_by_hash.find(hash); it != cell_by_hash.end()) {
      auto const* cell = it->second;
      hf.home = dag_scope(cell->home_modes);
      container::svector<Index> uses_modes;
      for (auto const& occ : cell->occurrences)
        for (auto const& [mode, range] : occ.ectx) uses_modes.push_back(mode);
      hf.uses = dag_scope(uses_modes);
    }
    report.home_fidelity.push_back(std::move(hf));
  }

  std::sort(report.home_fidelity.begin(), report.home_fidelity.end(),
            [](HomeFidelity const& a, HomeFidelity const& b) {
              return a.builds > b.builds;
            });

  return report;
}

///
/// \brief Runs the policy-selected executor (forest descent, whole-scope, or
/// ordered, per \p policy.scheduler) over \p forest through the
/// DryRun sizing backend, metering the replay with a fresh, \c PeakMonitor
/// -wired, build-tallying cache, and returns the assembled \c MeterReport.
///
/// Mirrors MPQC's wet dispatch: this drives the SAME Task-6 coexistence entry
/// point (\c sequant::evaluate(Nodes const&, BatchPolicy const&, layout, F,
/// CacheManager&, mode_order, ScopeGuardFactory), \c scope_executor.hpp) a
/// real solve would use under \p policy -- all three executors are selected
/// by the SAME \c policy.scheduler, not independently maintained code paths --
/// so the metered replay is exactly the run \p policy describes, not a
/// hand-rolled proxy of it. Non-throwing wrapper (if desired) is the
/// caller's responsibility; an exception from the replay propagates out of
/// this call, but the RAII logger-state guard still restores
/// \c Logger::instance().eval on the way out.
///
/// \param forest the evaluation forest (a range of \c EvalNodeDryRun).
/// \param policy the batch policy driving the coexistence entry -- in
///        particular \c scheduler (executor selection) and
///        \c batch_target_size (the batch-partition source; also the source
///        of the \c block_of function this call builds its OWN \c rich
///        schedule with, for \c assemble_report -- the coexistence entry
///        builds an independent, internal \c RichSchedule of its own from
///        the SAME \p policy.batch_target_size to drive the executor).
/// \param regime the size regime supplying the DryRun \c CostModel.
/// \param cfg cache configuration (footprint gate, min repeats, volatility)
///        for the metered cache, built exactly as \c build_dryrun_cache does
///        (same footprint arithmetic, same is_volatile default) -- NOT via
///        that builder directly, since its is_volatile default (substituted
///        for an empty \c cfg.is_volatile) is internal to it and would
///        otherwise be invisible to \c assemble_report below, which also
///        needs a callable predicate (an empty \c cfg.is_volatile passed to
///        it directly throws \c std::bad_function_call from
///        \c compute_volatility). The SAME locally-defaulted predicate is
///        used for both.
/// \param router optional placement override, installed on the metered cache
///        when non-null.
/// \param trace optional sink for the eval trace; when non-null,
///        \c Logger::instance().eval.stream is redirected there for the
///        duration of the call (restored on exit, along with the elevated
///        \c eval.level and the installed \c eval.node_meta).
/// \return the assembled \c MeterReport (peak, persistent/volatile
///         FLOPs+time, build-vs-home fidelity), stamped with
///         \p policy.scheduler.
///
inline MeterReport meter(
    std::vector<EvalNodeDryRun> const& forest, BatchPolicy const& policy,
    SizeRegime const& regime, CacheConfig const& cfg,
    PlacementRouter<EvalNodeDryRun> const* router = nullptr,
    std::ostream* trace = nullptr) {
  auto cm = std::make_shared<CostModel const>(regime);
  DryRunLeafEvaluator yield{cm};

  // Default is_volatile the SAME way build_dryrun_cache does (an empty
  // cfg.is_volatile means nothing is volatile) -- but keep the defaulted
  // function LOCAL rather than routing through that builder, so the exact
  // same predicate can also be threaded to assemble_report below.
  std::function<bool(EvalNodeDryRun const&)> const is_volatile =
      cfg.is_volatile ? cfg.is_volatile
                      : std::function<bool(EvalNodeDryRun const&)>(
                            [](EvalNodeDryRun const&) { return false; });

  // Footprint (bytes) of a node's RESULT, identical to build_dryrun_cache's
  // footprint_of (cost_profile.hpp): the moment-aware memsize counter over
  // canon_indices(), scaled to bytes, so cfg.max_footprint gates like-for-like.
  auto memsize = sequant::opt::detail::memsize_counter(regime.idx_to_extent(),
                                                       regime.inner_pow_fn());
  auto footprint_of =
      [memsize = std::move(memsize)](EvalNodeDryRun const& n) -> double {
    std::vector<Index> const result(n->canon_indices().begin(),
                                    n->canon_indices().end());
    return memsize(std::vector<Index>{}, std::vector<Index>{}, result) * 8.0;
  };

  auto cache =
      sequant::cache_manager(forest, is_volatile, cfg.min_repeats,
                             std::move(footprint_of), cfg.max_footprint);
  cache.set_recompute_tally_enabled(true);
  if (router) cache.set_placement_router(router);

  PeakMonitor mon;
  cache.set_peak_monitor(&mon);

  // Backend array-ops for the dry-run backend: the batched executors build a
  // scatter destination / enumerate axis batches through this seam (the SAME
  // one the wet run uses), so the dry replay realizes the same scatter
  // footprint and batch count. Sourced from the shared CostModel -- no array in
  // the DAG is consulted. Must outlive the replay below (it is a local here).
  auto const dry_aops = make_dryrun_array_ops(cm);
  cache.set_array_ops(&dry_aops);

  // The SAME block_of source the coexistence entry itself derives from
  // policy.batch_target_size (scope_executor.hpp's evaluate(Nodes const&,
  // BatchPolicy const&, ...)) -- an empty batch_target_size means "no
  // batching", guarded identically so compute_dag_boulevard never invokes an
  // empty std::function.
  std::function<std::size_t(Index const&)> const block_of =
      policy.batch_target_size
          ? policy.batch_target_size
          : std::function<std::size_t(Index const&)>(
                [](Index const&) -> std::size_t { return 1; });
  RichSchedule const rich = compute_dag_boulevard(forest, *cm, block_of);

  // RAII save/restore of every Logger::eval field this call touches, so an
  // exception from the replay below still leaves the process-wide Singleton
  // exactly as this call found it.
  auto& logger = Logger::instance();
  struct LoggerStateGuard {
    Logger& l;
    std::size_t prev_level;
    std::ostream* prev_stream;
    std::function<std::string(std::size_t)> prev_node_meta;
    ~LoggerStateGuard() {
      l.eval.level = prev_level;
      l.eval.stream = prev_stream;
      l.eval.node_meta = std::move(prev_node_meta);
    }
  } guard{logger, logger.eval.level, logger.eval.stream, logger.eval.node_meta};

  // Ensure printing() so DryRunOps::prod records flops/exec (feeding
  // cache.tally_build) and note_working_set() actually observes the
  // PeakMonitor -- without raising the level any HIGHER than a caller who
  // already wants a louder trace.
  logger.eval.level = std::max<std::size_t>(logger.eval.level, 1);
  if (trace) logger.eval.stream = trace;
  logger.eval.node_meta = make_node_meta(rich);

  // Forest descent (BatchScheduler::forest_descent) needs the SAME batched
  // custom evaluator MPQC's wet forest path installs (cck.ipp's `else`
  // branch, `cache.set_custom_evaluator(sequant::make_evaluator(ctx.
  // batch_policy, yielder, make_scope_guard))`): without it, plain
  // sequant::evaluate(Nodes const&, ...) ignores every node_slice_mask()
  // stamp and runs an unbatched, no-schedule single pass -- an infidelity
  // vs. the wet run this meter is supposed to mirror. Installed ONLY on
  // this branch: evaluate_impl consults cache.custom_evaluator() on every
  // non-leaf node, so installing it unconditionally would also fire on the
  // whole-scope AND ordered paths' own evaluate_impl calls too. Both
  // whole_scope (evaluate_whole_scope) and ordered (evaluate_ordered_
  // schedule) drive their OWN executor via the coexistence entry
  // (sequant::evaluate(forest, policy, ...) below, which dispatches on
  // policy.scheduler) -- installing the forest custom evaluator for ordered
  // would silently reroute its builds through the forest evaluator instead
  // of run_ordered_contracted_block, diverging from what the wet ordered run
  // (which installs NO custom evaluator) actually does. That divergence was
  // a real dry-run/wet-run fidelity bug; restricting this install to forest
  // descent fixes it.
  if (policy.scheduler == BatchScheduler::forest_descent)
    cache.set_custom_evaluator(
        sequant::make_evaluator(policy, yield, sequant::make_no_scope_guard{}));

  (void)sequant::evaluate<Trace::On>(forest, policy, std::wstring{}, yield,
                                     cache, {}, sequant::make_no_scope_guard{});

  // INSTRUMENTATION (SEQUANT_RECOMPUTE_DUMP, analysis-only): genuine,
  // batching-aware recompute is a SINGLE (node, slice) built more than once
  // (a value tiled over distinct batch slices has count==1 per slice -- see
  // BuildRecord). Dump every such slice with its hash so the recompute can be
  // localized without the batch-slice / accumulation confound.
  if (std::getenv("SEQUANT_RECOMPUTE_DUMP")) {
    std::size_t genuine = 0, slices_gt1 = 0;
    for (auto const& [node, bt] : cache.recompute_tally())
      for (auto const& [slice, rec] : bt.slices)
        if (rec.count > 1) {
          ++slices_gt1;
          genuine += rec.count - 1;
          std::cerr << "RECOMPUTE count=" << rec.count
                    << " hash=" << node->hash_value() << " slice=[" << slice
                    << "]\n";
        }
    std::cerr << "GENUINE-RECOMPUTE: distinct (node,slice) built >1 = "
              << slices_gt1 << " ; total avoidable rebuilds = " << genuine
              << "\n";
  }

  return assemble_report(cache, mon, rich, forest, is_volatile,
                         policy.scheduler);
}

}  // namespace sequant::eval::dryrun

#endif  // SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_METER_HPP
