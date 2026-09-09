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
  /// Summed unweighted static contraction FLOPs over all internal nodes.
  /// NOT CSE-aware across summands: a cross-term shared intermediate is
  /// walked (and its FLOPs counted) once per occurrence, not once overall.
  double flops = 0;
  /// Summed roofline-projected execution cost over all internal nodes. Same
  /// per-occurrence (not CSE-deduplicated) accounting caveat as \c flops.
  double exec_cost = 0;
  /// Number of internal (contraction) nodes across the forest.
  std::size_t n_ops = 0;
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
/// \par FLOPs / exec_cost accounting
/// \c CostProfile::flops and \c CostProfile::exec_cost are accumulated by a
/// static walk that sums a contribution per BINARIZED internal node of the
/// forest; they are NOT CSE-aware across summands. A shared intermediate
/// that recurs across multiple summand trees (or multiple times within one)
/// is counted once per occurrence, not once overall -- unlike \c peak_bytes,
/// which is driven by the gated cache and so does reflect cross-term reuse.
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
        profile.n_ops += 1;
        double const node_flops =
            flops_of(n.left()->canon_indices(), n.right()->canon_indices(),
                     n->canon_indices());
        profile.flops += node_flops;
        container::svector<Index> const left(n.left()->canon_indices().begin(),
                                             n.left()->canon_indices().end());
        profile.exec_cost += cm->exec_cost(node_flops, cm->memsize(left), 4096);
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
