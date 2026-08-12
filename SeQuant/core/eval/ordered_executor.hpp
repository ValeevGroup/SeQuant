#ifndef SEQUANT_EVAL_ORDERED_EXECUTOR_HPP
#define SEQUANT_EVAL_ORDERED_EXECUTOR_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/forest_combine.hpp>
#include <SeQuant/core/eval/ordered_schedule.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/eval/value_node_map.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <cstddef>
#include <functional>
#include <type_traits>
#include <unordered_map>
#include <variant>

namespace sequant::eval {

///
/// \brief SP3 Task 1 of the ordered-scope batched-eval design: the ORDERED
/// executor skeleton. Walks \c ordered.root.steps in sequence, building one
/// value per \c BuildStep via \c evaluate_impl (which cache-probes every
/// operand at its own \c Stage::Enter, so a value already \c cache.store'd
/// by an earlier step short-circuits its own re-descent -- see \c
/// evaluate_impl's doc comment), then combines the forest roots' results
/// into the final \c ResultPtr exactly as \c evaluate_whole_scope's shared
/// root-combine loop does (permute-to-\p layout + cross-root \c add_inplace,
/// with the identical Term/Permute/SumInplace trace bookkeeping).
///
/// \details Task 1 handles ONLY a schedule whose \c ordered.root contains no
/// nested \c ScopeBlock step (i.e. no batch loop was realized -- every value
/// is a plain root-level \c BuildStep, the shape \c build_ordered_schedule
/// produces when no index is batchable under the \c BatchPolicy that built
/// \p rich / \p legality). A root-level child \c ScopeBlock step -- SP3
/// Task 2's job -- trips a loud \c SEQUANT_ASSERT rather than silently
/// mis-evaluating: no caller of the Task-1 gating flag exercises batching
/// yet (see \c BatchPolicy::ordered_schedule_execution), so this is an inert
/// tripwire today, not a live restriction.
///
/// \param forest The forest whose per-root results are evaluated and summed;
///        same requirement as \c evaluate_whole_scope's \p forest.
/// \param ordered The \c OrderedSchedule (\c build_ordered_schedule) for
///        \p forest, built from \p rich / a \c LegalitySchedule over it.
/// \param rich The \c RichSchedule that produced \p ordered (\c
///        compute_dag_boulevard) -- used to resolve a \c BuildStep's
///        \c value_id back to a forest node via its \c ValueCell::hash and
///        \c build_value_node_map, and to resolve each forest root's own
///        \c value_id for the final combine.
/// \param layout The layout each root's result is permuted to before being
///        summed; same meaning as \c evaluate_whole_scope's \p layout.
/// \param leaf_evaluator The leaf evaluator, as in \c sequant::evaluate.
/// \param cache The cache for common sub-expression elimination, as in
///        \c sequant::evaluate; a value repeated across \c BuildStep's
///        operand needs is deduped exactly as \c evaluate_impl's Checked
///        cache-probe already does for forest descent.
/// \param target Per-index batch partition size (elements); unused by the
///        Task-1 (unbatched) walk, threaded for interface symmetry with
///        later tasks that recurse into a batch loop. Kept as a named
///        parameter (not `[[maybe_unused]]`-dropped) so Task 2 can extend
///        this signature's BODY without an interface break.
/// \param make_scope_guard Backend scope-guard factory; unused by the
///        Task-1 walk (no batch loop is realized here), threaded for the
///        same interface-symmetry reason as \p target.
/// \return The summed, per-root-permuted result, as \c evaluate_whole_scope
///         (and, for an unbatched forest, forest descent itself) would
///         produce for the same \p forest.
///
template <Trace EvalTrace = Trace::Default, meta::can_evaluate_range Nodes,
          typename F, typename N, bool FHC,
          typename ScopeGuardFactory = ::sequant::make_no_scope_guard>
  requires meta::leaf_node_evaluator<std::ranges::range_value_t<Nodes>, F>
ResultPtr evaluate_ordered_schedule(
    Nodes const& forest, OrderedSchedule const& ordered,
    RichSchedule const& rich, auto const& layout, F const& leaf_evaluator,
    CacheManager<N, FHC>& cache,
    [[maybe_unused]] std::function<std::size_t(Index const&)> const& target,
    [[maybe_unused]] ScopeGuardFactory const& make_scope_guard = {}) {
  using node_t = std::ranges::range_value_t<Nodes>;
  static_assert(std::is_same_v<node_t, N>,
                "the forest's node type and the cache's node type must match");

  // hash -> node, resolving a BuildStep's value_id (via rich.cells[vid].hash)
  // to the forest node evaluate_impl builds.
  auto const vmap = build_value_node_map(forest);

  // Per-value results, indexed by value_id (== a ValueCell's own slot in
  // rich.cells -- see peak_profile.hpp's ValueCell::value_id doc comment),
  // so the root-combine step below reads a forest root's own build directly
  // rather than re-resolving it through the cache.
  container::vector<ResultPtr> value_results(ordered.num_values);

  // -------- Task 1: walk the root block's own BuildSteps, in order. --------
  // A root-level child ScopeBlock (a realized batch loop) is SP3 Task 2's
  // job; see the function doc comment for why this is a loud tripwire, not a
  // silent gap, on the schedules this task's gating flag can produce today.
  for (Step const& step : ordered.root.steps) {
    if (auto const* build = std::get_if<BuildStep>(&step.value)) {
      std::size_t const vid = build->value_id;
      SEQUANT_ASSERT(vid < rich.cells.size());
      auto const hash = rich.cells[vid].hash;
      auto const it = vmap.find(hash);
      SEQUANT_ASSERT(it != vmap.end() &&
                     "evaluate_ordered_schedule: BuildStep value not found "
                     "in the forest's value-node map");
      value_results[vid] =
          evaluate_impl<EvalTrace>(it->second, leaf_evaluator, cache);
    } else {
      SEQUANT_ASSERT(
          false &&
          "evaluate_ordered_schedule: a root-level child ScopeBlock (a "
          "realized batch loop) is not yet supported -- SP3 Task 1 handles "
          "only root BuildSteps; Task 2 adds loop-block interpretation.");
    }
  }

  // hash -> value_id, to resolve each forest root's own build above into the
  // shared combine step's per-root pre_results below.
  std::unordered_map<std::size_t, std::size_t> hash_to_vid;
  hash_to_vid.reserve(rich.cells.size());
  for (auto const& c : rich.cells) hash_to_vid.emplace(c.hash, c.value_id);

  container::svector<node_t> roots;
  for (auto&& n : forest) roots.push_back(n);

  container::svector<ResultPtr> pre_results(roots.size());
  for (std::size_t i = 0; i != roots.size(); ++i) {
    auto const vid_it = hash_to_vid.find(roots[i]->hash_value());
    SEQUANT_ASSERT(vid_it != hash_to_vid.end() &&
                   "evaluate_ordered_schedule: forest root not found in the "
                   "schedule's value map");
    pre_results[i] = value_results[vid_it->second];
    SEQUANT_ASSERT(pre_results[i] &&
                   "evaluate_ordered_schedule: forest root was never "
                   "produced by a root BuildStep");
  }

  // -------- Shared combine: permute each root to layout and sum. --------
  // combine_forest_roots (forest_combine.hpp) is the SAME helper
  // evaluate_whole_scope's tail calls, so the two executors emit
  // byte-identical Term/Permute/SumInplace trace bookkeeping without a
  // hand-synced second copy.
  return combine_forest_roots<EvalTrace>(roots, pre_results, layout, cache);
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_ORDERED_EXECUTOR_HPP
