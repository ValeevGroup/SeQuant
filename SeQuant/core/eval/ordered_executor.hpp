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
#include <optional>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

namespace sequant::eval {

namespace detail {

///
/// \brief \c find_leaf_carrying (eval.hpp), but resolving \p axis against
/// whichever of \p block's OWN production sites (its direct \c BuildStep
/// values first, then its own \c outputs values, then -- if neither carries
/// it -- recursing into its nested child \c ScopeBlock's own production
/// sites) actually carries it below. Used to source \c mode_batches: an
/// \c AccumulateSum output's own node does NOT carry the axis free (it is
/// summed away AT that node -- see \c LoopRole::Reduction), but \c
/// find_leaf_carrying still finds a carrying LEAF further down its subtree,
/// so trying every production site in turn (rather than requiring one
/// specific site) is what makes this work uniformly for both \c BuildStep
/// (LoopLocal, carries the axis at its own node) and \c outputs
/// (Reduction/AccumulateSum, carries it only below).
///
template <typename node_t>
[[nodiscard]] std::optional<std::pair<node_t, std::size_t>> ordered_axis_leaf(
    ScopeBlock const& block, Index const& axis,
    std::unordered_map<std::size_t, node_t> const& vmap,
    RichSchedule const& rich) {
  auto const resolve = [&](std::size_t vid) -> node_t const* {
    auto const hash = rich.cells[vid].hash;
    auto const it = vmap.find(hash);
    return it == vmap.end() ? nullptr : &it->second;
  };
  for (Step const& step : block.steps) {
    if (auto const* build = std::get_if<BuildStep>(&step.value)) {
      if (auto const* n = resolve(build->value_id))
        if (auto found = find_leaf_carrying(*n, axis)) return found;
    }
  }
  for (auto const& [vid, kind] : block.outputs) {
    (void)kind;
    if (auto const* n = resolve(vid))
      if (auto found = find_leaf_carrying(*n, axis)) return found;
  }
  for (Step const& step : block.steps) {
    if (auto const* child = std::get_if<ScopeBlock>(&step.value)) {
      if (auto found = ordered_axis_leaf<node_t>(*child, axis, vmap, rich))
        return found;
    }
  }
  return std::nullopt;
}

///
/// \brief SP3 Task 2: realize one Contracted \c ScopeBlock's batch loop
/// against \p parent_cache -- the ordered-schedule counterpart of \c
/// scope_executor.hpp's \c detail::walk_scope, simplified to the shape the
/// \c OrderedSchedule IR already gives: \p block's own \c steps are a
/// topologically ORDERED interleaving of \c BuildStep's (values homed AT
/// this block -- \c LoopRole::LoopLocal, "Transient" per \c
/// ordered_schedule.hpp's \c OutputKind doc comment: built fresh every batch
/// on the per-batch scratch and simply dropped by the next \c reset(), never
/// stored anywhere) and nested child \c ScopeBlock steps (realized
/// recursively, in full, once per iteration of THIS loop -- an outer
/// Contracted loop with an inner child re-runs the WHOLE inner loop every
/// outer batch, exactly as \c walk_scope's nested case does); \p block's own
/// \c outputs list the values that escape it on close (SP3 Task 3's
/// AccumulateScatter aside, only \c AccumulateSum here): summed via \c
/// add_inplace across every batch, then stored at \p parent_cache -- the
/// scope one level OUT of this loop, i.e. where the \c ScopeBlock step
/// itself sits in ITS OWN enclosing block's \c steps -- so a later step at
/// that level reads it whole via the ordinary Checked cache probe, and
/// mirrored into \p value_results (keyed by the SAME global \c value_id
/// space \c evaluate_ordered_schedule's root walk uses) so the final
/// per-root combine can resolve a forest root produced INSIDE a loop block
/// exactly as it resolves one produced by a plain root \c BuildStep.
///
/// \details Per batch: \c bs.cache.reset() (drops every Transient/LoopLocal
/// value from the PRIOR batch), the batch context is extended by this
/// block's own axis over the batch's element range, then \p block's own \c
/// steps run in schedule order (a \c BuildStep evaluated on the scratch; a
/// nested \c ScopeBlock realized recursively, its own full loop nested
/// inside this single batch), and finally each \c AccumulateSum output is
/// (re)built on the SAME scratch (its operands -- this block's own \c
/// BuildStep values and/or a nested block's already-summed output, both
/// alive on \c bs.cache for this batch -- are already resolved) and summed
/// into the running accumulator. \c make_batched_scratch's shared-scratch
/// CSE dedups any sub-intermediate repeated across this block's own
/// \c BuildStep/\c outputs production sites (co-evaluated as one \p members
/// group, mirroring \c walk_scope's single-aux-loop \c group), the same
/// mechanism \c walk_scope's Task-3 leaf case relies on.
///
/// \note Every value inside \p block is keyed off \p block's own \c axis (a
/// single canonical \c Index, not a per-value physical remap) -- unlike \c
/// walk_scope's \c member_contracted_axis, which maps the schedule's
/// canonical mode to EACH member's own physical label because \c
/// ScopeSchedule's members are independently-labeled forest ROOTS. \c
/// OrderedSchedule's own bucketing (\c build_ordered_schedule's step 2)
/// instead groups by axis TYPE across the WHOLE forest, so a schedule built
/// from a forest whose cells name that TYPE under more than one physical
/// \c Index would need the same remap \c walk_scope applies; no fixture
/// exercises that yet (a single-physical-label forest, e.g. every operand
/// naming the same literal \c x_1), so this is not handled here -- a later
/// task's job if it becomes live.
///
template <Trace EvalTrace, typename node_t, typename F, typename N, bool FHC>
void run_ordered_contracted_block(
    ScopeBlock const& block,
    std::unordered_map<std::size_t, node_t> const& vmap,
    RichSchedule const& rich, F const& leaf_evaluator,
    CacheManager<N, FHC>& parent_cache,
    std::function<std::size_t(Index const&)> const& target,
    typename CacheManager<N, FHC>::BatchContext const& ectx,
    container::vector<ResultPtr>& value_results) {
  using Cache = CacheManager<N, FHC>;
  using BatchContext = typename Cache::BatchContext;
  using member_t = std::pair<node_t const*, Index>;

  SEQUANT_ASSERT(
      block.kind == BatchModeType::Contracted &&
      "evaluate_ordered_schedule: an External loop block (AccumulateScatter "
      "on close) is SP3 Task 3's job -- Task 2 handles Contracted only.");

  auto const resolve = [&](std::size_t vid) -> node_t const& {
    auto const hash = rich.cells[vid].hash;
    auto const it = vmap.find(hash);
    SEQUANT_ASSERT(it != vmap.end() &&
                   "evaluate_ordered_schedule: a loop-block value_id was not "
                   "found in the forest's value-node map");
    return it->second;
  };

  // members co-evaluated within one batch pass: this block's own direct
  // BuildStep values plus its own AccumulateSum output values -- one shared
  // scratch, so a sub-intermediate repeated between/within them is built
  // once per batch (mirrors walk_scope's Task-3 single-aux-loop `group`).
  std::vector<member_t> members;
  for (Step const& step : block.steps)
    if (auto const* build = std::get_if<BuildStep>(&step.value))
      members.push_back({&resolve(build->value_id), block.axis});
  for (auto const& [vid, kind] : block.outputs) {
    SEQUANT_ASSERT(kind == OutputKind::AccumulateSum &&
                   "evaluate_ordered_schedule: an AccumulateScatter output "
                   "is SP3 Task 3's job -- Task 2 handles AccumulateSum "
                   "(reduction) only.");
    members.push_back({&resolve(vid), block.axis});
  }

  auto bs = sequant::detail::make_batched_scratch(members, parent_cache);
  for (auto const* s : bs.seeds)
    (void)bs.cache.store(*s, parent_cache.access(*s));
  bs.cache.set_parent(&parent_cache);

  auto const lf = ordered_axis_leaf<node_t>(block, block.axis, vmap, rich);
  SEQUANT_ASSERT(lf &&
                 "evaluate_ordered_schedule: no leaf below this loop block "
                 "carries its own axis -- cannot source mode_batches");
  auto const batches =
      leaf_evaluator(lf->first)->mode_batches(lf->second, target(block.axis));

  container::vector<ResultPtr> acc(block.outputs.size());
  for (auto const& [e_lo, e_hi] : batches) {
    if (e_lo == e_hi) continue;
    bs.cache.reset();
    BatchContext ctx = ectx;
    ctx.push_back({block.axis, {e_lo, e_hi}});
    bs.cache.set_batch_context(ctx);

    for (Step const& step : block.steps) {
      if (auto const* build = std::get_if<BuildStep>(&step.value)) {
        (void)evaluate_impl<EvalTrace>(resolve(build->value_id), leaf_evaluator,
                                       bs.cache);
      } else {
        auto const& child = std::get<ScopeBlock>(step.value);
        run_ordered_contracted_block<EvalTrace>(child, vmap, rich,
                                                leaf_evaluator, bs.cache,
                                                target, ctx, value_results);
      }
    }

    for (std::size_t k = 0; k != block.outputs.size(); ++k) {
      auto const vid = block.outputs[k].first;
      ResultPtr part =
          evaluate_impl<EvalTrace>(resolve(vid), leaf_evaluator, bs.cache);
      if (!acc[k])
        acc[k] = std::move(part);
      else
        acc[k]->add_inplace(*part);
    }
  }

  for (std::size_t k = 0; k != block.outputs.size(); ++k) {
    SEQUANT_ASSERT(acc[k] &&
                   "evaluate_ordered_schedule: a loop block realized zero "
                   "batches for an AccumulateSum output");
    auto const vid = block.outputs[k].first;
    value_results[vid] = acc[k];
    (void)parent_cache.store(resolve(vid), std::move(acc[k]));
  }
}

}  // namespace detail

///
/// \brief SP3 Tasks 1-2 of the ordered-scope batched-eval design: the
/// ORDERED executor. Walks \c ordered.root.steps in sequence, building one
/// value per root-level \c BuildStep via \c evaluate_impl (which cache-probes
/// every operand at its own \c Stage::Enter, so a value already \c
/// cache.store'd by an earlier step short-circuits its own re-descent -- see
/// \c evaluate_impl's doc comment) and one value per root-level Contracted
/// \c ScopeBlock via \c detail::run_ordered_contracted_block (Task 2: a
/// realized batch loop, its own \c BuildStep's/nested child blocks run per
/// batch on a scratch cache, its \c AccumulateSum outputs summed across
/// batches and stored at the level the block itself sits in), then combines
/// the forest roots' results into the final \c ResultPtr exactly as \c
/// evaluate_whole_scope's shared root-combine loop does (permute-to-\p
/// layout + cross-root \c add_inplace, with the identical Term/Permute/
/// SumInplace trace bookkeeping).
///
/// \details A root-level External \c ScopeBlock (\c AccumulateScatter on
/// close) is SP3 Task 3's job -- \c detail::run_ordered_contracted_block
/// trips a loud \c SEQUANT_ASSERT on \c BatchModeType::External or an
/// \c AccumulateScatter output rather than silently mis-evaluating; no
/// caller of the Task 1/2 gating flag exercises External batching yet (see
/// \c BatchPolicy::ordered_schedule_execution), so this is an inert
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
/// \param target Per-index batch partition size (elements): the source of
///        each realized loop block's batch partition (\c
///        detail::run_ordered_contracted_block's \c mode_batches call).
/// \param make_scope_guard Backend scope-guard factory; unused by Tasks 1-2
///        (no backend screening relaxation is threaded into the loop-block
///        walk yet), threaded for interface symmetry with later tasks.
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
    std::function<std::size_t(Index const&)> const& target,
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

  // -------- Tasks 1-2: walk the root block's own steps, in order. --------
  // A BuildStep is built directly (Task 1); a child ScopeBlock (a realized
  // batch loop) is run via detail::run_ordered_contracted_block (Task 2),
  // which mirrors its own AccumulateSum outputs into value_results so the
  // root-combine below resolves a forest root produced inside a loop block
  // exactly like one produced by a plain root BuildStep. An External child
  // block (AccumulateScatter) is SP3 Task 3's job; see the function doc
  // comment for why that is a loud tripwire, not a silent gap, today.
  typename CacheManager<N, FHC>::BatchContext const root_ectx;
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
      auto const& block = std::get<ScopeBlock>(step.value);
      detail::run_ordered_contracted_block<EvalTrace>(
          block, vmap, rich, leaf_evaluator, cache, target, root_ectx,
          value_results);
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
