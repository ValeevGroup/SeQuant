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
/// \brief Task 2 of the eager-home-release plan: a TYPE-keyed (\c
/// IndexSpace::base_key()) block-count function over \p ordered's whole
/// \c ScopeBlock tree -- the \c OrderedSchedule counterpart of \c
/// walk_scope's \c nblocks_by_type (scope_executor.hpp). Walks every nested
/// \c ScopeBlock reachable from \p ordered.root.steps (root's own axis is
/// the sentinel default, never a real loop, so it is not itself registered);
/// for each block's axis TYPE not yet seen, sources a carrying leaf via \c
/// ordered_axis_leaf (the same leaf \c run_ordered_contracted_block reads
/// \c mode_batches from) and memoizes that leaf's own \c mode_batches(...)
/// .size() under the TYPE. The returned lambda looks up an \c Index's own
/// TYPE in that memo, defaulting to 1 (unbatched) for any type the walk
/// never saw as a block axis.
///
/// \note Mirrors the brief's reference shape exactly but for one thing: a
/// TYPE realized as more than one nested loop (an outer block and an inner
/// child block sharing the same axis TYPE) would still only memoize the
/// FIRST (outermost) one's \c mode_batches().size() here, undercounting the
/// true nested block count (which multiplies). No fixture in this repo
/// nests the SAME axis type today, so this is not implemented -- a later
/// task's job if/when a nested-same-type fixture goes live (the Task 5
/// validate assert is expected to catch an undercount if this is wrong).
///
template <typename node_t, typename F>
[[nodiscard]] std::function<std::size_t(Index const&)> ordered_n_blocks(
    OrderedSchedule const& ordered, RichSchedule const& rich,
    std::unordered_map<std::size_t, node_t> const& vmap,
    F const& leaf_evaluator,
    std::function<std::size_t(Index const&)> const& target) {
  auto by_type =
      std::make_shared<std::unordered_map<std::wstring, std::size_t>>();
  auto const add = [&](auto&& self, ScopeBlock const& b) -> void {
    std::wstring const bk(b.axis.space().base_key());
    if (!by_type->count(bk)) {
      if (auto const lf = ordered_axis_leaf<node_t>(b, b.axis, vmap, rich)) {
        by_type->emplace(bk, leaf_evaluator(lf->first)
                                 ->mode_batches(lf->second, target(b.axis))
                                 .size());
      }
    }
    for (Step const& s : b.steps)
      if (auto const* child = std::get_if<ScopeBlock>(&s.value))
        self(self, *child);
  };
  // Root's own axis is the sentinel default (never a real loop -- see
  // ScopeBlock::axis's doc comment), so only its CHILD blocks are walked.
  for (Step const& s : ordered.root.steps)
    if (auto const* child = std::get_if<ScopeBlock>(&s.value)) add(add, *child);
  return [by_type](Index const& m) -> std::size_t {
    auto const it = by_type->find(std::wstring(m.space().base_key()));
    return it == by_type->end() ? std::size_t{1} : it->second;
  };
}

///
/// \brief Defensive check for the single-physical-label simplification \c
/// run_ordered_contracted_block's own doc comment (\c \\note) documents: \p
/// node's subtree carries \p axis's TYPE (\c IndexSpace::base_key()) at some
/// leaf, but \p node's own subtree never carries the EXACT \p axis \c Index
/// (identity: space + ordinal + proto-indices) anywhere below it.
///
/// \details \c run_ordered_contracted_block slices/accumulates every value
/// in its block by \p axis's exact identity (\c ctx.push_back({block.axis,
/// ...})); \c evaluate_impl's Enter-stage \c slice_to_use then looks up that
/// EXACT \c Index on each fetched node (\c index_position(nd, block.axis)).
/// A member that instead carries a DIFFERENT physical label of the SAME
/// TYPE (e.g. a forest whose cells name the aux space under both \c x_1 and
/// \c x_2) would silently fail that lookup -- \c index_position returns
/// \c nullopt, so the node is never sliced, builds FULL every batch, and is
/// \c add_inplace'd once per batch: an N-fold-inflated reduction with NO
/// diagnostic. This predicate lets the caller turn that into a loud \c
/// SEQUANT_ASSERT instead, mirroring the AccumulateScatter/External asserts
/// already in this file.
///
template <typename node_t>
[[nodiscard]] bool ordered_axis_label_mismatch(node_t const& node,
                                               Index const& axis) {
  auto const base = axis.space().base_key();
  bool carries_type = false;
  bool carries_exact = false;
  auto const walk = [&](auto&& self, node_t const& n) -> void {
    if (carries_exact) return;  // already proven safe; no need to keep going
    if (n.leaf()) {
      for (Index const& ix : n->canon_indices()) {
        if (ix.space().base_key() != base) continue;
        carries_type = true;
        if (ix == axis) {
          carries_exact = true;
          return;
        }
      }
      return;
    }
    self(self, n.left());
    self(self, n.right());
  };
  walk(walk, node);
  return carries_type && !carries_exact;
}

///
/// \brief SP3 Tasks 2-3: realize one \c ScopeBlock's batch loop against
/// \p parent_cache -- the ordered-schedule counterpart of \c
/// scope_executor.hpp's \c detail::walk_scope, simplified to the shape the
/// \c OrderedSchedule IR already gives: \p block's own \c steps are a
/// topologically ORDERED interleaving of \c BuildStep's (values homed AT
/// this block -- \c LoopRole::LoopLocal, "Transient" per \c
/// ordered_schedule.hpp's \c OutputKind doc comment: built fresh every batch
/// on the per-batch scratch and simply dropped by the next \c reset(), never
/// stored anywhere) and nested child \c ScopeBlock steps (realized
/// recursively, in full, once per iteration of THIS loop -- an outer loop
/// with an inner child re-runs the WHOLE inner loop every outer batch,
/// exactly as \c walk_scope's nested case does); \p block's own \c outputs
/// list the values that escape it on close, either \c AccumulateSum
/// (reduction: summed via \c add_inplace across every batch -- Task 2) or \c
/// AccumulateScatter (loop-carried: written into a disjoint slice of a
/// pre-sized destination every batch -- Task 3, mirroring \c walk_scope's
/// External branch's \c pre_sized_zeros_over_mode / \c write_into_slice
/// pair). A block may mix both kinds of output freely (and \c BuildStep
/// Transients alongside them): the batch loop is the SAME either way -- only
/// how each output's per-batch \c part is folded into its own running result
/// differs by its \c OutputKind. Every output, once closed, is stored at \p
/// parent_cache -- the scope one level OUT of this loop, i.e. where the \c
/// ScopeBlock step itself sits in ITS OWN enclosing block's \c steps -- so a
/// later step at that level reads it whole via the ordinary Checked cache
/// probe, and mirrored into \p value_results (keyed by the SAME global \c
/// value_id space \c evaluate_ordered_schedule's root walk uses) so the
/// final per-root combine can resolve a forest root produced INSIDE a loop
/// block exactly as it resolves one produced by a plain root \c BuildStep.
///
/// \details Per batch: \c bs.cache.reset() (drops every Transient/LoopLocal
/// value from the PRIOR batch), the batch context is extended by this
/// block's own axis over the batch's element range, then \p block's own \c
/// steps run in schedule order (a \c BuildStep evaluated on the scratch; a
/// nested \c ScopeBlock realized recursively, its own full loop nested
/// inside this single batch), and finally each output is (re)built on the
/// SAME scratch (its operands -- this block's own \c BuildStep values and/or
/// a nested block's already-closed output, both alive on \c bs.cache for
/// this batch -- are already resolved): an \c AccumulateSum output is summed
/// into its running accumulator; an \c AccumulateScatter output is \c
/// write_into_slice'd into its running destination (built once, on the FIRST
/// nonempty batch, via \c pre_sized_zeros_over_mode against the SAME carrier
/// leaf \c mode_batches itself was sourced from -- sound because of the \c
/// \note below: every value in this block, escape output included, shares
/// ONE physical axis identity, hence one tiling). \c make_batched_scratch's
/// shared-scratch CSE dedups any sub-intermediate repeated across this
/// block's own \c BuildStep/\c outputs production sites (co-evaluated as one
/// \p members group, mirroring \c walk_scope's single-aux-loop \c group),
/// the same mechanism \c walk_scope's Task-3 leaf case relies on -- and,
/// because that CSE and the batch partition are shared across BOTH output
/// kinds here (unlike \c walk_scope's External branch, which cannot share a
/// scratch or a batch partition across members that may bind independently-
/// labeled physical axes -- see the \note below for why \c OrderedSchedule's
/// members can), a mixed-kind block realizes its \c AccumulateSum and \c
/// AccumulateScatter outputs in the SAME single pass over \c batches, not
/// two.
///
/// \par Forced-split producer/consumer passes (SP2 Task 4)
/// A forced split realizes its axis as TWO sibling \c ScopeBlock \c Step's
/// at the SAME nesting level (ordinal 0 the producer, ordinal 1 the
/// consumer) rather than one nested inside the other -- see \c
/// build_ordered_schedule's step 2b. No special-casing is needed here for
/// that shape: \c evaluate_ordered_schedule's root walk (and this function's
/// own \c steps loop, for a split at a non-root level) already runs sibling
/// \c Step's in \p block.steps SEQUENTIALLY, and \c
/// ordered_schedule_topo_sort_steps already guarantees the consumer pass
/// sorts AFTER the producer pass (the consumer's \c requires_ names the
/// producer's escaped -- now \c AccumulateScatter'd to FULL -- outputs, a
/// real dependency edge). Since each pass, on closing, \c stores its outputs
/// at THIS level's shared \p parent_cache (or \c cache at the root) before
/// the next sibling step begins, the consumer pass's own inner \c BuildStep
/// probes find the producer's completed value there via the ordinary
/// Checked cache lookup -- reading it WHOLE, exactly as design intends. This
/// is the identical mechanism Task 2 already relies on for a plain nested
/// child block reading an enclosing block's homed value; a forced split
/// changes only which axis realizes two blocks instead of one, not how
/// cross-step visibility works.
///
/// \note Every value inside \p block is keyed off \p block's own \c axis (a
/// single canonical \c Index, not a per-value physical remap) -- unlike \c
/// walk_scope's \c member_contracted_axis / \c member_external_axis, which
/// map the schedule's canonical mode to EACH member's own physical label
/// because \c ScopeSchedule's members are independently-labeled forest
/// ROOTS. \c OrderedSchedule's own bucketing (\c build_ordered_schedule's
/// step 2) instead groups by axis TYPE across the WHOLE forest, so a
/// schedule built from a forest whose cells name that TYPE under more than
/// one physical \c Index would need the same remap \c walk_scope applies;
/// \c ordered_axis_label_mismatch defensively asserts against that case (see
/// its own doc comment) rather than silently mis-evaluating -- no fixture
/// exercises the multi-label case yet (every current fixture, including
/// External/\c AccumulateScatter ones, names its batch axis under a single
/// literal \c Index throughout), so the remap itself is not implemented
/// here -- a later task's job if it becomes live.
///
/// \note \p is_volatile (trailing, defaulted) is the NODE-level lift of \c
/// BatchPolicy::is_volatile_leaf -- see \c evaluate_ordered_schedule's own
/// doc comment for where it is produced. Threaded through this function's
/// own recursion (the nested \c ScopeBlock call below) so it reaches every
/// level, but not yet CONSULTED here -- a later task's job.
///
template <Trace EvalTrace, typename node_t, typename F, typename N, bool FHC>
void run_ordered_contracted_block(
    ScopeBlock const& block,
    std::unordered_map<std::size_t, node_t> const& vmap,
    RichSchedule const& rich, F const& leaf_evaluator,
    CacheManager<N, FHC>& parent_cache,
    std::function<std::size_t(Index const&)> const& target,
    typename CacheManager<N, FHC>::BatchContext const& ectx,
    container::vector<ResultPtr>& value_results, container::vector<char>& built,
    std::function<bool(node_t const&)> const& is_volatile = {},
    std::function<std::size_t(Index const&)> const& n_blocks = {},
    std::function<std::size_t(std::size_t)> const& home_reads = {}) {
  using Cache = CacheManager<N, FHC>;
  using BatchContext = typename Cache::BatchContext;
  using member_t = std::pair<node_t const*, Index>;

  // R4 (SP3 Task 5): loud guard on this block's batch-mode kind. The batch-loop
  // primitive below realizes Contracted and External blocks UNIFORMLY (their
  // difference is carried entirely by each escape output's OutputKind, handled
  // per-output below), so both are supported; any OTHER BatchModeType value is
  // a schedule this executor cannot interpret and must refuse loudly rather
  // than silently mis-run. (An assert -- not a switch/default -- so the 2-value
  // enum does not trip -Wcovered-switch-default.)
  SEQUANT_ASSERT((block.kind == BatchModeType::Contracted ||
                  block.kind == BatchModeType::External) &&
                 "evaluate_ordered_schedule: unsupported ScopeBlock batch-mode "
                 "kind -- only Contracted/External are realizable");

  auto const resolve = [&](std::size_t vid) -> node_t const& {
    auto const hash = rich.cells[vid].hash;
    auto const it = vmap.find(hash);
    SEQUANT_ASSERT(it != vmap.end() &&
                   "evaluate_ordered_schedule: a loop-block value_id was not "
                   "found in the forest's value-node map");
    return it->second;
  };

  // members co-evaluated within one batch pass: this block's own direct
  // BuildStep values plus its own escape output values (AccumulateSum AND
  // AccumulateScatter alike -- see this function's own doc comment for why
  // they share one scratch/one batch partition here, unlike walk_scope's
  // External branch). Each member is defensively checked against the
  // single-physical-label simplification this function's own doc comment
  // (\note) documents: see ordered_axis_label_mismatch's doc comment for why
  // a silent skip here would otherwise be a silently N-fold-inflated
  // reduction, not a crash.
  std::vector<member_t> members;
  for (Step const& step : block.steps)
    if (auto const* build = std::get_if<BuildStep>(&step.value)) {
      node_t const& nd = resolve(build->value_id);
      SEQUANT_ASSERT(
          !ordered_axis_label_mismatch(nd, block.axis) &&
          "evaluate_ordered_schedule: multi-physical-label axis in one "
          "block not yet supported");
      members.push_back({&nd, block.axis});
    }

  // dm[k]: for an AccumulateScatter output, the position of block.axis on
  // ITS OWN result (computed once here, outside the batch loop, since a
  // value's structural mode order is the same across every batch) -- the
  // scatter mode walk_scope's External branch resolves per member via
  // index_position(*m, axis). Left at nullopt (and unused) for an
  // AccumulateSum output, which reduces block.axis away at its own node and
  // so never carries it on its own result.
  container::vector<std::optional<std::size_t>> dm(block.outputs.size());
  for (std::size_t k = 0; k != block.outputs.size(); ++k) {
    auto const vid = block.outputs[k].first;
    auto const kind = block.outputs[k].second;
    SEQUANT_ASSERT(kind == OutputKind::AccumulateSum ||
                   kind == OutputKind::AccumulateScatter);
    node_t const& nd = resolve(vid);
    SEQUANT_ASSERT(!ordered_axis_label_mismatch(nd, block.axis) &&
                   "evaluate_ordered_schedule: multi-physical-label axis in "
                   "one block not yet supported");
    members.push_back({&nd, block.axis});
    if (kind == OutputKind::AccumulateScatter) {
      dm[k] = index_position(nd, block.axis);
      SEQUANT_ASSERT(dm[k] &&
                     "evaluate_ordered_schedule: an AccumulateScatter output "
                     "does not carry its own escape axis on its result");
    }
  }

  auto bs = sequant::detail::make_batched_scratch(members, parent_cache,
                                                  /*read_from_home=*/true);
  bs.cache.set_parent(&parent_cache);

  auto const lf = ordered_axis_leaf<node_t>(block, block.axis, vmap, rich);
  SEQUANT_ASSERT(lf &&
                 "evaluate_ordered_schedule: no leaf below this loop block "
                 "carries its own axis -- cannot source mode_batches");
  ResultPtr const carrier_full = leaf_evaluator(lf->first);
  auto const batches =
      carrier_full->mode_batches(lf->second, target(block.axis));

  // Loop-structure trace. The ordered executor previously emitted NO
  // batch-execution marker (unlike scope_executor.hpp's whole-scope BatchGroup
  // and eval.hpp's forest SCHEDULE_RUN_GROUP), so whether batching engaged was
  // invisible. log::printing()-gated first-class marker, plus a structured
  // SEQUANT_SCHED_DUMP line for scriptable confirmation that this path realizes
  // batch blocks (block axis + batch count).
  if (log::printing()) {
    BatchContext s = ectx;
    s.push_back({block.axis, {0, 0}});
    log::log("BatchGroup", "Begin",
             std::format("{} members co-evaluated over {} batches of {} {}",
                         members.size(), batches.size(),
                         toUtf8(block.axis.full_label()), log::scope_annot(s)));
  }
  {
    static bool const sched_dump = std::getenv("SEQUANT_SCHED_DUMP") != nullptr;
    if (sched_dump)
      std::cerr << "ORDERED_RUN_BLOCK {\"kind\":\""
                << (block.kind == BatchModeType::External ? "external"
                                                          : "contracted")
                << "\",\"mode\":\"" << toUtf8(block.axis.full_label())
                << "\",\"blocks\":" << batches.size()
                << ",\"members\":" << members.size() << "}\n";
  }

  container::vector<ResultPtr> acc(block.outputs.size());
  container::vector<ResultPtr> dest(block.outputs.size());

  // An escape output that is ALREADY RESIDENT at its home -- a PERSISTENT,
  // batch-invariant output (it contracts/reduces the batch axis, so its result
  // does not carry it) that was built in a PRIOR CC iteration and kept across
  // cache.reset() -- must be REUSED, not re-accumulated. Re-running the batch
  // loop for it would read its full homed value as EACH batch's "partial" (the
  // per-batch evaluate_impl below resolves the output node from its home) and
  // sum it N_batches-fold, corrupting the persistent entry into the next
  // iteration. So skip such outputs here and reuse the resident value at close;
  // the block still runs for the volatile (and first-iteration) outputs and for
  // its loop-local Transients. On the first iteration nothing is resident yet,
  // so every output is accumulated normally.
  container::vector<char> reuse(block.outputs.size(), 0);
  for (std::size_t k = 0; k != block.outputs.size(); ++k)
    reuse[k] =
        parent_cache.resident_in_chain(resolve(block.outputs[k].first)) ? 1 : 0;

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
        // R3: record this loop-local Transient as produced (it is built fresh
        // every batch on the scratch and never lands in value_results, so the
        // built ledger is the only faithful "was built" slot for it).
        built[build->value_id] = 1;
      } else if (auto const* child = std::get_if<ScopeBlock>(&step.value)) {
        run_ordered_contracted_block<EvalTrace>(
            *child, vmap, rich, leaf_evaluator, bs.cache, target, ctx,
            value_results, built, is_volatile, n_blocks, home_reads);
      } else {
        // R4: the Step variant has exactly BuildStep/ScopeBlock alternatives;
        // a valueless-by-exception or future third alternative is a schedule
        // this executor cannot interpret.
        SEQUANT_ASSERT(false &&
                       "evaluate_ordered_schedule: unsupported Step variant");
      }
    }

    for (std::size_t k = 0; k != block.outputs.size(); ++k) {
      if (reuse[k])
        continue;  // resident persistent output: reuse, don't re-sum
      auto const vid = block.outputs[k].first;
      auto const kind = block.outputs[k].second;
      ResultPtr part =
          evaluate_impl<EvalTrace>(resolve(vid), leaf_evaluator, bs.cache);
      if (kind == OutputKind::AccumulateSum) {
        if (!acc[k])
          acc[k] = std::move(part);
        else
          acc[k]->add_inplace(*part);
      } else if (kind == OutputKind::AccumulateScatter) {
        if (!dest[k])
          dest[k] = part->pre_sized_zeros_over_mode(*dm[k], *carrier_full,
                                                    lf->second);
        dest[k]->write_into_slice(*part, *dm[k], e_lo, e_hi);
      } else {
        // R4: an escape output is AccumulateSum or AccumulateScatter (a
        // Transient never appears in outputs -- see OutputKind's doc comment);
        // any other kind is a construct this batch fold cannot realize.
        SEQUANT_ASSERT(false &&
                       "evaluate_ordered_schedule: unsupported escape "
                       "OutputKind in batch fold");
      }
    }
  }

  for (std::size_t k = 0; k != block.outputs.size(); ++k) {
    auto const vid = block.outputs[k].first;
    auto const kind = block.outputs[k].second;
    // Resident persistent output (see the reuse note above): reuse its home
    // value untouched -- it survives reset() and is invariant across CC
    // iterations -- rather than re-home/re-store it. It was NOT re-accumulated
    // in the loop, so acc[k]/dest[k] are null here.
    if (reuse[k]) {
      value_results[vid] = parent_cache.access(resolve(vid));
      SEQUANT_ASSERT(value_results[vid] &&
                     "evaluate_ordered_schedule: a resident output vanished "
                     "from its home before block close");
      built[vid] = 1;
      continue;
    }
    // R4: exhaustive on OutputKind before selecting the running result to home.
    SEQUANT_ASSERT((kind == OutputKind::AccumulateSum ||
                    kind == OutputKind::AccumulateScatter) &&
                   "evaluate_ordered_schedule: unsupported escape OutputKind "
                   "at block close");
    ResultPtr& out = (kind == OutputKind::AccumulateSum) ? acc[k] : dest[k];
    SEQUANT_ASSERT(out &&
                   "evaluate_ordered_schedule: a loop block realized zero "
                   "batches for an escape output");
    value_results[vid] = out;
    built[vid] = 1;  // R3: this escape output is produced (homed below).
    // Home the closed output at the scope one level OUT with a RESIDENT slot
    // so a later step there reads it WHOLE (via the ordinary Checked probe)
    // instead of re-descending and rebuilding it -- for an AccumulateSum
    // reduction that means rebuilding the FULL un-batched contraction (the
    // very peak this loop batched away). store() is a no-op on an unregistered
    // key, so without homing the output first, a consumer used only once (not
    // a CSE candidate) would silently recompute it.
    node_t const& out_node = resolve(vid);
    // Task 4: classify this homed escape output volatile-vs-persistent and set
    // its cache life. A subtree carrying a volatile leaf (is_volatile) is
    // NON-persistent -- released at its genuine last use (home_reads: the exact
    // number of times its home entry is accessed under read-from-home, from the
    // ordered schedule's realized scopes) rather than pinned resident; a
    // persistent (invariant) composite survives reset() and is built once, its
    // life ignored by the overload. With no volatility policy (is_volatile
    // empty) fall back to the unconditional resident pin, the pre-Task-4
    // behavior.
    if (!out_node.leaf()) {
      if (is_volatile) {
        bool const vol = subtree_any(out_node, is_volatile);
        std::size_t const life = home_reads ? home_reads(vid) : 1;
        parent_cache.ensure_home_slot(out_node, life, /*persistent=*/!vol);
      } else {
        parent_cache.ensure_home_slot(out_node);
      }
    }
    (void)parent_cache.store(out_node, std::move(out));
  }
}

}  // namespace detail

///
/// \brief SP3 Tasks 1-3 of the ordered-scope batched-eval design: the
/// ORDERED executor. Walks \c ordered.root.steps in sequence, building one
/// value per root-level \c BuildStep via \c evaluate_impl (which cache-probes
/// every operand at its own \c Stage::Enter, so a value already \c
/// cache.store'd by an earlier step short-circuits its own re-descent -- see
/// \c evaluate_impl's doc comment) and one value per root-level \c ScopeBlock
/// (Contracted or External alike) via \c detail::run_ordered_contracted_block
/// (Task 2: a realized batch loop, its own \c BuildStep's/nested child blocks
/// run per batch on a scratch cache, its \c AccumulateSum outputs summed
/// across batches; Task 3: its \c AccumulateScatter outputs written into a
/// disjoint slice of a pre-sized destination each batch -- both kinds stored
/// at the level the block itself sits in on close, including a forced
/// split's two ordinal-ordered producer/consumer sibling blocks, run in
/// schedule order like any other pair of steps), then combines the forest
/// roots' results into the final \c ResultPtr exactly as \c
/// evaluate_whole_scope's shared root-combine loop does (permute-to-\p
/// layout + cross-root \c add_inplace, with the identical Term/Permute/
/// SumInplace trace bookkeeping).
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
/// \param is_volatile NODE-level volatility predicate (empty means
///        "never volatile"), the lift of \c BatchPolicy::is_volatile_leaf
///        exactly as \c make_evaluator's own \c is_volatile_node lift
///        (eval.hpp) computes it -- threaded down through \c
///        detail::run_ordered_contracted_block's recursion but not yet
///        CONSULTED here (a later task's job: classifying a home value
///        volatile-vs-persistent via \c subtree_any at the homing sites).
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
    [[maybe_unused]] ScopeGuardFactory const& make_scope_guard = {},
    std::function<bool(std::ranges::range_value_t<Nodes> const&)> const&
        is_volatile = {},
    std::function<std::size_t(std::size_t)> const& home_life_override = {}) {
  using node_t = std::ranges::range_value_t<Nodes>;
  static_assert(std::is_same_v<node_t, N>,
                "the forest's node type and the cache's node type must match");

  // hash -> node, resolving a BuildStep's value_id (via rich.cells[vid].hash)
  // to the forest node evaluate_impl builds.
  auto const vmap = build_value_node_map(forest);

  // Task 4: the block-count function over the WHOLE ScopeBlock tree, built
  // ONCE here (Task 2's ordered_n_blocks) and threaded through the homing
  // sites -- both the root-composite site below and, via
  // run_ordered_contracted_block's recursion, every escape-output site. It
  // feeds home_reads (the per-batch loop factors) below. Cheap and
  // side-effect-free (it just sizes each realized loop axis's mode_batches
  // once).
  std::function<std::size_t(Index const&)> const n_blocks =
      detail::ordered_n_blocks<node_t>(ordered, rich, vmap, leaf_evaluator,
                                       target);

  // Task 4: the EXACT home use-count for each homed value under read-from-home,
  // computed ONCE from the ordered schedule's realized scopes + the DAG (with
  // multiplicity) -- see detail::ordered_home_reads. This is the non-persistent
  // life of a volatile homed value: it frees at its genuine last home access
  // instead of being pinned resident. A test may inject home_life_override (a
  // non-evicting life) to MEASURE the actual home reads and assert them equal
  // to home_reads -- the "predicted == measured" gate.
  std::function<std::size_t(std::size_t)> const home_reads =
      home_life_override
          ? home_life_override
          : detail::ordered_home_reads<node_t>(ordered, rich, vmap, n_blocks);

  // Per-value results, indexed by value_id (== a ValueCell's own slot in
  // rich.cells -- see peak_profile.hpp's ValueCell::value_id doc comment),
  // so the root-combine step below reads a forest root's own build directly
  // rather than re-resolving it through the cache.
  container::vector<ResultPtr> value_results(ordered.num_values);

  // R3 (SP3 Task 5): executor-side run-completeness ledger. Set to 1 at the
  // EXACT site each scheduled value is actually built -- a root-level BuildStep
  // below, a loop-local Transient BuildStep, or a block escape output's close
  // (see run_ordered_contracted_block). This is the true "was built" slot for
  // the completeness assert at the tail (value_results holds only root-scope
  // escaping values, so a loop-local Transient never lands there).
  container::vector<char> built(ordered.num_values, 0);

  // -------- Tasks 1-3: walk the root block's own steps, in order. --------
  // A BuildStep is built directly (Task 1); a child ScopeBlock (a realized
  // batch loop, Contracted or External alike) is run via
  // detail::run_ordered_contracted_block (Task 2 AccumulateSum, Task 3
  // AccumulateScatter), which mirrors its own outputs into value_results so
  // the root-combine below resolves a forest root produced inside a loop
  // block exactly like one produced by a plain root BuildStep.
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
      // Home this root-scope value at the root cache with a RESIDENT slot
      // (unbounded life until the per-term reset), so it is built ONCE and
      // read whole by every later consumer -- a root-level ancestor BuildStep
      // that would otherwise re-descend and rebuild it (the plain per-forest
      // CSE life, if any, is drained by its unbatched use count), and a
      // block-internal consumer reaching it up the scope chain (which reads a
      // root-homed invariant once per batch block). Without this, each such
      // consumer recomputes the composite; ensure_home_slot is what makes a
      // Κ-free home={} composite (I(i,i;a,a)) build exactly once. Leaves need
      // no slot (the leaf evaluator just hands back a precomputed input), so
      // this targets only the internal-node BuildSteps the schedule homes here.
      // Task 4: classify this root-homed composite volatile-vs-persistent and
      // seed its cache life (see the escape-output site in
      // run_ordered_contracted_block for the same pattern). A subtree carrying
      // a volatile leaf is NON-persistent (released at its genuine last use);
      // an invariant composite is persistent (built once, survives reset). With
      // no volatility policy (is_volatile empty) keep the unconditional
      // resident pin -- the pre-Task-4 behavior every non-volatility-aware
      // caller relies on.
      if (!it->second.leaf()) {
        if (is_volatile) {
          bool const vol = subtree_any(it->second, is_volatile);
          std::size_t const life = home_reads(vid);
          cache.ensure_home_slot(it->second, life, /*persistent=*/!vol);
        } else {
          cache.ensure_home_slot(it->second);
        }
      }
      value_results[vid] =
          evaluate_impl<EvalTrace>(it->second, leaf_evaluator, cache);
      built[vid] = 1;  // R3: this root-scope value is produced.
    } else if (auto const* block = std::get_if<ScopeBlock>(&step.value)) {
      detail::run_ordered_contracted_block<EvalTrace>(
          *block, vmap, rich, leaf_evaluator, cache, target, root_ectx,
          value_results, built, is_volatile, n_blocks, home_reads);
    } else {
      // R4: the Step variant has exactly BuildStep/ScopeBlock alternatives; any
      // other state is a schedule this executor cannot interpret.
      SEQUANT_ASSERT(false &&
                     "evaluate_ordered_schedule: unsupported Step variant at "
                     "root scope");
    }
  }

  // R3: run-completeness. Every value_id the schedule PROMISES to produce --
  // every BuildStep (root-level or loop-local Transient) and every block escape
  // output, enumerated by collect_production_ids -- must have been actually
  // built by the walk above. Complements well_formed's STATIC single-producer
  // check, which checks no DUPLICATE production but NOT completeness (see
  // ordered_schedule.hpp well_formed's \note). A gap here means the executor
  // skipped a scheduled value: a silent under-execution to refuse, not ignore.
  {
    container::vector<std::size_t> production_ids;
    detail::collect_production_ids(ordered.root, production_ids);
    for (std::size_t const vid : production_ids) {
      SEQUANT_ASSERT(vid < ordered.num_values &&
                     "evaluate_ordered_schedule: production value_id out of "
                     "range");
      SEQUANT_ASSERT(built[vid] &&
                     "evaluate_ordered_schedule: a scheduled value_id was "
                     "never produced during the run (incomplete execution)");
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
