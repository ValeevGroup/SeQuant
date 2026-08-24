#ifndef SEQUANT_EVAL_ORDERED_EXECUTOR_HPP
#define SEQUANT_EVAL_ORDERED_EXECUTOR_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backend_array_ops.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/forest_combine.hpp>
#include <SeQuant/core/eval/member_axis.hpp>
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
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

namespace sequant::eval {

namespace detail {

///
/// \brief Task 2 of the eager-home-release plan: a TYPE-keyed (\c
/// IndexSpace::base_key()) block-count function over \p ordered's whole
/// \c ScopeBlock tree -- the \c OrderedSchedule counterpart of \c
/// walk_scope's \c nblocks_by_type (scope_executor.hpp). Walks every nested
/// \c ScopeBlock reachable from \p ordered.root.steps (root's own axis is
/// the sentinel default, never a real loop, so it is not itself registered);
/// for each block's axis TYPE not yet seen, memoizes the count of the SAME
/// backend \c axis_batches \c run_ordered_contracted_block iterates over that
/// axis. The returned lambda looks up an \c Index's own TYPE in that memo,
/// defaulting to 1 (unbatched) for any type the walk never saw as a block
/// axis.
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
    std::function<std::size_t(Index const&)> const& target,
    BackendArrayOps const* aops) {
  auto by_type =
      std::make_shared<std::unordered_map<std::wstring, std::size_t>>();
  auto const add = [&](auto&& self, ScopeBlock const& b) -> void {
    std::wstring const bk(b.axis.space().base_key());
    if (!by_type->count(bk)) {
      // Single source of truth with the actual batch loop: the block COUNT is
      // the size of the SAME axis_batches the executor iterates (per-space).
      SEQUANT_ASSERT(aops &&
                     "ordered_n_blocks: batched schedule requires backend "
                     "array-ops (CacheManager::set_array_ops)");
      by_type->emplace(bk, aops->axis_batches(b.axis, target(b.axis)).size());
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
    std::function<std::size_t(std::size_t)> const& home_reads = {},
    container::set<std::size_t> const* needed_build = nullptr) {
  using Cache = CacheManager<N, FHC>;
  using BatchContext = typename Cache::BatchContext;
  using member_t = std::pair<node_t const*, Index>;
  // Backend array-ops (zero destination + axis chunking), sourced from the
  // cache chain (the backend -- mpqc's registries or a test's leaf source --
  // wires the root cache). A batched block cannot be realized without it.
  BackendArrayOps const* const aops = parent_cache.array_ops();
  SEQUANT_ASSERT(aops &&
                 "evaluate_ordered_schedule: batched eval requires backend "
                 "array-ops (CacheManager::set_array_ops)");

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
  // Map the block's DAG-scope axis (a SPACE) to each member's OWN physical mode
  // (member_axis.hpp). Index labels are tree-local, so a member is sliced /
  // scattered on the physical label IT carries for the space, never the block's
  // canonical representative reused across members. A Contracted member (and an
  // AccumulateSum output, which reduces the space away) reads it off its own
  // batched_here()/carrying leaf (member_contracted_axis, leaf-side); an
  // External member (and an AccumulateScatter output, which carries the space
  // FREE on its result) reads it off its result (member_external_axis). Falling
  // back to block.axis when the member carries no mode of the space leaves a
  // space-invariant member unsliced (correct). This replaces the former
  // ordered_axis_label_mismatch assert: the multi-label case it guarded against
  // is now HANDLED (per-member remap), not rejected.
  auto const base = std::wstring(block.axis.space().base_key());
  bool const block_external = block.kind == BatchModeType::External;
  auto const member_axis = [&](node_t const& nd, bool external) -> Index {
    auto const ax = external ? member_external_axis(nd, base)
                             : member_contracted_axis(nd, base);
    return ax.value_or(block.axis);
  };

  std::vector<member_t> members;
  for (Step const& step : block.steps)
    if (auto const* build = std::get_if<BuildStep>(&step.value)) {
      node_t const& nd = resolve(build->value_id);
      members.push_back({&nd, member_axis(nd, block_external)});
    }

  // dm[k]: for an AccumulateScatter output, the position of its OWN escape axis
  // on ITS OWN result (computed once here, outside the batch loop, since a
  // value's structural mode order is the same across every batch). Left at
  // nullopt (and unused) for an AccumulateSum output, which reduces its space
  // away at its own node and so never carries it on its own result.
  container::vector<std::optional<std::size_t>> dm(block.outputs.size());
  for (std::size_t k = 0; k != block.outputs.size(); ++k) {
    auto const vid = block.outputs[k].first;
    auto const kind = block.outputs[k].second;
    SEQUANT_ASSERT(kind == OutputKind::AccumulateSum ||
                   kind == OutputKind::AccumulateScatter);
    node_t const& nd = resolve(vid);
    bool const scatter = kind == OutputKind::AccumulateScatter;
    Index const mx = member_axis(nd, scatter);
    members.push_back({&nd, mx});
    if (scatter) {
      dm[k] = index_position(nd, mx);
      SEQUANT_ASSERT(dm[k] &&
                     "evaluate_ordered_schedule: an AccumulateScatter output "
                     "does not carry its own escape axis on its result");
    }
  }

  auto bs = [&]() {
    PhaseTimer::Scope _pt("A.make_scratch");
    return sequant::detail::make_batched_scratch(members, parent_cache,
                                                 /*read_from_home=*/true);
  }();
  bs.cache.set_parent(&parent_cache);
  // This block co-evaluates a whole type-bucketed group under ONE canonical
  // block.axis; a member binding the axis space under a different physical
  // label must be sliced on its own mode (member_axis.hpp). Ask slice-on-use to
  // SPACE-map its axes on this scratch (and every nested child block's, which
  // inherit this flag through their own make_batched_scratch off this cache).
  bs.cache.set_space_mapped_slicing(true);

  // Batch chunks over the loop axis, sourced per-space by the backend (no
  // carrier array is consulted).
  container::svector<std::pair<std::size_t, std::size_t>> const batches =
      aops->axis_batches(block.axis, target(block.axis));

  // Loop-structure trace. The ordered executor previously emitted NO
  // batch-execution marker (unlike scope_executor.hpp's whole-scope BatchGroup
  // and eval.hpp's forest SCHEDULE_RUN_GROUP), so whether batching engaged was
  // invisible. log::printing()-gated first-class marker, plus a structured
  // SEQUANT_SCHED_DUMP line for scriptable confirmation that this path realizes
  // batch blocks (block axis + batch count).
  if (log::printing()) {
    BatchContext s = ectx;
    s.push_back({block.axis, block.level, {0, 0}, std::nullopt});
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
    {
      PhaseTimer::Scope _pt("C.scratch_reset");
      bs.cache.reset();
    }
    BatchContext ctx = ectx;
    ctx.push_back({block.axis, block.level, {e_lo, e_hi}, std::nullopt});
    bs.cache.set_batch_context(ctx);

    for (Step const& step : block.steps) {
      if (auto const* build = std::get_if<BuildStep>(&step.value)) {
        // Cache-halt: skip a dead loop-local Transient -- one whose value is
        // not read this iteration because its only (transitive) consumers are
        // persistent nodes already resident in the cache. This is the SAME gate
        // the top-level scope applies to its BuildSteps (see
        // run_ordered_schedule_pre_results): the needed_build set is the forest
        // BFS that halts at cache-alive nodes, so a resident composite's
        // batch-loop prerequisites are absent from it after iteration 1 and are
        // not re-formed. built[] is marked so the R3 completeness check
        // accounts for the intentional skip. Without a gate (needed_build ==
        // nullptr) the step runs unconditionally -- byte-identical to before.
        if (needed_build &&
            !needed_build->count(rich.cells[build->value_id].hash)) {
          built[build->value_id] = 1;
          continue;
        }
        if (std::getenv("SEQUANT_UT_BLOCK_DIAG"))
          std::cerr << "[BLOCK] axis=" << toUtf8(block.axis.full_label())
                    << " batch=[" << e_lo << "," << e_hi
                    << ") BUILD vid=" << build->value_id << std::endl;
        // PILLAR 2: mark the use-site currently fetching values so
        // slice_to_use can disambiguate a shared symmetric value's free mode.
        node_t const& build_node = resolve(build->value_id);
        auto const prev_consumer = bs.cache.current_consumer();
        bs.cache.set_current_consumer(build_node->hash_value());
        (void)evaluate_impl<EvalTrace>(build_node, leaf_evaluator, bs.cache);
        bs.cache.set_current_consumer(prev_consumer);
        // R3: record this loop-local Transient as produced (it is built fresh
        // every batch on the scratch and never lands in value_results, so the
        // built ledger is the only faithful "was built" slot for it).
        built[build->value_id] = 1;
      } else if (auto const* child = std::get_if<ScopeBlock>(&step.value)) {
        run_ordered_contracted_block<EvalTrace>(
            *child, vmap, rich, leaf_evaluator, bs.cache, target, ctx,
            value_results, built, is_volatile, n_blocks, home_reads,
            needed_build);
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
      if (std::getenv("SEQUANT_UT_BLOCK_DIAG"))
        std::cerr << "[BLOCK] axis=" << toUtf8(block.axis.full_label())
                  << " batch=[" << e_lo << "," << e_hi << ") OUTPUT vid=" << vid
                  << " kind="
                  << (kind == OutputKind::AccumulateSum       ? "SUM"
                      : kind == OutputKind::AccumulateScatter ? "SCATTER"
                                                              : "?")
                  << std::endl;
      // PILLAR 2: this output is the use-site fetching (and slicing) shared
      // values this iteration -- name it as the current consumer so
      // slice_to_use binds a symmetric value's free mode to THIS output's own
      // mode (design sec.2). RAII-restored so a sibling output's fetch is not
      // mis-attributed.
      node_t const& out_eval_node = resolve(vid);
      auto const prev_consumer = bs.cache.current_consumer();
      bs.cache.set_current_consumer(out_eval_node->hash_value());
      ResultPtr part =
          evaluate_impl<EvalTrace>(out_eval_node, leaf_evaluator, bs.cache);
      bs.cache.set_current_consumer(prev_consumer);
      if (kind == OutputKind::AccumulateSum) {
        if (!acc[k])
          acc[k] = std::move(part);
        else
          acc[k]->add_inplace(*part);
      } else if (kind == OutputKind::AccumulateScatter) {
        // The full-extent zero destination is shaped by the node's own
        // (unsliced) index list; the backend realizes it from the spaces alone
        // -- no carrier, no Result-type reconciliation.
        if (!dest[k]) dest[k] = aops->make_zeros(resolve(vid)->canon_indices());
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

///
/// \brief Task 4 (multi-root single-DAG eval): the shared CORE every ordered
/// whole-forest entry point (\c evaluate_ordered_schedule's forest-wide SUM
/// and \c evaluate_ordered_multiroot's per-root MAP alike) delegates to --
/// walks \p ordered.root.steps exactly as \c evaluate_ordered_schedule always
/// has (Tasks 1-3: a root-level \c BuildStep built directly, a root-level \c
/// ScopeBlock realized via \c run_ordered_contracted_block) and returns each
/// forest root's own UNPERMUTED, already-built \c value_result, aligned
/// index-for-index with \p forest -- i.e. exactly the \c pre_results
/// \c combine_forest_roots expects, computed but NOT yet consumed by it. Pure
/// extraction of \c evaluate_ordered_schedule's original body (SP3 Tasks 1-4
/// through the original's own \c pre_results loop): no upstream logic
/// (schedule walk, \c ordered_n_blocks, \c ordered_home_reads, the
/// run-completeness assert) is touched -- see this file's own \note on why a
/// concatenated multi-root forest gets cross-root CSE for free from
/// \c compute_dag_boulevard's hash-keyed \c ValueCell bucketing (built
/// upstream of this function, in \p rich) with NO new dedup logic needed
/// here: a value shared across two \e independent root trees is just another
/// repeated hash, indistinguishable from a value shared across two summands
/// of one root's own forest, which this same walk already builds once.
///
/// \param forest Same requirement as \c evaluate_ordered_schedule's \p
///        forest: the roots whose results are computed by this call -- either
///        the summand terms of ONE equation (the pre-Task-4 caller) or
///        several INDEPENDENT equations' own root trees (Task 4's multi-root
///        caller); this function does not care which, since it produces one
///        unpermuted, unsummed result per element of \p forest either way.
/// \return Each element of \p forest's own already-built result, unpermuted,
///         same order and length as \p forest.
///
template <Trace EvalTrace = Trace::Default, meta::can_evaluate_range Nodes,
          typename F, typename N, bool FHC,
          typename ScopeGuardFactory = ::sequant::make_no_scope_guard>
  requires meta::leaf_node_evaluator<std::ranges::range_value_t<Nodes>, F>
[[nodiscard]] container::svector<ResultPtr> run_ordered_schedule_pre_results(
    Nodes const& forest, OrderedSchedule const& ordered,
    RichSchedule const& rich, F const& leaf_evaluator,
    CacheManager<N, FHC>& cache,
    std::function<std::size_t(Index const&)> const& target,
    [[maybe_unused]] ScopeGuardFactory const& make_scope_guard = {},
    std::function<bool(std::ranges::range_value_t<Nodes> const&)> const&
        is_volatile = {},
    std::function<std::size_t(std::size_t)> const& home_life_override = {}) {
  using node_t = std::ranges::range_value_t<Nodes>;
  static_assert(std::is_same_v<node_t, N>,
                "the forest's node type and the cache's node type must match");

  // Task 7 (dag-scope-runtime-slicing plan): self-wire the node->ModeToLevel
  // cache seam HERE, not only at the eval::evaluate dispatch wrapper
  // (scope_executor.hpp), which wires an equivalent map before calling
  // evaluate_ordered_schedule but is bypassed by any caller that invokes
  // evaluate_ordered_schedule / evaluate_ordered_multiroot directly (this
  // function is the shared core both delegate to) -- notably many unit tests
  // exercising the ordered executor at this level. slice_to_use's ordered-
  // path resolution (eval.hpp) now REQUIRES this map for the ordered push
  // style (BatchContextEntry::exact_axis == nullopt): with the old space-
  // fallback removed from that arm, an unwired seam would silently leave
  // such an entry unsliced rather than mis-slice it, so this is a
  // correctness dependency, not an optional convenience. Built from
  // `ordered`/`rich` -- already in hand -- via the same
  // detail::compute_cell_mode_to_level core populate_cell_mode_to_level uses
  // (that function additionally mirrors the result into each
  // ValueCell::mode_to_level for direct inspection by its own callers/tests;
  // this call site only needs the hash-keyed seam map). RAII-restored on
  // every exit (incl. exceptions), mirroring the wrapper's own
  // CellModeToLevelGuard so a caller's persistent, cross-iteration `cache`
  // never keeps a stale pointer into a map local to this call.
  container::svector<ModeToLevel> const per_cell_mode_to_level =
      compute_cell_mode_to_level(ordered, rich);
  std::unordered_map<std::size_t, ModeToLevel> cell_mode_to_level_map;
  cell_mode_to_level_map.reserve(rich.cells.size());
  for (std::size_t i = 0; i < rich.cells.size(); ++i)
    cell_mode_to_level_map.emplace(rich.cells[i].hash,
                                   per_cell_mode_to_level[i]);
  struct CellModeToLevelGuard {
    CacheManager<N, FHC>& c;
    std::unordered_map<std::size_t, ModeToLevel> const* prev;
    ~CellModeToLevelGuard() { c.set_cell_mode_to_level(prev); }
  } cell_mode_to_level_guard{cache, cache.cell_mode_to_level()};
  cache.set_cell_mode_to_level(&cell_mode_to_level_map);

  // Task 7 (sliced-value canonical-layout / loop-coloring design): the
  // loop-colored slice seam slice_to_use (eval.hpp) reads to resolve a
  // fetched value's physical slice mode off the loop-colored canonical layout
  // -- the SUCCESSOR to the mode_to_level seam above as the ordered path's
  // slice-mode source (sec.4). Built HERE, alongside that seam and from the
  // SAME already-in-hand `ordered`/`rich`, so every direct caller of the
  // shared ordered core (not only the eval::evaluate dispatch wrapper) gets
  // it. The per-(value,sliced-mode)->loop facts come from
  // compute_sliced_mode_assignment (value_id-keyed); projected here onto the
  // value hash slice_to_use has in hand at a fetch. RAII-restored on every
  // exit, mirroring the mode_to_level guard.
  SlicedModeAssignment const sliced_mode_assignment =
      compute_sliced_mode_assignment(ordered, rich);
  LoopColoredSliceSeam loop_colored_slice_seam;
  loop_colored_slice_seam.levels = sliced_mode_assignment.levels;
  loop_colored_slice_seam.by_hash.reserve(rich.cells.size());
  for (ValueCell const& c : rich.cells) {
    auto const it = sliced_mode_assignment.by_value.find(c.value_id);
    if (it != sliced_mode_assignment.by_value.end())
      loop_colored_slice_seam.by_hash.emplace(c.hash, it->second);
  }
  // PILLAR 2: project the per-occurrence (value, mode, loop, consumer) facts
  // onto the hash-keyed consumer-disambiguation map. value_id -> hash for both
  // the fetched value and the consumer use-site: at runtime slice_to_use has
  // the fetched node's hash (nd->hash_value()) and the current consumer's node
  // hash (CacheManager::current_consumer, set by this executor around each
  // evaluate_impl -- rich.cells[consumer_vid].hash IS that same node's hash for
  // the w8 case, where the consumer use-site is itself the evaluated output).
  // Only values with >1 mode under one loop ever consult this at runtime, so
  // recording every fact (single-mode values included) is harmless.
  for (auto const& [vid, mode, loop, consumer_vid] :
       sliced_mode_assignment.occ_facts) {
    SEQUANT_ASSERT(vid < rich.cells.size() && consumer_vid < rich.cells.size());
    loop_colored_slice_seam.by_hash_consumer[rich.cells[vid].hash].push_back(
        std::make_tuple(mode, loop, rich.cells[consumer_vid].hash));
  }
  struct LoopColoredSliceSeamGuard {
    CacheManager<N, FHC>& c;
    LoopColoredSliceSeam const* prev;
    ~LoopColoredSliceSeamGuard() { c.set_loop_colored_slice_seam(prev); }
  } loop_colored_slice_seam_guard{cache, cache.loop_colored_slice_seam()};
  cache.set_loop_colored_slice_seam(&loop_colored_slice_seam);

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
  std::function<std::size_t(Index const&)> const n_blocks = [&]() {
    PhaseTimer::Scope _pt("B.sched_setup");
    return ordered_n_blocks<node_t>(ordered, rich, vmap, leaf_evaluator, target,
                                    cache.array_ops());
  }();

  // Task 4: the EXACT home use-count for each homed value under read-from-home,
  // computed ONCE from the ordered schedule's realized scopes + the DAG (with
  // multiplicity) -- see detail::ordered_home_reads. This is the non-persistent
  // life of a volatile homed value: it frees at its genuine last use instead
  // of being pinned resident. A test may inject home_life_override (a
  // non-evicting life) to MEASURE the actual home reads and assert them equal
  // to home_reads -- the "predicted == measured" gate.
  std::function<std::size_t(std::size_t)> const home_reads =
      home_life_override ? home_life_override : [&]() {
        PhaseTimer::Scope _pt("B.sched_setup");
        return ordered_home_reads<node_t>(ordered, rich, vmap, n_blocks);
      }();

  // Per-value results, indexed by value_id (== a ValueCell's own slot in
  // rich.cells -- see peak_profile.hpp's ValueCell::value_id doc comment),
  // so the pre_results resolution below reads a forest root's own build
  // directly rather than re-resolving it through the cache.
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
  // the pre_results resolution below resolves a forest root produced inside a
  // loop block exactly like one produced by a plain root BuildStep.
  typename CacheManager<N, FHC>::BatchContext const root_ectx;
  // A forest root's ONLY reader is the pre_results combine below -- it has zero
  // DAG consumers (else it would not be a root). That combine read is a real
  // read of the value from its home, so it is charged here as +1 on the root's
  // home life: with it, evaluate_impl's own production store-access drains the
  // root to a still-live cache entry, and the combine's cache read drains it to
  // zero -- one holder (the cache), exact use-count, roots no longer special.
  // Without it (life == 1) the store-access drains the root immediately and the
  // cache holds nothing, which is why the executor used to keep a PARALLEL
  // value_results reference to every produced value -- a ref living PAST the
  // cache node removal (it held every intermediate to end-of-iteration even
  // after the cache drained it at its genuine last use, absent in forest
  // descent). We no longer retain it: non-root consumers read from the cache
  // (home_reads life), and roots read from the cache too (this +1). Escape
  // outputs produced inside a batch block are still mirrored into value_results
  // by run_ordered_contracted_block (their home is one scope out); pre_results
  // prefers that mirror when present and falls back to the cache otherwise.
  container::set<std::size_t> forest_root_hashes;
  for (auto&& n : forest) forest_root_hashes.insert(n->hash_value());

  // "Needed this iteration" gate. The schedule lists every BuildStep, but a
  // non-persistent intermediate whose consumers are ALL persistent-and-cached
  // is read only in iteration 1 (when those persistent consumers are built);
  // thereafter the consumers are cache hits and never re-read it, so rebuilding
  // it each iteration is wasted work (the reason the DAG used to over-persist
  // it via !vol). Mirror forest descent's "stop at cache hits": BFS from the
  // volatile roots, descending only through nodes NOT currently alive in the
  // cache. An alive (persistent, still-holding) node is read but not rebuilt
  // and does NOT propagate need to its children. Computed ONCE per call against
  // the cache state at entry (persistent survivors alive, everything else
  // drained by the per-term reset), so in iteration 1 (nothing cached) it
  // admits every node, and in later iterations it prunes the
  // persistent-shadowed subtrees.
  container::set<std::size_t> needed_build;
  {
    container::svector<node_t> stack;
    container::set<std::size_t> visited;
    for (auto&& n : forest) stack.push_back(n);
    while (!stack.empty()) {
      node_t const n = stack.back();
      stack.pop_back();
      if (n.leaf()) continue;
      if (!visited.insert(n->hash_value()).second) continue;
      if (cache.alive(n)) continue;  // cache hit: read, not rebuilt, no descend
      needed_build.insert(n->hash_value());
      stack.push_back(n.left());
      stack.push_back(n.right());
    }
  }

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
      bool const is_root = forest_root_hashes.count(hash) != 0;
      // Skip a BuildStep the "needed" gate pruned: its value is not read this
      // iteration (all consumers are persistent cache hits), so rebuilding it
      // is wasted work. Roots are always needed (volatile, never cached). A
      // skipped value is intentionally not produced this iteration -- mark it
      // accounted for so the completeness check does not mistake the prune for
      // a gap.
      if (!it->second.leaf() && !needed_build.count(hash)) {
        built[vid] = 1;
        continue;
      }
      if (!it->second.leaf()) {
        if (is_volatile) {
          // Persistence is the shared cache's classification (non-volatile AND
          // has a volatile DIRECT consumer), which sequant::cache_manager
          // already computed and stamped on the entry -- NOT the local !vol,
          // which over-enrolls a non-volatile value with no volatile consumer.
          // With the "needed" gate above, un-persisting such a value no longer
          // forces a rebuild: it is built once (iteration 1) and pruned after.
          bool const persistent = cache.entry_is_persistent(it->second);
          // +1 for the pre_results combine read of a forest root (see the block
          // comment above); non-roots carry only their DAG consumers' reads.
          std::size_t const life = home_reads(vid) + (is_root ? 1u : 0u);
          cache.ensure_home_slot(it->second, life, persistent);
        } else {
          // No volatility policy: the root-homed value is pinned resident until
          // reset() anyway, so it is trivially still in the cache for the
          // combine read; no life bump is needed.
          cache.ensure_home_slot(it->second);
        }
      }
      // Build the value; it self-stores into its cache home (evaluate_impl's
      // finish_phase_b). Discard the returned pointer: the cache is the single
      // holder now -- consumers (and, for a root, the combine) read it back
      // from the cache. Retaining it in value_results is exactly the dead
      // reference that kept every intermediate alive past its cache drain.
      (void)evaluate_impl<EvalTrace>(it->second, leaf_evaluator, cache);
      built[vid] = 1;  // R3: this root-scope value is produced.
    } else if (auto const* block = std::get_if<ScopeBlock>(&step.value)) {
      run_ordered_contracted_block<EvalTrace>(
          *block, vmap, rich, leaf_evaluator, cache, target, root_ectx,
          value_results, built, is_volatile, n_blocks, home_reads,
          &needed_build);
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
    collect_production_ids(ordered.root, production_ids);
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
  // per-root pre_results below.
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
    std::size_t const vid = vid_it->second;
    if (value_results[vid]) {
      // Produced inside a batch block and mirrored out by
      // run_ordered_contracted_block (its home is one scope out of this level).
      pre_results[i] = std::move(value_results[vid]);
    } else {
      // Produced by a plain root BuildStep: it lives in the cache (self-stored
      // by evaluate_impl, kept for this combine read by the +1 on its home
      // life). access_at drains that last read and hands back the CANONICAL
      // stored value; orient it to this root's phase, matching evaluate_impl's
      // own canonical->orientation return convention (apply_phase).
      auto ptr = cache.access_at(roots[i]).ptr;
      SEQUANT_ASSERT(ptr &&
                     "evaluate_ordered_schedule: forest root not resident in "
                     "the cache at the combine read");
      auto const ph = roots[i]->canon_phase();
      pre_results[i] = (ph == 1) ? std::move(ptr) : ptr->mult_by_phase(ph);
    }
    SEQUANT_ASSERT(pre_results[i] &&
                   "evaluate_ordered_schedule: forest root was never "
                   "produced");
  }

  return pre_results;
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

  // Task 4: the schedule walk itself (Tasks 1-3, plus the block-count/
  // home-reads setup and the run-completeness assert) is factored into
  // detail::run_ordered_schedule_pre_results, shared byte-for-byte with
  // evaluate_ordered_multiroot below -- see that function's own doc comment.
  // Nothing about the walk changes here; only what happens to its per-root
  // pre_results output differs (summed here, mapped there).
  container::svector<ResultPtr> pre_results =
      detail::run_ordered_schedule_pre_results<EvalTrace>(
          forest, ordered, rich, leaf_evaluator, cache, target,
          make_scope_guard, is_volatile, home_life_override);

  container::svector<node_t> roots;
  for (auto&& n : forest) roots.push_back(n);

  // -------- Shared combine: permute each root to layout and sum. --------
  // combine_forest_roots (forest_combine.hpp) is the SAME helper
  // evaluate_whole_scope's tail calls, so the two executors emit
  // byte-identical Term/Permute/SumInplace trace bookkeeping without a
  // hand-synced second copy.
  return combine_forest_roots<EvalTrace>(roots, pre_results, layout, cache);
}

///
/// \brief Task 4 of the multi-root single-DAG eval plan: the MULTI-ROOT
/// ordered entry point. Identical inputs to \c evaluate_ordered_schedule
/// (same \p roots/\p ordered/\p rich/\p leaf_evaluator/\p cache/\p target
/// contract -- \p ordered and \p rich must already have been built from the
/// SAME concatenated \p roots, e.g. via \c compute_dag_boulevard(roots, ...)
/// + \c build_ordered_schedule, exactly as any \c evaluate_ordered_schedule
/// caller already does for its own forest), EXCEPT \p layout widens to \p
/// layouts, one entry per root (Task 7 of the plan): unlike a single summed
/// forest, independent roots need not share a layout (e.g. distinct CC
/// residual annotations R1 `{a;i}` vs R2 `{ab;ij}`). Returns a MAP instead
/// of a SUM: each root's own (permuted to its own \p layouts entry, but NOT
/// cross-root-accumulated) result, aligned index-for-index with \p roots.
///
/// \details Runs the identical schedule walk \c evaluate_ordered_schedule
/// runs (via the same \c detail::run_ordered_schedule_pre_results this
/// function and \c evaluate_ordered_schedule both delegate to -- \c
/// ordered_home_reads, \c ordered_n_blocks, and the run-completeness assert
/// are UNCHANGED, exercised over the combined schedule exactly as they
/// already are for a multi-summand forest), so a value shared across two
/// \e independent roots is built exactly once for the identical reason a
/// value shared across two summands of one root's own forest already is:
/// \c compute_dag_boulevard buckets by structural hash into one \c ValueCell
/// regardless of which root(s) reference it, and \c well_formed's static
/// single-producer invariant (plus this walk's own run-completeness assert)
/// guarantees the schedule realizes exactly one \c BuildStep for it. The
/// ONLY change from \c evaluate_ordered_schedule is the tail: instead of one
/// forest-wide \c combine_forest_roots call (which sums every root into a
/// single \c ResultPtr), this calls \c combine_forest_roots once PER ROOT,
/// each with a singleton \c {roots[i]}/{pre_results[i]} pair -- reusing the
/// IDENTICAL per-root Term/Permute trace bookkeeping \c combine_forest_roots
/// already emits for each root of a summed forest, but never reaching its
/// cross-root \c add_inplace (a singleton call never has a second root to
/// sum against), so no roots are summed together. For a single-root \p roots
/// input this reduces to a singleton \c combine_forest_roots call byte-
/// identical to what \c evaluate_ordered_schedule would run for that same
/// one-root forest -- the regression anchor: \c evaluate_ordered_multiroot(
/// {r}, ...) == { evaluate_ordered_schedule({r}, ...) }.
///
/// \return One \c ResultPtr per element of \p roots, in \p roots' own order,
///         each permuted to its own \p layouts entry -- NOT summed together.
///
template <Trace EvalTrace = Trace::Default, meta::can_evaluate_range Nodes,
          typename F, typename N, bool FHC,
          typename ScopeGuardFactory = ::sequant::make_no_scope_guard>
  requires meta::leaf_node_evaluator<std::ranges::range_value_t<Nodes>, F>
[[nodiscard]] container::svector<ResultPtr> evaluate_ordered_multiroot(
    Nodes const& roots, OrderedSchedule const& ordered,
    RichSchedule const& rich, container::svector<std::string> const& layouts,
    F const& leaf_evaluator, CacheManager<N, FHC>& cache,
    std::function<std::size_t(Index const&)> const& target,
    [[maybe_unused]] ScopeGuardFactory const& make_scope_guard = {},
    std::function<bool(std::ranges::range_value_t<Nodes> const&)> const&
        is_volatile = {},
    std::function<std::size_t(std::size_t)> const& home_life_override = {}) {
  using node_t = std::ranges::range_value_t<Nodes>;

  container::svector<ResultPtr> pre_results =
      detail::run_ordered_schedule_pre_results<EvalTrace>(
          roots, ordered, rich, leaf_evaluator, cache, target, make_scope_guard,
          is_volatile, home_life_override);

  container::svector<node_t> root_nodes;
  for (auto&& n : roots) root_nodes.push_back(n);
  SEQUANT_ASSERT(pre_results.size() == root_nodes.size());
  SEQUANT_ASSERT(layouts.size() == root_nodes.size());

  // Per-root combine: reuses combine_forest_roots' Term/Permute bookkeeping
  // one root at a time (a singleton call never reaches its cross-root
  // add_inplace), so each root is permuted to its OWN layouts[i] exactly as
  // it would be as part of a summed forest, but no two roots are ever
  // accumulated together -- the map, not the sum.
  container::svector<ResultPtr> results(root_nodes.size());
  for (std::size_t i = 0; i != root_nodes.size(); ++i) {
    container::svector<node_t> one_root{root_nodes[i]};
    container::svector<ResultPtr> one_pre{std::move(pre_results[i])};
    results[i] =
        combine_forest_roots<EvalTrace>(one_root, one_pre, layouts[i], cache);
  }
  return results;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_ORDERED_EXECUTOR_HPP
