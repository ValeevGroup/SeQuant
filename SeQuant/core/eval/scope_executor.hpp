#ifndef SEQUANT_EVAL_SCOPE_EXECUTOR_HPP
#define SEQUANT_EVAL_SCOPE_EXECUTOR_HPP

#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/eval/scope_schedule.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <any>
#include <array>
#include <cstddef>
#include <functional>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sequant::eval {

///
/// \brief The value_id -> forest-node bridge (design integration point 1).
///
/// \details \c ScopeNode::homed_values are \c ValueCell::value_id's; a value's
/// eval node is recovered through the \c ValueCell::hash it carries -- exactly
/// the \c EvalExpr::hash_value() identity \c CacheManager dedups by. This maps
/// every distinct forest node by that hash, so a \c value_id resolves as
/// `map[rich.cells[value_id].hash]`. Under perfect CSE many occurrences share
/// one hash; a single representative node (the first visited) is kept, which is
/// all a homed value needs (every occurrence is the same value). Pure lookup
/// construction -- no execution.
///
template <meta::eval_node_range R>
[[nodiscard]] std::unordered_map<std::size_t, std::ranges::range_value_t<R>>
build_value_node_map(R const& forest) {
  using node_t = std::ranges::range_value_t<R>;
  std::unordered_map<std::size_t, node_t> out;
  auto visit = [&out](auto&& self, node_t const& n) -> void {
    out.emplace(n->hash_value(), n);
    if (!n.leaf()) {
      self(self, n.left());
      self(self, n.right());
    }
  };
  for (auto const& t : forest) visit(visit, t);
  return out;
}

///
/// \brief Task 2 of the whole-scope batched DAG execution design (see
/// `doc/dev/specs/2026-08-10-whole-scope-batched-dag-execution-design.md`):
/// the executor SKELETON -- a scope-tree walk that, for a forest with NO
/// batch loops (the scope tree is root-only), evaluates every forest root and
/// sums the results. Provably equivalent to the existing per-tree forest
/// descent (\c sequant::evaluate(Nodes const&, ...), eval.hpp) for that case.
/// Tasks 3-4 extend this to recurse into \c sched.root.children (one
/// accumulate/scatter block per realized batch loop); NONE of that batch-block
/// machinery lives here yet.
///
/// \details \c ScopeSchedule (Task 1) assigns every distinct forest VALUE --
/// not just forest roots, every hash-distinct node of the whole forest -- to
/// a \c ScopeNode by \c value_id, an index into the \c RichSchedule that built
/// it (see \c peak_profile.hpp's \c ValueCell::occurrences, which is the
/// value_id -> forest-node bridge in general). This function's signature,
/// though, is handed only the projected \c ScopeSchedule (no \c RichSchedule /
/// \c ValueCell in scope), so it cannot walk \c homed_values back to a forest
/// node in general here. For the UNBATCHED case that general bridge would
/// resolve to anyway is exactly "the forest's own top-level trees" (a root
/// occurrence, identified by \c OccurrenceRec::point == \c consumer_point, is
/// precisely a tree of the forest -- see \c scope_schedule.hpp's \c
/// detail::mode_is_external for the same identification), so this walk drives
/// directly off \p forest, per the Task-2 brief's ambiguity-resolution note.
/// \p sched is still consulted -- its shape gates dispatch (asserted below) --
/// so the entry point already has the scope-tree-walk SHAPE Tasks 3-4 extend
/// by recursing into \c sched.root.children; only the general value_id ->
/// node resolution is deferred to whichever later task plumbs \c
/// RichSchedule through this signature.
///
/// Reuses \c evaluate_impl (the single-node recursive engine, eval.hpp) for
/// every per-value build; no per-op evaluation is reimplemented here. The
/// per-root permute-to-\p layout and cross-root \c add_inplace accumulation
/// mirror \c sequant::evaluate(Node const&, layout, ...) and \c
/// sequant::evaluate(Nodes const&, layout, ...) respectively, inlined (rather
/// than delegated to those overloads) so that Tasks 3-4 have a natural seam
/// to inject batch-context management (slicing on entry, accumulate/scatter
/// on exit) around each per-value \c evaluate_impl call. The trace/hwmark
/// bookkeeping (\c log::term boundaries, the Permute \c EvalStat, the
/// SumInplace \c EvalStat) is reproduced line-for-line from those two
/// overloads too -- \c EvalTrace-gated, so it costs nothing when off, but
/// present so Tasks 3-4's batch-block trace/hwmark plumbing has a faithful
/// root-only baseline to extend rather than a silent gap to rediscover.
///
/// \param forest The forest whose per-root results are evaluated and summed;
///        same requirement as \c sequant::evaluate(Nodes const&, ...).
/// \param sched The scope-tree schedule (\c build_scope_schedule, Task 1)
///        for \p forest. Task 2 requires \c sched.root.children to be empty
///        (no realized batch loop): a non-empty scope tree needs the Task 3+
///        walk this function does not yet implement.
/// \param layout The layout each root's result is permuted to before being
///        summed; same meaning as \c sequant::evaluate(Node const&, layout,
///        ...)'s \p layout.
/// \param leaf_evaluator The leaf evaluator, as in \c sequant::evaluate.
/// \param cache The cache for common sub-expression elimination, as in \c
///        sequant::evaluate.
/// \param target_batch_size Per-index batch partition size (elements), the
///        source of the K partition for a realized batch loop; unused for the
///        root-only case. Same meaning as \c make_batched_custom_evaluator's
///        \p target_batch_size.
/// \param make_scope_guard Backend scope-guard factory, called with the batch
///        count and its RAII result HELD for the entire K-block loop -- exactly
///        as \c make_batched_custom_evaluator's contracted path does
///        (eval.hpp:2157). A screening backend uses it to relax block-sparse
///        screening scaled by the batch count so a contribution significant
///        over the FULL K range is not dropped within an individual batch (see
///        \c no_scope_guard). Defaults to \c make_no_scope_guard (a no-op for
///        the dense / DryRun backends), keeping every existing caller
///        byte-identical.
/// \return The summed, per-root-permuted result, as \c sequant::evaluate(
///         Nodes const&, layout, ...) would produce for the same \p forest.
///
template <Trace EvalTrace = Trace::Default, meta::can_evaluate_range Nodes,
          typename F, typename N, bool FHC,
          typename ScopeGuardFactory = ::sequant::make_no_scope_guard>
  requires meta::leaf_node_evaluator<std::ranges::range_value_t<Nodes>, F>
ResultPtr evaluate_whole_scope(
    Nodes const& forest, ScopeSchedule const& sched, auto const& layout,
    F const& leaf_evaluator, CacheManager<N, FHC>& cache,
    std::function<std::size_t(Index const&)> target_batch_size = {},
    ScopeGuardFactory make_scope_guard = {}) {
  using node_t = std::ranges::range_value_t<Nodes>;
  static_assert(std::is_same_v<node_t, N>,
                "the forest's node type and the cache's node type must match");

  // Task 3 handles the ROOT-ONLY tree (Task 2) plus ONE realized batch loop --
  // a single aux/contracted child scope. A deeper / branching scope tree needs
  // the Task 4+ recursion, not implemented here.
  SEQUANT_ASSERT(
      sched.root.children.size() <= 1 &&
      "evaluate_whole_scope (Task 3) implements the top scope plus at most ONE "
      "batch loop (one sched.root.children entry); a deeper/branching scope "
      "tree needs the Task 4+ batch-block recursion, not implemented here.");

  // if the layout is not the default constructed value need to permute --
  // mirrors sequant::evaluate(Node const&, layout, ...)'s identical check.
  bool const perm = layout != decltype(layout){};

  // The per-root, UNPERMUTED result, in forest order. The root-only case fills
  // each entry with a direct evaluate_impl; the single-loop case fills the
  // K-carrying roots from the accumulated batch loop (below) and the
  // K-invariant roots directly. The final permute-to-layout + cross-root
  // add_inplace combine (which reproduces sequant::evaluate(Nodes const&,
  // layout, ...)'s trace/hwmark bookkeeping line-for-line) is SHARED by both
  // cases, so Tasks 2-3 emit identical Term/SumInplace records.
  container::svector<node_t> roots;
  for (auto&& n : forest) roots.push_back(n);
  container::svector<ResultPtr> pre_results(roots.size());

  if (sched.root.children.empty()) {
    // -------- Task 2 top-scope walk: every root built directly. --------
    for (std::size_t i = 0; i != roots.size(); ++i)
      pre_results[i] =
          evaluate_impl<EvalTrace>(roots[i], leaf_evaluator, cache);
  } else {
    // -------- Task 3 single batch loop (aux Κ, contracted). --------
    //
    // The child scope is ONE realized contracted batch loop. Its mode K is
    // summed away below every forest root that carries it, so those roots are
    // homed at the ROOT scope (K not on their result) but built INCREMENTALLY:
    // per K-block partials that ACCUMULATE. Every value HOMED at the K scope
    // (the shared K-carrying gC composites) is built once per block on ONE
    // scratch shared by ALL K-carrying roots -- so a composite shared across
    // trees is built once per block, not once per consumer group. This is the
    // exhaustive, whole-forest generalization of the OPPORTUNISTIC
    // co-evaluation group in make_batched_custom_evaluator (eval.hpp): its
    // trigger-seeded group becomes the explicit whole-forest member set the
    // scope tree names, and the reused primitives (make_batched_scratch,
    // per-block evaluate_impl with Enter-stage slice-on-use, and cross-block
    // add_inplace) are DRIVEN from this walk rather than re-implemented.
    ScopeNode const& kscope = sched.root.children.front();
    SEQUANT_ASSERT(
        kscope.children.empty() &&
        "evaluate_whole_scope (Task 3): the single batch loop must be a LEAF "
        "scope (no nested loops); a nested loop nest is Task 4+.");
    SEQUANT_ASSERT(
        kscope.kind == BatchModeType::Contracted &&
        "evaluate_whole_scope (Task 3) realizes a CONTRACTED (accumulate) aux "
        "loop; an External (scatter) loop is a later increment.");
    SEQUANT_ASSERT(target_batch_size &&
                   "evaluate_whole_scope with a realized batch loop needs a "
                   "target_batch_size (the batch partition source).");
    Index const K = kscope.mode;

    using member_t = std::pair<node_t const*, Index>;

    // Partition the roots: those that CARRY K (their result contracts K below
    // them -- the co-evaluated members) vs those invariant to K (built once,
    // like the root-only path).
    container::svector<std::size_t> k_root_idx;
    container::svector<std::size_t> plain_root_idx;
    for (std::size_t i = 0; i != roots.size(); ++i) {
      if (sequant::find_leaf_carrying(roots[i], K))
        k_root_idx.push_back(i);
      else
        plain_root_idx.push_back(i);
    }

    // K-invariant roots: built directly on the root cache, exactly as the
    // top-scope path builds a root.
    for (auto i : plain_root_idx)
      pre_results[i] =
          evaluate_impl<EvalTrace>(roots[i], leaf_evaluator, cache);

    if (!k_root_idx.empty()) {
      // The K partition: identical (whole-tile element ranges) across every
      // member since K is one global aux mode. Read it once from any carrier
      // leaf, exactly as pick_sliceable does in make_batched_custom_evaluator.
      container::svector<std::pair<std::size_t, std::size_t>> batches;
      {
        auto const lf =
            sequant::find_leaf_carrying(roots[k_root_idx.front()], K);
        SEQUANT_ASSERT(lf);
        batches = leaf_evaluator(lf->first)->mode_batches(lf->second,
                                                          target_batch_size(K));
      }

      // The whole-forest co-evaluation group: every K-carrying root, paired
      // with K. ONE scratch dedups the shared K-homed composites across ALL
      // members (make_batched_scratch registers every repeated, consistently
      // sliced internal subnode), so each is built once per block.
      std::vector<member_t> group;
      group.reserve(k_root_idx.size());
      for (auto i : k_root_idx) group.emplace_back(&roots[i], K);

      auto bs = sequant::detail::make_batched_scratch(group, cache);
      for (auto const* s : bs.seeds) (void)bs.cache.store(*s, cache.access(*s));
      // Chain the scratch under the root cache so an above-homed (K-invariant)
      // value referenced mid-loop materializes lazily at its home and is sliced
      // on use, and so per-build tallies route to the root cache.
      bs.cache.set_parent(&cache);

      // RAII scope for the batched partial contractions, HELD for the ENTIRE
      // K-block loop -- mirroring make_batched_custom_evaluator's contracted
      // path (eval.hpp:2157). A screening backend's factory relaxes
      // block-sparse screening scaled by the batch count, so a contribution
      // whose norm clears the threshold over the FULL K range (but not within
      // one batch) is not dropped in every batch and lost from the K-sum. A
      // no-op for the dense / DryRun backends (make_no_scope_guard), so this is
      // byte-identical there.
      auto const scope_guard = make_scope_guard(batches.size());
      (void)scope_guard;

      // Accumulate: each block's per-root partial adds into the root's running
      // K-sum. `sum_K = sum_blocks sum_{K in block}`, exact.
      std::vector<ResultPtr> acc(group.size());
      for (auto const& [e_lo, e_hi] : batches) {
        if (e_lo == e_hi) continue;
        bs.cache.reset();
        // The enclosing batch context is just THIS block (the loop is at the
        // top scope). Setting it makes evaluate_impl's Enter-stage
        // slice-on-use slice every K-carrying leaf/value to the block.
        typename CacheManager<N, FHC>::BatchContext ctx;
        ctx.push_back({K, {e_lo, e_hi}});
        bs.cache.set_batch_context(std::move(ctx));
        for (std::size_t m = 0; m != group.size(); ++m) {
          ResultPtr part = evaluate_impl<EvalTrace>(*group[m].first,
                                                    leaf_evaluator, bs.cache);
          if (!acc[m])
            acc[m] = std::move(part);
          else
            acc[m]->add_inplace(*part);
        }
      }
      for (std::size_t m = 0; m != group.size(); ++m) {
        SEQUANT_ASSERT(acc[m]);
        pre_results[k_root_idx[m]] = std::move(acc[m]);
      }
      // Free the K scratch (bs goes out of scope here).
    }
  }

  // -------- Shared combine: permute each root to layout and sum. --------
  ResultPtr result;
  for (std::size_t i = 0; i != roots.size(); ++i) {
    node_t const& n = roots[i];
    // Per-root term boundary -- mirrors sequant::evaluate(Node const&,
    // layout, ...)'s identical log::term(Begin/End) bracket around its
    // evaluate_impl + permute.
    std::string xpr;
    if constexpr (sequant::detail::trace(EvalTrace)) {
      xpr = toUtf8(io::serialization::to_string(to_expr(n)));
      log::term(log::TermMode::Begin, xpr);
    }

    ResultPtr const pre = std::move(pre_results[i]);
    SEQUANT_ASSERT(pre);

    ResultPtr post;
    auto const permute_time = sequant::detail::timed_eval_inplace([&]() {
      post = perm ? pre->permute(std::array<std::any, 2>{n->annot(), layout})
                  : pre;
    });

    SEQUANT_ASSERT(post);

    // Permute EvalStat + term-End -- byte-for-byte the logging
    // sequant::evaluate(Node const&, layout, ...) emits around the identical
    // permute() call above.
    if constexpr (sequant::detail::trace(EvalTrace)) {
      if (perm) {
        size_t hwmark = log::bytes(cache, post).value;
        if (!cache.chain_holds(pre)) hwmark += log::bytes(pre).value;
        hwmark += cache.parent() ? cache.parent()->chain_residency() : 0;
        auto const stat =
            log::EvalStat{.mode = log::EvalMode::Permute,
                          .time = permute_time,
                          .mem_result = log::bytes(post),
                          .mem_alloc = log::bytes(post),
                          .mem_hwmark = {cache.note_working_set(hwmark)}};
        log::eval(stat, n->label());
      }
      log::term(log::TermMode::End, xpr);
    }

    if (!result) {
      result = post;
      continue;
    }

    // Cross-root accumulation -- mirrors sequant::evaluate(Nodes const&,
    // layout, ...)'s identical timed add_inplace() + SumInplace EvalStat
    // (there, `post` above is that overload's per-node `pre`).
    auto const sum_time = sequant::detail::timed_eval_inplace(
        [&]() { result->add_inplace(*post); });

    if constexpr (sequant::detail::trace(EvalTrace)) {
      size_t hwmark = log::bytes(cache, result).value;
      if (!cache.chain_holds(post)) hwmark += log::bytes(post).value;
      hwmark += cache.parent() ? cache.parent()->chain_residency() : 0;
      auto const stat =
          log::EvalStat{.mode = log::EvalMode::SumInplace,
                        .time = sum_time,
                        .mem_result = log::bytes(result),
                        .mem_alloc = {0},
                        .mem_hwmark = {cache.note_working_set(hwmark)}};
      log::eval(stat, n->label());
    }
  }

  return result;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_SCOPE_EXECUTOR_HPP
