#ifndef SEQUANT_EVAL_SCOPE_EXECUTOR_HPP
#define SEQUANT_EVAL_SCOPE_EXECUTOR_HPP

#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/eval/scope_schedule.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <any>
#include <array>

namespace sequant::eval {

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
/// on exit) around each per-value \c evaluate_impl call.
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
/// \return The summed, per-root-permuted result, as \c sequant::evaluate(
///         Nodes const&, layout, ...) would produce for the same \p forest.
///
template <Trace EvalTrace = Trace::Default, meta::can_evaluate_range Nodes,
          typename F, typename N, bool FHC>
  requires meta::leaf_node_evaluator<std::ranges::range_value_t<Nodes>, F>
ResultPtr evaluate_whole_scope(Nodes const& forest, ScopeSchedule const& sched,
                               auto const& layout, F const& leaf_evaluator,
                               CacheManager<N, FHC>& cache) {
  SEQUANT_ASSERT(
      sched.root.children.empty() &&
      "evaluate_whole_scope (Task 2) only implements the top-scope walk: a "
      "scope tree with realized batch loops (non-empty sched.root.children) "
      "needs the Task 3+ batch-block recursion, not implemented here.");

  // if the layout is not the default constructed value need to permute --
  // mirrors sequant::evaluate(Node const&, layout, ...)'s identical check.
  bool const perm = layout != decltype(layout){};

  ResultPtr result;
  for (auto&& n : forest) {
    ResultPtr pre = evaluate_impl<EvalTrace>(n, leaf_evaluator, cache);
    ResultPtr post =
        perm ? pre->permute(std::array<std::any, 2>{n->annot(), layout}) : pre;

    if (!result) {
      result = post;
      continue;
    }
    result->add_inplace(*post);
  }

  return result;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_SCOPE_EXECUTOR_HPP
