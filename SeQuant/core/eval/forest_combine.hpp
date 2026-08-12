#ifndef SEQUANT_EVAL_FOREST_COMBINE_HPP
#define SEQUANT_EVAL_FOREST_COMBINE_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <any>
#include <array>
#include <cstddef>
#include <string>

namespace sequant::eval {

///
/// \brief The shared per-root COMBINE step every whole-forest executor ends
/// with: permute each root's own (already-built) result to \p layout, then
/// sum across roots -- the identical Term Begin/End boundary, Permute
/// EvalStat, and cross-root \c add_inplace + SumInplace EvalStat trace/hwmark
/// bookkeeping \c sequant::evaluate(Nodes const&, layout, ...) (eval.hpp)
/// emits.
///
/// \details Factored out of \c evaluate_whole_scope's tail (its original, and
/// still primary, caller -- scope_executor.hpp) so a second whole-forest
/// executor -- \c evaluate_ordered_schedule (ordered_executor.hpp) -- can
/// reuse the IDENTICAL bookkeeping rather than hand-syncing a second copy
/// that would silently drift as later tasks extend either executor. Lives in
/// its own header (not inline in \c scope_executor.hpp, where the code
/// originated) for the same reason \c build_value_node_map does
/// (value_node_map.hpp's own doc comment): \c scope_executor.hpp's
/// coexistence driver dispatches INTO \c ordered_executor.hpp, so the reverse
/// direction (\c ordered_executor.hpp including \c scope_executor.hpp for
/// just this helper) would cycle. Both headers include this one instead.
///
/// \param roots The forest's own top-level trees, in forest order.
/// \param pre_results Each root's own UNPERMUTED, already-built result,
///        aligned index-for-index with \p roots. Moved out of on use (each
///        entry is consumed), so the caller's own copy is left empty after
///        this call returns.
/// \param layout The layout each root's result is permuted to before being
///        summed; same meaning as \c sequant::evaluate(Node const&, layout,
///        ...)'s \p layout.
/// \param cache The cache consulted for the trace/hwmark bookkeeping only
///        (byte counts, chain residency) -- no value is built or stored by
///        this function.
/// \return The summed, per-root-permuted result.
///
template <Trace EvalTrace, meta::can_evaluate node_t, typename N, bool FHC>
[[nodiscard]] ResultPtr combine_forest_roots(
    container::svector<node_t> const& roots,
    container::svector<ResultPtr>& pre_results, auto const& layout,
    CacheManager<N, FHC>& cache) {
  SEQUANT_ASSERT(pre_results.size() == roots.size());

  // if the layout is not the default constructed value need to permute --
  // mirrors sequant::evaluate(Node const&, layout, ...)'s identical check.
  bool const perm = layout != decltype(layout){};

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
        auto const stat = log::EvalStat{
            .mode = log::EvalMode::Permute,
            .time = permute_time,
            .mem_result = log::bytes(post),
            .mem_alloc = log::bytes(post),
            .mem_hwmark = {cache.note_working_set(hwmark, n->hash_value())}};
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
      auto const stat = log::EvalStat{
          .mode = log::EvalMode::SumInplace,
          .time = sum_time,
          .mem_result = log::bytes(result),
          .mem_alloc = {0},
          .mem_hwmark = {cache.note_working_set(hwmark, n->hash_value())}};
      log::eval(stat, n->label());
    }
  }

  return result;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_FOREST_COMBINE_HPP
