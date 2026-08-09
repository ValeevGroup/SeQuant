#ifndef SEQUANT_EVAL_PLACEMENT_REMAT_HPP
#define SEQUANT_EVAL_PLACEMENT_REMAT_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/occurrence_key.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/eval/placement_router.hpp>
#include <SeQuant/core/eval/slicing_signature.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <cstddef>
#include <optional>
#include <ranges>
#include <unordered_map>
#include <utility>

namespace sequant::eval {

///
/// \brief The rematerialization placement pass (spec item "O2" of design doc
/// `doc/dev/specs/2026-08-04-phase4-threshold-o2-cutover-design.md`): a
/// whole-forest greedy pass that relaxes the peak-maximal perfect-CSE seed to
/// fit a peak budget by REMATERIALIZING values -- rebuilding them in smaller
/// (sliced) or shorter-lived pieces rather than holding them full -- trading
/// recompute for peak. In register-allocation terms it is rematerialization
/// (https://en.wikipedia.org/wiki/Rematerialization), not spilling: no value is
/// moved to slower storage; pressure drops because a value is recomputed nearer
/// its use. STANDALONE in Phase 4a (zero production
/// callers). It works on the rich per-value \c ValueCell (from \c
/// compute_dag_boulevard); \c home_modes is the field a remat move mutates (a
/// SHRINK adds a demoted carried batch mode, block-sizing it -- the zero-cost
/// remat case). This header lands the flat input/projection (\c RematInput, \c
/// to_schedule), the candidate query (\c shrink_candidates), the SHRINK move,
/// and the greedy loop (\c rematerialize_to_budget).
///

///
/// \brief The flat input to the remat pass: every \c ValueCell of a linearized
/// forest, plus the timeline length needed to re-sweep it.
///
struct RematInput {
  container::svector<ValueCell> cells;
  std::size_t num_points = 0;
};

///
/// \brief Seed the remat pass's working cells from an eval \p forest: the
/// \c ValueCell records of \c compute_dag_boulevard, plus \c num_points.
///
/// \p block_of is any `Index -> std::size_t` callable giving the block
/// (sliced) element count for a mode; forwarded to \c compute_dag_boulevard.
///
template <meta::eval_node_range R>
RematInput remat_cells(R const& forest, dryrun::CostModel const& cm,
                       auto const& block_of) {
  RichSchedule rich = compute_dag_boulevard(forest, cm, block_of);
  RematInput out;
  out.cells.assign(std::make_move_iterator(rich.cells.begin()),
                   std::make_move_iterator(rich.cells.end()));
  out.num_points = rich.num_points;
  return out;
}

///
/// \brief Project the remat pass's working \p cells back to a flat \c Schedule,
/// ready for \c peak_profile_sweep.
///
/// \details Per cell: `Cell{value_id, home_depth, cell_footprint(carried,
/// home_modes, cm, block_of), first_use, last_use}` -- the same projection
/// \c compute_dag_path applies to a \c RichSchedule, just over the remat pass's
/// (possibly SHRINK-mutated) \c home_modes instead of the as-seeded ones.
///
/// \p block_of is any `Index -> std::size_t` callable giving the block
/// (sliced) element count for a mode; forwarded to \c detail::cell_footprint.
///
inline Schedule to_schedule(container::svector<ValueCell> const& cells,
                            dryrun::CostModel const& cm, auto const& block_of,
                            std::size_t num_points) {
  Schedule out;
  out.num_points = num_points;
  out.cells.reserve(cells.size());
  for (auto const& c : cells) {
    Cell fc;
    fc.value_id = c.value_id;
    fc.home_depth = c.home_depth;
    fc.footprint = detail::cell_footprint(c.carried, c.home_modes, cm, block_of,
                                          c.divergent_modes);
    fc.first_use = c.first_use;
    fc.last_use = c.last_use;
    out.cells.push_back(fc);
  }
  return out;
}

///
/// \brief The demoted carried batch modes of \p c: modes that EVER enclose
/// \p c (\c enclosing_modes), are among \p c's own carried indices (\c
/// carried), but are NOT (yet) part of its home (\c home_modes).
///
/// \details `(enclosing_modes INTERSECT carried) MINUS home_modes`. Each
/// returned mode is a SHRINK candidate: adding it to \c c.home_modes (T2)
/// block-sizes it in \c c's footprint instead of leaving it FULL.
///
inline container::svector<Index> shrink_candidates(ValueCell const& c) {
  container::svector<Index> out;
  for (auto const& m : c.enclosing_modes) {
    bool const in_carried =
        std::find(c.carried.begin(), c.carried.end(), m) != c.carried.end();
    bool const in_home = std::find(c.home_modes.begin(), c.home_modes.end(),
                                   m) != c.home_modes.end();
    if (in_carried && !in_home) out.push_back(m);
  }
  return out;
}

///
/// \brief Apply a SHRINK move to \p c: add candidate mode \p m to its home,
/// so \c to_schedule then BLOCK-sizes \p m instead of leaving it FULL.
///
/// \details The demoted-giant lever (design section 3a): a cell homed ABOVE a
/// batch loop it carries holds that mode full; adding the mode to \c home_modes
/// re-homes the cell INSIDE the loop, dropping its footprint by the block
/// factor at ZERO recompute (result-mode slicing is free per the \c W model).
/// The liveness interval is unchanged. \p m MUST be a current shrink candidate
/// of \p c (asserted).
///
inline void apply_shrink(ValueCell& c, Index const& m) {
  auto const cands = shrink_candidates(c);
  SEQUANT_ASSERT(std::find(cands.begin(), cands.end(), m) != cands.end());
  c.home_modes.push_back(m);
}

///
/// \brief Why \c rematerialize_to_budget stopped.
///
/// \details \c Feasible: `peak <= threshold`. The two infeasibilities (design
/// section 4) are NOT repairable within the remat pass and feed Phase 5: \c
/// RebatchNeeded (every cell live at the binding point has no shrink candidate
/// -- a giant needing a batch loop the per-term factorizer never added) and \c
/// FactorizationInherent (a live cell had a candidate, but no single shrink
/// reduced the binding peak -- the factorization itself must change).
///
enum class RematStatus { Feasible, FactorizationInherent, RebatchNeeded };

///
/// \brief The result of \c rematerialize_to_budget: the final (possibly shrunk)
/// cells, the resulting \c PeakProfile, and why it stopped.
///
struct RematResult {
  container::svector<ValueCell> cells;
  PeakProfile profile;
  RematStatus status;
};

namespace detail {

/// Classify an infeasible \c rematerialize_to_budget stop (design section 4):
/// \c RebatchNeeded if EVERY cell live at the binding point has an EMPTY \c
/// shrink_candidates; else \c FactorizationInherent (some live cell had a
/// candidate, but no shrink lowered the binding peak).
inline RematStatus classify_infeasible(
    container::svector<ValueCell> const& cells, PeakProfile const& pp) {
  for (auto ci : pp.live_at_binding)
    if (!shrink_candidates(cells[ci]).empty())
      return RematStatus::FactorizationInherent;
  return RematStatus::RebatchNeeded;
}

}  // namespace detail

///
/// \brief the remat pass's greedy loop: shrink demoted giants until the modeled
/// peak fits \p threshold, or report an infeasibility (design sections 2-4).
///
/// \details At each step it sweeps the current placement (\c to_schedule +
/// \c peak_profile_sweep); if `peak <= threshold` it is \c Feasible. Otherwise
/// it considers every SHRINK of every cell live at the binding peak point and
/// applies the one with the largest \c DeltaPeak (all shrinks are
/// zero-recompute this phase, so the metric is pure \c DeltaPeak); if no shrink
/// reduces the peak it stops infeasible. Termination: each accepted shrink
/// strictly appends a mode to some cell's \c home_modes from that cell's finite
/// (and strictly shrinking) \c shrink_candidates set, so the move supply is
/// finite. Un-hoist / split (the evict moves) are DEFERRED (design section 3a).
///
/// \p threshold a byte budget; `+infinity` makes remat a no-op (the raw
/// perfect-CSE seed).
///
template <typename BlockOfFn>
RematResult rematerialize_to_budget(container::svector<ValueCell> cells,
                                    dryrun::CostModel const& cm,
                                    BlockOfFn const& block_of,
                                    std::size_t num_points, double threshold) {
  for (;;) {
    Schedule const s = to_schedule(cells, cm, block_of, num_points);
    PeakProfile pp = peak_profile_sweep(s);
    if (pp.peak_bytes <= threshold)
      return {std::move(cells), std::move(pp), RematStatus::Feasible};

    // Among the cells live at the binding peak point, pick the shrink with the
    // largest DeltaPeak (re-sweeping a trial placement per candidate).
    double best_drop = 0;
    std::optional<std::pair<std::size_t, Index>> best;
    for (auto ci : pp.live_at_binding) {
      for (auto const& m : shrink_candidates(cells[ci])) {
        auto trial = cells;
        apply_shrink(trial[ci], m);
        double const p2 =
            peak_profile_sweep(to_schedule(trial, cm, block_of, num_points))
                .peak_bytes;
        double const drop = pp.peak_bytes - p2;
        if (drop > best_drop) {
          best_drop = drop;
          best = std::pair<std::size_t, Index>{ci, m};
        }
      }
    }

    if (!best) {  // no placement move reduces the binding point: infeasible
      RematStatus const status = detail::classify_infeasible(cells, pp);
      return {std::move(cells), std::move(pp), status};
    }

    apply_shrink(cells[best->first], best->second);  // full re-sweep next iter
  }
}

///
/// \brief Emit \c PlacementRouter overrides for the MOVED cells of a remat
/// result (Phase 4b-2). Additive: produces a populated router; nothing is
/// wired into runtime (that is Phase 4b-3).
///
/// \details A remat cell is one-per-VALUE (its \c ValueCell::hash, the CSE
/// identity). The router is keyed by the batched-slot-AWARE OCCURRENCE KEY
/// (\c occurrence_key). So a single MOVED value emits one override per DISTINCT
/// occurrence key of that value, all pointing to the SAME \c HomeTarget (its
/// meet home). Only MOVED cells (final \c home_modes != seed \c home_modes) are
/// emitted: an unmoved cell's seed home equals the runtime's default
/// derivation, so an override would be redundant (Phase-2 empty-router-is-inert
/// seam). The value HASH links a \c RematResult cell back to its forest nodes;
/// the (expensive, bliss) occurrence keys are computed LAZILY here, for
/// moved-cell occurrences only -- \c compute_dag_boulevard / \c peak_profile
/// pay no per-node bliss cost.
///
/// \param seed_cells the seed placement (e.g. \c remat_cells(forest).cells);
///        \c home_modes is each value's seed \c home_scope.
/// \param final_cells the post-spill placement (\c
///        rematerialize_to_budget(...).cells); \c home_modes is the final home.
/// \param forest the eval forest \p seed_cells / \p final_cells were built
/// from.
/// \return a \c PlacementRouter with one override per (moved value, occurrence
///         key). Empty iff nothing moved.
///
template <meta::eval_node_range R>
[[nodiscard]] PlacementRouter<std::ranges::range_value_t<R>> remat_to_router(
    container::svector<ValueCell> const& seed_cells,
    container::svector<ValueCell> const& final_cells, R const& forest) {
  using Node = std::ranges::range_value_t<R>;

  // Seed home_modes by value_id (robust to any reordering between seed and
  // result), then the MOVED map: value hash -> final home_modes.
  std::unordered_map<std::size_t, container::svector<Index> const*> seed_home;
  for (auto const& c : seed_cells) seed_home.emplace(c.value_id, &c.home_modes);
  // MOVED map: value hash -> the DAG-global HomeTarget applied at EVERY
  // occurrence key of that value. Its DAG-scope is computed ONCE, label-free,
  // from the CELL's `home_modes`: each home mode contributes its SPACE, ordered
  // by the mode's position in the enclosing NEST (`enclosing_modes`, folded
  // outer->inner across every occurrence). It is a nest prefix, not a
  // value-relative coordinate: a home mode need not appear on another
  // occurrence's own slots (the g.C legs' i_3 vs i_4), and a home mode the
  // value is invariant to (homed WITHIN but does not carry) is just another
  // DAG-scope space -- no free-mode special case. home_depth resolves each
  // space to THIS occurrence's physical index of that space per use.
  std::unordered_map<std::size_t, HomeTarget> moved;
  for (auto const& c : final_cells) {
    auto const sh = seed_home.find(c.value_id);
    if (sh == seed_home.end() || c.home_modes == *sh->second) continue;
    // Order the home modes by nest depth (their index in enclosing_modes;
    // modes absent from it sort last, keeping their home_modes order).
    auto const nest_pos = [&](Index const& m) -> std::size_t {
      auto const it =
          std::find(c.enclosing_modes.begin(), c.enclosing_modes.end(), m);
      return it == c.enclosing_modes.end()
                 ? c.enclosing_modes.size()
                 : static_cast<std::size_t>(it - c.enclosing_modes.begin());
    };
    container::svector<Index> ordered(c.home_modes.begin(), c.home_modes.end());
    std::stable_sort(ordered.begin(), ordered.end(),
                     [&](Index const& a, Index const& b) {
                       return nest_pos(a) < nest_pos(b);
                     });
    HomeTarget ht;
    for (auto const& m : ordered) ht.dag_scope.push_back(m.space());
    moved.emplace(c.hash, std::move(ht));
  }

  PlacementRouter<Node> router;
  if (moved.empty()) return router;  // nothing moved => inert seed seam
  // Context-invariant "this value is demoted" set, so an OUTER-scope hoist
  // (whose occurrence-key query would miss) still learns not to build a moved
  // value full at the root (see place_at_this_level in eval.hpp).
  for (auto const& [h, ht] : moved) router.mark_moved(h);

  // Re-walk the forest (same ectx accumulation as compute_dag_boulevard). For
  // each node whose value is moved, emit the precomputed overlay keyed by its
  // occurrence key at the ambient (enclosing) batch context.
  auto visit = [&](auto&& self, Node const& n,
                   container::svector<Index> ctx_modes) -> void {
    if (auto const it = moved.find(n->hash_value()); it != moved.end()) {
      auto key = occurrence_key(n, ctx_modes);
      router.set_override(std::move(key), it->second);
    }
    if (n.leaf()) return;
    // Children see this node's OWN realized loops on top of ctx_modes; the node
    // itself is keyed with the enclosing ctx_modes (excludes its own loops).
    container::svector<Index> child_ctx = ctx_modes;
    for (auto const& [ix, kind] : n->batched_here()) {
      container::svector<Index> expanded;
      sequant::detail::proto_expand_into(expanded, ix);
      for (auto const& m : expanded) child_ctx.push_back(m);
    }
    self(self, n.left(), child_ctx);
    self(self, n.right(), child_ctx);
  };
  for (auto const& tree : forest)
    visit(visit, tree, container::svector<Index>{});
  return router;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_PLACEMENT_REMAT_HPP
