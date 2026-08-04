#ifndef SEQUANT_EVAL_PLACEMENT_O2_HPP
#define SEQUANT_EVAL_PLACEMENT_O2_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/index.hpp>

#include <algorithm>
#include <cstddef>

namespace sequant::eval {

///
/// \brief O2, a whole-forest greedy spill pass (design doc
/// `doc/dev/specs/2026-08-04-phase4-threshold-o2-cutover-design.md`), is
/// STANDALONE in Phase 4a -- zero production callers. This header lands its
/// working representation: a mutable per-value cell (\c O2Cell), the flat
/// input/output it is built from and projected back to (\c O2Input, \c
/// to_schedule), and the candidate-selection query a SHRINK move needs (\c
/// shrink_candidates).
///

///
/// \brief O2's mutable per-value working record.
///
/// \details Identical in shape to \c ValueCell (Phase 3b/4a's rich
/// per-value fold from \c linearize_rich) -- reused directly rather than
/// duplicated, since O2 mutates exactly the fields \c ValueCell already
/// carries. \c home_modes is the field a later SHRINK move (Phase 4a T2)
/// mutates in place (adding a demoted carried batch mode, block-sizing it);
/// the other fields stay fixed once seeded.
///
using O2Cell = ValueCell;

///
/// \brief The flat input to O2: every \c O2Cell of a linearized forest, plus
/// the timeline length needed to re-sweep it.
///
struct O2Input {
  container::svector<O2Cell> cells;
  std::size_t num_points = 0;
};

///
/// \brief Seed O2's working cells from an eval \p forest: \c linearize_rich
/// mapped 1:1 to \c O2Cell (which IS \c ValueCell), plus \c num_points.
///
/// \p block_of is any `Index -> std::size_t` callable giving the block
/// (sliced) element count for a mode; forwarded to \c linearize_rich.
///
template <meta::eval_node_range R>
O2Input o2_cells(R const& forest, dryrun::CostModel const& cm,
                 auto const& block_of) {
  RichSchedule rich = linearize_rich(forest, cm, block_of);
  O2Input out;
  out.cells.assign(std::make_move_iterator(rich.cells.begin()),
                   std::make_move_iterator(rich.cells.end()));
  out.num_points = rich.num_points;
  return out;
}

///
/// \brief Project O2's working \p cells back to a flat \c Schedule, ready
/// for \c peak_profile_sweep.
///
/// \details Per cell: `Cell{value_id, home_depth, cell_footprint(carried,
/// home_modes, cm, block_of), first_use, last_use}` -- the same projection
/// \c linearize applies to a \c RichSchedule, just over O2's (possibly
/// SHRINK-mutated) \c home_modes instead of the as-seeded ones.
///
/// \p block_of is any `Index -> std::size_t` callable giving the block
/// (sliced) element count for a mode; forwarded to \c detail::cell_footprint.
///
inline Schedule to_schedule(container::svector<O2Cell> const& cells,
                            dryrun::CostModel const& cm, auto const& block_of,
                            std::size_t num_points) {
  Schedule out;
  out.num_points = num_points;
  out.cells.reserve(cells.size());
  for (auto const& c : cells) {
    Cell fc;
    fc.value_id = c.value_id;
    fc.home_depth = c.home_depth;
    fc.footprint =
        detail::cell_footprint(c.carried, c.home_modes, cm, block_of);
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
inline container::svector<Index> shrink_candidates(O2Cell const& c) {
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

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_PLACEMENT_O2_HPP
