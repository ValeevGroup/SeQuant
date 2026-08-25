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
/// carried), are NOT (yet) part of its home (\c home_modes), and are
/// SHARED-LABEL (not \c divergent_modes).
///
/// \details `((enclosing_modes INTERSECT carried) MINUS home_modes) MINUS
/// divergent_modes`. Each returned mode is an IN-PLACE SHRINK candidate: adding
/// it to \c c.home_modes (via \c apply_shrink) block-sizes it in \c c's
/// footprint instead of leaving it FULL, at ZERO recompute. A RELABELED
/// (divergent) mode is EXCLUDED here -- slicing it cannot be shared across
/// occurrences, so it is offered only as a SPLIT (see \c split_candidates /
/// \c apply_split). For a cell with no divergent modes (every seed cell of a
/// forest whose CSE folds no relabeled value) this is byte-identical to the
/// pre-CSE-aware behavior.
///
inline container::svector<Index> shrink_candidates(ValueCell const& c) {
  container::svector<Index> out;
  for (auto const& m : c.enclosing_modes) {
    bool const in_carried =
        std::find(c.carried.begin(), c.carried.end(), m) != c.carried.end();
    bool const in_home = std::find(c.home_modes.begin(), c.home_modes.end(),
                                   m) != c.home_modes.end();
    bool const divergent =
        std::find(c.divergent_modes.begin(), c.divergent_modes.end(), m) !=
        c.divergent_modes.end();
    if (in_carried && !in_home && !divergent) out.push_back(m);
  }
  return out;
}

///
/// \brief The RELABELED (divergent) batch modes of \p c that a SPLIT can slice:
/// modes in \c divergent_modes that EVER enclose \p c (\c enclosing_modes).
///
/// \details These are the modes \c shrink_candidates deliberately excludes: a
/// mode whose PHYSICAL label differs across \p c's occurrences (the g.C legs'
/// \c i_3 vs \c i_4) cannot be sliced while the value stays CSE-folded -- one
/// sliced copy would feed the wrong slice to the relabeling occurrence (design
/// spec, "The verified correctness hole"). It is offered only as a SPLIT: \c
/// apply_split un-folds \p c into one cell per physical binding, each of which
/// slices the mode under its OWN label (\c apply_split). Membership in \c
/// enclosing_modes mirrors \c shrink_candidates' gate: only a relabeled mode
/// that actually appears as an enclosing loop somewhere is sliceable at all.
///
inline container::svector<Index> split_candidates(ValueCell const& c) {
  container::svector<Index> out;
  for (auto const& m : c.divergent_modes) {
    bool const in_enclosing =
        std::find(c.enclosing_modes.begin(), c.enclosing_modes.end(), m) !=
        c.enclosing_modes.end();
    if (in_enclosing) out.push_back(m);
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

namespace detail {

/// The REPLICATION FACTOR of a split cell \p c: the product, over every
/// enclosing loop \p c is homed-WITHIN but does NOT carry (\c enclosing_modes
/// MINUS \c carried), of that loop's block COUNT (number of iterations). A cell
/// homed at a nest prefix is rebuilt once per iteration of each such loop (it
/// is INVARIANT to that loop -- its result has no slot for it -- so it cannot
/// be sliced there and is re-materialized full every block). This is the
/// DAG-scope-only term the value-relative \c home_modes coordinate cannot name
/// (leg B's \c i_3 is absent from B's \c carried); see the design spec, Design
/// section 3. The canonical case: leg A carries its enclosing loop -> factor 1;
/// leg B is enclosed by an extra loop it does not carry -> factor = that loop's
/// block count N3, so the two split cells recompute `(1 + N3) x` a single
/// build, equal to a flat 2x only at N3 = 1.
///
/// The block COUNT is `ceil(full extent / block size)`: \p block_of returns the
/// block SIZE (sliced element count, as \c cell_footprint consumes it), so the
/// iteration count is the full regime extent (\p cm) divided by it -- NOT
/// \c block_of itself.
template <typename BlockOfFn>
[[nodiscard]] inline std::size_t replication_factor(ValueCell const& c,
                                                    dryrun::CostModel const& cm,
                                                    BlockOfFn const& block_of) {
  std::size_t factor = 1;
  for (auto const& m : c.enclosing_modes)
    if (std::find(c.carried.begin(), c.carried.end(), m) == c.carried.end()) {
      std::size_t const bsz = block_of(m);
      std::size_t const ext = cm.regime().extent(m);
      factor *= (bsz == 0) ? 1 : (ext + bsz - 1) / bsz;  // ceil-div block count
    }
  return factor;
}

}  // namespace detail

///
/// \brief Apply a SPLIT move: un-fold the divergent cell at index \p ci of
/// \p cells along relabeled mode \p m into TWO-OR-MORE same-hash cells, one per
/// physical binding of \p m across the cell's occurrences, and RETURN the
/// forecast replication recompute of the split (report-only; see below).
///
/// \details The CSE-aware counterpart of \c apply_shrink. \p m must be a
/// current \c split_candidate of the cell (a relabeled/divergent enclosing
/// mode; asserted). The cell's occurrences are partitioned by the PHYSICAL
/// index each binds at \p m's canonical result slot (leg A -> \c i_3, leg B ->
/// \c i_4); each partition becomes a NON-divergent \c ValueCell of the SAME
/// hash with:
///   - subset-local \c carried (the partition's occurrences agree on it),
///   - \c home_modes = the shared home PLUS this subset's physical binding of
///     \p m (the relabeled mode folded into the home under its own label),
///   - subset-local \c enclosing_modes and liveness (\c first_use / \c
///     last_use) folded over just this subset's occurrences,
///   - \c divergent_modes = {} (each split cell is single-binding),
///   - a FRESH \c value_id (so \c peak_profile_sweep prices the two cells'
///     co-residency by their liveness overlap, keyed on \c value_id -- two
///     cells of one hash need no sweep change).
/// The one original cell is REPLACED by the partition cells (the container
/// grows). Peak is left to the sweep; this returns the SEPARATE, report-only
/// RECOMPUTE forecast: the sum over the split cells of
/// `cell_footprint(cell) x replication_factor(cell)` (see \c
/// detail::replication_factor). remat's selection stays PEAK-driven -- the
/// recompute is recorded for the dry-run cost forecast, NOT added to the
/// objective (design spec, "Does this recompute term change the schedule?").
///
template <typename BlockOfFn>
inline std::size_t apply_split(container::svector<ValueCell>& cells,
                               std::size_t ci, Index const& m,
                               dryrun::CostModel const& cm,
                               BlockOfFn const& block_of) {
  SEQUANT_ASSERT(ci < cells.size());
  ValueCell const orig = cells[ci];
  {
    auto const cands = split_candidates(orig);
    SEQUANT_ASSERT(std::find(cands.begin(), cands.end(), m) != cands.end());
  }

  // The canonical result slot of the relabeled mode: its position in an
  // occurrence that carries it. Canonically-equal occurrences (one hash) place
  // it at the SAME position, so every occurrence's binding is `carried[slot]`.
  std::size_t slot = orig.carried.size();  // sentinel "no such slot"
  for (auto const& occ : orig.occurrences) {
    auto const it = std::find(occ.carried.begin(), occ.carried.end(), m);
    if (it != occ.carried.end()) {
      slot = static_cast<std::size_t>(it - occ.carried.begin());
      break;
    }
  }

  // Partition occurrences by their physical binding at `slot` (order-preserving
  // so the first-seen binding leads). Occurrences whose result lacks that slot
  // (the mode is genuinely absent, not merely relabeled) collect in a separate
  // absent-binding bucket.
  container::svector<Index> bindings;  // distinct binding labels, in order
  container::svector<container::svector<std::size_t>> groups;  // occ indices
  container::svector<std::size_t> absent_group;
  for (std::size_t oi = 0; oi < orig.occurrences.size(); ++oi) {
    auto const& oc = orig.occurrences[oi].carried;
    if (slot >= oc.size()) {
      absent_group.push_back(oi);
      continue;
    }
    Index const& bind = oc[slot];
    auto const it = std::find(bindings.begin(), bindings.end(), bind);
    std::size_t const b = static_cast<std::size_t>(it - bindings.begin());
    if (it == bindings.end()) {
      bindings.push_back(bind);
      groups.emplace_back();
    }
    groups[b].push_back(oi);
  }

  // Fresh value ids: max existing + 1, bumped per new cell.
  std::size_t next_id = 0;
  for (auto const& c : cells) next_id = std::max(next_id, c.value_id + 1);

  // Build one split cell per binding group (and one for the absent bucket).
  auto build_cell = [&](container::svector<std::size_t> const& group,
                        Index const* binding) -> ValueCell {
    ValueCell nc;
    nc.hash = orig.hash;
    nc.value_id = next_id++;
    // subset-local carried: occurrences of one group agree; read off the first.
    nc.carried = orig.occurrences[group.front()].carried;
    // home = shared home + this subset's physical binding of the relabeled
    // mode.
    nc.home_modes = orig.home_modes;
    if (binding && std::find(nc.home_modes.begin(), nc.home_modes.end(),
                             *binding) == nc.home_modes.end())
      nc.home_modes.push_back(*binding);
    // subset-local enclosing_modes + liveness.
    std::size_t fu = orig.occurrences[group.front()].point;
    std::size_t lu = orig.occurrences[group.front()].consumer_point;
    for (auto oi : group) {
      auto const& occ = orig.occurrences[oi];
      for (auto const& e : occ.ectx)
        if (std::find(nc.enclosing_modes.begin(), nc.enclosing_modes.end(),
                      e.first) == nc.enclosing_modes.end())
          nc.enclosing_modes.push_back(e.first);
      fu = std::min(fu, occ.point);
      lu = std::max(lu, occ.consumer_point);
      nc.occurrences.push_back(occ);
    }
    nc.first_use = fu;
    nc.last_use = lu;
    nc.home_depth = detail::home_depth_of(nc.home_modes,
                                          orig.occurrences[group.front()].ectx);
    // divergent_modes stays empty: each split cell is single-binding.
    return nc;
  };

  container::svector<ValueCell> split_cells;
  for (std::size_t b = 0; b < groups.size(); ++b)
    split_cells.push_back(build_cell(groups[b], &bindings[b]));
  if (!absent_group.empty())
    split_cells.push_back(build_cell(absent_group, nullptr));

  // Report-only replication recompute over the split cells.
  std::size_t recompute = 0;
  for (auto const& nc : split_cells)
    recompute +=
        detail::cell_footprint(nc.carried, nc.home_modes, cm, block_of) *
        detail::replication_factor(nc, cm, block_of);

  // Replace cells[ci] with the split cells (grow the container).
  container::svector<ValueCell> rebuilt;
  rebuilt.reserve(cells.size() + split_cells.size());
  for (std::size_t i = 0; i < cells.size(); ++i) {
    if (i == ci)
      for (auto& nc : split_cells) rebuilt.push_back(std::move(nc));
    else
      rebuilt.push_back(std::move(cells[i]));
  }
  cells = std::move(rebuilt);
  return recompute;
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
  //!< Report-only FORECAST of the replication recompute incurred by the SPLIT
  //!< moves this pass committed (sum of each split cell's
  //!< `cell_footprint x replication_factor`; see \c apply_split). ZERO when no
  //!< divergent value was split (every existing shrink-only case). It does NOT
  //!< enter the selection objective -- that stays pure \c DeltaPeak (design
  //!< spec, "Does this recompute term change the schedule?").
  std::size_t modeled_recompute = 0;
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
/// it considers, over every cell live at the binding peak point, both every
/// in-place SHRINK (a shared-label demoted mode, zero recompute) and every
/// SPLIT (un-fold a divergent value along a relabeled mode; see \c apply_split)
/// and applies the one with the largest \c DeltaPeak. The selection objective
/// is PEAK-ONLY: a split's replication recompute is RECORDED (accumulated into
/// \c RematResult::modeled_recompute) for the dry-run cost forecast but does
/// NOT enter the objective -- making recompute a secondary objective is a
/// separate, deliberate decision (design spec, "Does this recompute term change
/// the schedule?"). If no move reduces the peak it stops infeasible.
///
/// Termination: a SHRINK strictly appends a mode to some cell's \c home_modes
/// from that cell's finite (and strictly shrinking) \c shrink_candidates set; a
/// SPLIT strictly reduces the forest's total divergent-mode count (each split
/// cell carries an EMPTY \c divergent_modes, so its members can never be split
/// again), and the divergent-mode count is finite. Both move supplies are
/// therefore finite. Un-hoist (the other evict move) is DEFERRED (design
/// section 3a).
///
/// \p threshold a byte budget; `+infinity` makes remat a no-op (the raw
/// perfect-CSE seed).
///
template <typename BlockOfFn>
RematResult rematerialize_to_budget(container::svector<ValueCell> cells,
                                    dryrun::CostModel const& cm,
                                    BlockOfFn const& block_of,
                                    std::size_t num_points, double threshold) {
  std::size_t total_recompute = 0;
  for (;;) {
    Schedule const s = to_schedule(cells, cm, block_of, num_points);
    PeakProfile pp = peak_profile_sweep(s);
    if (pp.peak_bytes <= threshold)
      return {std::move(cells), std::move(pp), RematStatus::Feasible,
              total_recompute};

    // Among the cells live at the binding peak point, pick the move (shrink or
    // split) with the largest DeltaPeak (re-sweeping a trial placement per
    // candidate). A split grows the container, so a split candidate is trialed
    // by copying `cells` and applying the split to the copy.
    enum class Kind { Shrink, Split };
    double best_drop = 0;
    std::optional<std::tuple<std::size_t, Index, Kind>> best;
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
          best = std::tuple{ci, m, Kind::Shrink};
        }
      }
      for (auto const& m : split_candidates(cells[ci])) {
        auto trial = cells;
        (void)apply_split(trial, ci, m, cm, block_of);
        double const p2 =
            peak_profile_sweep(to_schedule(trial, cm, block_of, num_points))
                .peak_bytes;
        double const drop = pp.peak_bytes - p2;
        if (drop > best_drop) {
          best_drop = drop;
          best = std::tuple{ci, m, Kind::Split};
        }
      }
    }

    if (!best) {  // no placement move reduces the binding point: infeasible
      RematStatus const status = detail::classify_infeasible(cells, pp);
      return {std::move(cells), std::move(pp), status, total_recompute};
    }

    if (std::get<2>(*best) == Kind::Shrink)
      apply_shrink(cells[std::get<0>(*best)], std::get<1>(*best));
    else
      total_recompute +=
          apply_split(cells, std::get<0>(*best), std::get<1>(*best), cm,
                      block_of);  // full re-sweep next iter
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

  // Moved-detection is HASH-keyed, not value_id-keyed: a SPLIT (see
  // apply_split) replaces one seed cell with two FINAL cells of the same hash
  // but FRESH value_ids, which a value_id-keyed seed map would miss and drop as
  // "unmoved". Keyed by hash: a value is moved iff any of its final cells has a
  // home_modes differing from the seed's, OR the split changed how many cells
  // the hash has (one seed -> two final). Seed cells are one-per-hash (grouped
  // by compute_dag_boulevard), so `seed_home` is well-defined; a value never
  // split at seed has exactly one final cell whose home_modes an unmoved value
  // leaves equal to seed.
  std::unordered_map<std::size_t, container::svector<Index> const*> seed_home;
  std::unordered_map<std::size_t, std::size_t> seed_count, final_count;
  for (auto const& c : seed_cells) {
    seed_home.emplace(c.hash, &c.home_modes);
    ++seed_count[c.hash];
  }
  for (auto const& c : final_cells) ++final_count[c.hash];
  // MOVED map: value hash -> the DAG-global HomeTarget applied at EVERY
  // occurrence key of that value. Its DAG-scope is computed ONCE, label-free,
  // from the CELL's `home_modes`: each home mode contributes its SPACE, ordered
  // by the mode's position in the enclosing NEST (`enclosing_modes`, folded
  // outer->inner across every occurrence). It is a nest prefix, not a
  // value-relative coordinate: a home mode need not appear on another
  // occurrence's own slots (the g.C legs' i_3 vs i_4), and a home mode the
  // value is invariant to (homed WITHIN but does not carry) is just another
  // DAG-scope space -- no free-mode special case. home_depth resolves each
  // space to THIS occurrence's physical index of that space per use. Both split
  // cells of one value are homed on the SAME space (e.g. both on `occ`, under
  // labels i_3 vs i_4), so they yield the SAME dag_scope: we emit ONE overlay
  // per hash and ASSERT the two agree (rather than emplace-first-wins).
  std::unordered_map<std::size_t, HomeTarget> moved;
  for (auto const& c : final_cells) {
    auto const sh = seed_home.find(c.hash);
    bool const home_changed =
        sh == seed_home.end() || c.home_modes != *sh->second;
    bool const count_changed = final_count[c.hash] != seed_count[c.hash];
    if (!home_changed && !count_changed) continue;
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
    container::svector<IndexSpace> scope;
    for (auto const& m : ordered) scope.push_back(m.space());
    auto const it = moved.find(c.hash);
    if (it == moved.end())
      moved.emplace(c.hash, HomeTarget{std::move(scope)});
    else  // another split cell of this hash: its DAG-scope must AGREE.
      SEQUANT_ASSERT(it->second.dag_scope == scope);
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
    for (auto const& [ix, kind] : n->node_slice_mask()) {
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
