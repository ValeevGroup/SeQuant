#ifndef SEQUANT_EVAL_PEAK_PROFILE_HPP
#define SEQUANT_EVAL_PEAK_PROFILE_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/lifetime_mask.hpp>
#include <SeQuant/core/index.hpp>

#include <algorithm>
#include <cstddef>
#include <unordered_map>
#include <utility>

namespace sequant::eval {

///
/// \brief NEW, purely-additive static peak-profile analysis (spec section 9,
/// `doc/dev/specs/2026-08-03-meet-based-home-scope-phase3-design.md`) that
/// consumes the Phase-3a `home_scope` seed. Phase 3b T1 landed the two sizing
/// primitives (`detail::home_depth_of`, `detail::cell_footprint`); T2 adds the
/// forest linearization (`compute_dag_path`) and the interval-event sweep
/// (`peak_profile_sweep`) that turns an eval forest into a `PeakProfile`. Still
/// ZERO production callers -- runtime untouched; O2 (a later task) consumes the
/// `PeakProfile`.
///
namespace detail {

///
/// \brief The enclosing-batch-context type consulted below.
///
/// \details An ordered stack (outermost-first) of the enclosing realized
/// batch loops, one entry per loop, `{axis mode, {block_lo, block_hi}}`
/// (element range). This is verbatim the same underlying type as \c
/// CacheManager<TreeNode, force_hash_collisions>::BatchContext (see \c
/// cache_manager.hpp:95) -- that alias does not itself depend on \c
/// TreeNode, so this is not a look-alike but literally the same type, and
/// the real runtime \c BatchContext can be passed to \c home_depth_of /
/// \c cell_footprint as-is once a later task wires them in.
///
using BatchContext =
    container::svector<std::pair<Index, std::pair<std::size_t, std::size_t>>>;

///
/// \brief Resolve a residency mode-set to an enclosing-batch-context loop
/// depth.
///
/// \details Returns the deepest (innermost) level \p i of \p ectx whose loop
/// mode `ectx[i].first` is a member of \p home_modes; -1 if no level matches
/// (the cell is invariant to the whole nest -- the chain root). Mirrors the
/// runtime rl-walk at `eval.hpp:1776-1782` verbatim (same innermost-to-
/// outermost scan, same membership test), just against a caller-supplied
/// mode-set instead of the sliced/contracted-modes union computed there.
///
[[nodiscard]] inline int home_depth_of(
    container::svector<Index> const& home_modes,
    BatchContext const& ectx) noexcept {
  for (int i = static_cast<int>(ectx.size()) - 1; i >= 0; --i)
    if (std::find(home_modes.begin(), home_modes.end(), ectx[i].first) !=
        home_modes.end())
      return i;
  return -1;
}

///
/// \brief Home-relative footprint of a cell via the existing \c
/// dryrun::CostModel::memsize.
///
/// \details A cell homed at \p home_modes (depth `d = home_depth_of(
/// home_modes, ectx)`) sizes each CARRIED mode `m` at BLOCK extent (via \p
/// block_of) if `m`'s loop ENCLOSES the home -- i.e. `m` appears in \p ectx
/// at a level `<= d` -- else at FULL (nominal regime) extent. Builds the \c
/// dryrun::ExtentOverrides mapping exactly those block-sized modes to their
/// block extent, then delegates the actual extent-product / composite-moment
/// math to \p cm.memsize() verbatim: no sizing logic is reimplemented here.
///
/// \p block_of is any `Index -> std::size_t` callable giving the block
/// (sliced) element count for a mode.
///
/// \p divergent_modes are the value's RELABELED modes (carried by some
/// occurrences but not all -- see \c ValueCell::divergent_modes). Slicing one
/// of them cannot be shared across occurrences (they bind different physical
/// labels to that canonical slot), so the runtime SPLITS the value into
/// per-occurrence copies (see \c place_at_this_level). When the home slices a
/// divergent mode, the co-resident footprint is therefore that of TWO sliced
/// copies, not one -- priced here as 2x. (A pairwise model: the dominant
/// divergent-CSE case is two legs of one contraction, e.g. the g.C legs.)
template <typename BlockOfFn>
[[nodiscard]] inline std::size_t cell_footprint(
    container::svector<Index> const& carried,
    container::svector<Index> const& home_modes, dryrun::CostModel const& cm,
    BlockOfFn const& block_of,
    container::svector<Index> const& divergent_modes = {}) {
  dryrun::ExtentOverrides ov;
  // block iff in the meet-home. Overrides are POSITIONAL against `carried`:
  // map each home mode to its position there. memsize() re-expands the position
  // to that Index and applies the block extent wherever the Index recurs
  // (including as a composite's outer proto), so composite slicing is
  // preserved.
  for (auto const& m : home_modes) {
    auto const it = std::find(carried.begin(), carried.end(), m);
    if (it != carried.end())
      ov[static_cast<std::size_t>(it - carried.begin())] = block_of(m);
  }
  std::size_t base = cm.memsize(carried, ov);  // non-meet carried modes FULL
  bool split = false;
  for (auto const& m : home_modes)
    if (std::find(divergent_modes.begin(), divergent_modes.end(), m) !=
        divergent_modes.end()) {
      split = true;
      break;
    }
  return split ? 2 * base : base;  // split -> two co-resident sliced copies
}

}  // namespace detail

///
/// \brief One value group in a linearized schedule: a single logical value
/// (all its perfect-CSE occurrences folded together) with its home-relative
/// footprint and its inclusive static-point liveness range.
///
struct Cell {
  std::size_t value_id;   //!< stable index of the value group (== its slot in
                          //!< \c Schedule::cells)
  int home_depth;         //!< \c home_depth_of(home_scope, ectx) at any
                          //!< occurrence (-1 == above the whole nest)
  std::size_t footprint;  //!< home-relative size in BYTES
  std::size_t first_use;  //!< earliest static point the value is live at
  std::size_t last_use;   //!< latest static point the value is live at
                          //!< (its last consumer), inclusive
};

///
/// \brief A whole forest linearized to a flat list of value cells over a
/// single monotone static-point timeline.
///
struct Schedule {
  container::svector<Cell> cells;
  std::size_t num_points = 0;  //!< one past the last static point
};

///
/// \brief The result of the interval-event sweep over a \c Schedule: the peak
/// live footprint, the (lowest) point achieving it, and which cells are live
/// there.
///
struct PeakProfile {
  double peak_bytes = 0;
  std::size_t binding_point = 0;
  container::svector<std::size_t>
      live_at_binding;  //!< indices into \c Schedule::cells live at \c
                        //!< binding_point
};

///
/// \brief Sweep the per-cell liveness intervals of \p s to find the peak live
/// footprint, the lowest static point achieving it, and the set of cells live
/// there.
///
/// \details A textbook +delta/-delta interval-event sweep: each cell deposits
/// +footprint at \c first_use and -footprint just past \c last_use, then a
/// single left-to-right prefix scan tracks the running live total. The strict
/// `>` comparison keeps the FIRST (lowest) point among equal-height maxima.
///
inline PeakProfile peak_profile_sweep(Schedule const& s) {
  container::svector<double> delta(s.num_points + 1, 0.0);
  for (auto const& c : s.cells) {
    delta[c.first_use] += double(c.footprint);
    delta[c.last_use + 1] -= double(c.footprint);
  }
  double run = 0, peak = 0;
  std::size_t arg = 0;
  for (std::size_t p = 0; p < s.num_points; ++p) {
    run += delta[p];
    if (run > peak) {
      peak = run;
      arg = p;
    }  // strict > => lowest-point tie-break
  }
  PeakProfile out;
  out.peak_bytes = peak;
  out.binding_point = arg;
  for (std::size_t i = 0; i < s.cells.size(); ++i)
    if (s.cells[i].first_use <= arg && arg <= s.cells[i].last_use)
      out.live_at_binding.push_back(i);
  return out;
}

///
/// \brief Independent REPLAY ORACLE for \c peak_profile_sweep: the peak live
/// footprint of \p s computed by an explicit per-static-point live-set sum.
///
/// \details Deliberately a DIFFERENT algorithm from \c peak_profile_sweep's
/// +delta/-delta interval-event difference array: for every static point it
/// re-scans all cells and sums the footprints of those whose inclusive
/// `[first_use, last_use]` range covers the point, taking the running max.
/// Because it shares nothing with the sweep's interval bookkeeping, exact
/// agreement between the two on a given \c Schedule is a real cross-check of
/// the sweep's interval logic (validation step 9.6). O(points * cells); used
/// only in tests / analysis, never on the runtime path.
///
inline double peak_profile_replay(Schedule const& s) {
  double peak = 0;
  for (std::size_t p = 0; p < s.num_points; ++p) {
    double sum = 0;
    for (auto const& c : s.cells)
      if (c.first_use <= p && p <= c.last_use) sum += double(c.footprint);
    peak = std::max(peak, sum);
  }
  return peak;
}

///
/// \brief One value group in a RICH linearized schedule (Phase 4a O2 working
/// representation): the same per-value fold as \c Cell, but keeping the
/// pieces \c Cell already collapses into \c footprint -- \c carried and
/// \c home_modes separately -- plus the one genuinely NEW field, \c
/// enclosing_modes, so a later spill pass (O2) has enough to consider
/// demoting a carried mode INTO the home.
///
struct ValueCell {
  std::size_t value_id;  //!< stable index of the value group (== its slot in
                         //!< \c RichSchedule::cells)
  std::size_t hash;      //!< the value's \c EvalExpr::hash_value() -- the CSE
                         //!< identity that folds occurrences into this cell.
                         //!< Batched-slot-BLIND; DISTINCT from the
                         //!< batched-slot-aware occurrence key the router uses
                         //!< (see remat_to_router). Links a cell back to its
                         //!< forest nodes.
  int home_depth;        //!< \c home_depth_of(home_scope, ectx) at the
                         //!< FIRST occurrence -- informational, as \c
                         //!< Cell::home_depth
  container::svector<Index> carried;     //!< canon_indices (same across
                                         //!< occurrences)
  container::svector<Index> home_modes;  //!< the Phase-3b footprint home:
                                         //!< \c r.home MINUS
                                         //!< own_modes_union[hash], read off
                                         //!< the FIRST occurrence
  container::svector<Index>
      enclosing_modes;  //!< NEW: union, over ALL occurrences, of every
                        //!< loop mode that EVER encloses this value (\c
                        //!< ectx[i].first for each level of each
                        //!< occurrence's ectx)
  container::svector<Index>
      divergent_modes;    //!< RELABELED modes: carried by SOME occurrences but
                          //!< not all (union MINUS intersection of the
                          //!< occurrences' canon_indices). Slicing one cannot
                          //!< be shared -- the runtime SPLITS the value, so
                          //!< cell_footprint prices a divergent-mode home 2x.
  std::size_t first_use;  //!< earliest static point the value is live at
  std::size_t last_use;   //!< latest static point the value is live at (its
                          //!< last consumer), inclusive
};

///
/// \brief A whole forest linearized to a flat list of RICH value cells over a
/// single monotone static-point timeline. The \c Schedule consumed by \c
/// peak_profile_sweep is a pure PROJECTION of this (see \c compute_dag_path).
///
struct RichSchedule {
  container::svector<ValueCell> cells;
  std::size_t num_points = 0;  //!< one past the last static point
};

///
/// \brief Linearize an eval \p forest into a \c RichSchedule of RICH value
/// cells over a single post-order static-point timeline.
///
/// \details Stamps the lifetime masks first (so \c home_scope, which reads
/// \c sliced_modes, is populated), then walks every tree in post-order
/// (children before parent),
/// assigning each visited node -- leaves included -- a monotone static point.
/// On descent each node's \c batched_here() loops (proto-expanded, matching
/// \c home_scope's proto expansion) are pushed onto an enclosing-batch-context
/// stack visible to that node's CHILDREN, and popped before the node itself is
/// recorded: a node's own realized loop encloses its operands but not its own
/// (loop-result) value.
///
/// Nodes are then grouped by \c hash_value() -- the same value identity \c
/// CacheManager uses -- into one \c ValueCell per distinct value. Under
/// perfect CSE the group's \c first_use is its single (earliest) production
/// point and its \c last_use is the latest structural consumer (the max
/// parent point over the group; a root with no parent contributes its own
/// point). \c home_depth, \c carried, and \c home_modes are read off the
/// FIRST occurrence (the seed-residency meet is identical across occurrences
/// of a hoisted value); \c enclosing_modes accumulates across EVERY
/// occurrence, since a demoted mode may only enclose the value at SOME of
/// its occurrences.
///
/// \p block_of is any `Index -> std::size_t` callable giving the block
/// (sliced) element count for a mode; only used while walking (to size the
/// enclosing-batch-context entries pushed on descent) -- \p cm is accepted
/// for signature symmetry with \c compute_dag_path / \c detail::cell_footprint
/// but not otherwise used here (no footprint is computed at this stage).
///
template <meta::eval_node_range R, typename BlockOfFn>
RichSchedule compute_dag_boulevard(R const& forest,
                                   [[maybe_unused]] dryrun::CostModel const& cm,
                                   BlockOfFn const& block_of) {
  using Node = std::ranges::range_value_t<R>;

  // Populate home_scope (EvalExpr::sliced_modes) on every internal node.
  stamp_lifetime_masks(forest);

  // Per-occurrence record captured during the post-order walk. consumer_point
  // is the parent's point (its structural consumer); a root keeps its own
  // point (set at construction, overwritten by the parent if any).
  struct NodeRec {
    std::size_t hash;
    std::size_t point;
    std::size_t consumer_point;
    container::svector<Index> home;     // home_scope (proto-expanded)
    container::svector<Index> carried;  // canon_indices
    detail::BatchContext ectx;  // ENCLOSING context (excludes own loops)
    container::svector<Index> own_modes;  // THIS occurrence's OWN realized
                                          // loop modes (proto-expanded) --
                                          // see the own_modes_union note
                                          // below.
  };

  container::svector<NodeRec> recs;
  std::size_t counter = 0;

  auto visit = [&](auto&& self, Node const& n,
                   detail::BatchContext ectx) -> std::size_t {
    // Children see this node's own realized loops on top of the enclosing
    // context; the node itself does NOT (it is recorded with `ectx`).
    detail::BatchContext child_ectx = ectx;
    container::svector<Index> own_modes;
    for (auto const& [ix, kind] : n->batched_here()) {
      container::svector<Index> expanded;
      sequant::detail::proto_expand_into(expanded, ix);
      for (auto const& m : expanded) {
        child_ectx.push_back({m, {std::size_t{0}, block_of(m)}});
        own_modes.push_back(m);
      }
    }

    container::svector<std::size_t> child_recs;
    if (!n.leaf()) {
      child_recs.push_back(self(self, n.left(), child_ectx));
      child_recs.push_back(self(self, n.right(), child_ectx));
    }

    std::size_t const point = counter++;
    NodeRec r;
    r.hash = n->hash_value();
    r.point = point;
    r.consumer_point = point;  // root default; overwritten by parent below
    r.home = home_scope(n);    // proto-expanded sliced_modes (empty on leaves)
    r.carried.assign(n->canon_indices().begin(), n->canon_indices().end());
    r.ectx = std::move(ectx);
    r.own_modes = std::move(own_modes);
    std::size_t const idx = recs.size();
    recs.push_back(std::move(r));
    for (auto ci : child_recs) recs[ci].consumer_point = point;
    return idx;
  };

  for (auto const& tree : forest) visit(visit, tree, detail::BatchContext{});

  // Cross-occurrence union, per canonical value (hash), of the modes that
  // value EVER realizes as its OWN loop (\c batched_here() at the value's own
  // node, not an ancestor's). \c home_scope already folds a node's own
  // batched-here contribution into its meet (stamp_lifetime_masks: "the node
  // AND all its ancestors"), which is right for OTHER consumers of that
  // meet, but wrong for THIS cell's own footprint: a mode a value realizes as
  // its OWN loop slices that value's OPERANDS on the way down, never the
  // value's own result -- the node is the loop's full accumulated output
  // (see the post-order walk's doc comment above: "the node itself does
  // NOT" see its own loop). Subtracting this union keeps the exclusion
  // occurrence-independent by the same reasoning as the meet itself (built
  // from ALL occurrences, not just the one that seeds the cell).
  std::unordered_map<std::size_t, container::svector<Index>> own_modes_union;
  for (auto const& r : recs) {
    auto& acc = own_modes_union[r.hash];
    for (auto const& m : r.own_modes)
      if (std::find(acc.begin(), acc.end(), m) == acc.end()) acc.push_back(m);
  }

  // Per-value RELABELED modes: carried by SOME occurrences but not all -- the
  // UNION minus the INTERSECTION of the occurrences' canon_indices (r.carried).
  // A mode in the intersection appears (at the same canonical slot) in every
  // occurrence, so slicing it is shareable; one present only in the union binds
  // a different physical label in some occurrence, so slicing it forces a SPLIT
  // (see ValueCell::divergent_modes / cell_footprint's 2x pricing).
  std::unordered_map<std::size_t, container::svector<Index>> carried_union,
      carried_isect;
  std::unordered_map<std::size_t, bool> carried_seeded;
  for (auto const& r : recs) {
    auto& u = carried_union[r.hash];
    for (auto const& m : r.carried)
      if (std::find(u.begin(), u.end(), m) == u.end()) u.push_back(m);
    auto& is = carried_isect[r.hash];
    if (!carried_seeded[r.hash]) {
      carried_seeded[r.hash] = true;
      is.assign(r.carried.begin(), r.carried.end());
    } else {
      container::svector<Index> keep;
      for (auto const& m : is)
        if (std::find(r.carried.begin(), r.carried.end(), m) != r.carried.end())
          keep.push_back(m);
      is = std::move(keep);
    }
  }

  // Group occurrences by value identity (hash_value): one ValueCell per
  // group.
  RichSchedule out;
  out.num_points = counter;
  std::unordered_map<std::size_t, std::size_t> hash_to_cell;
  for (auto const& r : recs) {
    // Fold this occurrence's ectx into the running enclosing_modes union
    // (every occurrence contributes, unlike home_modes/carried which are
    // read off the FIRST occurrence only).
    auto fold_enclosing = [&](container::svector<Index>& enclosing_modes) {
      for (auto const& e : r.ectx)
        if (std::find(enclosing_modes.begin(), enclosing_modes.end(),
                      e.first) == enclosing_modes.end())
          enclosing_modes.push_back(e.first);
    };

    auto const it = hash_to_cell.find(r.hash);
    if (it == hash_to_cell.end()) {
      ValueCell c;
      c.value_id = out.cells.size();
      c.hash = r.hash;  // the value's hash_value() (== the hash_to_cell key)
      c.first_use = r.point;
      c.last_use = r.consumer_point;
      c.home_depth = detail::home_depth_of(r.home, r.ectx);
      auto const& self_modes = own_modes_union[r.hash];
      container::svector<Index> home_modes;
      for (auto const& m : r.home)
        if (std::find(self_modes.begin(), self_modes.end(), m) ==
            self_modes.end())
          home_modes.push_back(m);
      c.carried = r.carried;
      c.home_modes = std::move(home_modes);
      {
        auto const& u = carried_union[r.hash];
        auto const& is = carried_isect[r.hash];
        for (auto const& m : u)
          if (std::find(is.begin(), is.end(), m) == is.end())
            c.divergent_modes.push_back(m);
      }
      fold_enclosing(c.enclosing_modes);
      hash_to_cell.emplace(r.hash, c.value_id);
      out.cells.push_back(std::move(c));
    } else {
      ValueCell& c = out.cells[it->second];
      c.first_use = std::min(c.first_use, r.point);
      c.last_use = std::max(c.last_use, r.consumer_point);
      fold_enclosing(c.enclosing_modes);
    }
  }
  return out;
}

///
/// \brief Linearize an eval \p forest into a \c Schedule of value cells over a
/// single post-order static-point timeline, ready for \c peak_profile_sweep.
///
/// \details A thin PROJECTION of \c compute_dag_boulevard: runs the one
/// post-order walk there, then collapses each \c ValueCell's \c carried / \c
/// home_modes down to a single \c Cell::footprint via \c
/// detail::cell_footprint. The returned \c Schedule is BYTE-IDENTICAL to what
/// the (pre-Phase-4a) inline-walk version of \c compute_dag_path produced -- \c
/// enclosing_modes is the only new information \c compute_dag_boulevard
/// computes, and it does not reach this projection.
///
/// \p block_of is any `Index -> std::size_t` callable giving the block
/// (sliced) element count for a mode; forwarded to \c detail::cell_footprint.
///
template <meta::eval_node_range R, typename BlockOfFn>
Schedule compute_dag_path(R const& forest, dryrun::CostModel const& cm,
                          BlockOfFn const& block_of) {
  RichSchedule const rich = compute_dag_boulevard(forest, cm, block_of);

  Schedule out;
  out.num_points = rich.num_points;
  out.cells.reserve(rich.cells.size());
  for (auto const& vc : rich.cells) {
    Cell c;
    c.value_id = vc.value_id;
    c.home_depth = vc.home_depth;
    c.footprint = detail::cell_footprint(vc.carried, vc.home_modes, cm,
                                         block_of, vc.divergent_modes);
    c.first_use = vc.first_use;
    c.last_use = vc.last_use;
    out.cells.push_back(c);
  }
  return out;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_PEAK_PROFILE_HPP
