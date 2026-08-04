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
/// forest linearization (`linearize`) and the interval-event sweep
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
template <typename BlockOfFn>
[[nodiscard]] inline std::size_t cell_footprint(
    container::svector<Index> const& carried,
    container::svector<Index> const& home_modes, BatchContext const& ectx,
    dryrun::CostModel const& cm, BlockOfFn const& block_of) {
  int const d = home_depth_of(home_modes, ectx);
  dryrun::ExtentOverrides ov;
  for (int i = 0; i <= d; ++i) {  // loops at level <= d ENCLOSE the home
    Index const& m = ectx[i].first;
    if (std::find(carried.begin(), carried.end(), m) != carried.end())
      ov[m] = block_of(m);  // sliced inside an enclosing loop
  }
  return cm.memsize(carried, ov);  // non-overridden carried modes size FULL
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
/// \brief Linearize an eval \p forest into a \c Schedule of value cells over a
/// single post-order static-point timeline, ready for \c peak_profile_sweep.
///
/// \details Stamps the Phase-3a seed residency first (so \c home_scope is
/// populated), then walks every tree in post-order (children before parent),
/// assigning each visited node -- leaves included -- a monotone static point.
/// On descent each node's \c batched_here() loops (proto-expanded, matching
/// \c home_scope's proto expansion) are pushed onto an enclosing-batch-context
/// stack visible to that node's CHILDREN, and popped before the node itself is
/// recorded: a node's own realized loop encloses its operands but not its own
/// (loop-result) value.
///
/// Nodes are then grouped by \c hash_value() -- the same value identity \c
/// CacheManager uses -- into one \c Cell per distinct value. Under perfect CSE
/// the group's \c first_use is its single (earliest) production point and its
/// \c last_use is the latest structural consumer (the max parent point over
/// the group; a root with no parent contributes its own point). The home depth
/// and footprint are read off any occurrence (the seed-residency meet is
/// identical across occurrences of a hoisted value).
///
/// \p block_of is any `Index -> std::size_t` callable giving the block
/// (sliced) element count for a mode; forwarded to \c detail::cell_footprint.
///
template <meta::eval_node_range R, typename BlockOfFn>
Schedule linearize(R const& forest, dryrun::CostModel const& cm,
                   BlockOfFn const& block_of) {
  using Node = std::ranges::range_value_t<R>;

  // Populate home_scope (EvalExpr::seed_residency) on every internal node.
  stamp_seed_residency(forest);

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
  };

  container::svector<NodeRec> recs;
  std::size_t counter = 0;

  auto visit = [&](auto&& self, Node const& n,
                   detail::BatchContext ectx) -> std::size_t {
    // Children see this node's own realized loops on top of the enclosing
    // context; the node itself does NOT (it is recorded with `ectx`).
    detail::BatchContext child_ectx = ectx;
    for (auto const& [ix, kind] : n->batched_here()) {
      container::svector<Index> expanded;
      sequant::detail::proto_expand_into(expanded, ix);
      for (auto const& m : expanded)
        child_ectx.push_back({m, {std::size_t{0}, block_of(m)}});
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
    r.home = home_scope(n);  // proto-expanded seed residency (empty on leaves)
    r.carried.assign(n->canon_indices().begin(), n->canon_indices().end());
    r.ectx = std::move(ectx);
    std::size_t const idx = recs.size();
    recs.push_back(std::move(r));
    for (auto ci : child_recs) recs[ci].consumer_point = point;
    return idx;
  };

  for (auto const& tree : forest) visit(visit, tree, detail::BatchContext{});

  // Group occurrences by value identity (hash_value): one Cell per group.
  Schedule out;
  out.num_points = counter;
  std::unordered_map<std::size_t, std::size_t> hash_to_cell;
  for (auto const& r : recs) {
    auto const it = hash_to_cell.find(r.hash);
    if (it == hash_to_cell.end()) {
      Cell c;
      c.value_id = out.cells.size();
      c.first_use = r.point;
      c.last_use = r.consumer_point;
      c.home_depth = detail::home_depth_of(r.home, r.ectx);
      c.footprint =
          detail::cell_footprint(r.carried, r.home, r.ectx, cm, block_of);
      hash_to_cell.emplace(r.hash, c.value_id);
      out.cells.push_back(std::move(c));
    } else {
      Cell& c = out.cells[it->second];
      c.first_use = std::min(c.first_use, r.point);
      c.last_use = std::max(c.last_use, r.consumer_point);
    }
  }
  return out;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_PEAK_PROFILE_HPP
