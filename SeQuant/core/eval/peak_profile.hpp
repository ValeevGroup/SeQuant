#ifndef SEQUANT_EVAL_PEAK_PROFILE_HPP
#define SEQUANT_EVAL_PEAK_PROFILE_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/lifetime_mask.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <map>
#include <optional>
#include <set>
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
/// \p divergent_modes is now INFORMATIONAL only. It once triggered a flat 2x
/// pricing fudge here (the placeholder `46b495eba` shipped): a home that sliced
/// a RELABELED mode was priced as TWO co-resident copies. That flat 2x captured
/// neither of the split's two real costs (peak co-residency and the replication
/// recompute) and silently dropped the dominant, mis-priced recompute term. It
/// is DELETED: a divergent value that remat homes at a sub-scope is UN-FOLDED
/// into two real, non-divergent \c ValueCell s (see \c apply_split in
/// placement_remat.hpp), each priced ONCE here at its own home; peak
/// co-residency is priced structurally by \c peak_profile_sweep over the two
/// cells' liveness intervals (which keys on \c value_id, so two cells of one
/// hash need no sweep change), and the replication recompute is a SEPARATE,
/// report-only term (\c apply_split's return / \c
/// RematResult::modeled_recompute). The parameter is retained for signature
/// compatibility with existing callers.
template <typename BlockOfFn>
[[nodiscard]] inline std::size_t cell_footprint(
    container::svector<Index> const& carried,
    container::svector<Index> const& home_modes, dryrun::CostModel const& cm,
    BlockOfFn const& block_of,
    [[maybe_unused]] container::svector<Index> const& divergent_modes = {}) {
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
  return cm.memsize(carried, ov);  // non-meet carried modes FULL
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
/// \brief One USE-SITE of a value, kept alongside the folded \c ValueCell so a
/// CSE-aware remat SPLIT (see \c apply_split) can partition the value's
/// occurrences by their PHYSICAL binding of a relabeled mode and re-derive each
/// split cell's subset-local \c carried / \c home / liveness / enclosing nest.
///
/// \details These are the per-\c NodeRec fields \c compute_dag_boulevard
/// computes during its post-order walk and once DISCARDED at grouping (keeping
/// only the first occurrence's home/carried + the union/min/max). They are now
/// retained: (a) two occurrences that bind a relabeled mode to different
/// physical labels (the g.C legs' \c i_3 vs \c i_4) are told apart by their \c
/// carried; (b) each split cell's replication factor is a product over the
/// levels it is homed-within-but-does-not-carry, read from its subset-local \c
/// ectx (enclosing loops) minus \c carried.
///
struct OccurrenceRec {
  std::size_t point;           //!< this occurrence's production static point
  std::size_t consumer_point;  //!< its structural consumer (parent) point
  container::svector<Index> carried;  //!< this occurrence's canon_indices
  container::svector<Index> home;     //!< home_scope (plain modes)
  detail::BatchContext ectx;  //!< ENCLOSING loops (excludes this node's own)
  //!< The static point of the occurrence (in the same tree) that OPENED each
  //!< \c ectx entry, parallel to \c ectx: names the loop INSTANCE each entry
  //!< is, so the boulevard can read the DP's nesting order off the
  //!< occurrences (RichSchedule::loop_order).
  container::svector<std::size_t> ectx_opener_point;
  //!< Static points of this occurrence's operand occurrences, left then
  //!< right (empty on a leaf). An operand's LEG is its index here; the
  //!< sliced-mode seam attributes its facts per leg, so one value read on
  //!< both legs of a node under different labels (a self-product of a shared
  //!< intermediate) gets two distinct slicings.
  container::svector<std::size_t> operand_points;
  //!< The loops this occurrence's node OPENS, with their kind (Contracted:
  //!< a mode this node contracts in batches; External: a carried mode whose
  //!< physical loop is introduced here). Carried into
  //!< RichSchedule::loop_kind once the loop slots are numbered.
  container::svector<std::pair<Index, BatchModeType>> opens;
  //!< Task 2 (loop identity): per \c carried position, the \c loop_slot of the
  //!< batch loop that slices it (which MEMBER of its same-space group), or -1
  //!< where the position is not a batched (loop-sliced) mode. Assigned by the
  //!< union-find over producer->consumer slot connectivity in \c
  //!< compute_dag_boulevard (spec 2026-08-28 sec.4-5). Parallel to \c carried.
  container::svector<int> loop_slot;
  //!< Loop identity for the modes this occurrence's value CONTRACTS (reduces)
  //!< in batches at its own node: assigned either by uniting a producing
  //!< operand's HOME-SLICED carried-mode node with a synthetic reduction node
  //!< (the reduction loop and the operand's slice loop are one physical loop
  //!< and must share \c loop_slot), or -- when no operand is home-sliced on
  //!< the mode -- by seeding that synthetic node directly (see \c
  //!< contracted_batched below). A reduced mode has no \c carried position,
  //!< so its slot is recorded here as (mode, loop_slot). Read by \c
  //!< ordered_schedule's \c fusion_slot when it places a Reduction escape; a
  //!< Reduction mode that still resolves to no slot here is a hard error
  //!< there (\c build_ordered_schedule throws in its escape-placement loop),
  //!< not a slot-0 default.
  container::svector<std::pair<Index, int>> reduced_slot;
  //!< Modes THIS occurrence's value contracts IN BATCHES at its own node
  //!< (the legality \c build_site_of CONTRACTED test, mirrored -- see \c
  //!< NodeRec::contracted_batched). A mode here owns a loop identity even
  //!< when no operand of the contraction is itself home-sliced on it (an
  //!< all-input reduction): the union-find below seeds a component for every
  //!< entry here that \c classify_axis would actually call \c Reduction (no
  //!< carried index of the SAME SPACE as the mode), instead of relying
  //!< solely on a home-sliced child to create one; an entry beside a
  //!< same-space carried index is \c classify_axis LoopLocal/LoopCarried, not
  //!< Reduction, and is left to the ordinary carried-position path.
  container::svector<Index> contracted_batched;
};

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
  bool is_leaf = false;  //!< the value is a forest LEAF (an input fetched on
                         //!< demand), not a computed intermediate -- so it is
                         //!< never scheduled as a BuildStep.
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
                          //!< be shared -- remat SPLITS the value into two
                          //!< non-divergent cells (see \c apply_split); the
                          //!< split cells then carry an EMPTY divergent_modes.
  std::size_t first_use;  //!< earliest static point the value is live at
  std::size_t last_use;   //!< latest static point the value is live at (its
                          //!< last consumer), inclusive
  container::svector<OccurrenceRec>
      occurrences;  //!< every use-site of this value (retained so a remat
                    //!< SPLIT can partition them by physical binding and
                    //!< re-derive each split cell's subset-local records).
                    //!< A split cell (one occurrence subset) keeps only its
                    //!< subset here.
};

///
/// \brief A whole forest linearized to a flat list of RICH value cells over a
/// single monotone static-point timeline. The \c Schedule consumed by \c
/// peak_profile_sweep is a pure PROJECTION of this (see \c compute_dag_path).
///
struct RichSchedule {
  container::svector<ValueCell> cells;
  std::size_t num_points = 0;  //!< one past the last static point
  //!< The kind of every numbered loop instance, keyed by (space base_key,
  //!< loop_slot): Contracted when the open that created the instance
  //!< contracts the mode in batches at its node, External when it introduces
  //!< a carried mode's physical loop. A space may hold instances of both kinds
  //!< (an occupied pair contracted in batches beside an occupied external
  //!< pair), so the kind is a property of the INSTANCE, not of the space; the
  //!< ordered schedule builder reads a block's kind here.
  std::map<std::pair<std::wstring, int>, BatchModeType> loop_kind;
  //!< Loop NESTING constraints read off the DP's realization: (outer, inner)
  //!< pairs of loop instances, each (space base_key, loop_slot), such that
  //!< some occurrence sits inside `outer` and `inner` is opened inside it
  //!< (a consecutive pair of its enclosing context, or its enclosing context
  //!< and a loop it opens itself). The ordered schedule builder nests the
  //!< realized chain to satisfy every pair (a contradiction is a builder
  //!< error: the loop identity fused two loops that nest in opposite orders).
  //!< The mapped value is a witness: the value id of the first occurrence
  //!< that produced the pair (diagnostics only).
  std::map<
      std::pair<std::pair<std::wstring, int>, std::pair<std::wstring, int>>,
      std::size_t>
      loop_order;
};

///
/// \brief Linearize an eval \p forest into a \c RichSchedule of RICH value
/// cells over a single post-order static-point timeline.
///
/// \details Stamps the lifetime masks first (so \c home_scope, which reads
/// \c sliced_modes, is populated), then walks every tree in post-order
/// (children before parent),
/// assigning each visited node -- leaves included -- a monotone static point.
/// On descent each node's \c node_slice_mask() loops are pushed onto an
/// enclosing-batch-context
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
    bool is_leaf;
    std::size_t point;
    std::size_t consumer_point;
    container::svector<Index> home;     // home_scope
    container::svector<Index> carried;  // canon_indices
    detail::BatchContext ectx;  // ENCLOSING context (excludes own loops)
    // Pre-order id of the node that OPENED each ectx entry (parallel to
    // ectx); resolved to that node's static point on the occurrence record.
    container::svector<std::size_t> ectx_opener;
    // Static points of this node's operands, left then right (empty on a
    // leaf): names each operand OCCURRENCE by its leg, so a value that is
    // both operands of one node under different labels keeps two identities.
    container::svector<std::size_t> operand_points;
    // The batch loops this node OPENS (batch_loops_opened_here), with the
    // kind of each: the source of every loop instance's kind (see
    // RichSchedule::loop_kind).
    container::svector<std::pair<Index, BatchModeType>> opens;
    container::svector<Index> own_modes;  // THIS occurrence's OWN realized
                                          // loop modes -- see the
                                          // own_modes_union note below.
    container::svector<Index>
        contracted_batched;  // modes THIS node contracts in batches at its
                             // OWN node (mirrors legality::build_site_of's
                             // CONTRACTED test -- contracted_indices(n)
                             // intersected with a Contracted-kind
                             // node_slice_mask() stamp -- inlined here rather
                             // than called, since legality.hpp already
                             // depends on this header). A value that reduces
                             // one of these modes owns a loop identity for
                             // it even when no operand is itself home-sliced
                             // on the mode (see the union-find pass below).
  };

  container::svector<NodeRec> recs;
  std::size_t counter = 0;
  // Pre-order ids name a node BEFORE its subtree is walked (its post-order
  // point is not known yet when the children record it as their opener).
  std::size_t pre_counter = 0;
  std::unordered_map<std::size_t, std::size_t> pre_to_point;

  auto visit = [&](auto&& self, Node const& n, detail::BatchContext ectx,
                   container::svector<std::size_t> ectx_opener) -> std::size_t {
    std::size_t const pre = pre_counter++;
    // Children see this node's own realized loops on top of the enclosing
    // context; the node itself does NOT (it is recorded with `ectx`).
    detail::BatchContext child_ectx = ectx;
    container::svector<std::size_t> child_opener = ectx_opener;
    container::svector<Index> own_modes;
    // Enclosing-loop context is built from the loops OPENED at this node
    // (batch_loops_opened_here), NOT the per-node sliced mask
    // (node_slice_mask). The DP stamps an external mode's sliced mask on EVERY
    // carrying node, so reading node_slice_mask here counted one physical loop
    // once-per-carrying-node
    // -- ectx piled up duplicates of the same occ index (i i i ...) and no
    // longer matched the DAG scope's one-loop-one-level de-duplication. Opens
    // name each physical loop exactly once, so child_ectx is a true loop nest.
    // own_modes likewise: a value's OWN loop is the one it OPENS ("its own
    // node, not an ancestor's" -- own_modes doc), which the sliced mask
    // over-reported for every carried descendant.
    // Opened loops are over plain occ/aux modes (never a composite result
    // slot), so each is taken as itself.
    for (auto const& [ix, kind] : n->batch_loops_opened_here()) {
      child_ectx.push_back({ix, {std::size_t{0}, block_of(ix)}});
      child_opener.push_back(pre);
      own_modes.push_back(ix);
    }

    // TEMP instrumentation (Task 2 verify): dump this node's opened_here (the
    // group nest structure the factorizer emits) -- which same-space modes open
    // at ONE node (a multi-loop group) vs at different nodes (separate groups).
    // Guarded by SEQUANT_DUMP_OPENS.
    if (std::getenv("SEQUANT_DUMP_OPENS") &&
        !n->batch_loops_opened_here().empty()) {
      std::wcerr << L"[opens] hash=" << n->hash_value() << L" opened_here={";
      for (auto const& [ix, kind] : n->batch_loops_opened_here())
        std::wcerr << ix.full_label() << L":" << ix.space().base_key() << L":"
                   << (kind == BatchModeType::Contracted ? L"C" : L"E") << L" ";
      std::wcerr << L"} carried={";
      for (auto const& c : n->canon_indices())
        std::wcerr << c.full_label() << L" ";
      std::wcerr << L"}\n";
    }

    container::svector<std::size_t> child_recs;
    if (!n.leaf()) {
      child_recs.push_back(self(self, n.left(), child_ectx, child_opener));
      child_recs.push_back(self(self, n.right(), child_ectx, child_opener));
    }

    std::size_t const point = counter++;
    NodeRec r;
    r.hash = n->hash_value();
    r.is_leaf = n.leaf();
    r.point = point;
    r.consumer_point = point;  // root default; overwritten by parent below
    r.home = home_scope(n);    // sliced_modes (empty on leaves)
    r.carried.assign(n->canon_indices().begin(), n->canon_indices().end());
    // Modes THIS node contracts in batches at its own node -- mirrors
    // legality::build_site_of's CONTRACTED test (eval.hpp's
    // contracted_indices(n), intersected with a Contracted-kind
    // node_slice_mask() stamp) verbatim, inlined rather than shared: this
    // header cannot include legality.hpp / eval.hpp (legality.hpp already
    // depends on peak_profile.hpp).
    if (!n.leaf() && n->is_product()) {
      auto const& l = n.left()->canon_indices();
      auto const& rr = n.right()->canon_indices();
      auto const& c = n->canon_indices();
      auto const contains = [](auto const& vec, Index const& ix) {
        return std::find(vec.begin(), vec.end(), ix) != vec.end();
      };
      auto const& stamps = n->node_slice_mask();
      for (Index const& ix : l) {
        if (!contains(rr, ix) || contains(c, ix))
          continue;  // not contracted at this node
        bool const batched =
            std::any_of(stamps.begin(), stamps.end(), [&](auto const& p) {
              return p.second == BatchModeType::Contracted && p.first == ix;
            });
        if (batched && !contains(r.contracted_batched, ix))
          r.contracted_batched.push_back(ix);
      }
    }
    r.ectx = std::move(ectx);
    r.ectx_opener = std::move(ectx_opener);
    pre_to_point[pre] = point;
    r.own_modes = std::move(own_modes);
    for (auto const& [ix, kind] : n->batch_loops_opened_here())
      r.opens.push_back({ix, kind});
    for (auto ci : child_recs) r.operand_points.push_back(recs[ci].point);
    std::size_t const idx = recs.size();
    recs.push_back(std::move(r));
    for (auto ci : child_recs) recs[ci].consumer_point = point;
    return idx;
  };

  for (auto const& tree : forest)
    visit(visit, tree, detail::BatchContext{}, {});

  // Cross-occurrence union, per canonical value (hash), of the modes that
  // value EVER realizes as its OWN loop (\c node_slice_mask() at the value's
  // own node, not an ancestor's). \c home_scope already folds a node's own
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

    // The per-occurrence record retained on the cell for a CSE-aware split.
    auto make_occ = [&]() -> OccurrenceRec {
      OccurrenceRec o;
      o.point = r.point;
      o.consumer_point = r.consumer_point;
      o.carried = r.carried;
      o.home = r.home;
      o.ectx = r.ectx;
      for (std::size_t pre : r.ectx_opener)
        o.ectx_opener_point.push_back(pre_to_point.at(pre));
      o.contracted_batched = r.contracted_batched;
      o.opens = r.opens;
      o.operand_points = r.operand_points;
      return o;
    };

    auto const it = hash_to_cell.find(r.hash);
    if (it == hash_to_cell.end()) {
      ValueCell c;
      c.value_id = out.cells.size();
      c.is_leaf = r.is_leaf;
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
      c.occurrences.push_back(make_occ());
      hash_to_cell.emplace(r.hash, c.value_id);
      out.cells.push_back(std::move(c));
    } else {
      ValueCell& c = out.cells[it->second];
      c.first_use = std::min(c.first_use, r.point);
      c.last_use = std::max(c.last_use, r.consumer_point);
      fold_enclosing(c.enclosing_modes);
      c.occurrences.push_back(make_occ());
    }
  }

  // Task 2 (loop identity, spec 2026-08-28 sec.4-5): assign a per-occurrence
  // `loop_slot` -- which MEMBER of a same-space loop group slices each batched
  // carried mode.
  //
  // A loop's identity comes from PHYSICAL producer->consumer connectivity,
  // never from the loop-colored value-id (that would be circular -- the
  // value-id needs `loop_slot`). Loops are grouped by connected components; a
  // component is one physical loop and its `loop_slot` numbers the components
  // of each space.
  //
  // A component is a set of (value_id, STRUCTURAL slot) nodes, where the slot
  // is a mode's position in `canon_indices`. That order is invariant across a
  // value's occurrences (same structure => same canonical order; only the
  // LABELS differ between terms), so keying by position folds the occurrences
  // across trees. Edges connect a child value's slot to its parent's slot for
  // the SAME physical mode, matched by `Index ==` in the PARENT OCCURRENCE's
  // own frame -- the child and its parent are in one tree, so labels agree
  // there (eval.hpp `contracted_indices` relies on exactly this). A same-space
  // transposition between two occurrences (spec sec.4(b)) surfaces as a
  // rejected union (below), not as a relabeling of the position key; recording
  // it is deferred.
  //
  // Symmetry stays OUT of this pass: a symmetric tensor's two modes are still
  // two distinct physical loops; folding `A(_,i)` with `A(i,_)` is a downstream
  // VALUE-ID decision (Task 4) via the existing loop-colored occurrence_key.
  // Where connectivity is contradictory (spec sec.4 C/D) or symmetry leaves the
  // mapping free, the safe direction is to KEEP loops distinct (reject the
  // merge): an over-split loop is a missed fusion (still correct), a collapse
  // is the crash.
  {
    std::unordered_map<std::size_t, std::size_t>
        point_value;  // point->value_id
    std::unordered_map<std::size_t, OccurrenceRec const*> point_occ;
    for (ValueCell const& c : out.cells)
      for (OccurrenceRec const& o : c.occurrences) {
        point_value[o.point] = c.value_id;
        point_occ[o.point] = &o;
      }

    std::size_t constexpr POS_BITS = 20;  // a value's index count is tiny
    auto const encode = [](std::size_t vid, std::size_t pos) -> std::size_t {
      return (vid << POS_BITS) | pos;
    };
    auto const dec_vid = [](std::size_t n) { return n >> POS_BITS; };
    auto const dec_pos = [](std::size_t n) {
      return n & ((std::size_t{1} << POS_BITS) - 1);
    };
    // Reduction-mode nodes live in the UPPER half of a value's position space
    // (>= CONTRACTED_BASE) so they never collide with a carried position (the
    // lower half). A value's index count is tiny, so both halves are ample. A
    // (value, reduced-mode) pair maps to one stable synthetic position, keyed
    // by the mode's canonical label in the PARENT's own frame (the same frame
    // the carried-mode edge match uses), so two operands reducing the SAME
    // physical mode at one parent resolve to ONE reduction node.
    std::size_t constexpr CONTRACTED_BASE = std::size_t{1} << (POS_BITS - 1);
    std::map<std::pair<std::size_t, std::wstring>, std::size_t> red_pos;
    std::size_t red_next = CONTRACTED_BASE;
    auto const reduction_node = [&](std::size_t par_vid,
                                    Index const& m) -> std::size_t {
      auto const key = std::make_pair(par_vid, std::wstring{m.full_label()});
      auto const it = red_pos.find(key);
      if (it != red_pos.end()) return encode(par_vid, it->second);
      std::size_t const pos = red_next++;
      SEQUANT_ASSERT(red_next < (std::size_t{1} << POS_BITS));
      red_pos[key] = pos;
      return encode(par_vid, pos);
    };
    // (parent value_id, reduced mode) pairs to stamp onto the parent's
    // occurrences once components are numbered.
    container::svector<std::pair<std::size_t, Index>> reduction_stamps;

    std::unordered_map<std::size_t, std::size_t>
        uf;  // node -> parent (self=root)
    // CONFLICT-AWARE union-find: a component is one physical loop, so it must
    // never hold two DISTINCT slots of one value. members[root] maps value_id
    // -> the single slot that value contributes; a union that would violate
    // this is REJECTED (keep loops distinct -- the safe direction) rather than
    // merged, which otherwise cascades a single transposition into a per-space
    // collapse.
    std::unordered_map<std::size_t, std::map<std::size_t, std::size_t>> members;
    auto find = [&](std::size_t x) -> std::size_t {
      auto it = uf.find(x);
      if (it == uf.end()) {
        uf.emplace(x, x);
        members[x][dec_vid(x)] = dec_pos(x);
        return x;
      }
      std::size_t root = x;
      while (uf[root] != root) root = uf[root];
      while (uf[x] != root) {  // path-compress
        std::size_t const nxt = uf[x];
        uf[x] = root;
        x = nxt;
      }
      return root;
    };
    auto try_unite = [&](std::size_t a, std::size_t b) -> bool {
      std::size_t const ra = find(a), rb = find(b);
      if (ra == rb) return true;
      auto& ma = members[ra];
      auto& mb = members[rb];
      for (auto const& [vid, pos] : ma) {
        auto const jt = mb.find(vid);
        if (jt != mb.end() && jt->second != pos) return false;  // conflict
      }
      for (auto const& [vid, pos] : ma) mb[vid] = pos;  // merge ma -> mb
      members.erase(ra);
      uf[ra] = rb;
      return true;
    };

    // A value's mode gets a loop component ONLY where the value is HOME-SLICED
    // on it (materialized slice-by-slice), NOT merely USED inside a loop
    // (occ.ectx). A full-homed value materializes nothing slice-by-slice, so it
    // must contribute NO loop -- otherwise its varying use-sites (bound to i_1
    // at one consumer, i_2 at another) force the union-find to invent a
    // spurious extra member. Its per-consumer USE-site slices are a separate,
    // derived fact.
    auto const is_batched = [](OccurrenceRec const& occ,
                               Index const& m) -> bool {
      for (auto const& h : occ.home)
        if (h == m) return true;
      return false;
    };

    bool const conflict_dump = std::getenv("SEQUANT_DUMP_LOOP_SLOT") != nullptr;
    // Edges + node creation, keyed by CANONICAL slot (occurrence-invariant).
    for (ValueCell const& c : out.cells)
      for (OccurrenceRec const& occ : c.occurrences) {
        OccurrenceRec const* par = nullptr;
        std::size_t par_vid = 0;
        if (occ.consumer_point != occ.point) {
          auto const pit = point_occ.find(occ.consumer_point);
          auto const vit = point_value.find(occ.consumer_point);
          if (pit != point_occ.end() && vit != point_value.end()) {
            par = pit->second;
            par_vid = vit->second;
          }
        }
        for (std::size_t pV = 0; pV < occ.carried.size(); ++pV) {
          Index const& m = occ.carried[pV];
          if (!is_batched(occ, m)) continue;
          (void)find(encode(c.value_id, pV));  // seed the node
          if (!par) continue;
          // Match the SAME physical mode in the parent OCCURRENCE's own frame
          // (same tree as the child, so labels agree -- contracted_indices
          // relies on this). NOT cell.carried, whose labels come from a
          // possibly- different tree. Node keys use STRUCTURAL position
          // (canon_indices order is invariant across a value's occurrences;
          // only the labels differ), which folds the value's occurrences across
          // trees.
          auto const pj =
              std::find(par->carried.begin(), par->carried.end(), m);
          if (pj == par->carried.end()) {
            // CONTRACTED at parent (m is home-sliced by this child but absent
            // from the parent's result): the parent REDUCES m, and its
            // reduction loop is the SAME physical loop as this child's slice
            // loop -- they must share loop_slot. A reduced mode has no carried
            // position at the parent, so unite the child's carried-mode node
            // with a synthetic reduction node for (parent, m); the parent's
            // reduction escape then inherits this slot instead of landing in
            // a different same-space nest than the operand it reads (the
            // eviction this fixes) or, absent any slot at all, hitting
            // build_ordered_schedule's hard-error throw. try_unite stays
            // conflict-aware, so a genuine transposition is still kept
            // distinct.
            std::size_t const rn = reduction_node(par_vid, m);
            if (try_unite(encode(c.value_id, pV), rn))
              reduction_stamps.push_back({par_vid, m});
            continue;
          }
          std::size_t const pC =
              static_cast<std::size_t>(pj - par->carried.begin());
          if (!try_unite(encode(c.value_id, pV), encode(par_vid, pC)) &&
              conflict_dump) {
            std::wcerr << L"[loop_slot] REJECTED edge (value " << c.value_id
                       << L" pos " << pV << L") ~ (value " << par_vid
                       << L" pos " << pC << L") mode " << m.full_label()
                       << L" -- would collapse two members; kept distinct\n";
            auto lbls = [](container::svector<Index> const& v) {
              std::wstring s;
              for (auto const& ix : v) {
                s += std::wstring(ix.full_label());
                s += L" ";
              }
              return s;
            };
            std::wcerr << L"[loop_slot]   child hash=" << (c.hash % 100000u)
                       << L" cell.carried=[" << lbls(c.carried)
                       << L"] occ.carried=[" << lbls(occ.carried) << L"]\n";
            std::wcerr << L"[loop_slot]   parent hash="
                       << (out.cells[par_vid].hash % 100000u)
                       << L" cell.carried=[" << lbls(out.cells[par_vid].carried)
                       << L"] occ.carried=[" << lbls(par->carried) << L"]\n";
          }
        }
      }

    // A value's OWN contracted-in-batches modes (the legality build_site_of
    // CONTRACTED test, mirrored as NodeRec::contracted_batched /
    // OccurrenceRec::contracted_batched above) own a loop identity even when
    // every operand of the contraction is an INPUT -- a leaf, or a
    // root-resident value that is not itself home-sliced on the mode. With no
    // home-sliced child there is no edge for the loop above to unite through,
    // so it never creates a component for (value, mode) and the value's
    // reduction escape is later placed with no loop identity at all. Seed the
    // reduction node here for every occurrence's own contracted_batched mode.
    // The child-driven union above still fires whenever a home-sliced operand
    // exists and unites its slice loop into this same component
    // (reduction_node is idempotent -- same (value, mode) always encodes to
    // the same synthetic node).
    // Seeded for EVERY mode contracted in batches at the node, whether or not
    // the value also carries an index of that space: legality classifies the
    // contracted axis per axis (classify_axis decides by the axis's identity,
    // not its space), so such a mode is a Reduction on its own loop instance
    // and needs a loop identity the escape placement can resolve. (An earlier
    // guard skipped the same-space-carried case to mirror the per-space
    // classification that is now gone.)
    for (ValueCell const& c : out.cells)
      for (OccurrenceRec const& occ : c.occurrences)
        for (Index const& m : occ.contracted_batched) {
          (void)find(reduction_node(c.value_id, m));
          reduction_stamps.push_back({c.value_id, m});
        }

    // Number the components: one loop_slot per component, ranked per SPACE in
    // first-seen order over each value's canonical frame. A component is
    // single-space (edges join only identical physical modes).
    std::unordered_map<std::size_t, int> root_slot;  // root -> loop_slot
    std::map<std::wstring, int> next_slot;           // space -> next slot #
    for (ValueCell const& c : out.cells)
      for (OccurrenceRec const& occ : c.occurrences)
        for (std::size_t pV = 0; pV < occ.carried.size(); ++pV) {
          if (!is_batched(occ, occ.carried[pV])) continue;
          std::size_t const root = find(encode(c.value_id, pV));
          if (root_slot.find(root) != root_slot.end()) continue;
          std::wstring const sp{occ.carried[pV].space().base_key()};
          root_slot.emplace(root, next_slot[sp]++);
        }
    // A reduction-only component -- one seeded above by contracted_batched
    // that no home-sliced operand ever united into -- never surfaces in the
    // carried-position walk just above (it walks occ.carried positions, and
    // a reduced mode has none), so it would otherwise stay unnumbered and
    // the final stamping loop below would silently skip it. Number it here;
    // find() already resolves to the SAME root as above wherever a union did
    // happen, so root_slot.find de-duplicates that case for free.
    for (auto const& [par_vid, m] : reduction_stamps) {
      std::size_t const root = find(reduction_node(par_vid, m));
      if (root_slot.find(root) != root_slot.end()) continue;
      std::wstring const sp{m.space().base_key()};
      root_slot.emplace(root, next_slot[sp]++);
    }

    // Component-membership dump: per (space, slot), the distinct member values
    // (hash:pos:mode). Two same-slot fragments of one physical loop appear as
    // DIFFERENT slots here; this shows which values fragmented apart. Guarded
    // by SEQUANT_DUMP_LOOP_SLOT.
    if (conflict_dump) {
      std::map<std::pair<std::wstring, int>, std::map<std::wstring, int>> comp;
      for (ValueCell const& c : out.cells)
        for (OccurrenceRec const& occ : c.occurrences)
          for (std::size_t pV = 0; pV < occ.carried.size(); ++pV) {
            if (!is_batched(occ, occ.carried[pV])) continue;
            std::size_t const root = find(encode(c.value_id, pV));
            std::wstring const sp{occ.carried[pV].space().base_key()};
            std::wstring const mem =
                std::to_wstring(c.hash % 100000u) + L"(" +
                std::wstring(occ.carried[pV].full_label()) + L")";
            comp[{sp, root_slot.at(root)}][mem] = 1;
          }
      for (auto const& [key, members] : comp) {
        std::wcerr << L"[comp] " << key.first << L"#slot" << key.second
                   << L" members={ ";
        for (auto const& [mem, _] : members) std::wcerr << mem << L" ";
        std::wcerr << L"}\n";
      }
    }

    // Stamp each occurrence's per-position loop_slot (structural position keys
    // fold occurrences of one value across trees).
    bool const dump = conflict_dump;
    for (ValueCell& c : out.cells)
      for (OccurrenceRec& occ : c.occurrences) {
        occ.loop_slot.assign(occ.carried.size(), -1);
        for (std::size_t pV = 0; pV < occ.carried.size(); ++pV) {
          if (!is_batched(occ, occ.carried[pV])) continue;
          occ.loop_slot[pV] = root_slot.at(find(encode(c.value_id, pV)));
        }
        if (dump) {
          std::wcerr << L"[loop_slot] value_id=" << c.value_id << L" hash="
                     << (c.hash % 100000u) << L" point=" << occ.point
                     << L" carried={";
          for (std::size_t k = 0; k < occ.carried.size(); ++k)
            std::wcerr << occ.carried[k].full_label() << L":"
                       << occ.loop_slot[k] << L" ";
          std::wcerr << L"} home={";
          for (auto const& h : occ.home) std::wcerr << h.full_label() << L" ";
          std::wcerr << L"} ectx={";
          for (auto const& e : occ.ectx)
            std::wcerr << e.first.full_label() << L" ";
          std::wcerr << L"}\n";
        }
      }

    // Stamp reduced-mode slots on the reducing parents: for each
    // (parent, reduced-mode) whose reduction node was united with an operand's
    // slice loop, record (mode, slot) on EVERY occurrence of that parent (the
    // value reduces the mode consistently across its occurrences). Deduplicated
    // so two operands reducing the same mode stamp it once.
    std::set<std::pair<std::size_t, std::wstring>> stamped;
    for (auto const& [par_vid, m] : reduction_stamps) {
      auto const key = std::make_pair(par_vid, std::wstring{m.full_label()});
      if (!stamped.insert(key).second) continue;
      auto const rit = root_slot.find(find(reduction_node(par_vid, m)));
      if (rit == root_slot.end()) continue;  // component never numbered
      for (OccurrenceRec& occ : out.cells[par_vid].occurrences)
        occ.reduced_slot.push_back({m, rit->second});
    }
    // Per-instance loop KIND and NESTING: every open names one physical loop
    // -- a Contracted open is the reduction node of (value, mode), an
    // External open is the carried position of the mode at the opening node
    // -- and that component's slot takes the open's kind. Two opens of one
    // instance must agree (an instance is either reduced or carried),
    // asserted. The same resolution names the loop instance of every ectx
    // entry (through its opener occurrence), which yields the DP's nesting
    // order as (outer, inner) constraints.
    auto const instance_slot =
        [&](std::size_t vid, OccurrenceRec const& occ, Index const& ix,
            BatchModeType kind) -> std::optional<std::pair<std::wstring, int>> {
      std::optional<std::size_t> root;
      if (kind == BatchModeType::Contracted) {
        root = find(reduction_node(vid, ix));
      } else {
        for (std::size_t pV = 0; pV < occ.carried.size(); ++pV)
          if (occ.carried[pV] == ix) {
            root = find(encode(vid, pV));
            break;
          }
      }
      if (!root) return std::nullopt;
      auto const rit = root_slot.find(*root);
      if (rit == root_slot.end()) return std::nullopt;  // never numbered
      return std::make_pair(std::wstring{ix.space().base_key()}, rit->second);
    };
    for (ValueCell const& c : out.cells)
      for (OccurrenceRec const& occ : c.occurrences)
        for (auto const& [ix, kind] : occ.opens) {
          auto const key = instance_slot(c.value_id, occ, ix, kind);
          if (!key) continue;
          auto const [kit, inserted] = out.loop_kind.emplace(*key, kind);
          SEQUANT_ASSERT(kit->second == kind &&
                         "compute_dag_boulevard: one loop instance opened "
                         "with two kinds (contracted and external)");
        }
    // Two EXTERNAL loops carry no data dependence between them (a value
    // sliced on both is produced per pair of batches either way, and the
    // escape placement is order-free), so trees may open them in either
    // order and the fused nest picks one: no constraint. A pair with a
    // contracted loop is a real constraint: a reduction completes inside the
    // loops that enclose it, so its consumers' loops must nest outside.
    auto const constrains = [&](std::pair<std::wstring, int> const& a,
                                std::pair<std::wstring, int> const& b) {
      auto const ka = out.loop_kind.find(a);
      auto const kb = out.loop_kind.find(b);
      bool const a_ext =
          ka != out.loop_kind.end() && ka->second == BatchModeType::External;
      bool const b_ext =
          kb != out.loop_kind.end() && kb->second == BatchModeType::External;
      return !(a_ext && b_ext);
    };
    for (ValueCell const& c : out.cells)
      for (OccurrenceRec const& occ : c.occurrences) {
        // Enclosing chain, outermost first, each entry resolved through its
        // opener occurrence (same tree, so the opened Index matches verbatim).
        container::svector<std::pair<std::wstring, int>> chain;
        for (std::size_t k = 0;
             k < occ.ectx.size() && k < occ.ectx_opener_point.size(); ++k) {
          auto const pit = point_occ.find(occ.ectx_opener_point[k]);
          auto const vit = point_value.find(occ.ectx_opener_point[k]);
          if (pit == point_occ.end() || vit == point_value.end()) continue;
          OccurrenceRec const& opener = *pit->second;
          for (auto const& [oix, okind] : opener.opens)
            if (oix == occ.ectx[k].first) {
              if (auto const key =
                      instance_slot(vit->second, opener, oix, okind))
                chain.push_back(*key);
              break;
            }
        }
        for (auto const& [ix, kind] : occ.opens)
          if (auto const key = instance_slot(c.value_id, occ, ix, kind))
            if (!chain.empty() && chain.back() != *key &&
                constrains(chain.back(), *key))
              out.loop_order.emplace(std::make_pair(chain.back(), *key),
                                     c.value_id);
        for (std::size_t k = 1; k < chain.size(); ++k)
          if (chain[k - 1] != chain[k] && constrains(chain[k - 1], chain[k]))
            out.loop_order.emplace(std::make_pair(chain[k - 1], chain[k]),
                                   c.value_id);
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
