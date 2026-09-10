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
  std::size_t hash;      //!< the value's NODE id (\c EvalExpr::hash_value(),
                         //!< the canonical colored graph): links a cell back
                         //!< to its forest nodes; batched-slot-BLIND.
  std::size_t key = 0;   //!< the VALUE id (\c value_key: node id combined
                         //!< with the home-sliced canonical positions,
                         //!< explicit-cells design section 11) -- the
                         //!< identity that folds occurrences into this cell.
                         //!< 0 (a hand-built cell) means "== hash"; read it
                         //!< through \c value_key_of(ValueCell const&).
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
/// The value id of \p c: its \c key, or its node id for a cell built without
/// one (every hand-built fixture; a value with nothing home-sliced has the
/// two equal anyway).
[[nodiscard]] inline std::size_t value_key_of(ValueCell const& c) noexcept {
  return c.key ? c.key : c.hash;
}

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
  using Data = typename Node::value_type;

  // Populate the forest path's residency meet (EvalExpr::sliced_modes) and,
  // for THIS engine, every occurrence's own home (EvalExpr::occurrence_home,
  // what home_scope / value_key_of read).
  stamp_lifetime_masks(forest);
  stamp_occurrence_homes(forest);

  // Per-occurrence record captured during the post-order walk. consumer_point
  // is the parent's point (its structural consumer); a root keeps its own
  // point (set at construction, overwritten by the parent if any).
  struct NodeRec {
    std::size_t hash;  // node id
    std::size_t key;   // value id (value_key_of(n))
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
    container::svector<std::size_t> child_recs;  // rec indices, left, right
    container::svector<int> loop_slot;           // per carried position
    container::svector<std::pair<Index, int>> reduced_slot;
    Node const* node = nullptr;  // the forest node, to stamp its value key
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
    r.key = n->hash_value();  // the FINAL value key is assigned below
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
    r.child_recs = child_recs;
    r.node = &n;
    std::size_t const idx = recs.size();
    recs.push_back(std::move(r));
    for (auto ci : child_recs) recs[ci].consumer_point = point;
    return idx;
  };

  for (auto const& tree : forest)
    visit(visit, tree, detail::BatchContext{}, {});

  // ---------------------------------------------------------------------
  // Loop identity FIRST, over OCCURRENCES (explicit-cells design section
  // 11, as amended): a loop instance is a connected component of
  // (occurrence, position) nodes joined by producer->consumer edges within a
  // tree and by conflict-aware FOLDS across trees; value identity is defined
  // AFTER the components are numbered -- node id + (position, loop slot) of
  // every home-sliced position + the operands' keys -- so one node sliced at
  // one position by two different loops in two terms (the residual's two
  // external loops, say) is two values, while occurrences one physical loop
  // does slice fold into one. The old order (value ids by hash, then loops
  // keyed by (value, position)) forced the two loops through the shared
  // node and either collapsed them or mis-stamped one family's slots.
  // ---------------------------------------------------------------------
  std::size_t const nrec = recs.size();
  std::unordered_map<std::size_t, std::size_t> rec_of_point;
  for (std::size_t i = 0; i < nrec; ++i) rec_of_point[recs[i].point] = i;

  // The GROUP key an occurrence folds under: node id + home-sliced canonical
  // positions. Occurrences of one group are candidates for one loop
  // instance per position; the fold below decides.
  auto const home_positions = [](NodeRec const& r) {
    container::svector<std::size_t> pos;
    for (Index const& m : r.home)
      for (std::size_t p = 0; p < r.carried.size(); ++p)
        if (r.carried[p] == m) {
          pos.push_back(p);
          break;
        }
    std::sort(pos.begin(), pos.end());
    return pos;
  };
  container::svector<std::size_t> group_key(nrec);
  for (std::size_t i = 0; i < nrec; ++i)
    group_key[i] = value_key(recs[i].hash, home_positions(recs[i]));

  std::size_t constexpr POS_BITS = 20;
  auto const encode = [](std::size_t idx, std::size_t pos) -> std::size_t {
    return (idx << POS_BITS) | pos;
  };
  auto const dec_idx = [](std::size_t n) { return n >> POS_BITS; };
  auto const dec_pos = [](std::size_t n) {
    return n & ((std::size_t{1} << POS_BITS) - 1);
  };
  // Reduction-mode nodes live in the UPPER half of an occurrence's position
  // space, one per (occurrence, reduced mode label in its own frame).
  std::size_t constexpr CONTRACTED_BASE = std::size_t{1} << (POS_BITS - 1);
  std::map<std::pair<std::size_t, std::wstring>, std::size_t> red_pos;
  std::size_t red_next = CONTRACTED_BASE;
  auto const reduction_node = [&](std::size_t idx,
                                  Index const& m) -> std::size_t {
    auto const key = std::make_pair(idx, std::wstring{m.full_label()});
    auto const it = red_pos.find(key);
    if (it != red_pos.end()) return encode(idx, it->second);
    std::size_t const pos = red_next++;
    SEQUANT_ASSERT(red_next < (std::size_t{1} << POS_BITS));
    red_pos[key] = pos;
    return encode(idx, pos);
  };
  // (occurrence, reduced mode) pairs to stamp once components are numbered.
  container::svector<std::pair<std::size_t, Index>> reduction_stamps;

  std::unordered_map<std::size_t, std::size_t> uf;  // node -> parent
  // CONFLICT-AWARE union-find: a component is one physical loop, and one
  // batch loop slices ONE mode of any array OCCURRENCE -- so a component must
  // never hold two distinct positions of one occurrence. That is the only
  // physical constraint: one loop may slice different positions of one NODE
  // in different occurrences (a term contracting an index that sits on
  // position 2 of one operand's subtree and position 0 of the other's, both
  // the same intermediate), and a term's two external loops may pair with
  // another term's in either assignment (a value the mirror term slices by
  // the other external loop folds into that loop; the terms' roots keep the
  // two loops apart). members[root]: occurrence index -> its single position.
  std::unordered_map<std::size_t, std::map<std::size_t, std::size_t>> members;
  // Components that must NEVER be united, keyed by root (symmetric): the
  // reduction loop of a value and every loop enclosing (or opened by) an
  // occurrence that reads that value COMPLETE. Were they one physical loop,
  // the reader would sit inside the loop the value sums over and read a
  // partial sum -- the illegal fusion the pass split and rule 4 used to
  // paper over by re-running the nest and materializing every loop-local
  // value whole. Seeded after the within-tree edges, enforced by every
  // fold: the reader's loop stays its own instance, in its own nest,
  // sequenced after the reduction's, and its sliced operands are its own
  // productions (the recompute the forest performs and the per-term
  // optimizer costs).
  std::unordered_map<std::size_t, std::set<std::size_t>> forbid;
  // The KIND of a component's loop (Contracted: a batched reduction at its
  // opener; External: a carried mode's physical loop), keyed by root: two
  // components of different kinds never unite -- a term's external occupied
  // loop and another term's contracted occupied loop may slice one mode of a
  // shared value, but they are two loops (one scatters into the result, the
  // other sums), and the builder realizes a block as one kind.
  std::unordered_map<std::size_t, BatchModeType> comp_kind;
  auto find = [&](std::size_t x) -> std::size_t {
    auto it = uf.find(x);
    if (it == uf.end()) {
      uf.emplace(x, x);
      members[x][dec_idx(x)] = dec_pos(x);
      return x;
    }
    std::size_t root = x;
    while (uf[root] != root) root = uf[root];
    while (uf[x] != root) {
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
    for (auto const& [o, pos] : ma) {
      auto const jt = mb.find(o);
      if (jt != mb.end() && jt->second != pos) return false;  // conflict
    }
    if (auto const fa = forbid.find(ra);
        fa != forbid.end() && fa->second.count(rb))
      return false;  // a reader's loop and the reduction it reads complete
    if (auto const ka = comp_kind.find(ra), kb = comp_kind.find(rb);
        ka != comp_kind.end() && kb != comp_kind.end() &&
        ka->second != kb->second)
      return false;  // an external loop and a contracted loop stay distinct
    for (auto const& [o, pos] : ma) mb[o] = pos;
    members.erase(ra);
    uf[ra] = rb;
    if (auto const ka = comp_kind.find(ra); ka != comp_kind.end()) {
      comp_kind[rb] = ka->second;
      comp_kind.erase(ka);
    }
    if (auto const fa = forbid.find(ra); fa != forbid.end()) {
      auto moved = std::move(fa->second);
      forbid.erase(fa);
      for (std::size_t x : moved) {
        forbid[x].erase(ra);
        forbid[x].insert(rb);
        forbid[rb].insert(x);
      }
    }
    return true;
  };
  auto const is_batched = [](NodeRec const& r, Index const& m) -> bool {
    for (auto const& h : r.home)
      if (h == m) return true;
    return false;
  };
  bool const conflict_dump = std::getenv("SEQUANT_DUMP_LOOP_SLOT") != nullptr;

  // Edges within a tree: a home-sliced carried position of an occurrence to
  // the SAME physical mode at its parent occurrence (label match in the
  // parent's frame -- same tree, labels agree), or, where the parent reduces
  // the mode, to the parent's reduction node for it.
  for (std::size_t i = 0; i < nrec; ++i) {
    NodeRec const& r = recs[i];
    NodeRec const* par = nullptr;
    std::size_t par_idx = 0;
    if (r.consumer_point != r.point) {
      auto const pit = rec_of_point.find(r.consumer_point);
      if (pit != rec_of_point.end()) {
        par_idx = pit->second;
        par = &recs[par_idx];
      }
    }
    for (std::size_t pV = 0; pV < r.carried.size(); ++pV) {
      Index const& m = r.carried[pV];
      if (!is_batched(r, m)) continue;
      (void)find(encode(i, pV));
      if (!par) continue;
      auto const pj = std::find(par->carried.begin(), par->carried.end(), m);
      if (pj == par->carried.end()) {
        std::size_t const rn = reduction_node(par_idx, m);
        if (try_unite(encode(i, pV), rn)) {
          reduction_stamps.push_back({par_idx, m});
          if (conflict_dump)
            std::wcerr << L"[loop_slot] edge occ" << i << L"@" << pV << L"("
                       << m.full_label() << L") ~ red(occ" << par_idx << L","
                       << m.full_label() << L")\n";
        }
        continue;
      }
      std::size_t const pC =
          static_cast<std::size_t>(pj - par->carried.begin());
      bool const united = try_unite(encode(i, pV), encode(par_idx, pC));
      if (conflict_dump)
        std::wcerr << L"[loop_slot] " << (united ? L"edge" : L"REJECTED edge")
                   << L" occ" << i << L"@" << pV << L"(" << m.full_label()
                   << L") ~ occ" << par_idx << L"@" << pC << L"\n";
    }
  }
  // Every mode an occurrence contracts in batches at its own node owns a
  // loop identity even when every operand is an input.
  for (std::size_t i = 0; i < nrec; ++i)
    for (Index const& m : recs[i].contracted_batched) {
      (void)find(reduction_node(i, m));
      reduction_stamps.push_back({i, m});
    }

  // Component kinds (see `comp_kind`), from every open, before any fold.
  for (std::size_t i = 0; i < nrec; ++i)
    for (auto const& [ix, kind] : recs[i].opens) {
      std::optional<std::size_t> node;
      if (kind == BatchModeType::Contracted) {
        node = reduction_node(i, ix);
      } else {
        for (std::size_t pV = 0; pV < recs[i].carried.size(); ++pV)
          if (recs[i].carried[pV] == ix) {
            node = encode(i, pV);
            break;
          }
      }
      if (!node) continue;
      std::size_t const r = find(*node);
      auto const it = comp_kind.find(r);
      SEQUANT_ASSERT((it == comp_kind.end() || it->second == kind) &&
                     "compute_dag_boulevard: one loop instance opened with "
                     "two kinds within a tree");
      comp_kind[r] = kind;
    }

  // Reader-versus-reduction constraints (see `forbid`): for every occurrence
  // O and every operand occurrence c that reduces modes in batches, c's
  // reduction components must stay apart from every loop instance enclosing
  // O (through its openers) and every loop O opens itself.
  {
    auto const instance_node_of =
        [&](std::size_t opener_idx,
            Index const& ix) -> std::optional<std::size_t> {
      NodeRec const& op = recs[opener_idx];
      for (auto const& [oix, okind] : op.opens)
        if (oix == ix) {
          if (okind == BatchModeType::Contracted)
            return reduction_node(opener_idx, ix);
          for (std::size_t pV = 0; pV < op.carried.size(); ++pV)
            if (op.carried[pV] == ix) return encode(opener_idx, pV);
          return std::nullopt;
        }
      return std::nullopt;
    };
    for (std::size_t i = 0; i < nrec; ++i) {
      NodeRec const& o = recs[i];
      container::svector<std::size_t> reader_loops;
      for (std::size_t k = 0; k < o.ectx.size() && k < o.ectx_opener.size();
           ++k) {
        auto const pit = pre_to_point.find(o.ectx_opener[k]);
        if (pit == pre_to_point.end()) continue;
        auto const rit = rec_of_point.find(pit->second);
        if (rit == rec_of_point.end()) continue;
        if (auto const nd = instance_node_of(rit->second, o.ectx[k].first))
          reader_loops.push_back(*nd);
      }
      for (auto const& [ix, kind] : o.opens)
        if (auto const nd = instance_node_of(i, ix))
          reader_loops.push_back(*nd);
      if (reader_loops.empty()) continue;
      for (std::size_t ci : o.child_recs)
        for (Index const& m : recs[ci].contracted_batched) {
          std::size_t const r = find(reduction_node(ci, m));
          for (std::size_t ln : reader_loops) {
            std::size_t const x = find(ln);
            if (x == r) continue;
            forbid[r].insert(x);
            forbid[x].insert(r);
          }
        }
    }
  }

  // FOLD across trees: occurrences of one group are one value where one
  // physical loop slices them. Attempt the union position by position (and
  // reduction by reduction, by index in the node's contracted order, which
  // is canonical); a rejected fold (the two trees' loops would collapse two
  // distinct positions of some group) leaves the occurrences in different
  // instances, and the final key below tells them apart.
  {
    std::unordered_map<std::size_t, std::size_t> first_of_group;
    for (std::size_t i = 0; i < nrec; ++i) {
      auto const [fit, inserted] = first_of_group.emplace(group_key[i], i);
      if (inserted) continue;
      std::size_t const f = fit->second;
      NodeRec const& rf = recs[f];
      NodeRec const& ro = recs[i];
      for (std::size_t p = 0; p < ro.carried.size() && p < rf.carried.size();
           ++p) {
        if (!is_batched(ro, ro.carried[p]) || !is_batched(rf, rf.carried[p]))
          continue;
        bool const united = try_unite(encode(f, p), encode(i, p));
        if (conflict_dump && !united)
          std::wcerr << L"[loop_slot] REJECTED fold occ" << f << L"@" << p
                     << L" ~ occ" << i << L"@" << p << L" (hash "
                     << (ro.hash % 100000u) << L")\n";
      }
      for (std::size_t j = 0;
           j < ro.contracted_batched.size() && j < rf.contracted_batched.size();
           ++j)
        (void)try_unite(reduction_node(f, rf.contracted_batched[j]),
                        reduction_node(i, ro.contracted_batched[j]));
    }
  }

  // Number the components: one loop_slot per component, per SPACE in
  // first-seen order.
  std::unordered_map<std::size_t, int> root_slot;
  std::map<std::wstring, int> next_slot;
  for (std::size_t i = 0; i < nrec; ++i)
    for (std::size_t pV = 0; pV < recs[i].carried.size(); ++pV) {
      if (!is_batched(recs[i], recs[i].carried[pV])) continue;
      std::size_t const root = find(encode(i, pV));
      if (root_slot.find(root) != root_slot.end()) continue;
      std::wstring const sp{recs[i].carried[pV].space().base_key()};
      root_slot.emplace(root, next_slot[sp]++);
    }
  for (auto const& [idx, m] : reduction_stamps) {
    std::size_t const root = find(reduction_node(idx, m));
    if (root_slot.find(root) != root_slot.end()) continue;
    std::wstring const sp{m.space().base_key()};
    root_slot.emplace(root, next_slot[sp]++);
  }

  // Stamp each occurrence's per-position loop_slot and reduced_slot.
  for (std::size_t i = 0; i < nrec; ++i) {
    NodeRec& r = recs[i];
    r.loop_slot.assign(r.carried.size(), -1);
    for (std::size_t pV = 0; pV < r.carried.size(); ++pV)
      if (is_batched(r, r.carried[pV]))
        r.loop_slot[pV] = root_slot.at(find(encode(i, pV)));
  }
  {
    std::set<std::pair<std::size_t, std::wstring>> stamped;
    for (auto const& [idx, m] : reduction_stamps) {
      if (!stamped.insert({idx, std::wstring{m.full_label()}}).second) continue;
      auto const rit = root_slot.find(find(reduction_node(idx, m)));
      if (rit == root_slot.end()) continue;
      recs[idx].reduced_slot.push_back({m, rit->second});
    }
  }

  // FINAL value key, bottom-up (recs are in post-order: operands precede
  // their consumer): node id combined with (position, slot) of every
  // home-sliced position, (index, slot) of every mode reduced in batches,
  // and the operands' keys -- the node id alone when nothing below is
  // sliced. Stamped on the node so value_key_of(node) agrees everywhere.
  container::svector<std::size_t> final_key(nrec);
  for (std::size_t i = 0; i < nrec; ++i) {
    NodeRec const& r = recs[i];
    bool sliced = false;
    std::size_t h = r.hash;
    hash::combine(h, std::size_t{0x5eed});
    for (std::size_t pV = 0; pV < r.carried.size(); ++pV)
      if (r.loop_slot[pV] >= 0) {
        sliced = true;
        hash::combine(h, pV);
        hash::combine(h, static_cast<std::size_t>(r.loop_slot[pV]));
      }
    for (std::size_t j = 0; j < r.reduced_slot.size(); ++j) {
      sliced = true;
      hash::combine(h, CONTRACTED_BASE + j);
      hash::combine(h, static_cast<std::size_t>(r.reduced_slot[j].second));
    }
    for (std::size_t ci : r.child_recs) {
      if (final_key[ci] != recs[ci].hash) sliced = true;
      hash::combine(h, final_key[ci]);
    }
    final_key[i] = sliced ? h : r.hash;
    recs[i].key = final_key[i];
    if (r.node) const_cast<Data&>(**r.node).set_value_key(final_key[i]);
  }

  // Cross-occurrence union, per VALUE, of the modes that value EVER realizes
  // as its OWN loop (\c node_slice_mask() at the value's own node, not an
  // ancestor's): a mode a value realizes as its OWN loop slices that value's
  // OPERANDS on the way down, never the value's own result, so it is
  // excluded from the cell's footprint home.
  std::unordered_map<std::size_t, container::svector<Index>> own_modes_union;
  for (auto const& r : recs) {
    auto& acc = own_modes_union[r.key];
    for (auto const& m : r.own_modes)
      if (std::find(acc.begin(), acc.end(), m) == acc.end()) acc.push_back(m);
  }

  // Per-value RELABELED modes: carried by SOME occurrences but not all -- the
  // UNION minus the INTERSECTION of the occurrences' canon_indices.
  std::unordered_map<std::size_t, container::svector<Index>> carried_union,
      carried_isect;
  std::unordered_map<std::size_t, bool> carried_seeded;
  for (auto const& r : recs) {
    auto& u = carried_union[r.key];
    for (auto const& m : r.carried)
      if (std::find(u.begin(), u.end(), m) == u.end()) u.push_back(m);
    auto& is = carried_isect[r.key];
    if (!carried_seeded[r.key]) {
      carried_seeded[r.key] = true;
      is.assign(r.carried.begin(), r.carried.end());
    } else {
      container::svector<Index> keep;
      for (auto const& m : is)
        if (std::find(r.carried.begin(), r.carried.end(), m) != r.carried.end())
          keep.push_back(m);
      is = std::move(keep);
    }
  }

  // Group occurrences by value identity (the final key): one ValueCell per
  // group.
  RichSchedule out;
  out.num_points = counter;
  std::unordered_map<std::size_t, std::size_t> hash_to_cell;
  for (auto const& r : recs) {
    auto fold_enclosing = [&](container::svector<Index>& enclosing_modes) {
      for (auto const& e : r.ectx)
        if (std::find(enclosing_modes.begin(), enclosing_modes.end(),
                      e.first) == enclosing_modes.end())
          enclosing_modes.push_back(e.first);
    };
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
      o.loop_slot = r.loop_slot;
      o.reduced_slot = r.reduced_slot;
      return o;
    };
    auto const it = hash_to_cell.find(r.key);
    if (it == hash_to_cell.end()) {
      ValueCell c;
      c.value_id = out.cells.size();
      c.is_leaf = r.is_leaf;
      c.hash = r.hash;  // the node id
      c.key = r.key;    // the value id (== the hash_to_cell key)
      c.first_use = r.point;
      c.last_use = r.consumer_point;
      c.home_depth = detail::home_depth_of(r.home, r.ectx);
      auto const& self_modes = own_modes_union[r.key];
      container::svector<Index> home_modes;
      for (auto const& m : r.home)
        if (std::find(self_modes.begin(), self_modes.end(), m) ==
            self_modes.end())
          home_modes.push_back(m);
      c.carried = r.carried;
      c.home_modes = std::move(home_modes);
      {
        auto const& u = carried_union[r.key];
        auto const& is = carried_isect[r.key];
        for (auto const& m : u)
          if (std::find(is.begin(), is.end(), m) == is.end())
            c.divergent_modes.push_back(m);
      }
      fold_enclosing(c.enclosing_modes);
      c.occurrences.push_back(make_occ());
      hash_to_cell.emplace(r.key, c.value_id);
      out.cells.push_back(std::move(c));
    } else {
      ValueCell& c = out.cells[it->second];
      c.first_use = std::min(c.first_use, r.point);
      c.last_use = std::max(c.last_use, r.consumer_point);
      fold_enclosing(c.enclosing_modes);
      c.occurrences.push_back(make_occ());
    }
  }

  // Component-membership dump: per (space, slot), the member values
  // (value:hash(label@position)). Guarded by SEQUANT_DUMP_LOOP_SLOT.
  if (conflict_dump) {
    std::map<std::pair<std::wstring, int>, std::map<std::wstring, int>> comp;
    for (ValueCell const& c : out.cells)
      for (OccurrenceRec const& occ : c.occurrences)
        for (std::size_t pV = 0; pV < occ.carried.size(); ++pV) {
          if (pV >= occ.loop_slot.size() || occ.loop_slot[pV] < 0) continue;
          std::wstring const sp{occ.carried[pV].space().base_key()};
          std::wstring const mem = L"v" + std::to_wstring(c.value_id) + L":" +
                                   std::to_wstring(c.hash % 100000u) + L"(" +
                                   std::wstring(occ.carried[pV].full_label()) +
                                   L"@" + std::to_wstring(pV) + L")";
          comp[{sp, occ.loop_slot[pV]}][mem] = 1;
        }
    for (auto const& [key, mems] : comp) {
      std::wcerr << L"[comp] " << key.first << L"#slot" << key.second
                 << L" members={ ";
      for (auto const& [mem, _] : mems) std::wcerr << mem << L" ";
      std::wcerr << L"}\n";
    }
  }

  // Per-instance loop KIND and NESTING (unchanged in substance): every open
  // names one physical loop -- a Contracted open is the reduction node of
  // (occurrence, mode), an External open the carried position of the mode at
  // the opening occurrence -- and that component's slot takes the open's
  // kind. The same resolution names the loop instance of every ectx entry
  // (through its opener occurrence), which yields the DP's nesting order.
  {
    std::unordered_map<std::size_t, OccurrenceRec const*> point_occ;
    for (ValueCell const& c : out.cells)
      for (OccurrenceRec const& o : c.occurrences) point_occ[o.point] = &o;
    auto const instance_slot =
        [&](OccurrenceRec const& occ, Index const& ix,
            BatchModeType kind) -> std::optional<std::pair<std::wstring, int>> {
      auto const rit_idx = rec_of_point.find(occ.point);
      if (rit_idx == rec_of_point.end()) return std::nullopt;
      std::size_t const idx = rit_idx->second;
      std::optional<std::size_t> root;
      if (kind == BatchModeType::Contracted) {
        root = find(reduction_node(idx, ix));
      } else {
        for (std::size_t pV = 0; pV < occ.carried.size(); ++pV)
          if (occ.carried[pV] == ix) {
            root = find(encode(idx, pV));
            break;
          }
      }
      if (!root) return std::nullopt;
      auto const rit = root_slot.find(*root);
      if (rit == root_slot.end()) return std::nullopt;
      return std::make_pair(std::wstring{ix.space().base_key()}, rit->second);
    };
    for (ValueCell const& c : out.cells)
      for (OccurrenceRec const& occ : c.occurrences)
        for (auto const& [ix, kind] : occ.opens) {
          auto const key = instance_slot(occ, ix, kind);
          if (!key) continue;
          auto const [kit, inserted] = out.loop_kind.emplace(*key, kind);
          SEQUANT_ASSERT(kit->second == kind &&
                         "compute_dag_boulevard: one loop instance opened "
                         "with two kinds (contracted and external)");
        }
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
        container::svector<std::pair<std::wstring, int>> chain;
        container::svector<std::size_t> chain_opener;
        for (std::size_t k = 0;
             k < occ.ectx.size() && k < occ.ectx_opener_point.size(); ++k) {
          auto const pit = point_occ.find(occ.ectx_opener_point[k]);
          if (pit == point_occ.end()) continue;
          OccurrenceRec const& opener = *pit->second;
          for (auto const& [oix, okind] : opener.opens)
            if (oix == occ.ectx[k].first) {
              if (auto const key = instance_slot(opener, oix, okind)) {
                chain.push_back(*key);
                chain_opener.push_back(occ.ectx_opener_point[k]);
              }
              break;
            }
        }
        for (auto const& [ix, kind] : occ.opens)
          if (auto const key = instance_slot(occ, ix, kind))
            if (!chain.empty() && chain.back() != *key &&
                constrains(chain.back(), *key))
              out.loop_order.emplace(std::make_pair(chain.back(), *key),
                                     c.value_id);
        for (std::size_t k = 1; k < chain.size(); ++k) {
          if (chain[k - 1] == chain[k]) continue;
          if (chain_opener[k - 1] == chain_opener[k]) continue;
          if (!constrains(chain[k - 1], chain[k])) continue;
          out.loop_order.emplace(std::make_pair(chain[k - 1], chain[k]),
                                 c.value_id);
        }
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
