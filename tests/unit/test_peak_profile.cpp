// Phase 3b T1: tests for the static peak-profile sizing primitives
// (SeQuant/core/eval/peak_profile.hpp). Two free functions are pinned here:
//   - home_depth_of: resolve a residency mode-set to an enclosing-batch-
//     context loop depth (mirrors the runtime rl-walk at eval.hpp:1776-1782).
//   - cell_footprint: home-relative footprint of a carried-index set via the
//     existing dryrun::CostModel::memsize (no sizing logic reimplemented).
// Neither function has a production caller yet (T2/T3 wire them in); this
// task only pins the two primitives' contracts.

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/backends/dryrun/meter.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <set>
#include <string_view>
#include <vector>

namespace {

using sequant::BatchModeType;
using sequant::BatchPolicy;
using sequant::binarize;
using sequant::deserialize;
using sequant::EvalExpr;
using sequant::EvalNode;
using sequant::ExprPtr;
using sequant::Index;
using sequant::container::svector;
using sequant::eval::Cell;
using sequant::eval::compute_dag_boulevard;
using sequant::eval::peak_profile_replay;
using sequant::eval::peak_profile_sweep;
using sequant::eval::PeakProfile;
using sequant::eval::RichSchedule;
using sequant::eval::Schedule;
using sequant::eval::ValueCell;
using sequant::eval::detail::BatchContext;
using sequant::eval::detail::cell_footprint;
using sequant::eval::detail::home_depth_of;
using sequant::eval::dryrun::CacheConfig;
using sequant::eval::dryrun::CostModel;
using sequant::eval::dryrun::EvalExprDryRun;
using sequant::eval::dryrun::EvalNodeDryRun;
using sequant::eval::dryrun::ExtentOverrides;
using sequant::eval::dryrun::SizeRegime;

// Build an EvalExpr from a single-tensor spec (e.g. "R{i_1;a_5}"); its
// canon_indices are exactly the tensor's bra+ket slots. Mirrors the helper in
// test_lifetime_mask.cpp.
EvalExpr eval_tensor(std::string_view tensor) {
  auto expr = sequant::deserialize<ExprPtr>(std::string(tensor));
  REQUIRE(static_cast<bool>(expr));
  return EvalExpr{expr->as<sequant::Tensor>()};
}

// A leaf eval node carrying the given tensor's slots.
EvalNode<EvalExpr> leaf(std::string_view tensor) {
  return EvalNode<EvalExpr>{eval_tensor(tensor)};
}

// An internal eval node whose OWN result slots are the given tensor's slots,
// with the two supplied child subtrees.
EvalNode<EvalExpr> inode(std::string_view result, EvalNode<EvalExpr> l,
                         EvalNode<EvalExpr> r) {
  return EvalNode<EvalExpr>{eval_tensor(result), std::move(l), std::move(r)};
}

// Stamp a single External batch loop mode at a node.
void stamp_ext(EvalNode<EvalExpr>& n, Index ix) {
  // Realizes an External loop AT n: stamp both the per-node sliced mask AND the
  // loop-OPEN annotation the (opens-based) enclosing-context walk
  // (OccurrenceRec::ectx) reads to reconstruct the loop nest.
  n->set_node_slice_mask({{ix, BatchModeType::External}});
  n->set_batch_loops_opened_here({{std::move(ix), BatchModeType::External}});
}

// Local stand-in for the retired \c eval::compute_dag_path, which had no
// production caller and was deleted with the placement-router cleanup: the
// same thin PROJECTION of \c compute_dag_boulevard onto flat \c Cell s that
// it performed -- one Cell per ValueCell, its carried/home_modes collapsed to
// a single footprint by \c detail::cell_footprint. Kept HERE (not restored to
// the header) so the sweep/replay equalities below still run against REAL
// linearized forests -- \c compute_dag_boulevard, \c cell_footprint, \c
// peak_profile_sweep and \c peak_profile_replay are all live production code
// -- instead of only against hand-built Schedules.
template <typename R, typename BlockOfFn>
Schedule dag_path(R const& forest, CostModel const& cm,
                  BlockOfFn const& block_of) {
  RichSchedule const rich = compute_dag_boulevard(forest, cm, block_of);

  Schedule out;
  out.num_points = rich.num_points;
  out.cells.reserve(rich.cells.size());
  for (auto const& vc : rich.cells) {
    Cell c;
    c.value_id = vc.value_id;
    c.home_depth = vc.home_depth;
    c.footprint = cell_footprint(vc.carried, vc.home_modes, cm, block_of,
                                 vc.divergent_modes);
    c.first_use = vc.first_use;
    c.last_use = vc.last_use;
    out.cells.push_back(c);
  }
  return out;
}

}  // namespace

TEST_CASE("home_depth_of resolves the deepest enclosing loop in home_modes",
          "[peak_profile]") {
  // Three enclosing loops, outermost-first (matches BatchContext's documented
  // order): o (level 0), i (level 1), a (level 2, innermost).
  Index const o1{L"o_1"}, i1{L"i_1"}, a1{L"a_1"}, x1{L"x_1"};
  BatchContext const ectx{
      {o1, {0, 1}},
      {i1, {0, 1}},
      {a1, {0, 1}},
  };

  // Deepest (innermost) match.
  CHECK(home_depth_of(svector<Index>{a1}, ectx) == 2);
  // Outermost match.
  CHECK(home_depth_of(svector<Index>{o1}, ectx) == 0);
  // No match at all -> chain root.
  CHECK(home_depth_of(svector<Index>{}, ectx) == -1);
  // Mixed set: one member present (i1, level 1), one absent (x1) -> the
  // present member's level, not disturbed by the absent one.
  CHECK(home_depth_of(svector<Index>{i1, x1}, ectx) == 1);
  // Two members present at DIFFERENT levels (o1 @ 0, a1 @ 2) -> the DEEPEST
  // (innermost) match wins, proving the innermost-to-outermost scan returns
  // the first (deepest) hit rather than the shallowest.
  CHECK(home_depth_of(svector<Index>{o1, a1}, ectx) == 2);
}

TEST_CASE("cell_footprint sizes enclosing-home modes at BLOCK, others FULL",
          "[peak_profile]") {
  // full(i) = full(a) = 10; block(i) = 2 (a fake per-mode slice extent).
  SizeRegime r;
  r.space_extent = {{L"i", 10}, {L"a", 10}};
  CostModel const cm{r};

  Index const p{L"i_1"};  // the mode with an enclosing loop
  Index const q{L"a_1"};  // unbatched: no enclosing loop, no override
  svector<Index> const carried{p, q};

  auto const block_of = [](Index const& ix) -> std::size_t {
    return ix == Index{L"i_1"} ? 2 : 10;
  };

  // Cell homed AT p (home_modes = {p}) -> p is in the meet -> p sized at
  // BLOCK (2), q (not in the meet) sized at FULL (10).
  {
    svector<Index> const home_modes{p};
    auto const got = cell_footprint(carried, home_modes, cm, block_of);
    ExtentOverrides ov;
    ov[0] = 2;  // mode 0 (p) at block extent
    CHECK(got == cm.memsize(carried, ov));
    // memsize() reports BYTES (elems * CostModel's numeric_size_ = 8.0), not
    // a raw element count: 2 (block p) * 10 (full q) * 8 bytes/elem.
    CHECK(got == 160);
  }

  // Same cell homed ABOVE p's loop (home_modes = {}) -> neither mode is in
  // the meet -> both modes sized FULL.
  {
    svector<Index> const home_modes{};
    auto const got = cell_footprint(carried, home_modes, cm, block_of);
    CHECK(got == cm.memsize(carried));
    // 10 (full p) * 10 (full q) * 8 bytes/elem.
    CHECK(got == 800);
  }
}

TEST_CASE(
    "cell_footprint no longer doubles a divergent-mode home (the 2x fudge is "
    "gone; divergent_modes is informational)",
    "[peak_profile]") {
  // The flat 2x pricing fudge (placeholder 46b495eba) is DELETED: a divergent
  // value homed at a sub-scope is UN-FOLDED into two real, non-divergent
  // ValueCells, each priced ONCE here at its own home; peak co-residency is
  // priced structurally by
  // peak_profile_sweep over the two cells' liveness intervals, and the
  // replication recompute is a separate report-only term. So cell_footprint
  // IGNORES divergent_modes -- passing it never changes the size.
  SizeRegime r;
  r.space_extent = {{L"i", 10}, {L"a", 10}};
  CostModel const cm{r};
  Index const p{L"i_1"}, q{L"a_1"};
  svector<Index> const carried{p, q};
  auto const block_of = [](Index const& ix) -> std::size_t {
    return ix == Index{L"i_1"} ? 2 : 10;
  };
  svector<Index> const home_modes{p};

  auto const shared = cell_footprint(carried, home_modes, cm, block_of, {});
  CHECK(shared == 160);  // one sliced copy: block(p)=2 * full(q)=10 * 8
  // A divergent mode that IS sliced (in the home) no longer doubles.
  CHECK(cell_footprint(carried, home_modes, cm, block_of, svector<Index>{p}) ==
        shared);
  // A divergent mode that is NOT sliced (not in the home) also does not change.
  CHECK(cell_footprint(carried, home_modes, cm, block_of, svector<Index>{q}) ==
        shared);
}

// ---------------------------------------------------------------------
// T2: the interval-event sweep over a hand-built Schedule (bypasses
// compute_dag_path entirely; Cells are constructed directly).
// ---------------------------------------------------------------------

TEST_CASE("peak_profile_sweep finds the peak, its point, and the live set",
          "[peak_profile]") {
  // Three overlapping lifetimes on a 5-point timeline:
  //   A = [0,3] fp=100, B = [1,2] fp=40, C = [2,4] fp=10.
  // Point 2 is the sole point where all three are live => 100+40+10 = 150.
  Schedule s;
  s.num_points = 5;
  s.cells.push_back(Cell{/*value_id=*/0, /*home_depth=*/-1, /*footprint=*/100,
                         /*first_use=*/0, /*last_use=*/3});
  s.cells.push_back(Cell{/*value_id=*/1, /*home_depth=*/-1, /*footprint=*/40,
                         /*first_use=*/1, /*last_use=*/2});
  s.cells.push_back(Cell{/*value_id=*/2, /*home_depth=*/-1, /*footprint=*/10,
                         /*first_use=*/2, /*last_use=*/4});

  auto const p = peak_profile_sweep(s);
  CHECK(p.peak_bytes == 150.0);
  CHECK(p.binding_point == 2);
  REQUIRE(p.live_at_binding.size() == 3);
  CHECK(p.live_at_binding[0] == 0);
  CHECK(p.live_at_binding[1] == 1);
  CHECK(p.live_at_binding[2] == 2);
}

TEST_CASE("peak_profile_sweep breaks a peak tie toward the LOWER point",
          "[peak_profile]") {
  // Two disjoint plateaus of EQUAL height 100 on a 4-point timeline:
  //   A = [0,1] fp=100 (live at points 0,1), B = [2,3] fp=100 (points 2,3).
  // The running peak is 100 at point 0 and again at points 2,3; the strict
  // `>` comparison keeps the FIRST (lowest) argmax.
  Schedule s;
  s.num_points = 4;
  s.cells.push_back(Cell{0, -1, 100, 0, 1});
  s.cells.push_back(Cell{1, -1, 100, 2, 3});

  auto const p = peak_profile_sweep(s);
  CHECK(p.peak_bytes == 100.0);
  CHECK(p.binding_point == 0);  // lowest-point tie-break
  // Only A is live at point 0.
  REQUIRE(p.live_at_binding.size() == 1);
  CHECK(p.live_at_binding[0] == 0);
}

TEST_CASE("compute_dag_boulevard threads the value hash onto each ValueCell",
          "[peak_profile]") {
  // The rich cell must carry the value's hash_value() -- the CSE identity
  // that links a cell back to its forest nodes (the cell table's own
  // value-to-node lookups key by it). Same CSE forest as above: the shared V
  // folds into ONE cell whose hash must equal the shared node's
  // hash_value().
  auto make_V = [] { return leaf("V{i_1;i_2}"); };
  auto const v_node = make_V();
  std::size_t const v_hash = v_node->hash_value();

  auto P1 = inode("P{i_1;a_3}", make_V(), leaf("W{a_1;a_3}"));
  auto P2 = inode("P{i_1;a_3}", make_V(), leaf("W{a_1;a_3}"));

  SizeRegime r;
  r.space_extent = {{L"i", 5}, {L"a", 10}};
  CostModel const cm{r};
  auto const block_of = [](Index const&) -> std::size_t { return 1; };
  std::vector<EvalNode<EvalExpr>> forest{P1, P2};

  RichSchedule const rich = compute_dag_boulevard(forest, cm, block_of);
  REQUIRE(rich.cells.size() == 3);  // V, W, P fold pairwise

  // EXACTLY one cell carries V's hash, and it IS the V value (carried
  // {i_1,i_2}).
  std::size_t v_cells = 0;
  ValueCell const* v = nullptr;
  for (auto const& c : rich.cells)
    if (c.hash == v_hash) {
      ++v_cells;
      v = &c;
    }
  REQUIRE(v_cells == 1);
  REQUIRE(v != nullptr);
  std::set<Index> const carried_set(v->carried.begin(), v->carried.end());
  std::set<Index> const v_slots(v_node->canon_indices().begin(),
                                v_node->canon_indices().end());
  CHECK(carried_set == v_slots);

  // Every cell's hash is a real value hash present in the forest (no cell has a
  // zero/garbage hash), and the three cells' hashes are pairwise distinct.
  std::set<std::size_t> hashes;
  for (auto const& c : rich.cells) {
    CHECK(c.hash != 0);
    hashes.insert(c.hash);
  }
  CHECK(hashes.size() == 3);
}

// =====================================================================
// T3: the independent REPLAY ORACLE (peak_profile_replay) -- oracle ==
// sweep on hand-built Schedules (design section 9.6).
// =====================================================================

TEST_CASE("peak_profile_replay agrees with the sweep on a hand-built Schedule",
          "[peak_profile]") {
  // The same three-lifetime Schedule the T2 sweep test pins:
  //   A = [0,3] fp=100, B = [1,2] fp=40, C = [2,4] fp=10, 5 points => 150.
  // peak_profile_replay is a DELIBERATELY DIFFERENT algorithm (explicit
  // per-point live-set sum, not the +delta/-delta difference array), so an
  // exact match is a real cross-check of the sweep's interval logic.
  {
    Schedule s;
    s.num_points = 5;
    s.cells.push_back(Cell{0, -1, 100, 0, 3});
    s.cells.push_back(Cell{1, -1, 40, 1, 2});
    s.cells.push_back(Cell{2, -1, 10, 2, 4});
    CHECK(peak_profile_replay(s) == 150.0);
    CHECK(peak_profile_replay(s) == peak_profile_sweep(s).peak_bytes);
  }

  // The T2 tie-break Schedule: two disjoint equal-height plateaus.
  {
    Schedule s;
    s.num_points = 4;
    s.cells.push_back(Cell{0, -1, 100, 0, 1});
    s.cells.push_back(Cell{1, -1, 100, 2, 3});
    CHECK(peak_profile_replay(s) == 100.0);
    CHECK(peak_profile_replay(s) == peak_profile_sweep(s).peak_bytes);
  }
}

// ---- Step B: oracle == sweep on real linearized forests -------------
//
// These run on forests linearized by the LIVE compute_dag_boulevard (via the
// local `dag_path` projection above), not on hand-built Schedules, so they
// cover cell_footprint's real inputs and the sweep's real interval shapes.

TEST_CASE(
    "peak_profile_replay == sweep on linearized forests (incl. a CSE fold)",
    "[peak_profile]") {
  SizeRegime r;
  r.space_extent = {{L"i", 5}, {L"a", 10}};
  CostModel const cm{r};
  auto const block_of = [](Index const&) -> std::size_t { return 1; };

  // (1) A single contraction of two leaves: A*B -> R.
  {
    auto R = inode("R{i_1;a_5}", leaf("A{i_1;a_3}"), leaf("B{a_3;a_5}"));
    std::vector<EvalNode<EvalExpr>> forest{R};
    auto const s = dag_path(forest, cm, block_of);
    CHECK(peak_profile_replay(s) == peak_profile_sweep(s).peak_bytes);
  }

  // (2) The CSE forest (a shared V consumed by two parents): the fold widens
  // V's interval across occurrences, exercising min/max over the grouped
  // cell -- the sweep and the replay must still agree exactly.
  {
    auto make_V = [] { return leaf("V{i_1;i_2}"); };
    auto P1 = inode("P{i_1;a_3}", make_V(), leaf("W{a_1;a_3}"));
    auto P2 = inode("P{i_1;a_3}", make_V(), leaf("W{a_1;a_3}"));
    std::vector<EvalNode<EvalExpr>> forest{P1, P2};
    auto const s = dag_path(forest, cm, block_of);
    CHECK(peak_profile_replay(s) == peak_profile_sweep(s).peak_bytes);
  }
}

TEST_CASE(
    "peak_profile_replay == sweep with a DEMOTED value (empty cross-occ meet)",
    "[peak_profile]") {
  // The core O3b demoted case: an internal value V occurs TWICE.
  //   - In tree 1, V sits under a parent P that realizes an External i_1 loop.
  //     stamp_lifetime_masks accumulates i_1 down to V, and i_1 is one of V's
  //     own slots, so V's occurrence-1 residency contribution is {i_1}.
  //   - In tree 2, the structurally-identical V sits under a parent Q with NO
  //     batch loop, so its occurrence-2 contribution is {} (empty).
  // The cross-occurrence MEET is {i_1} INTERSECT {} = {} => V's home_scope is
  // empty => it homes at the chain root (home_depth == -1) and is sized FULL,
  // EVEN THOUGH in tree 1 it is carried inside the i_1 loop's ectx. That is a
  // demotion: a value carried full above a loop it slices in another
  // occurrence. Both static algorithms read the SAME resulting Schedule, so
  // their equality must hold regardless of the demotion.
  auto make_V = [] {
    return inode("V{i_1;i_2}", leaf("V1{i_1;x_1}"), leaf("V2{x_1;i_2}"));
  };
  auto P = inode("P{i_2;a_1}", make_V(), leaf("W{i_1;a_1}"));
  stamp_ext(P, Index{L"i_1"});  // realized i_1 loop above V in tree 1
  auto Q = inode("Q{i_1;i_2}", make_V(), leaf("U{i_2;i_1}"));  // no loop

  SizeRegime r;
  r.space_extent = {{L"i", 10}, {L"a", 10}, {L"x", 10}};
  CostModel const cm{r};
  auto const block_of = [](Index const& ix) -> std::size_t {
    return ix == Index{L"i_1"} ? 2 : 10;
  };

  std::vector<EvalNode<EvalExpr>> forest{P, Q};
  auto const s = dag_path(forest, cm, block_of);

  // V (carried {i_1,i_2}) must be the demoted cell: home_depth == -1 and sized
  // FULL (10*10*8 = 800), NOT block-narrowed to i_1=2 despite tree 1's loop.
  Cell const* v = nullptr;
  for (auto const& c : s.cells)
    if (c.footprint == 800 && c.home_depth == -1) v = &c;
  REQUIRE(v != nullptr);  // the demotion actually happened
  CHECK(v->footprint == 800);
  CHECK(v->home_depth == -1);

  // THE Step-B equality on the demoted forest.
  CHECK(peak_profile_replay(s) == peak_profile_sweep(s).peak_bytes);
}

// ---- Step C: the Phase-1 anchor (non-demoted forest) ----------------

// A small order-of-magnitude regime (occ i=10, virt a=20), matching
// backend_test_regime() in test_eval_dryrun.cpp: same extents => the static
// footprints and the metered replay's sizes are drawn from the identical
// memsize model.
namespace {
SizeRegime anchor_regime() {
  SizeRegime r;
  r.space_extent = {{L"i", 10}, {L"a", 20}};
  double const pno = 4.0;
  for (std::size_t k = 0; k <= 4; ++k) r.csv_pno_moment[k] = std::pow(pno, k);
  r.csv_osv_moment = r.csv_pno_moment;
  return r;
}
}  // namespace

TEST_CASE("peak-profile anchor: static sweep vs metered replay co-resident sum",
          "[peak_profile]") {
  // WHY this is gated to a NON-DEMOTED forest (design section 9.6): the static
  // seed placement and the runtime External-only heuristic (sliced_modes)
  // diverge on a demoted value (empty cross-occurrence meet homes it at the
  // root while the heuristic would slice it in the occurrence that carries the
  // loop). The anchor is only a valid equality where the two COINCIDE -- i.e.
  // a forest with NO batching at all, where every home_scope is empty on both
  // sides and every value sizes FULL. This forest has zero node_slice_mask
  // loops.
  auto ctx = sequant::get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  auto const regime = anchor_regime();

  // A plain two-contraction forest, no batch loops: X = (g * t) is one product;
  // the outer contracts X with another leaf. Deliberately CSE-free so the
  // runtime cache holds a single monotone working set.
  auto expr =
      deserialize<ExprPtr>("(g{i_1,i_2;a_1,a_2} * t{a_1,i_3;}) * u{a_2;i_3}");
  REQUIRE(static_cast<bool>(expr));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = binarize<EvalExprDryRun>(expr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  std::vector<EvalNodeDryRun> forest{node};

  // ---- runtime measurement: metered-replay co-resident-sum peak ----
  BatchPolicy policy;  // no batchable indices => no batching engages
  CacheConfig cfg;
  cfg.max_footprint = 1e11;
  cfg.min_repeats = 1;
  cfg.is_volatile = [](EvalNodeDryRun const&) { return false; };
  double const replay_peak_bytes =
      sequant::eval::dryrun::meter(forest, policy, regime, cfg).peak_bytes;

  // ---- static measurement: dag_path + sweep, SAME regime ----
  CostModel const cm{regime};
  auto const block_of = [](Index const&) -> std::size_t { return 1; };
  auto const sched = dag_path(forest, cm, block_of);
  auto const sweep = peak_profile_sweep(sched);

  // Both static algorithms agree (Step B holds here too).
  CHECK(peak_profile_replay(sched) == sweep.peak_bytes);

  std::wcerr << L"\n[peak_profile-anchor] static sweep peak_bytes="
             << sweep.peak_bytes << L"  metered replay peak_bytes="
             << replay_peak_bytes << L"\n";

  // THE anchor: the static continuous-liveness sweep over seed cells must match
  // the runtime co-resident-sum measurement on this non-demoted forest.
  CHECK(sweep.peak_bytes == replay_peak_bytes);
}
