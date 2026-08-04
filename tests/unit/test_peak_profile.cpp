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
#include <SeQuant/core/eval/backends/dryrun/cost_profile.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
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
using sequant::eval::linearize;
using sequant::eval::peak_profile_replay;
using sequant::eval::peak_profile_sweep;
using sequant::eval::PeakProfile;
using sequant::eval::Schedule;
using sequant::eval::detail::BatchContext;
using sequant::eval::detail::cell_footprint;
using sequant::eval::detail::home_depth_of;
using sequant::eval::dryrun::CacheConfig;
using sequant::eval::dryrun::cost_profile;
using sequant::eval::dryrun::CostModel;
using sequant::eval::dryrun::CostProfile;
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
  n->set_batched_here({{std::move(ix), BatchModeType::External}});
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

  // p's own loop is the sole enclosing context; block(p) narrows it to 2.
  BatchContext const ectx{{p, {0, 2}}};
  auto const block_of = [](Index const& ix) -> std::size_t {
    return ix == Index{L"i_1"} ? 2 : 10;
  };

  // Cell homed AT p (home_modes = {p}) -> d = 0 -> p's loop (level 0 <= d)
  // encloses the home -> p sized at BLOCK (2), q sized at FULL (10).
  {
    svector<Index> const home_modes{p};
    auto const got = cell_footprint(carried, home_modes, ectx, cm, block_of);
    ExtentOverrides ov;
    ov[p] = 2;
    CHECK(got == cm.memsize(carried, ov));
    // memsize() reports BYTES (elems * CostModel's numeric_size_ = 8.0), not
    // a raw element count: 2 (block p) * 10 (full q) * 8 bytes/elem.
    CHECK(got == 160);
  }

  // Same cell homed ABOVE p's loop (home_modes = {}) -> d = -1 -> no level
  // encloses the home -> both modes sized FULL.
  {
    svector<Index> const home_modes{};
    auto const got = cell_footprint(carried, home_modes, ectx, cm, block_of);
    CHECK(got == cm.memsize(carried));
    // 10 (full p) * 10 (full q) * 8 bytes/elem.
    CHECK(got == 800);
  }
}

// ---------------------------------------------------------------------
// T2: the interval-event sweep over a hand-built Schedule (bypasses
// linearize entirely; Cells are constructed directly).
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

// ---------------------------------------------------------------------
// T2: linearize -- the real forest walk (CSE grouping by hash_value +
// post-order interval assignment + home-relative footprint).
// ---------------------------------------------------------------------

TEST_CASE(
    "linearize groups a CSE value into ONE cell spanning to its LAST consumer",
    "[peak_profile]") {
  // A shared value V (a common sub-value) is consumed by TWO parents P1, P2,
  // one per tree. The two structurally-identical copies of V carry the SAME
  // hash_value(), so linearize must fold them into exactly ONE cell whose
  // first_use is V's (earliest) production point and whose last_use is the
  // LATER of the two parents' points.
  //
  // Sizes are chosen so every distinct value has a DISTINCT footprint, making
  // the shared value's cell unambiguous: with full(i)=5, full(a)=10,
  //   V{i_1;i_2} = 5*5*8   = 200  (the shared value -- the only 200-byte cell)
  //   W{a_1;a_3} = 10*10*8 = 800  (the shared sibling)
  //   P{i_1;a_3} = 5*10*8  = 400  (the shared root)
  auto make_V = [] { return leaf("V{i_1;i_2}"); };
  auto P1 = inode("P{i_1;a_3}", make_V(), leaf("W{a_1;a_3}"));
  auto P2 = inode("P{i_1;a_3}", make_V(), leaf("W{a_1;a_3}"));

  SizeRegime r;
  r.space_extent = {{L"i", 5}, {L"a", 10}};
  CostModel const cm{r};
  auto const block_of = [](Index const&) -> std::size_t { return 1; };

  std::vector<EvalNode<EvalExpr>> forest{P1, P2};
  auto const sched = linearize(forest, cm, block_of);

  // Post-order per tree: V(0),W(1),P(2) and V(3),W(4),P(5) => 6 points. The
  // pairwise-identical V, W, and P values each fold => exactly 3 cells.
  CHECK(sched.num_points == 6);
  CHECK(sched.cells.size() == 3);

  // Find V's cell by its unique 200-byte footprint.
  std::size_t v_cells = 0;
  Cell const* v = nullptr;
  for (auto const& c : sched.cells)
    if (c.footprint == 200) {
      ++v_cells;
      v = &c;
    }
  REQUIRE(v_cells == 1);  // exactly ONE cell for the shared value
  REQUIRE(v != nullptr);
  CHECK(v->home_depth == -1);  // no enclosing loop
  // first_use is the FIRST tree's V production (point 0); last_use is the
  // SECOND tree's parent P (point 5, the global last point) -- the LATER
  // consumer, proving the fold widened the interval across occurrences.
  CHECK(v->first_use == 0);
  CHECK(v->last_use == 5);
  CHECK(v->last_use == sched.num_points - 1);
}

TEST_CASE(
    "linearize sizes a loop-carried value at BLOCK and its loop result at FULL",
    "[peak_profile]") {
  // R = (C * D) with an External i_1 loop realized AT R. R's own result slot
  // carries i_1, and so does its left child C. stamp_seed_residency (invoked
  // inside linearize) therefore stamps BOTH with home_scope = {i_1}.
  //
  // During the post-order walk R's i_1 loop is pushed onto the enclosing
  // context for its CHILDREN (so C sees ectx=[i_1] and is sized per-block) but
  // popped before R itself is recorded (so R sees ectx=[] -- it is the loop
  // RESULT, sized FULL over the whole i_1 range).
  auto C = inode("C{i_1;a_3}", leaf("C1{i_1;a_1}"), leaf("C2{a_1;a_3}"));
  auto D = leaf("D1{a_3;a_5}");
  auto R = inode("R{i_1;a_5}", std::move(C), std::move(D));
  stamp_ext(R, Index{L"i_1"});

  // full(i)=full(a)=10; block(i_1)=2.
  SizeRegime r;
  r.space_extent = {{L"i", 10}, {L"a", 10}};
  CostModel const cm{r};
  auto const block_of = [](Index const& ix) -> std::size_t {
    return ix == Index{L"i_1"} ? 2 : 10;
  };

  std::vector<EvalNode<EvalExpr>> forest{R};
  auto const sched = linearize(forest, cm, block_of);

  // Locate C's and R's cells by their carried-slot footprints. C carries
  // {i_1,a_3}: homed AT the i_1 loop (depth 0) => i_1 sized BLOCK(2), a_3
  // FULL(10) => 2*10*8 = 160, and home_depth == 0. R carries {i_1,a_5}: it is
  // the loop RESULT (ectx empty at R) => i_1 FULL(10), a_5 FULL(10) =>
  // 10*10*8 = 800, and home_depth == -1.
  Cell const* c_cell = nullptr;
  Cell const* r_cell = nullptr;
  for (auto const& c : sched.cells) {
    if (c.footprint == 160 && c.home_depth == 0) c_cell = &c;
    if (c.footprint == 800 && c.home_depth == -1) r_cell = &c;
  }
  REQUIRE(c_cell != nullptr);  // loop-carried value sized at BLOCK
  REQUIRE(r_cell != nullptr);  // loop result sized at FULL
  CHECK(c_cell->footprint == 160);
  CHECK(r_cell->footprint == 800);
}

// =====================================================================
// T3: the independent REPLAY ORACLE (peak_profile_replay) and the two-
// part validation from design section 9.6.
//   Step A -- oracle == sweep on hand-built Schedules (isolated).
//   Step B -- oracle == sweep on real linearized forests, INCLUDING a
//             demoted case (both sides consume the SAME Schedule, so the
//             equality validates the sweep's interval algorithm).
//   Step C -- the Phase-1 anchor: static sweep vs runtime cost_profile()
//             co-resident-sum peak on a NON-demoted forest.
// =====================================================================

// ---- Step A: oracle == sweep on hand-built Schedules ----------------

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
    auto const s = linearize(forest, cm, block_of);
    CHECK(peak_profile_replay(s) == peak_profile_sweep(s).peak_bytes);
  }

  // (2) The T2 CSE forest (a shared V consumed by two parents): the fold
  // widens V's interval across occurrences, exercising min/max over the
  // grouped cell -- the sweep and the replay must still agree exactly.
  {
    auto make_V = [] { return leaf("V{i_1;i_2}"); };
    auto P1 = inode("P{i_1;a_3}", make_V(), leaf("W{a_1;a_3}"));
    auto P2 = inode("P{i_1;a_3}", make_V(), leaf("W{a_1;a_3}"));
    std::vector<EvalNode<EvalExpr>> forest{P1, P2};
    auto const s = linearize(forest, cm, block_of);
    CHECK(peak_profile_replay(s) == peak_profile_sweep(s).peak_bytes);
  }
}

TEST_CASE(
    "peak_profile_replay == sweep with a DEMOTED value (empty cross-occ meet)",
    "[peak_profile]") {
  // The core O3b demoted case: an internal value V occurs TWICE.
  //   - In tree 1, V sits under a parent P that realizes an External i_1 loop.
  //     stamp_seed_residency accumulates i_1 down to V, and i_1 is one of V's
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
  auto const s = linearize(forest, cm, block_of);

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
// linearize footprints and the runtime cost_profile sizes are drawn from the
// identical memsize model.
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

TEST_CASE(
    "peak-profile anchor: static sweep vs runtime cost_profile co-resident sum",
    "[peak_profile]") {
  // WHY this is gated to a NON-DEMOTED forest (design section 9.6): the static
  // seed placement and the runtime External-only heuristic (sliced_modes)
  // diverge on a demoted value (empty cross-occurrence meet homes it at the
  // root while the heuristic would slice it in the occurrence that carries the
  // loop). The anchor is only a valid equality where the two COINCIDE -- i.e.
  // a forest with NO batching at all, where every home_scope is empty on both
  // sides and every value sizes FULL. This forest has zero batched_here loops.
  auto ctx = sequant::get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  auto const regime = anchor_regime();

  // A plain two-contraction forest, no batch loops: X = (g * t) is one product;
  // the outer contracts X with another leaf. Deliberately CSE-free so the
  // runtime cache holds a single monotone working set.
  auto expr =
      deserialize<ExprPtr>("(g{i_1,i_2;a_1,a_2} * t{a_1;i_3}) * u{a_2;i_3}");
  REQUIRE(static_cast<bool>(expr));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = binarize<EvalExprDryRun>(expr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  std::vector<EvalNodeDryRun> forest{node};

  // ---- runtime measurement: cost_profile co-resident-sum peak ----
  BatchPolicy policy;  // no batchable indices => no batching engages
  CacheConfig cfg;
  cfg.max_footprint = 1e11;
  cfg.min_repeats = 1;
  cfg.is_volatile = [](EvalNodeDryRun const&) { return false; };
  CostProfile const cp = cost_profile(forest, policy, cfg, regime,
                                      /*trace=*/nullptr);

  // ---- static measurement: linearize + sweep, SAME regime ----
  CostModel const cm{regime};
  auto const block_of = [](Index const&) -> std::size_t { return 1; };
  auto const sched = linearize(forest, cm, block_of);
  auto const sweep = peak_profile_sweep(sched);

  // Both static algorithms agree (Step B holds here too).
  CHECK(peak_profile_replay(sched) == sweep.peak_bytes);

  std::wcerr << L"\n[peak_profile-anchor] static sweep peak_bytes="
             << sweep.peak_bytes << L"  runtime cost_profile peak_bytes="
             << cp.peak_bytes << L"\n";

  // THE anchor: the static continuous-liveness sweep over seed cells must match
  // the runtime co-resident-sum measurement on this non-demoted forest.
  CHECK(sweep.peak_bytes == cp.peak_bytes);
}
