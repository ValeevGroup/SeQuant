// Phase 3b T1: tests for the static peak-profile sizing primitives
// (SeQuant/core/eval/peak_profile.hpp). Two free functions are pinned here:
//   - home_depth_of: resolve a residency mode-set to an enclosing-batch-
//     context loop depth (mirrors the runtime rl-walk at eval.hpp:1776-1782).
//   - cell_footprint: home-relative footprint of a carried-index set via the
//     existing dryrun::CostModel::memsize (no sizing logic reimplemented).
// Neither function has a production caller yet (T2/T3 wire them in); this
// task only pins the two primitives' contracts.

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
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
using sequant::EvalExpr;
using sequant::EvalNode;
using sequant::ExprPtr;
using sequant::Index;
using sequant::container::svector;
using sequant::eval::Cell;
using sequant::eval::linearize;
using sequant::eval::peak_profile_sweep;
using sequant::eval::PeakProfile;
using sequant::eval::Schedule;
using sequant::eval::detail::BatchContext;
using sequant::eval::detail::cell_footprint;
using sequant::eval::detail::home_depth_of;
using sequant::eval::dryrun::CostModel;
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
