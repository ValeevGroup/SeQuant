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
#include <SeQuant/core/eval/legality.hpp>
#include <SeQuant/core/eval/ordered_schedule.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/eval/value_node_map.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <cstdlib>
#include <set>
#include <string_view>
#include <variant>
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

// An EvalExpr carrying \p result's slots but a CALLER-CHOSEN node id (hash):
// the public EvalExpr ctor takes the hash outright, which is what lets a test
// force two structurally different values onto one 64-bit key (and, in the
// deep-spine test, mint thousands of distinct node ids without re-parsing).
EvalExpr eval_tensor_hashed(std::string_view result, std::size_t hash,
                            sequant::EvalOp op) {
  auto expr = sequant::deserialize<ExprPtr>(std::string(result));
  REQUIRE(static_cast<bool>(expr));
  EvalExpr const base{expr->as<sequant::Tensor>()};
  EvalExpr::index_vector ixs{base.canon_indices().begin(),
                             base.canon_indices().end()};
  return EvalExpr{op,          sequant::ResultType::Tensor,
                  expr,        std::move(ixs),
                  /*phase=*/1, hash,
                  nullptr};
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

// ---------------------------------------------------------------------
// Copilot review (PR #613): value-cell grouping must survive a value-key
// collision, and the linearizing prepass must survive a deep Sum spine.
// ---------------------------------------------------------------------

// Shared by the two collision tests: build the rich schedule for `forest`,
// then push it through the downstream resolutions the review flagged --
// ordered_schedule_dep_graph's `value_id_of`, analyze_legality's
// `CellLegality::hash`, build_ordered_schedule -- and assert that every cell
// keeps its OWN value id all the way through. All of those are first-wins
// `emplace`s on `value_key_of(ValueCell)`, so two cells sharing one key would
// silently collapse to the first one's value_id and the split the structural
// check made would not propagate at all.
void check_split_propagates(std::vector<EvalNode<EvalExpr>> const& forest,
                            std::size_t expected_cells,
                            std::size_t collide_hash) {
  SizeRegime r;
  r.space_extent = {{L"i", 5}, {L"a", 10}};
  CostModel const cm{r};
  auto const block_of = [](Index const&) -> std::size_t { return 1; };

  RichSchedule const rich = compute_dag_boulevard(forest, cm, block_of);
  REQUIRE(rich.cells.size() == expected_cells);

  // Two cells carry the colliding NODE id, one occurrence each ...
  std::size_t colliding_cells = 0;
  for (auto const& c : rich.cells)
    if (c.hash == collide_hash) {
      ++colliding_cells;
      CHECK(c.occurrences.size() == 1);
    }
  CHECK(colliding_cells == 2);

  // ... and every cell's VALUE id is distinct, which is what makes the split
  // visible downstream.
  std::set<std::size_t> keys;
  for (auto const& c : rich.cells) keys.insert(sequant::eval::value_key_of(c));
  CHECK(keys.size() == rich.cells.size());

  // Node-side agreement: compute_dag_boulevard re-stamps EvalExpr::value_key,
  // so the key->node maps the executor and legality join on resolve each cell
  // to one of its OWN nodes.
  std::set<std::size_t> node_keys;
  for (auto const& t : forest)
    t.visit([&](auto const& n) { node_keys.insert(sequant::value_key_of(n)); });
  CHECK(node_keys == keys);

  // The dep graph's value_id_of is the first-wins map: one entry per cell.
  auto const g = sequant::eval::detail::ordered_schedule_dep_graph(rich);
  CHECK(g.value_id_of.size() == rich.cells.size());

  sequant::BatchPolicy policy;
  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  REQUIRE(legality.cells.size() == rich.cells.size());
  std::set<std::size_t> legality_keys;
  for (auto const& cl : legality.cells) legality_keys.insert(cl.hash);
  CHECK(legality_keys == keys);

  sequant::eval::OrderedSchedule sched;
  REQUIRE_NOTHROW(sched = sequant::eval::build_ordered_schedule(rich, legality,
                                                                policy, {}));
  // NOT `sched.num_values == rich.cells.size()`: build_ordered_schedule
  // assigns exactly that, so the comparison cannot fail. What bites on a
  // first-wins collapse is the schedule's actual PRODUCTIONS -- every non-leaf
  // cell has to be built somewhere in it, under its OWN value id, so a
  // collapsed colliding pair would leave one of the two products unproduced
  // (and the other built once for both).
  std::multiset<std::size_t> built;
  auto collect = [&](auto&& self, sequant::eval::ScopeBlock const& b) -> void {
    for (auto const& st : b.steps) {
      if (auto const* bs = std::get_if<sequant::eval::BuildStep>(&st.value))
        built.insert(bs->value_id);
      else if (auto const* sb =
                   std::get_if<sequant::eval::ScopeBlock>(&st.value))
        self(self, *sb);
    }
  };
  collect(collect, sched.root);
  std::multiset<std::size_t> expected_built;
  for (auto const& c : rich.cells)
    if (!c.is_leaf) expected_built.insert(c.value_id);
  CHECK(built == expected_built);
}

TEST_CASE(
    "compute_dag_boulevard opens separate cells for two values that COLLIDE "
    "on the value key",
    "[peak_profile][value-cell]") {
  // The grouping used to be `hash_to_cell.find(r.key)` alone: a bucket hit on
  // the 64-bit value key WAS identity. Two distinct values colliding there
  // would merge into one ValueCell, and the schedule would then build one of
  // them and read it as the other -- silently the wrong value. A bucket hit is
  // now confirmed structurally (TreeNodeEqualityComparator through the
  // canonical child view, plus the home slicing the key folds in and the
  // operands' CELLS), and a mismatch opens a new cell whose value id is salted
  // so the split reaches the schedule.
  //
  // The collision is forced through the public EvalExpr ctor, which takes the
  // node id outright. Here the two products differ in their OPERANDS, so the
  // inductive child-cell comparison is what separates them.
  std::size_t const collide = 0xC0111DEULL;

  auto X = EvalNode<EvalExpr>{
      eval_tensor_hashed("I{i_1;a_1}", collide, sequant::EvalOp::Product),
      leaf("t{i_1;a_3}"), leaf("g{a_3;a_1}")};
  auto Y = EvalNode<EvalExpr>{
      eval_tensor_hashed("J{i_1;a_1}", collide, sequant::EvalOp::Product),
      leaf("u{i_1;a_4}"), leaf("h{a_4;a_1}")};
  REQUIRE(X->hash_value() == Y->hash_value());

  // 4 distinct leaves + the two colliding products = 6 cells. (Before the fix
  // this was 5: the products merged into one cell with two occurrences.)
  check_split_propagates({X, Y}, 6, collide);
}

TEST_CASE(
    "compute_dag_boulevard separates key-colliding values with IDENTICAL "
    "operands",
    "[peak_profile][value-cell]") {
  // The companion case: the two colliding products contract the SAME two
  // leaves, so they fold to the same child CELLS and the same (empty)
  // loop_slot / reduced_slot -- every cheap discriminator ties, and only
  // TreeNodeEqualityComparator on the nodes themselves can tell them apart
  // (here on the result tensor's block, I{i_1;a_1} vs J{i_1;a_1}, since
  // neither carries a connectivity graph). Without that last comparison the
  // two would merge.
  std::size_t const collide = 0xC0111DE2ULL;

  auto X = EvalNode<EvalExpr>{
      eval_tensor_hashed("I{i_1;a_1}", collide, sequant::EvalOp::Product),
      leaf("t{i_1;a_3}"), leaf("g{a_3;a_1}")};
  auto Y = EvalNode<EvalExpr>{
      eval_tensor_hashed("J{i_1;a_1}", collide, sequant::EvalOp::Product),
      leaf("t{i_1;a_3}"), leaf("g{a_3;a_1}")};
  REQUIRE(X->hash_value() == Y->hash_value());
  // The operands really are one value each (2 leaf cells, not 4).
  check_split_propagates({X, Y}, 4, collide);
}

TEST_CASE(
    "the ordered prepasses walk a deep Sum spine without overflowing the "
    "stack",
    "[.][peak_profile][stack-safety]") {
  // An equation's residual/energy reaches the ordered path as a SINGLE
  // in-place Sum tree with one node per summand, so its LEFT SPINE is as deep
  // as the number of terms -- thousands for a large equation. Every prepass
  // that walks the forest therefore has to unwind that spine iteratively (the
  // executor itself already does); a recursive descent overflows the call
  // stack long before the executor is reached.
  //
  // Hidden ([.]) only because 40000 nodes cost ~4 GB to hold; the walks
  // themselves are sub-second. Run it by name or by [stack-safety];
  // SEQUANT_UT_SPINE_N overrides the depth.
  //
  // build_value_node_map / build_value_key_node_map are covered here too (see
  // the end of the test). They used to hold each node BY VALUE, and Node's
  // copy constructor deep-copies the subtree, so one entry per node cost
  // O(nodes x subtree) memory -- measured on this very tree: 1.8 GB at N=500,
  // 6.9 GB at N=1000, 15.7 GB at N=2000, OOM past that, with
  // evaluate_ordered_schedule building one per run. They now hold non-owning
  // pointers into the forest (ValueNodeMap), which is what makes them
  // affordable at this depth at all.
  std::size_t N = 20000;
  if (char const* n = std::getenv("SEQUANT_UT_SPINE_N"))
    N = static_cast<std::size_t>(std::atoll(n));

  // Built by COPYING two parsed EvalExprs (no re-parsing per level), so the
  // test measures the walks, not the parser. Each spine node gets its own node
  // id, as the real per-summand accumulators do.
  EvalExpr const lf = eval_tensor("t{i_1;a_1}");
  EvalNode<EvalExpr> tree{lf};
  for (std::size_t k = 1; k < N; ++k)
    tree = EvalNode<EvalExpr>{
        eval_tensor_hashed("R{i_1;a_1}", k, sequant::EvalOp::Sum),
        std::move(tree), EvalNode<EvalExpr>{lf}};
  REQUIRE(tree.size() == 2 * N - 1);

  std::vector<EvalNode<EvalExpr>> forest{tree};

  // The full ordered prepass entry point: stamp_lifetime_masks +
  // stamp_occurrence_homes + the post-order linearizing walk, all three of
  // which were converted.
  SizeRegime r;
  r.space_extent = {{L"i", 5}, {L"a", 10}};
  CostModel const cm{r};
  auto const block_of = [](Index const&) -> std::size_t { return 1; };
  RichSchedule const rich = compute_dag_boulevard(forest, cm, block_of);
  CHECK(rich.num_points == 2 * N - 1);
  // N-1 distinct spine values plus the one folded leaf value.
  CHECK(rich.cells.size() == N);

  // Bottom-up value-key fold over the whole spine. Nothing is home-sliced
  // here, so every key collapses to the node id -- the point is that the fold
  // REACHES the bottom of the spine at all.
  CHECK(sequant::value_key_of(tree) == tree->hash_value());
  CHECK(sequant::value_key_of(tree.left()) == tree.left()->hash_value());

  // The remaining converted walks over the same forest: the two value->node
  // bridges and analyze_legality's node_of pre-order.
  sequant::BatchPolicy policy;
  sequant::eval::LegalitySchedule legality;
  REQUIRE_NOTHROW(legality =
                      sequant::eval::analyze_legality(rich, forest, policy));
  CHECK(legality.cells.size() == rich.cells.size());

  // The two value->node bridges over the same forest. Nothing here is
  // home-sliced, so every node's value key IS its node id: N-1 distinct spine
  // ids plus the one folded leaf id, i.e. one entry per cell in both maps.
  auto const vmap = sequant::eval::build_value_node_map(forest);
  auto const vkmap = sequant::eval::build_value_key_node_map(forest);
  CHECK(vmap.size() == rich.cells.size());
  CHECK(vkmap.size() == rich.cells.size());
  // ... and an entry is a VIEW of the forest node, not a deep copy of its
  // subtree: the root's entry is the root's own address. (This is the memory
  // property the numbers in the comment at the top measured; holding it by
  // value is what made a by-value map quadratic in the tree.)
  CHECK(vmap.at(sequant::value_key_of(forest.front())) == &forest.front());
  CHECK(vkmap.at(sequant::value_key_of(forest.front())) == &forest.front());
}
