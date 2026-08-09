// Phase 4a T1: tests for O2's standalone working representation
// (SeQuant/core/eval/placement_remat.hpp) plus the peak_profile.hpp refactor
// that exposes it (compute_dag_boulevard / ValueCell / RichSchedule).
//
// Two things are pinned here:
//   - shrink_candidates finds the demoted carried mode of a partially-homed
//     value (reusing test_peak_profile.cpp's PARTIAL-DEMOTION forest), and
//     to_schedule's projected footprint matches the Phase-3b value.
//   - GUARDRAIL: compute_dag_path's flat Schedule (built via the NEW
//     compute_dag_boulevard + thin projection) is BYTE-IDENTICAL, cell for
//     cell, to what remat_cells + to_schedule (the O2 round trip over the SAME
//     rich record) produce -- proof the rich split changed nothing observable.
//
// Zero production callers of the O2 symbols; this task is standalone.

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_profile.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/eval/placement_remat.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstddef>
#include <limits>
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
using sequant::eval::apply_split;
using sequant::eval::Cell;
using sequant::eval::compute_dag_path;
using sequant::eval::HomeTarget;
using sequant::eval::occurrence_key;
using sequant::eval::OccurrenceRec;
using sequant::eval::peak_profile_sweep;
using sequant::eval::PeakProfile;
using sequant::eval::PlacementRouter;
using sequant::eval::remat_cells;
using sequant::eval::remat_to_router;
using sequant::eval::rematerialize_to_budget;
using sequant::eval::RematInput;
using sequant::eval::RematResult;
using sequant::eval::RematStatus;
using sequant::eval::Schedule;
using sequant::eval::shrink_candidates;
using sequant::eval::split_candidates;
using sequant::eval::to_schedule;
using sequant::eval::ValueCell;
using sequant::eval::dryrun::CostModel;
using sequant::eval::dryrun::EvalExprDryRun;
using sequant::eval::dryrun::EvalNodeDryRun;
using sequant::eval::dryrun::SizeRegime;

// The DAG-scope a moved HomeTarget carries for a cell's \p home_modes --
// mirrors remat_to_router's construction: each home mode's SPACE, ordered by
// the mode's position in the cell's enclosing NEST (\p enclosing_modes, absent
// modes sort last). Suffixed to avoid a Unity-build anonymous-namespace clash.
svector<sequant::IndexSpace> expected_dag_scope_o2(
    svector<Index> const& home_modes, svector<Index> const& enclosing_modes) {
  auto const nest_pos = [&](Index const& m) -> std::size_t {
    auto const it =
        std::find(enclosing_modes.begin(), enclosing_modes.end(), m);
    return it == enclosing_modes.end()
               ? enclosing_modes.size()
               : static_cast<std::size_t>(it - enclosing_modes.begin());
  };
  svector<Index> ordered(home_modes.begin(), home_modes.end());
  std::stable_sort(ordered.begin(), ordered.end(),
                   [&](Index const& a, Index const& b) {
                     return nest_pos(a) < nest_pos(b);
                   });
  svector<sequant::IndexSpace> out;
  for (auto const& m : ordered) out.push_back(m.space());
  return out;
}

// Build an EvalExpr from a single-tensor spec (e.g. "R{i_1;a_5}"); its
// canon_indices are exactly the tensor's bra+ket slots. Mirrors the helper in
// test_peak_profile.cpp / test_lifetime_mask.cpp. Suffixed `_o2` (rather than
// reusing the exact same name) because a Unity build unifies this file with
// test_peak_profile.cpp into one translation unit, where an anonymous
// namespace does not prevent a same-named free function from colliding.
EvalExpr eval_tensor_remat(std::string_view tensor) {
  auto expr = sequant::deserialize<ExprPtr>(std::string(tensor));
  REQUIRE(static_cast<bool>(expr));
  return EvalExpr{expr->as<sequant::Tensor>()};
}

// A leaf eval node carrying the given tensor's slots.
EvalNode<EvalExpr> leaf_remat(std::string_view tensor) {
  return EvalNode<EvalExpr>{eval_tensor_remat(tensor)};
}

// An internal eval node whose OWN result slots are the given tensor's slots,
// with the two supplied child subtrees.
EvalNode<EvalExpr> inode_remat(std::string_view result, EvalNode<EvalExpr> l,
                               EvalNode<EvalExpr> r) {
  return EvalNode<EvalExpr>{eval_tensor_remat(result), std::move(l),
                            std::move(r)};
}

// Stamp a single External batch loop mode at a node.
void stamp_ext_remat(EvalNode<EvalExpr>& n, Index ix) {
  n->set_batched_here({{std::move(ix), BatchModeType::External}});
}

// Stamp TWO External batch loop modes at a node, outer first: realizes a
// two-level nest (ix_outer at level 0, ix_inner at level 1) in ONE node
// rather than needing a separate ancestor per level.
void stamp_ext_pair_remat(EvalNode<EvalExpr>& n, Index ix_outer,
                          Index ix_inner) {
  n->set_batched_here({{std::move(ix_outer), BatchModeType::External},
                       {std::move(ix_inner), BatchModeType::External}});
}

// Assert two flat Schedules are cell-for-cell identical: same cell count,
// same value_id / home_depth / footprint / first_use / last_use in the SAME
// order, and the same num_points. This is the byte-identity guardrail's
// actual check.
void require_schedules_identical(Schedule const& a, Schedule const& b) {
  REQUIRE(a.num_points == b.num_points);
  REQUIRE(a.cells.size() == b.cells.size());
  for (std::size_t i = 0; i < a.cells.size(); ++i) {
    Cell const& ca = a.cells[i];
    Cell const& cb = b.cells[i];
    CHECK(ca.value_id == cb.value_id);
    CHECK(ca.home_depth == cb.home_depth);
    CHECK(ca.footprint == cb.footprint);
    CHECK(ca.first_use == cb.first_use);
    CHECK(ca.last_use == cb.last_use);
  }
}

// The O2 round trip over the SAME forest: remat_cells (compute_dag_boulevard
// mapped to ValueCell) immediately projected back via to_schedule, with NO
// shrink move applied. Used only to prove the rich split is faithful
// (guardrail) -- no production caller does this.
template <typename Forest, typename BlockOfFn>
Schedule remat_round_trip(Forest const& forest, CostModel const& cm,
                          BlockOfFn const& block_of) {
  RematInput const in = remat_cells(forest, cm, block_of);
  return to_schedule(in.cells, cm, block_of, in.num_points);
}

// Build an ValueCell directly (bypassing compute_dag_path) so
// rematerialize_to_budget's greedy logic can be tested in isolation with
// precise control over footprints and co-residency. home_depth is set to 0 when
// a home is given, else -1 (only informational; rematerialize_to_budget does
// not read it).
ValueCell mk_cell(std::size_t id, std::vector<Index> carried,
                  std::vector<Index> home, std::vector<Index> enclosing,
                  std::size_t first_use, std::size_t last_use) {
  ValueCell c;
  c.value_id = id;
  c.home_depth = home.empty() ? -1 : 0;
  c.carried.assign(carried.begin(), carried.end());
  c.home_modes.assign(home.begin(), home.end());
  c.enclosing_modes.assign(enclosing.begin(), enclosing.end());
  c.first_use = first_use;
  c.last_use = last_use;
  return c;
}

// A single occurrence record for a hand-built divergent ValueCell.
OccurrenceRec mk_occ(std::size_t point, std::size_t consumer_point,
                     std::vector<Index> carried, std::vector<Index> enclosing) {
  OccurrenceRec o;
  o.point = point;
  o.consumer_point = consumer_point;
  o.carried.assign(carried.begin(), carried.end());
  // home stays empty (the divergent seed's cross-occurrence meet is empty).
  for (auto const& m : enclosing) o.ectx.push_back({m, {std::size_t{0}, 5}});
  return o;
}

// Build a DIVERGENT CSE ValueCell whose two occurrences bind a RELABELED mode
// to i_3 (outer occurrence, carried at slot 1) and i_4 (inner occurrence). The
// value's structural hash is the same across both (they are CSE-folded); only
// the physical label at the relabeled slot differs. Leg A is enclosed by i_3
// only (and carries it); leg B is enclosed by BOTH i_3 and i_4 but carries only
// i_4 -- so leg B is invariant to i_3 (homed within it, does not carry it) and
// must be rebuilt per i_3 block: the replication factor the split prices.
ValueCell mk_divergent_cell(std::size_t hash) {
  ValueCell c;
  c.value_id = 0;
  c.hash = hash;
  c.carried = svector<Index>{Index{L"x_1"}, Index{L"i_3"}};  // first-occ order
  c.home_modes = {};  // divergent meet is empty
  c.enclosing_modes = svector<Index>{Index{L"i_3"}, Index{L"i_4"}};
  c.divergent_modes = svector<Index>{Index{L"i_3"}, Index{L"i_4"}};
  c.first_use = 0;
  c.last_use = 2;
  c.home_depth = -1;
  // occ A (leg carrying i_3): enclosed by i_3 only.
  c.occurrences.push_back(mk_occ(/*point=*/0, /*consumer=*/2,
                                 {Index{L"x_1"}, Index{L"i_3"}},
                                 {Index{L"i_3"}}));
  // occ B (leg carrying i_4): enclosed by i_3 (outer) and i_4 (inner).
  c.occurrences.push_back(mk_occ(/*point=*/1, /*consumer=*/2,
                                 {Index{L"x_1"}, Index{L"i_4"}},
                                 {Index{L"i_3"}, Index{L"i_4"}}));
  return c;
}

}  // namespace

TEST_CASE(
    "remat splits along a relabeled mode, prices replication not 2x "
    "(Task 5)",
    "[placement_remat]") {
  // i-space extent 10, block SIZE 5 -> the i_3 loop has block COUNT
  // ceil(10/5) = 2 = N3, the replication factor of the inner cell. cell A homes
  // {i_3} (carries i_3) -> replication 1; cell B homes {i_4} but is invariant
  // to the enclosing i_3 loop it does not carry -> replication = i_3's block
  // COUNT (2), NOT its block size (5). Both cells' single build is identical
  // (i_3/i_4 same space, same block), so the split recomputes (1 + 2) = 3x that
  // build. The three models are distinguishable here: flat-2x fudge -> 2x; a
  // wrong block-SIZE replication -> (1 + 5) = 6x; correct block-COUNT -> 3x. So
  // asserting == 3x (and != 2x) pins the correct model against both.
  SizeRegime r;
  r.space_extent = {{L"i", 10}, {L"x", 10}};
  CostModel const cm{r};
  auto const block_of = [](Index const& ix) -> std::size_t {
    return ix.space() == Index{L"i_1"}.space() ? 5 : 2;
  };
  std::size_t const N3 = 2;  // block COUNT = ceil(extent 10 / block size 5)

  ValueCell const cell = mk_divergent_cell(/*hash=*/42);
  Index const i3{L"i_3"}, i4{L"i_4"};

  // The relabeled mode is offered ONLY as a split, never as an in-place shrink.
  auto const shrinks = shrink_candidates(cell);
  CHECK(std::find(shrinks.begin(), shrinks.end(), i3) == shrinks.end());
  CHECK(std::find(shrinks.begin(), shrinks.end(), i4) == shrinks.end());
  auto const splits = split_candidates(cell);
  CHECK(std::find(splits.begin(), splits.end(), i3) != splits.end());
  CHECK(std::find(splits.begin(), splits.end(), i4) != splits.end());

  // apply_split un-folds the one cell into TWO same-hash cells at DISTINCT
  // homes ({i_3} and {i_4}), each non-divergent, and prices the replication.
  sequant::container::svector<ValueCell> cells;
  cells.push_back(cell);
  std::size_t const recompute = apply_split(cells, /*ci=*/0, i3, cm, block_of);

  std::size_t n_hash42 = 0;
  svector<svector<Index>> homes;
  for (auto const& c : cells)
    if (c.hash == 42) {
      ++n_hash42;
      homes.push_back(c.home_modes);
      CHECK(c.divergent_modes.empty());  // each split cell is single-binding
    }
  CHECK(n_hash42 == 2);
  REQUIRE(homes.size() == 2);
  CHECK(homes[0] != homes[1]);  // distinct homes: {i_3} vs {i_4}

  // The priced recompute is (1 + N3) x a single build, NOT 2x.
  std::size_t const full_build = sequant::eval::detail::cell_footprint(
      svector<Index>{Index{L"x_1"}, i3}, svector<Index>{i3}, cm, block_of);
  CHECK(full_build == 5 * 10 * 8);  // block(i_3)=5 * full(x_1)=10 * 8 bytes
  CHECK(recompute == (1 + N3) * full_build);
  CHECK(recompute != 2 * full_build);  // the old flat-2x fudge is gone
}

TEST_CASE(
    "shrink_candidates finds the demoted carried mode of a "
    "PARTIALLY-demoted value",
    "[placement_remat]") {
  // Reuses test_peak_profile.cpp's exact PARTIAL-DEMOTION forest: a shared
  // value V carries TWO slots {o_1,i_1}. Occurrence 1 sits under a single
  // i_1 loop; occurrence 2 sits under BOTH an outer o_1 loop and an inner
  // i_1 loop. The cross-occurrence MEET home is {i_1} -- o_1 is batched in
  // only ONE occurrence, so it stays OUT of the home and must be a SHRINK
  // candidate: it is carried, and it DOES enclose V at occurrence 2, but it
  // is not (yet) part of V's home.
  auto make_V = [] {
    return inode_remat("V{o_1;i_1}", leaf_remat("V1{o_1;x_1}"),
                       leaf_remat("V2{x_1;i_1}"));
  };

  SizeRegime r;
  r.space_extent = {{L"o", 10}, {L"i", 10}, {L"x", 10}, {L"a", 10}};
  CostModel const cm{r};
  auto const block_of = [](Index const&) -> std::size_t { return 2; };
  // Same value the Phase-3b test pins: i_1 BLOCK (in the meet home), o_1
  // FULL (not in the meet) => 2 (block i_1) * 10 (full o_1) * 8 bytes/elem.
  std::size_t const expected_footprint = 160;

  auto P1 = inode_remat("P1{i_1;a_1}", make_V(), leaf_remat("W1{i_1;a_1}"));
  stamp_ext_remat(P1, Index{L"i_1"});  // occurrence 1: i_1 loop only
  auto P2 = inode_remat("P2{a_1;a_3}", make_V(), leaf_remat("W2{a_1;a_3}"));
  stamp_ext_pair_remat(
      P2, Index{L"o_1"},
      Index{L"i_1"});  // occurrence 2: o_1 (outer), i_1 (inner)

  std::vector<EvalNode<EvalExpr>> const forest{P1, P2};

  RematInput const in = remat_cells(forest, cm, block_of);

  // V is the unique cell with a non-trivial home_depth (same reasoning as
  // the Phase-3b test: P1/P2 are forest roots with empty ectx, every leaf
  // stays unstamped, only V sits under a parent loop).
  ValueCell const* v = nullptr;
  for (auto const& c : in.cells)
    if (c.home_depth != -1) v = &c;
  REQUIRE(v != nullptr);

  CHECK(v->home_modes == svector<Index>{Index{L"i_1"}});
  CHECK(shrink_candidates(*v) == svector<Index>{Index{L"o_1"}});

  Schedule const sched = to_schedule(in.cells, cm, block_of, in.num_points);
  CHECK(sched.cells[v->value_id].footprint == expected_footprint);
}

TEST_CASE(
    "GUARDRAIL: compute_dag_path's flat Schedule == remat_cells + "
    "to_schedule's "
    "round trip, cell for cell",
    "[placement_remat]") {
  SizeRegime r;
  r.space_extent = {{L"o", 10}, {L"i", 10}, {L"x", 10}, {L"a", 10}};
  CostModel const cm{r};
  auto const block_of = [](Index const&) -> std::size_t { return 2; };

  SECTION("plain two-contraction forest, no batching") {
    auto expr = deserialize<ExprPtr>(
        "(g{i_1,i_2;a_1,a_2} * t{a_1;i_3}) * "
        "u{a_2;i_3}");
    REQUIRE(static_cast<bool>(expr));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize<EvalExprDryRun>(expr);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    std::vector<EvalNodeDryRun> const forest{node};

    Schedule const legacy = compute_dag_path(forest, cm, block_of);
    Schedule const via_o2 = remat_round_trip(forest, cm, block_of);
    require_schedules_identical(legacy, via_o2);

    PeakProfile const legacy_peak = peak_profile_sweep(legacy);
    PeakProfile const via_o2_peak = peak_profile_sweep(via_o2);
    CHECK(legacy_peak.peak_bytes == via_o2_peak.peak_bytes);
    CHECK(legacy_peak.binding_point == via_o2_peak.binding_point);
  }

  SECTION("partially-demoted CSE forest (occ-1-first)") {
    auto make_V = [] {
      return inode_remat("V{o_1;i_1}", leaf_remat("V1{o_1;x_1}"),
                         leaf_remat("V2{x_1;i_1}"));
    };
    auto P1 = inode_remat("P1{i_1;a_1}", make_V(), leaf_remat("W1{i_1;a_1}"));
    stamp_ext_remat(P1, Index{L"i_1"});
    auto P2 = inode_remat("P2{a_1;a_3}", make_V(), leaf_remat("W2{a_1;a_3}"));
    stamp_ext_pair_remat(P2, Index{L"o_1"}, Index{L"i_1"});
    std::vector<EvalNode<EvalExpr>> const forest{P1, P2};

    Schedule const legacy = compute_dag_path(forest, cm, block_of);
    Schedule const via_o2 = remat_round_trip(forest, cm, block_of);
    require_schedules_identical(legacy, via_o2);

    PeakProfile const legacy_peak = peak_profile_sweep(legacy);
    PeakProfile const via_o2_peak = peak_profile_sweep(via_o2);
    CHECK(legacy_peak.peak_bytes == via_o2_peak.peak_bytes);
    CHECK(legacy_peak.binding_point == via_o2_peak.binding_point);
  }

  SECTION("partially-demoted CSE forest (occ-2-first)") {
    auto make_V = [] {
      return inode_remat("V{o_1;i_1}", leaf_remat("V1{o_1;x_1}"),
                         leaf_remat("V2{x_1;i_1}"));
    };
    auto P1 = inode_remat("P1{i_1;a_1}", make_V(), leaf_remat("W1{i_1;a_1}"));
    stamp_ext_remat(P1, Index{L"i_1"});
    auto P2 = inode_remat("P2{a_1;a_3}", make_V(), leaf_remat("W2{a_1;a_3}"));
    stamp_ext_pair_remat(P2, Index{L"o_1"}, Index{L"i_1"});
    std::vector<EvalNode<EvalExpr>> const forest{P2, P1};

    Schedule const legacy = compute_dag_path(forest, cm, block_of);
    Schedule const via_o2 = remat_round_trip(forest, cm, block_of);
    require_schedules_identical(legacy, via_o2);

    PeakProfile const legacy_peak = peak_profile_sweep(legacy);
    PeakProfile const via_o2_peak = peak_profile_sweep(via_o2);
    CHECK(legacy_peak.peak_bytes == via_o2_peak.peak_bytes);
    CHECK(legacy_peak.binding_point == via_o2_peak.binding_point);
  }
}

// All extents 10 (so a single-mode leaf sizes 10*8 = 80 bytes FULL).
namespace {
CostModel unit_regime_cm() {
  SizeRegime r;
  r.space_extent = {{L"o", 10}, {L"i", 10}, {L"x", 10}, {L"a", 10}};
  return CostModel{r};
}
// Block extent 2 (NOT 1): a fully-sliced tensor with all extents == 1 hits
// memsize's scalar special-case (product 1.0 => 0 bytes), which would mask the
// footprint arithmetic. With block 2 a mode sizes 80 FULL (10*8) and 16 BLOCK
// (2*8); a carried {o_1,i_1} cell homed at {i_1} => 2*10*8 = 160, shrink o_1
// => 2*2*8 = 32.
auto const unit_block_of = [](Index const&) -> std::size_t { return 2; };
}  // namespace

TEST_CASE(
    "rematerialize_to_budget shrinks a demoted giant under the threshold "
    "(Feasible)",
    "[placement_remat]") {
  CostModel const cm = unit_regime_cm();
  // One demoted giant: carried {o_1,i_1}, homed at {i_1} (o_1 demoted); o_1 is
  // a shrink candidate (carried + enclosing, not in home). Seed footprint 80.
  sequant::container::svector<ValueCell> cells;
  cells.push_back(mk_cell(0, {Index{L"o_1"}, Index{L"i_1"}}, {Index{L"i_1"}},
                          {Index{L"o_1"}, Index{L"i_1"}}, 0, 0));

  RematResult const res =
      rematerialize_to_budget(cells, cm, unit_block_of, /*num_points=*/1,
                              /*threshold=*/50.0);

  CHECK(res.status == RematStatus::Feasible);
  CHECK(res.profile.peak_bytes == 32.0);  // exact: o_1 shrunk => 2*2*8
  // the giant's home now carries o_1 (it was shrunk into the o_1 loop)
  auto const& hm = res.cells[0].home_modes;
  CHECK(std::find(hm.begin(), hm.end(), Index{L"o_1"}) != hm.end());
}

TEST_CASE(
    "rematerialize_to_budget at threshold = +infinity is a strict no-op (seed "
    "placement)",
    "[placement_remat]") {
  // Drive the seed from a REAL forest (the partial-demotion CSE case), then
  // assert O2 makes no move at an infinite budget.
  CostModel const cm = unit_regime_cm();
  auto make_V = [] {
    return inode_remat("V{o_1;i_1}", leaf_remat("V1{o_1;x_1}"),
                       leaf_remat("V2{x_1;i_1}"));
  };
  auto P1 = inode_remat("P1{i_1;a_1}", make_V(), leaf_remat("W1{i_1;a_1}"));
  stamp_ext_remat(P1, Index{L"i_1"});
  auto P2 = inode_remat("P2{a_1;a_3}", make_V(), leaf_remat("W2{a_1;a_3}"));
  stamp_ext_pair_remat(P2, Index{L"o_1"}, Index{L"i_1"});
  std::vector<EvalNode<EvalExpr>> const forest{P1, P2};

  RematInput const in = remat_cells(forest, cm, unit_block_of);
  double const seed_peak =
      peak_profile_sweep(
          to_schedule(in.cells, cm, unit_block_of, in.num_points))
          .peak_bytes;

  RematResult const res =
      rematerialize_to_budget(in.cells, cm, unit_block_of, in.num_points,
                              std::numeric_limits<double>::infinity());

  CHECK(res.status == RematStatus::Feasible);
  CHECK(res.profile.peak_bytes == seed_peak);
  REQUIRE(res.cells.size() == in.cells.size());
  for (std::size_t i = 0; i < in.cells.size(); ++i)
    CHECK(res.cells[i].home_modes == in.cells[i].home_modes);  // unchanged
}

TEST_CASE(
    "rematerialize_to_budget reports RebatchNeeded for an irreducible giant "
    "(no shrink "
    "candidate) and terminates",
    "[placement_remat]") {
  CostModel const cm = unit_regime_cm();
  // A leaf-like giant: carried {a_1}, no home, no enclosing loops => empty
  // shrink_candidates. Full footprint 80 > threshold, nothing to shrink into.
  sequant::container::svector<ValueCell> cells;
  cells.push_back(mk_cell(0, {Index{L"a_1"}}, {}, {}, 0, 0));

  RematResult const res =
      rematerialize_to_budget(cells, cm, unit_block_of, /*num_points=*/1,
                              /*threshold=*/50.0);

  CHECK(res.status == RematStatus::RebatchNeeded);  // returned, i.e. terminated
}

TEST_CASE(
    "rematerialize_to_budget reports FactorizationInherent when a shrinkable "
    "binding cell "
    "cannot lower a peak tied elsewhere",
    "[placement_remat]") {
  CostModel const cm = unit_regime_cm();
  // Point 0: a shrinkable giant G (160). Point 1: two irreducible leaves
  // La+Lo (80+80 = 160). The peak (160) is TIED across the two points; binding
  // is point 0 (lowest). Shrinking G lowers point 0 (160 -> 32) but leaves
  // point 1 at 160, so NO single move drops the GLOBAL peak (drop == 0). A live
  // cell (G) DID have a shrink candidate => FactorizationInherent (not
  // RebatchNeeded, which requires EVERY live cell to be un-shrinkable).
  sequant::container::svector<ValueCell> cells;
  cells.push_back(mk_cell(0, {Index{L"o_1"}, Index{L"i_1"}}, {Index{L"i_1"}},
                          {Index{L"o_1"}, Index{L"i_1"}}, 0, 0));  // G(160)@p0
  cells.push_back(mk_cell(1, {Index{L"a_1"}}, {}, {}, 1, 1));      // La(80)@p1
  cells.push_back(mk_cell(2, {Index{L"o_1"}}, {}, {}, 1, 1));      // Lo(80)@p1

  RematResult const res =
      rematerialize_to_budget(cells, cm, unit_block_of, /*num_points=*/2,
                              /*threshold=*/50.0);

  CHECK(res.status == RematStatus::FactorizationInherent);
}

// Phase 4a T3: END-TO-END validation on a REALISTIC linearized forest (built
// via the inode_remat/leaf_remat/stamp_ext_remat/stamp_ext_pair_remat
// helpers, not hand-built ValueCells) plus the runtime-witness
// non-regression proof (see the design doc's section 5). T1/T2 above pinned
// the working-cell projection and the greedy loop's four outcomes on
// hand-built cells or on a SMALL forest; this closes Phase 4a by driving the
// SAME remat_cells -> rematerialize_to_budget pipeline over a forest whose
// (demoted) giant dominates the WHOLE forest's peak, so the SHRINK move is
// exercised end to end rather than in isolation.
namespace {

// A cost model where the demoted mode (o) is a THOUSAND elements and every
// other space is pinned AT the block extent (2): any OTHER cell's shrink
// candidate (if it has one -- see below) has a ZERO byte delta, since its
// full extent already equals its block extent, so only the mode that is
// actually huge (o) can ever move the peak. This isolates the giant
// end-to-end without hand-built cells or a bespoke non-uniform block_of.
CostModel giant_regime_cm() {
  SizeRegime r;
  r.space_extent = {{L"o", 1000}, {L"i", 2}, {L"x", 2}, {L"a", 2}};
  return CostModel{r};
}

// Same topology as the partial-demotion CSE forest above (a shared value G
// carries {o_1,i_1}; occurrence 1 sits under a single i_1 loop, occurrence 2
// under BOTH an outer o_1 loop and an inner i_1 loop -- the cross-occurrence
// MEET home is {i_1}, so o_1 stays OUT of the home and is G's SHRINK
// candidate), but G's own two leaf children (A, B) carry ONLY a/x-space
// slots -- disjoint from {o_1,i_1} -- so they inherit the enclosing o_1/i_1
// context (like every descendant of P1/P2) WITHOUT it ever intersecting
// their own carried indices, i.e. they have NO shrink candidates. Combined
// with giant_regime_cm's huge o-space extent, G is the forest's unique
// giant: every other cell's footprint is a constant 32 bytes (2*2*8)
// regardless of any move, while G's is 16000 bytes (2*1000*8) FULL-o,
// dropping to 32 bytes (2*2*8) once o_1 is shrunk into its home -- a ~500x
// swing that dwarfs the rest of the forest at every point in the timeline.
std::vector<EvalNode<EvalExpr>> build_demoted_giant_forest_remat() {
  auto make_G = [] {
    return inode_remat("G{o_1;i_1}", leaf_remat("A{a_5;x_2}"),
                       leaf_remat("B{x_2;a_6}"));
  };
  auto P1 = inode_remat("P1{i_1;a_1}", make_G(), leaf_remat("W1{i_1;a_1}"));
  stamp_ext_remat(P1, Index{L"i_1"});  // occurrence 1: i_1 loop only
  auto P2 = inode_remat("P2{a_1;a_3}", make_G(), leaf_remat("W2{a_1;a_3}"));
  stamp_ext_pair_remat(
      P2, Index{L"o_1"},
      Index{L"i_1"});  // occurrence 2: o_1 (outer), i_1 (inner)
  return {P1, P2};
}

}  // namespace

TEST_CASE(
    "END-TO-END: remat lowers a realistic linearized forest's peak via "
    "SHRINK (demoted giant dominates the whole forest)",
    "[placement_remat]") {
  CostModel const cm = giant_regime_cm();
  std::vector<EvalNode<EvalExpr>> const forest =
      build_demoted_giant_forest_remat();

  RematInput const in = remat_cells(forest, cm, unit_block_of);
  double const seed_peak =
      peak_profile_sweep(
          to_schedule(in.cells, cm, unit_block_of, in.num_points))
          .peak_bytes;

  // G is the unique cell with a non-trivial home_depth (same reasoning as
  // the shrink_candidates test at the top of this file): P1/P2 are forest
  // roots with empty ancestor context, and leaves are never given a seed
  // residency, so only G sits under a parent's batch loop.
  std::size_t giant_id = in.cells.size();
  for (auto const& c : in.cells)
    if (c.home_depth != -1) giant_id = c.value_id;
  REQUIRE(giant_id < in.cells.size());
  CHECK(in.cells[giant_id].home_modes == svector<Index>{Index{L"i_1"}});
  CHECK(shrink_candidates(in.cells[giant_id]) == svector<Index>{Index{L"o_1"}});

  // block-sized-giant-peak (160, every cell block-sized) <= T < seed_peak
  // (16128, G held FULL at its o_1 occurrence): forces exactly the SHRINK
  // this test is pinning.
  double const T = 200.0;
  RematResult const res =
      rematerialize_to_budget(in.cells, cm, unit_block_of, in.num_points, T);

  CHECK(res.status == RematStatus::Feasible);
  CHECK(res.profile.peak_bytes <= T);
  CHECK(res.profile.peak_bytes < seed_peak);  // it actually dropped
  auto const& hm = res.cells[giant_id].home_modes;
  CHECK(std::find(hm.begin(), hm.end(), Index{L"o_1"}) !=
        hm.end());  // the giant was shrunk
}

TEST_CASE(
    "END-TO-END: rematerialize_to_budget at threshold = +infinity is a "
    "strict no-op on the realistic (demoted-giant) forest",
    "[placement_remat]") {
  // T2 pinned the INF no-op on the SMALL partial-demotion forest; this pins
  // it on the LARGER end-to-end forest above (same construction, so the
  // no-op holds even though the seed placement now has a much bigger peak
  // to leave alone).
  CostModel const cm = giant_regime_cm();
  std::vector<EvalNode<EvalExpr>> const forest =
      build_demoted_giant_forest_remat();

  RematInput const in = remat_cells(forest, cm, unit_block_of);
  double const seed_peak =
      peak_profile_sweep(
          to_schedule(in.cells, cm, unit_block_of, in.num_points))
          .peak_bytes;

  RematResult const res =
      rematerialize_to_budget(in.cells, cm, unit_block_of, in.num_points,
                              std::numeric_limits<double>::infinity());

  CHECK(res.status == RematStatus::Feasible);
  CHECK(res.profile.peak_bytes == seed_peak);
  REQUIRE(res.cells.size() == in.cells.size());
  for (std::size_t i = 0; i < in.cells.size(); ++i)
    CHECK(res.cells[i].home_modes == in.cells[i].home_modes);  // unchanged
}

// ---------------------------------------------------------------------
// Phase 4b-2: remat_to_router -- emit PlacementRouter overrides for the MOVED
// cells of a remat result, keyed by occurrence key.
// ---------------------------------------------------------------------

TEST_CASE(
    "remat_to_router emits an override for the moved giant keyed by its "
    "occurrence key, and nothing for unmoved cells (Phase 4b-2)",
    "[placement_remat]") {
  CostModel const cm = giant_regime_cm();
  std::vector<EvalNode<EvalExpr>> const forest =
      build_demoted_giant_forest_remat();

  RematInput const in = remat_cells(forest, cm, unit_block_of);
  RematResult const res = rematerialize_to_budget(in.cells, cm, unit_block_of,
                                                  in.num_points, /*T=*/200.0);
  REQUIRE(res.status == RematStatus::Feasible);

  auto const router = remat_to_router(in.cells, res.cells, forest);
  CHECK_FALSE(router.empty());  // the giant moved -> at least one override

  // The giant cell's final (shrunk) home_modes.
  std::size_t giant_id = in.cells.size();
  for (auto const& c : in.cells)
    if (c.home_depth != -1) giant_id = c.value_id;
  REQUIRE(giant_id < in.cells.size());
  auto const& giant_home = res.cells[giant_id].home_modes;

  // G is P1's left child; under P1 (stamped i_1) its ambient (enclosing) batch
  // context is {i_1}. Its occurrence key must route to a HomeTarget whose
  // DAG-scope == the spaces of the shrunk home_modes (i_1 kept + o_1 added by
  // the shrink), ordered by nest.
  auto const& giant_encl = res.cells[giant_id].enclosing_modes;
  auto const& G = forest[0].left();
  auto const g_key = occurrence_key(G, svector<Index>{Index{L"i_1"}});
  HomeTarget const* g_home = router.route(g_key);
  REQUIRE(g_home != nullptr);
  CHECK(g_home->dag_scope == expected_dag_scope_o2(giant_home, giant_encl));

  // G's SECOND occurrence (under P2, ambient context {o_1,i_1}) also routes to
  // the SAME HomeTarget (one value -> one meet home). NOTE: on THIS fixture the
  // two occurrence keys are IDENTICAL -- but only because the fixture is
  // SYNTHETIC: G's declared result slots {o_1;i_1} are INCONSISTENT with its
  // A{a_5;x_2}*B{x_2;a_6} subtree (which really contracts to {a_5;a_6}), a
  // shortcut that is fine for the mask/footprint machinery (it reads only
  // canon_indices + batched_here) but hides o_1/i_1 from occurrence_key, which
  // is built from the SUBTREE LEAVES. In a REAL eval tree an intermediate's
  // batched result modes are genuine subtree indices, so the two occurrences
  // would get DISTINCT keys -- that CONSISTENT case is the next test; here we
  // just pin that the second occurrence's key still resolves to the shared
  // home.
  auto const& G2 = forest[1].left();
  auto const g_key2 =
      occurrence_key(G2, svector<Index>{Index{L"o_1"}, Index{L"i_1"}});
  CHECK(sequant::eval::RouterKeyEqual{}(g_key, g_key2));  // same key here
  HomeTarget const* g_home2 = router.route(g_key2);
  REQUIRE(g_home2 != nullptr);
  CHECK(g_home2->dag_scope == expected_dag_scope_o2(giant_home, giant_encl));

  // An UNMOVED cell (W1 = P1's right leaf) has NO override.
  auto const& W1 = forest[0].right();
  auto const w_key = occurrence_key(W1, svector<Index>{Index{L"i_1"}});
  CHECK(router.route(w_key) == nullptr);
}

TEST_CASE(
    "remat_to_router returns an EMPTY router when nothing moved (inf "
    "threshold) -- the Phase-2 inert seed seam",
    "[placement_remat]") {
  CostModel const cm = giant_regime_cm();
  std::vector<EvalNode<EvalExpr>> const forest =
      build_demoted_giant_forest_remat();

  RematInput const in = remat_cells(forest, cm, unit_block_of);
  RematResult const res =
      rematerialize_to_budget(in.cells, cm, unit_block_of, in.num_points,
                              std::numeric_limits<double>::infinity());

  auto const router = remat_to_router(in.cells, res.cells, forest);
  CHECK(router.empty());
}

TEST_CASE(
    "remat_to_router: a moved value with TWO DISTINCT occurrence keys emits "
    "two overrides -> the SAME home (partial-overlap / meet-home case)",
    "[placement_remat]") {
  // A CONSISTENT fixture (unlike build_demoted_giant_forest_remat): V's batched
  // result modes i_1,i_2 are GENUINE subtree indices -- V{i_1;i_2} =
  // A{i_1;x_1} * B{x_1;i_2}, contracting x_1 -- so occurrence_key (built from
  // V's leaves) SEES them. V is sliced on {i_1,i_2} under P1 but only {i_1}
  // under P2, so the cross-occurrence MEET home is {i_1} and i_2 is demoted
  // (V's shrink candidate). The two occurrences therefore have DISTINCT
  // occurrence keys ({i_1,i_2} vs {i_1}); once i_2 is shrunk into V's home,
  // remat_to_router emits BOTH keys -- each -> the SAME HomeTarget.
  auto make_V = [] {
    return inode_remat("V{i_1;i_2}", leaf_remat("A{i_1;x_1}"),
                       leaf_remat("B{x_1;i_2}"));
  };
  auto P1 = inode_remat("P1{i_1;i_2}", make_V(), leaf_remat("W1{i_1;i_2}"));
  stamp_ext_pair_remat(P1, Index{L"i_1"},
                       Index{L"i_2"});  // occ1: slices {i_1,i_2}
  auto P2 = inode_remat("P2{i_1;a_1}", make_V(), leaf_remat("W2{i_1;a_1}"));
  stamp_ext_remat(P2, Index{L"i_1"});  // occ2: slices {i_1} only
  std::vector<EvalNode<EvalExpr>> const forest{P1, P2};

  SizeRegime r;
  r.space_extent = {{L"i", 10}, {L"x", 10}, {L"a", 10}};
  CostModel const cm{r};
  auto const block_of = [](Index const&) -> std::size_t { return 2; };

  RematInput const in = remat_cells(forest, cm, block_of);

  // Precondition: V is the unique cell whose meet home is {i_1}, with i_2 as
  // its (demoted) shrink candidate.
  std::size_t v_id = in.cells.size();
  for (auto const& c : in.cells)
    if (c.home_modes == svector<Index>{Index{L"i_1"}}) v_id = c.value_id;
  REQUIRE(v_id < in.cells.size());
  CHECK(shrink_candidates(in.cells[v_id]) == svector<Index>{Index{L"i_2"}});

  // Simulate the shrink of i_2 into V's home (final home = {i_1,i_2}); V moved.
  auto final_cells = in.cells;
  final_cells[v_id].home_modes.push_back(Index{L"i_2"});
  svector<Index> const v_home = final_cells[v_id].home_modes;  // {i_1,i_2}

  auto const router = remat_to_router(in.cells, final_cells, forest);

  // The two occurrences' keys are GENUINELY DISTINCT (the batched modes appear
  // in V's leaves here), yet BOTH route to the SAME home {i_1,i_2}.
  auto const& V1 = forest[0].left();
  auto const& V2 = forest[1].left();
  auto const k1 =
      occurrence_key(V1, svector<Index>{Index{L"i_1"}, Index{L"i_2"}});
  auto const k2 = occurrence_key(V2, svector<Index>{Index{L"i_1"}});
  CHECK_FALSE(sequant::eval::RouterKeyEqual{}(k1, k2));  // genuinely distinct
  HomeTarget const* h1 = router.route(k1);
  HomeTarget const* h2 = router.route(k2);
  REQUIRE(h1 != nullptr);
  REQUIRE(h2 != nullptr);
  auto const& v_encl = final_cells[v_id].enclosing_modes;
  CHECK(h1->dag_scope == expected_dag_scope_o2(v_home, v_encl));
  CHECK(h2->dag_scope == expected_dag_scope_o2(v_home, v_encl));
}

// ---------------------------------------------------------------------
// Phase 4b-3 T1: the router-vs-store-seam CONSISTENCY invariant. After the
// has_demoted_external veto is deleted, place_at_this_level is purely
// router-override-or-seed, so the whole cutover rests on the router (built by
// the MPQC pre-pass via remat_to_router) resolving a moved value's home
// EXACTLY the way place_at_this_level's runtime store seam does. This test
// pins that identity end to end.
// ---------------------------------------------------------------------

TEST_CASE(
    "INVARIANT: the remat router resolves a moved value's home exactly the "
    "way place_at_this_level's runtime store seam does (Phase 4b-3 T1)",
    "[placement_remat]") {
  // Design section 11.4: the router-override branch of place_at_this_level
  // (eval.hpp) derives home as
  //   key = eval::occurrence_key(d, ectx_modes);
  //   if (auto const* home = router->route(key))
  //     rl = router->home_depth(*home, ectx) - 1;
  // where ectx is the LIVE BatchContext and ectx_modes are its modes. This
  // test replicates that seam verbatim against a router the pre-pass built and
  // proves the two derive the SAME home. If they ever diverge the cutover
  // would silently misplace a moved value; this fails the moment they do.
  //
  // Reuse the CONSISTENT distinct-key fixture (V's batched result modes are
  // genuine subtree indices, so occurrence_key sees them): V{i_1;i_2} =
  // A{i_1;x_1} * B{x_1;i_2}, sliced {i_1,i_2} under P1 but only {i_1} under
  // P2, so the cross-occurrence meet home is {i_1} and i_2 is V's (demoted)
  // shrink candidate. A finite threshold shrinks i_2 into V's home -> V moves
  // -> the router emits one override per DISTINCT occurrence, all -> the same
  // shrunk home {i_1,i_2}.
  auto make_V = [] {
    return inode_remat("V{i_1;i_2}", leaf_remat("A{i_1;x_1}"),
                       leaf_remat("B{x_1;i_2}"));
  };
  auto P1 = inode_remat("P1{i_1;i_2}", make_V(), leaf_remat("W1{i_1;i_2}"));
  stamp_ext_pair_remat(P1, Index{L"i_1"}, Index{L"i_2"});  // occ1: {i_1,i_2}
  auto P2 = inode_remat("P2{i_1;a_1}", make_V(), leaf_remat("W2{i_1;a_1}"));
  stamp_ext_remat(P2, Index{L"i_1"});  // occ2: {i_1} only
  std::vector<EvalNode<EvalExpr>> const forest{P1, P2};

  SizeRegime r;
  r.space_extent = {{L"i", 10}, {L"x", 10}, {L"a", 10}};
  CostModel const cm{r};
  auto const block_of = [](Index const&) -> std::size_t { return 2; };

  // Build the router the way the pre-pass will (design section 11.4):
  // remat_cells (the seed) -> the shrink move (i_2 folded into V's home) ->
  // remat_to_router. The shrink is applied directly to V's cell (as the
  // companion distinct-key test does): on THIS consistent fixture V is not the
  // global binding cell, so rematerialize_to_budget reports RebatchNeeded and
  // never performs the move itself; the pre-pass consistency this test pins is
  // in remat_to_router + the store seam, not in the greedy loop's choice.
  RematInput const in = remat_cells(forest, cm, block_of);
  std::size_t v_id = in.cells.size();
  for (auto const& c : in.cells)
    if (c.home_modes == svector<Index>{Index{L"i_1"}}) v_id = c.value_id;
  REQUIRE(v_id < in.cells.size());
  CHECK(shrink_candidates(in.cells[v_id]) == svector<Index>{Index{L"i_2"}});
  auto final_cells = in.cells;
  final_cells[v_id].home_modes.push_back(
      Index{L"i_2"});  // shrink i_2 into home
  auto const router = remat_to_router(in.cells, final_cells, forest);
  REQUIRE_FALSE(router.empty());  // V moved -> at least one override

  using BatchContext =
      sequant::eval::PlacementRouter<EvalNode<EvalExpr>>::BatchContext;

  // --- Occurrence 1: V under P1, live loops {i_1 (outer), i_2 (inner)} ---
  // Replicate the runtime store seam: build the live BatchContext, derive its
  // modes, key the node, route, and resolve the home depth.
  BatchContext const ectx1{{Index{L"i_1"}, {0, 2}}, {Index{L"i_2"}, {0, 2}}};
  svector<Index> const ectx1_modes{Index{L"i_1"}, Index{L"i_2"}};
  auto const key1 = occurrence_key(forest[0].left(), ectx1_modes);  // V/P1
  HomeTarget const* home1 = router.route(key1);
  REQUIRE(home1 != nullptr);
  auto const& v_encl = final_cells[v_id].enclosing_modes;
  CHECK(home1->dag_scope ==
        expected_dag_scope_o2(svector<Index>{Index{L"i_1"}, Index{L"i_2"}},
                              v_encl));
  // V's post-shrink home is {i_1,i_2}: i_2 sits at level 1 of ectx1, so the
  // deepest residency match is level 1 => home_depth 2; the store seam sets
  // rl = home_depth - 1 = 1 (V homed at the inner i_2 loop, built once there).
  CHECK(router.home_depth(*home1, ectx1, key1) == 2);

  // --- Occurrence 2: V under P2, only the i_1 loop is live ---
  BatchContext const ectx2{{Index{L"i_1"}, {0, 2}}};
  svector<Index> const ectx2_modes{Index{L"i_1"}};
  auto const key2 = occurrence_key(forest[1].left(), ectx2_modes);  // V/P2
  HomeTarget const* home2 = router.route(key2);
  REQUIRE(home2 != nullptr);
  // Same value -> same home {i_1,i_2}; under P2 only i_1 is a live loop, so
  // the deepest residency match is level 0 => home_depth 1 (rl = 0).
  CHECK(router.home_depth(*home2, ectx2, key2) == 1);

  // The two occurrence keys are genuinely DISTINCT ({i_1,i_2} vs {i_1}) yet
  // both route to the SAME home -- the meet-home guarantee the runtime relies
  // on, resolved through the identical store-seam call sequence.
  CHECK_FALSE(sequant::eval::RouterKeyEqual{}(key1, key2));
  CHECK(home1->dag_scope == home2->dag_scope);
}
