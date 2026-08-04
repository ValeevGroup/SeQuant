// Phase 4a T1: tests for O2's standalone working representation
// (SeQuant/core/eval/placement_o2.hpp) plus the peak_profile.hpp refactor
// that exposes it (linearize_rich / ValueCell / RichSchedule).
//
// Two things are pinned here:
//   - shrink_candidates finds the demoted carried mode of a partially-homed
//     value (reusing test_peak_profile.cpp's PARTIAL-DEMOTION forest), and
//     to_schedule's projected footprint matches the Phase-3b value.
//   - GUARDRAIL: linearize's flat Schedule (built via the NEW
//     linearize_rich + thin projection) is BYTE-IDENTICAL, cell for cell, to
//     what o2_cells + to_schedule (the O2 round trip over the SAME rich
//     record) produce -- proof the rich split changed nothing observable.
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
#include <SeQuant/core/eval/placement_o2.hpp>
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
using sequant::eval::o2_cells;
using sequant::eval::O2Cell;
using sequant::eval::O2Input;
using sequant::eval::peak_profile_sweep;
using sequant::eval::PeakProfile;
using sequant::eval::Schedule;
using sequant::eval::shrink_candidates;
using sequant::eval::to_schedule;
using sequant::eval::dryrun::CostModel;
using sequant::eval::dryrun::EvalExprDryRun;
using sequant::eval::dryrun::EvalNodeDryRun;
using sequant::eval::dryrun::SizeRegime;

// Build an EvalExpr from a single-tensor spec (e.g. "R{i_1;a_5}"); its
// canon_indices are exactly the tensor's bra+ket slots. Mirrors the helper in
// test_peak_profile.cpp / test_lifetime_mask.cpp. Suffixed `_o2` (rather than
// reusing the exact same name) because a Unity build unifies this file with
// test_peak_profile.cpp into one translation unit, where an anonymous
// namespace does not prevent a same-named free function from colliding.
EvalExpr eval_tensor_o2(std::string_view tensor) {
  auto expr = sequant::deserialize<ExprPtr>(std::string(tensor));
  REQUIRE(static_cast<bool>(expr));
  return EvalExpr{expr->as<sequant::Tensor>()};
}

// A leaf eval node carrying the given tensor's slots.
EvalNode<EvalExpr> leaf_o2(std::string_view tensor) {
  return EvalNode<EvalExpr>{eval_tensor_o2(tensor)};
}

// An internal eval node whose OWN result slots are the given tensor's slots,
// with the two supplied child subtrees.
EvalNode<EvalExpr> inode_o2(std::string_view result, EvalNode<EvalExpr> l,
                            EvalNode<EvalExpr> r) {
  return EvalNode<EvalExpr>{eval_tensor_o2(result), std::move(l), std::move(r)};
}

// Stamp a single External batch loop mode at a node.
void stamp_ext_o2(EvalNode<EvalExpr>& n, Index ix) {
  n->set_batched_here({{std::move(ix), BatchModeType::External}});
}

// Stamp TWO External batch loop modes at a node, outer first: realizes a
// two-level nest (ix_outer at level 0, ix_inner at level 1) in ONE node
// rather than needing a separate ancestor per level.
void stamp_ext_pair_o2(EvalNode<EvalExpr>& n, Index ix_outer, Index ix_inner) {
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

// The O2 round trip over the SAME forest: o2_cells (linearize_rich mapped to
// O2Cell) immediately projected back via to_schedule, with NO shrink move
// applied. Used only to prove the rich split is faithful (guardrail) -- no
// production caller does this.
template <typename Forest, typename BlockOfFn>
Schedule o2_round_trip(Forest const& forest, CostModel const& cm,
                       BlockOfFn const& block_of) {
  O2Input const in = o2_cells(forest, cm, block_of);
  return to_schedule(in.cells, cm, block_of, in.num_points);
}

}  // namespace

TEST_CASE(
    "shrink_candidates finds the demoted carried mode of a "
    "PARTIALLY-demoted value",
    "[placement_o2]") {
  // Reuses test_peak_profile.cpp's exact PARTIAL-DEMOTION forest: a shared
  // value V carries TWO slots {o_1,i_1}. Occurrence 1 sits under a single
  // i_1 loop; occurrence 2 sits under BOTH an outer o_1 loop and an inner
  // i_1 loop. The cross-occurrence MEET home is {i_1} -- o_1 is batched in
  // only ONE occurrence, so it stays OUT of the home and must be a SHRINK
  // candidate: it is carried, and it DOES enclose V at occurrence 2, but it
  // is not (yet) part of V's home.
  auto make_V = [] {
    return inode_o2("V{o_1;i_1}", leaf_o2("V1{o_1;x_1}"),
                    leaf_o2("V2{x_1;i_1}"));
  };

  SizeRegime r;
  r.space_extent = {{L"o", 10}, {L"i", 10}, {L"x", 10}, {L"a", 10}};
  CostModel const cm{r};
  auto const block_of = [](Index const&) -> std::size_t { return 2; };
  // Same value the Phase-3b test pins: i_1 BLOCK (in the meet home), o_1
  // FULL (not in the meet) => 2 (block i_1) * 10 (full o_1) * 8 bytes/elem.
  std::size_t const expected_footprint = 160;

  auto P1 = inode_o2("P1{i_1;a_1}", make_V(), leaf_o2("W1{i_1;a_1}"));
  stamp_ext_o2(P1, Index{L"i_1"});  // occurrence 1: i_1 loop only
  auto P2 = inode_o2("P2{a_1;a_3}", make_V(), leaf_o2("W2{a_1;a_3}"));
  stamp_ext_pair_o2(P2, Index{L"o_1"},
                    Index{L"i_1"});  // occurrence 2: o_1 (outer), i_1 (inner)

  std::vector<EvalNode<EvalExpr>> const forest{P1, P2};

  O2Input const in = o2_cells(forest, cm, block_of);

  // V is the unique cell with a non-trivial home_depth (same reasoning as
  // the Phase-3b test: P1/P2 are forest roots with empty ectx, every leaf
  // stays unstamped, only V sits under a parent loop).
  O2Cell const* v = nullptr;
  for (auto const& c : in.cells)
    if (c.home_depth != -1) v = &c;
  REQUIRE(v != nullptr);

  CHECK(v->home_modes == svector<Index>{Index{L"i_1"}});
  CHECK(shrink_candidates(*v) == svector<Index>{Index{L"o_1"}});

  Schedule const sched = to_schedule(in.cells, cm, block_of, in.num_points);
  CHECK(sched.cells[v->value_id].footprint == expected_footprint);
}

TEST_CASE(
    "GUARDRAIL: linearize's flat Schedule == o2_cells + to_schedule's "
    "round trip, cell for cell",
    "[placement_o2]") {
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

    Schedule const legacy = linearize(forest, cm, block_of);
    Schedule const via_o2 = o2_round_trip(forest, cm, block_of);
    require_schedules_identical(legacy, via_o2);

    PeakProfile const legacy_peak = peak_profile_sweep(legacy);
    PeakProfile const via_o2_peak = peak_profile_sweep(via_o2);
    CHECK(legacy_peak.peak_bytes == via_o2_peak.peak_bytes);
    CHECK(legacy_peak.binding_point == via_o2_peak.binding_point);
  }

  SECTION("partially-demoted CSE forest (occ-1-first)") {
    auto make_V = [] {
      return inode_o2("V{o_1;i_1}", leaf_o2("V1{o_1;x_1}"),
                      leaf_o2("V2{x_1;i_1}"));
    };
    auto P1 = inode_o2("P1{i_1;a_1}", make_V(), leaf_o2("W1{i_1;a_1}"));
    stamp_ext_o2(P1, Index{L"i_1"});
    auto P2 = inode_o2("P2{a_1;a_3}", make_V(), leaf_o2("W2{a_1;a_3}"));
    stamp_ext_pair_o2(P2, Index{L"o_1"}, Index{L"i_1"});
    std::vector<EvalNode<EvalExpr>> const forest{P1, P2};

    Schedule const legacy = linearize(forest, cm, block_of);
    Schedule const via_o2 = o2_round_trip(forest, cm, block_of);
    require_schedules_identical(legacy, via_o2);

    PeakProfile const legacy_peak = peak_profile_sweep(legacy);
    PeakProfile const via_o2_peak = peak_profile_sweep(via_o2);
    CHECK(legacy_peak.peak_bytes == via_o2_peak.peak_bytes);
    CHECK(legacy_peak.binding_point == via_o2_peak.binding_point);
  }

  SECTION("partially-demoted CSE forest (occ-2-first)") {
    auto make_V = [] {
      return inode_o2("V{o_1;i_1}", leaf_o2("V1{o_1;x_1}"),
                      leaf_o2("V2{x_1;i_1}"));
    };
    auto P1 = inode_o2("P1{i_1;a_1}", make_V(), leaf_o2("W1{i_1;a_1}"));
    stamp_ext_o2(P1, Index{L"i_1"});
    auto P2 = inode_o2("P2{a_1;a_3}", make_V(), leaf_o2("W2{a_1;a_3}"));
    stamp_ext_pair_o2(P2, Index{L"o_1"}, Index{L"i_1"});
    std::vector<EvalNode<EvalExpr>> const forest{P2, P1};

    Schedule const legacy = linearize(forest, cm, block_of);
    Schedule const via_o2 = o2_round_trip(forest, cm, block_of);
    require_schedules_identical(legacy, via_o2);

    PeakProfile const legacy_peak = peak_profile_sweep(legacy);
    PeakProfile const via_o2_peak = peak_profile_sweep(via_o2);
    CHECK(legacy_peak.peak_bytes == via_o2_peak.peak_bytes);
    CHECK(legacy_peak.binding_point == via_o2_peak.binding_point);
  }
}
