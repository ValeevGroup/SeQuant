#include <SeQuant/core/basis.hpp>
#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>

#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/legality.hpp>
#include <SeQuant/core/eval/ordered_executor.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/optimize/options.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/space_qns.hpp>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <limits>
#include <unordered_map>
#include <utility>

TEST_CASE("build-site-contracted-axis", "[legality]") {
  using sequant::eval::dryrun::EvalExprDryRun;

  // g{i_1;a_3} * t{a_3;i_1}: contracts BOTH i_1 (occupied) and a_3
  // (virtual), fully -- so canon_indices() on the product node is empty and
  // contracted_indices() reports both i_1 and a_3. Mark only a_3 batchable
  // so the assertion distinguishes "in the build-site" from "in the
  // subtree at all".
  auto expr = sequant::deserialize<sequant::ExprPtr>("g{i_1;a_3} * t{a_3;i_1}");
  REQUIRE(static_cast<bool>(expr));

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = sequant::binarize<EvalExprDryRun>(expr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE_FALSE(node.leaf());

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const& ix) {
    return ix.basis().base_key() == L"a";
  };

  // build_site_of now reads the DP's ACTUAL per-node decision (a Contracted
  // node_slice_mask stamp == the node's chosen aprime), NOT policy.is_batchable
  // (legality.hpp) -- so stamp a_3 as the DP would when it batches this
  // contraction here. i_1 is left unstamped (not batched), so it stays out.
  node->set_node_slice_mask(
      {{sequant::Index{L"a_3"}, sequant::BatchModeType::Contracted}});

  auto const site = sequant::eval::build_site_of(node, policy);

  auto const contains = [&](std::wstring const& key) {
    return std::any_of(site.begin(), site.end(), [&](sequant::Index const& ix) {
      return ix.basis().base_key() == key;
    });
  };
  CHECK(contains(L"a"));
  CHECK_FALSE(contains(L"i"));
}

TEST_CASE("classify-axis-per-axis", "[legality]") {
  using sequant::Index;
  using sequant::eval::LoopRole;
  using sequant::eval::OccurrenceRec;
  // value V{i_1;a_1} = sum_{i_2} X{i_1,i_2;a_1} * Y{i_2}: carries i_1 (sliced,
  // lockstep with an enclosing i_1 loop) and contracts i_2 at its own node
  // (sliced there by the DP). Both indices are of the occupied space.
  Index const i1{L"i_1"}, i2{L"i_2"}, a1{L"a_1"};
  sequant::container::svector<Index> const carried{i1, a1};
  sequant::container::svector<Index> const contracted_below{i2};
  OccurrenceRec occ;
  occ.carried = carried;
  occ.home = {i1};  // this occurrence's own home (positional, same frame)
  occ.ectx = {{i1, std::pair<std::size_t, std::size_t>{0, 16}}};
  occ.ectx_opener_point = {1};
  occ.loop_slot = {0, -1};
  sequant::container::svector<OccurrenceRec> const occurrences{occ};
  // the enclosing i_1 loop is the instance (i, 0) that slices position 0
  auto const enclosing = [](OccurrenceRec const&) {
    return sequant::container::svector<sequant::eval::LoopInstance>{{L"i", 0}};
  };
  CHECK(sequant::eval::classify_axis(carried, contracted_below, i1, occurrences,
                                     enclosing) == LoopRole::LoopLocal);
  CHECK(sequant::eval::classify_axis(carried, contracted_below, i2, occurrences,
                                     enclosing) == LoopRole::Reduction);
  // and an axis the value neither carries nor contracts is invariant
  CHECK(sequant::eval::classify_axis(carried, contracted_below, Index{L"i_3"},
                                     occurrences,
                                     enclosing) == LoopRole::LoopInvariant);
}

// ===========================================================================
// classify_axis: a value carrying TWO same-space indices sliced by two loop
// instances, only one of which encloses the occurrence. The per-occurrence
// verdict is per axis and per recorded instance: the position whose instance
// encloses the read is LoopLocal, the other is LoopCarried. Direct unit test
// on classify_axis with a HAND-CONSTRUCTED OccurrenceRec -- no
// forest/optimize pipeline needed.
// ===========================================================================
TEST_CASE("classify-two-same-space", "[legality]") {
  using sequant::Index;
  using sequant::container::svector;
  using sequant::eval::classify_axis;
  using sequant::eval::LoopInstance;
  using sequant::eval::LoopRole;
  using sequant::eval::OccurrenceAxisRole;
  using sequant::eval::OccurrenceRec;
  using sequant::eval::OccurrenceRole;

  auto const i_1 = Index{L"i_1"};
  auto const i_2 = Index{L"i_2"};

  // B{i_1,i_2} = A{;i_1} * A{;i_2}: an outer product carrying two DISTINCT
  // occ indices into its own result, each sliced by its own loop instance
  // (slots 0 and 1), read under a realized loop that is instance 0 only.
  svector<Index> const carried{i_1, i_2};
  svector<Index> const contracted_below{};

  OccurrenceRec occ;
  occ.point = 7;
  occ.consumer_point = 8;
  occ.carried = carried;
  occ.home = {i_1, i_2};
  occ.loop_slot = {0, 1};
  occ.ectx.push_back({i_1, {0u, 4u}});
  occ.ectx_opener_point.push_back(9);
  svector<OccurrenceRec> const occurrences{occ};

  // The enclosing instances, as analyze_legality resolves them through the
  // openers: only (i, 0) is running at this occurrence.
  auto const enclosing = [](OccurrenceRec const&) {
    return svector<LoopInstance>{{L"i", 0}};
  };

  svector<OccurrenceAxisRole> per_occ;
  CHECK(classify_axis(carried, contracted_below, i_1, occurrences, enclosing,
                      &per_occ) == LoopRole::LoopLocal);
  REQUIRE(per_occ.size() == 1);
  CHECK(per_occ.front().point == 7);
  CHECK(per_occ.front().role == OccurrenceRole::Lockstep);

  per_occ.clear();
  CHECK(classify_axis(carried, contracted_below, i_2, occurrences, enclosing,
                      &per_occ) == LoopRole::LoopCarried);
  REQUIRE(per_occ.size() == 1);
  CHECK(per_occ.front().role == OccurrenceRole::Assembled);

  // A batched position with no recorded loop identity is a defect of the
  // record, not a case to guess at.
  OccurrenceRec unstamped = occ;
  unstamped.loop_slot = {0, -1};
  REQUIRE_THROWS_AS(classify_axis(carried, contracted_below, i_2,
                                  svector<OccurrenceRec>{unstamped}, enclosing),
                    sequant::Exception);
}

// ===========================================================================
// classify_axis: a value with two occurrences, one read inside the loop
// instance that slices it (lockstep) and one read where that instance is not
// running (the value opens the loop itself; its consumer reads it whole).
// The value is LoopLocal -- it has a lockstep reader -- and the assembled
// read is recorded per occurrence, so the scheduler also escapes it.
// (Issue #720: the outside occurrence used to be skipped, and the value was
// built inside the loop only.)
// ===========================================================================
TEST_CASE("classify-lockstep-assembled", "[legality]") {
  using sequant::Index;
  using sequant::container::svector;
  using sequant::eval::classify_axis;
  using sequant::eval::LoopInstance;
  using sequant::eval::LoopRole;
  using sequant::eval::OccurrenceAxisRole;
  using sequant::eval::OccurrenceRec;
  using sequant::eval::OccurrenceRole;

  auto const a_2 = Index{L"a_2"};
  svector<Index> const carried{Index{L"i_3"}, Index{L"i_1"}, a_2,
                               Index{L"a_3"}};
  svector<Index> const contracted_below{};

  OccurrenceRec in_loop;
  in_loop.point = 3;
  in_loop.consumer_point = 5;
  in_loop.carried = carried;
  in_loop.home = {a_2};
  in_loop.loop_slot = {-1, -1, 0, -1};
  in_loop.ectx.push_back({a_2, {0u, 4u}});
  in_loop.ectx_opener_point.push_back(6);

  OccurrenceRec at_root = in_loop;
  at_root.point = 11;
  at_root.consumer_point = 13;
  at_root.ectx.clear();
  at_root.ectx_opener_point.clear();
  at_root.opens.push_back({a_2, sequant::BatchModeType::External});

  svector<OccurrenceRec> const occurrences{in_loop, at_root};
  auto const enclosing = [](OccurrenceRec const& occ) {
    return occ.ectx.empty() ? svector<LoopInstance>{}
                            : svector<LoopInstance>{{L"a", 0}};
  };

  svector<OccurrenceAxisRole> per_occ;
  CHECK(classify_axis(carried, contracted_below, a_2, occurrences, enclosing,
                      &per_occ) == LoopRole::LoopLocal);
  REQUIRE(per_occ.size() == 2);
  CHECK(per_occ[0].point == 3);
  CHECK(per_occ[0].role == OccurrenceRole::Lockstep);
  CHECK(per_occ[1].point == 11);
  CHECK(per_occ[1].role == OccurrenceRole::Assembled);

  // Only the outside occurrence: the value is LoopCarried as a whole.
  per_occ.clear();
  CHECK(classify_axis(carried, contracted_below, a_2,
                      svector<OccurrenceRec>{at_root}, enclosing,
                      &per_occ) == LoopRole::LoopCarried);
  REQUIRE(per_occ.size() == 1);
  CHECK(per_occ.front().role == OccurrenceRole::Assembled);
}

// ===========================================================================
// enclosing_instances resolves each enclosing loop through the occurrence
// that opened it and that opener's stamps, and throws on every way the
// record can fail to name the instance -- never falling back to a label.
// ===========================================================================
TEST_CASE("enclosing-instances", "[legality]") {
  using sequant::Index;
  using sequant::container::svector;
  using sequant::eval::enclosing_instances;
  using sequant::eval::LoopInstance;
  using sequant::eval::OccurrenceRec;
  using sequant::eval::PointOccurrence;

  auto const i_1 = Index{L"i_1"}, k_1 = Index{L"a_1"};

  // The opener (point 1): carries i_1 sliced by (i, 2) and opens that loop as
  // External; also opens a Contracted loop over a_1 stamped (a, 0).
  OccurrenceRec opener;
  opener.point = 1;
  opener.consumer_point = 1;
  opener.carried = {i_1};
  opener.home = {i_1};
  opener.loop_slot = {2};
  opener.opens = {{i_1, sequant::BatchModeType::External},
                  {k_1, sequant::BatchModeType::Contracted}};
  opener.reduced_slot = {{k_1, 0}};
  std::unordered_map<std::size_t, PointOccurrence> points;
  points[1] = {0, &opener};

  // An occurrence under both loops names the opener for each.
  OccurrenceRec occ;
  occ.point = 5;
  occ.consumer_point = 1;
  occ.ectx = {{i_1, {0u, 4u}}, {k_1, {0u, 8u}}};
  occ.ectx_opener_point = {1, 1};
  auto const inst = enclosing_instances(occ, points);
  REQUIRE(inst.size() == 2);
  CHECK(inst[0] == LoopInstance{L"i", 2});
  CHECK(inst[1] == LoopInstance{L"a", 0});

  // Fewer opener points than enclosing loops: the record is inconsistent.
  OccurrenceRec short_rec = occ;
  short_rec.ectx_opener_point = {1};
  REQUIRE_THROWS_AS(enclosing_instances(short_rec, points), sequant::Exception);

  // An opener point that is no occurrence.
  OccurrenceRec lost = occ;
  lost.ectx_opener_point = {1, 7};
  REQUIRE_THROWS_AS(enclosing_instances(lost, points), sequant::Exception);

  // An opener that does not open the loop the occurrence says it does.
  OccurrenceRec mislabeled = occ;
  mislabeled.ectx = {{Index{L"i_2"}, {0u, 4u}}};
  mislabeled.ectx_opener_point = {1};
  REQUIRE_THROWS_AS(enclosing_instances(mislabeled, points),
                    sequant::Exception);

  // An opener that opens the loop but carries no stamp for it.
  OccurrenceRec unstamped_opener = opener;
  unstamped_opener.loop_slot = {-1};
  std::unordered_map<std::size_t, PointOccurrence> points2;
  points2[1] = {0, &unstamped_opener};
  OccurrenceRec under_i = occ;
  under_i.ectx = {{i_1, {0u, 4u}}};
  under_i.ectx_opener_point = {1};
  REQUIRE_THROWS_AS(enclosing_instances(under_i, points2), sequant::Exception);
}

// ===========================================================================
// classify_axis (spectator): a value carrying a batched, lockstep-bound index
// (i_1) BESIDE a same-space FULL spectator dimension (i_4 -- NOT batched, no
// enclosing loop) must classify the batched axis LoopLocal. The spectator is
// a retained full dimension, not a loop, so it must not drag the axis to
// LoopCarried. (Regression: w8 aux+occ value [i_4 mu mu i_1], sliced only on
// i_1, was wrongly escaped and evicted before an in-loop consumer read it.)
// ===========================================================================
TEST_CASE("classify-spectator", "[legality]") {
  using sequant::Index;
  using sequant::container::svector;
  using sequant::eval::classify_axis;
  using sequant::eval::LoopRole;
  using sequant::eval::OccurrenceRec;

  auto const i_1 = Index{L"i_1"};
  auto const i_4 = Index{L"i_4"};

  // Value carries i_4 (a full, un-batched occ dimension) and i_1 (batched).
  svector<Index> const carried{i_4, i_1};
  svector<Index> const contracted_below{};

  OccurrenceRec occ;
  occ.point = 0;
  occ.consumer_point = 1;
  occ.carried = carried;
  occ.home = {i_1};
  occ.loop_slot = {-1, 0};  // i_4 has no loop identity; i_1 is instance (i, 0)
  occ.ectx.push_back({i_1, {0u, 4u}});  // only i_1 has a realized loop
  occ.ectx_opener_point.push_back(2);

  svector<OccurrenceRec> const occurrences{occ};

  // Only i_1 is a batched (sliced) mode, by loop instance (i, 0); i_4 is a
  // full spectator dimension with no loop identity.
  auto const enclosing = [](OccurrenceRec const&) {
    return svector<sequant::eval::LoopInstance>{{L"i", 0}};
  };

  // i_1 is lockstep with its enclosing loop instance, so the axis is
  // LoopLocal -- the value is homed inside the i_1 loop, carrying full i_4
  // through. (i_4 is not a build-site axis and is never classified.)
  CHECK(classify_axis(carried, contracted_below, i_1, occurrences, enclosing) ==
        LoopRole::LoopLocal);
}

// ===========================================================================
// forced_split_types groups CellLegality::forced_split_axes by
// axis TYPE (IndexBasis::base_key()), collapsing multiple same-space Index
// INSTANCES (e.g. the cross-iteration fixture's i_3/i_4 outer product, both
// LoopCarried on occ) into the single loop (axis) they jointly force to
// split.
// ===========================================================================
TEST_CASE("forced-split-types-collapse", "[legality]") {
  using sequant::Index;
  using sequant::eval::CellLegality;
  using sequant::eval::forced_split_types;

  CellLegality cl;
  cl.hash = 0;
  // Two DISTINCT occ Index instances, both LoopCarried on the same space --
  // the per-instance record the synthetic B{i_3,i_4} = A{;i_3} * A{;i_4}
  // outer product produces.
  cl.forced_split_axes = {Index{L"i_3"}, Index{L"i_4"}};

  auto const types = forced_split_types(cl);
  REQUIRE(types.size() == 1);
  CHECK(types.front().basis().base_key() == L"i");

  // A cell with no forced splits collapses to an empty list.
  CellLegality const empty_cl;
  CHECK(forced_split_types(empty_cl).empty());

  // Two DIFFERENT-space LoopCarried axes (occ + virtual) must NOT collapse --
  // forced_split_types groups WITHIN a space type only. ("i"/"a" are present
  // in the default context without any scoped registry setup, unlike the DF
  // aux space "Κ", which needs a scoped registry.)
  CellLegality mixed_cl;
  mixed_cl.forced_split_axes = {Index{L"i_3"}, Index{L"i_4"}, Index{L"a_1"}};
  auto const mixed_types = forced_split_types(mixed_cl);
  REQUIRE(mixed_types.size() == 2);
}

TEST_CASE("forced_split_types keeps a named basis apart from its space",
          "[legality][basis]") {
  using namespace sequant;
  auto ibr = mbpt::make_min_sr_spaces();
  const IndexSpace a = ibr->retrieve_space(L"a"), i = ibr->retrieve_space(L"i");
  ibr->add(L"μ̃", a, std::numeric_limits<IndexBasis::instance_type>::max());
  auto ctx = set_scoped_default_context(
      Context({.index_basis_registry_shared_ptr = std::move(ibr),
               .vacuum = Vacuum::SingleProduct}));
  const Index pao(get_default_context().index_basis_registry()->retrieve(L"μ̃"),
                  1);
  const Index composite(a, 2, {Index(i, 1), Index(i, 2)});
  eval::CellLegality cell;
  cell.forced_split_axes = {composite, pao,
                            Index(a, 3, {Index(i, 1), Index(i, 3)})};
  const auto types = eval::forced_split_types(cell);
  REQUIRE(types.size() == 2);  // a (the composites) and μ̃
  CHECK(types[0] == composite);
  CHECK(types[1] == pao);
}
