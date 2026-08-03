// Phase 2 T2: unit tests for the PlacementRouter container --
// set_override/route/empty and the home_depth resolver. Pure additions with
// no runtime wiring (see occurrence_key.hpp / cache_manager.hpp).
//
// Phase 2 T3 adds the shadow-assert fixture below: SEQUANT_ROUTER_SHADOW
// compiles in a dev-only correctness check at the eval Enter-stage read seam
// (see eval.hpp), so this file must define it BEFORE eval.hpp's first
// (possibly transitive) inclusion in this translation unit -- header guards
// make a later #define a no-op. This file is excluded from the Unity build
// for exactly that reason (tests/unit/CMakeLists.txt): unifying it with
// test_eval_dryrun.cpp (which includes eval.hpp without the macro) would
// silently disable the very code path this file means to exercise.
#define SEQUANT_ROUTER_SHADOW 1

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/occurrence_key.hpp>
#include <SeQuant/core/eval/placement_router.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace {

using sequant::bra;
using sequant::CacheManager;
using sequant::ColumnSymmetry;
using sequant::ex;
using sequant::ExprPtr;
using sequant::Index;
using sequant::ket;
using sequant::ResultPtr;
using sequant::Symmetry;
using sequant::Tensor;
using sequant::container::svector;
using sequant::eval::HomeTarget;
using sequant::eval::occurrence_key;
using sequant::eval::PlacementRouter;
using sequant::eval::dryrun::CostModel;
using sequant::eval::dryrun::DryRunLeafEvaluator;
using sequant::eval::dryrun::EvalExprDryRun;
using sequant::eval::dryrun::EvalNodeDryRun;
using sequant::eval::dryrun::SizeRegime;

using router_type = PlacementRouter<EvalNodeDryRun>;
using ctx_type = router_type::BatchContext;

/// Builds a trivial (single-leaf) DryRun eval node wrapping \p t. Named
/// distinctly from test_occurrence_key.cpp's identical helper: both files
/// land in the same Unity translation unit, where two same-named functions
/// in file-local anonymous namespaces would collide.
EvalNodeDryRun router_test_leaf_node(ExprPtr const& t) {
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = sequant::binarize<EvalExprDryRun>(t);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  return node;
}

}  // namespace

TEST_CASE("placement_router: home_depth resolves the deepest matching scope",
          "[placement_router]") {
  // Outermost-first: [J, K, L] at ctx indices [0, 1, 2].
  Index J{L"i_1"}, K{L"i_2"}, Lx{L"i_3"}, M{L"i_9"};  // M absent from ctx

  ctx_type ctx{{J, {0, 1}}, {K, {0, 1}}, {Lx, {0, 1}}};

  router_type router;

  // K sits at ctx index 1 -> home_depth = 1 + 1 = 2.
  CHECK(router.home_depth(HomeTarget{svector<Index>{K}}, ctx) == 2);

  // Empty residency => invariant to the whole nest => chain root => 0.
  CHECK(router.home_depth(HomeTarget{}, ctx) == 0);

  // A residency mode absent from ctx => 0.
  CHECK(router.home_depth(HomeTarget{svector<Index>{M}}, ctx) == 0);

  // Residency naming more than one ctx mode resolves to the DEEPEST match.
  CHECK(router.home_depth(HomeTarget{svector<Index>{J, Lx}}, ctx) == 3);
}

TEST_CASE("placement_router: set_override/route/empty", "[placement_router]") {
  Index i1{L"i_1"}, i2{L"i_2"};
  auto t1 =
      ex<Tensor>(L"A", bra(svector<Index>{i1, i2}), ket{}, Symmetry::Nonsymm,
                 std::nullopt, ColumnSymmetry::Nonsymm);
  auto t2 =
      ex<Tensor>(L"B", bra(svector<Index>{i1, i2}), ket{}, Symmetry::Nonsymm,
                 std::nullopt, ColumnSymmetry::Nonsymm);

  auto n1 = router_test_leaf_node(t1);
  auto n2 = router_test_leaf_node(t2);

  svector<Index> ctx_modes{i1};
  auto key1 = occurrence_key(n1, ctx_modes);
  auto key2 = occurrence_key(n2, ctx_modes);

  router_type router;
  CHECK(router.empty());
  CHECK(router.route(key1) == nullptr);

  HomeTarget home{svector<Index>{i1}};
  router.set_override(key1, home);

  CHECK_FALSE(router.empty());
  auto const* routed = router.route(key1);
  REQUIRE(routed != nullptr);
  CHECK(routed->residency == home.residency);
  CHECK(routed->split_index == 0);

  // A different occurrence (distinct tensor label) stays unrouted.
  CHECK(router.route(key2) == nullptr);
}

// Phase 2 T3 shadow-assert fixture. Injects a NO-OP override -- the
// registered HomeTarget's residency (empty: invariant to the whole nest)
// resolves, via home_depth() against this (unbatched, standalone) cache's
// empty batch_context, to the value's ACTUAL current home: this same,
// only, cache level (hops == 0). This exercises eval.hpp's Enter-stage read
// seam end to end (via a real evaluate() call, not a direct CacheManager
// poke), proving access_at_hops() + home_depth() reproduce access_at()
// pointer-for-pointer BEFORE Phase 4 lands a real relocation.
//
// The proof is independent of SEQUANT_ASSERT's build-time behavior (THROW/
// ABORT/IGNORE -- see CMakeLists.txt's SEQUANT_ASSERT_BEHAVIOR, IGNORE by
// default outside Debug): the entry's life-count decrement is the
// discriminator. A registered entry's access() decays its life by exactly
// one per call. The router branch, on a hit, calls access_at_hops() once;
// the SEQUANT_ROUTER_SHADOW block then calls access_at() a SECOND time on
// the SAME entry. So a life decrement of 2 (not 0 or 1) on the second
// evaluate() call below is only possible if BOTH the routed fetch and the
// shadow comparison actually ran -- i.e. the router fired and the shadow
// landing executed (and, since the call did not throw/abort, its assertion
// held). Combined with the direct pointer-identity CHECK on evaluate()'s own
// return value (unaffected by the life bookkeeping), this pins the invariant
// regardless of how SEQUANT_ASSERT is configured to behave.
TEST_CASE(
    "placement_router: shadow-assert reproduces access_at for a no-op "
    "root-home override",
    "[placement_router]") {
  Index i1{L"i_1"};
  auto t = ex<Tensor>(L"A", bra(svector<Index>{i1}), ket{}, Symmetry::Nonsymm,
                      std::nullopt, ColumnSymmetry::Nonsymm);
  auto node = router_test_leaf_node(t);  // a single-leaf node IS a valid
                                         // cache key -- evaluate()'s Enter
                                         // stage checked-cache gate runs
                                         // ahead of the leaf/internal split.

  // Register `node` directly (bypassing the CSE-repeat-count factories) with
  // a generous life count, so the two evaluate() calls below (1 decay, then
  // 2 more) never drain it to a release.
  using Cache = CacheManager<EvalNodeDryRun, false>;
  std::vector<std::pair<EvalNodeDryRun, std::size_t>> decaying{{node, 10}};
  Cache cache{decaying};
  CHECK(cache.life(node) == 10);

  // The NO-OP override: empty residency resolves (home_depth, against the
  // empty batch_context this standalone cache always has) to depth 0 --
  // exactly where `node`, once stored, actually lives.
  router_type router;
  svector<Index> const ctx_modes{};  // no enclosing batch loops
  auto const key = occurrence_key(node, ctx_modes);
  router.set_override(key, HomeTarget{});
  CHECK_FALSE(router.empty());
  cache.set_placement_router(&router);

  SizeRegime regime;
  regime.space_extent = {{L"i", 4}};
  auto cm = std::make_shared<CostModel const>(regime);
  DryRunLeafEvaluator const leaf{cm};

  // First call: `node` is registered but not yet stored, so both the routed
  // fetch (a miss: the entry exists but holds no data) and the plain
  // access_at() fallback below it are misses -- misses never decay a life
  // count (CacheManager::entry::access() checks `!data_p` before decaying).
  // This call therefore builds `node` via the leaf evaluator and stores it --
  // CacheManager::store() itself IMPLICITLY ACCESSES the just-stored entry to
  // hand back a pointer (see its doc comment), so THIS is where the first
  // real decay happens: one, not zero.
  ResultPtr const r1 = sequant::evaluate(node, leaf, cache);
  REQUIRE(r1);
  CHECK(cache.life(node) == 9);

  // Second call: `node` is now stored, so the router's route() hit resolves
  // to a REAL hit this time (home_depth 0 == hops 0 against the still-empty
  // batch_context): access_at_hops() decays life once, then -- because
  // SEQUANT_ROUTER_SHADOW is defined for this TU -- the shadow landing calls
  // access_at() on the SAME entry, decaying it a second time, and asserts the
  // two pointers match.
  ResultPtr const r2 = sequant::evaluate(node, leaf, cache);
  REQUIRE(r2);
  CHECK(r2.get() == r1.get());
  CHECK(cache.life(node) == 7);
}
