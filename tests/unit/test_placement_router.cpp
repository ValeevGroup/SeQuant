// Phase 2 T2: unit tests for the PlacementRouter container --
// set_override/route/empty and the home_depth resolver. Pure additions with
// no runtime wiring (see occurrence_key.hpp / cache_manager.hpp); the eval
// read/store paths do not yet consult this router (Phase 2 T3).

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/occurrence_key.hpp>
#include <SeQuant/core/eval/placement_router.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <catch2/catch_test_macros.hpp>

#include <optional>

namespace {

using sequant::bra;
using sequant::ColumnSymmetry;
using sequant::ex;
using sequant::ExprPtr;
using sequant::Index;
using sequant::ket;
using sequant::Symmetry;
using sequant::Tensor;
using sequant::container::svector;
using sequant::eval::HomeTarget;
using sequant::eval::occurrence_key;
using sequant::eval::PlacementRouter;
using sequant::eval::dryrun::EvalExprDryRun;
using sequant::eval::dryrun::EvalNodeDryRun;

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
