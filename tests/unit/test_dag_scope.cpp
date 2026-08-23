// Unit tests for SeQuant/core/eval/dag_scope.hpp: DagScopeLevel and
// ModeToLevel::mode_of -- the mode<->DAG-scope-loop map that the runtime
// slicing rework keys the batched evaluator's placement decisions off of.

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/dag_scope.hpp>

#include <catch2/catch_test_macros.hpp>

#include <optional>

TEST_CASE("DagScopeLevel + ModeToLevel: mode_of round-trips", "[dag-scope]") {
  using sequant::DagScopeLevel;
  using sequant::ModeToLevel;
  using sequant::container::svector;

  DagScopeLevel level0{0, L"occ", 0};
  DagScopeLevel level1{1, L"aux", 0};
  DagScopeLevel level_absent{2, L"virt", 0};

  ModeToLevel m2l{svector<std::optional<DagScopeLevel>>{std::nullopt, level0,
                                                        std::nullopt, level1}};

  CHECK(m2l.mode_of(level0) == 1);
  CHECK(m2l.mode_of(level1) == 3);
  CHECK(m2l.mode_of(level_absent) == std::nullopt);

  // round-trip: the mode found for a level maps back to an equal level.
  REQUIRE(m2l.mode_of(level0).has_value());
  CHECK(m2l.by_mode[*m2l.mode_of(level0)] == level0);
  REQUIRE(m2l.mode_of(level1).has_value());
  CHECK(m2l.by_mode[*m2l.mode_of(level1)] == level1);
}
