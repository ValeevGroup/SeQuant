// Unit tests for SeQuant/core/eval/dag_scope.hpp: DagScopeLevel and
// ModeToLevel::mode_of -- the mode<->DAG-scope-loop map that the runtime
// slicing rework keys the batched evaluator's placement decisions off of.
// Also covers mode_to_level_from_signature (SeQuant/core/eval/
// slicing_signature.hpp): the single point turning "positions of loop axes
// on a node" (a slicing_signature) into "mode->level" (a ModeToLevel).

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/dag_scope.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/eval/slicing_signature.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <catch2/catch_test_macros.hpp>

#include <optional>

namespace {

// Distinctly named (Unity-build safe vs. test_slicing_signature.cpp's
// sig_test_leaf/sig_test_g): a single-tensor leaf eval node.
sequant::eval::dryrun::EvalNodeDryRun dag_leaf(sequant::ExprPtr const& t) {
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = sequant::binarize<sequant::eval::dryrun::EvalExprDryRun>(t);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  return node;
}

sequant::ExprPtr dag_g(
    sequant::container::svector<sequant::Index> const& idxs) {
  return sequant::ex<sequant::Tensor>(L"g", sequant::bra(idxs), sequant::ket{},
                                      sequant::Symmetry::Nonsymm, std::nullopt,
                                      sequant::ColumnSymmetry::Nonsymm);
}

}  // namespace

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

TEST_CASE(
    "mode_to_level_from_signature: zips a slicing_signature into a "
    "ModeToLevel",
    "[dag-scope]") {
  using sequant::DagScopeLevel;
  using sequant::Index;
  using sequant::ModeToLevel;
  using sequant::container::svector;

  // A rank-4 node whose canon_indices are, in order, [i_4, mu (a
  // stand-in for a PNO-like mu-tilde mode), i_1, i_2] -- i_2 (the chosen
  // loop axis) sits at position 3.
  Index i1{L"i_1"}, i2{L"i_2"}, i4{L"i_4"}, mu{L"i_5"};
  auto node = dag_leaf(dag_g({i4, mu, i1, i2}));

  std::size_t const rank = node->canon_indices().size();
  REQUIRE(rank == 4);

  auto sig = sequant::slicing_signature(node, svector<Index>{i2});
  REQUIRE(sig.size() == 1);
  REQUIRE(sig[0].has_value());
  CHECK(*sig[0] == 3);

  DagScopeLevel level_i2{1, L"aux", 0};
  auto m2l = sequant::mode_to_level_from_signature(
      rank, sig, svector<DagScopeLevel>{level_i2});

  REQUIRE(m2l.by_mode.size() == 4);
  CHECK_FALSE(m2l.by_mode[0].has_value());
  CHECK_FALSE(m2l.by_mode[1].has_value());
  CHECK_FALSE(m2l.by_mode[2].has_value());
  REQUIRE(m2l.by_mode[3].has_value());
  CHECK(*m2l.by_mode[3] == level_i2);
}
