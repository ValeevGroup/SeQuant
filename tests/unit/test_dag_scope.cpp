// Unit tests for SeQuant/core/eval/dag_scope.hpp: DagScopeLevel and
// ModeToLevel -- the mode<->DAG-scope-loop map that the runtime
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

#include <functional>
#include <optional>
#include <unordered_map>

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

TEST_CASE("LoopKey keys a map on the FULL loop identity, unpacked",
          "[dag-scope]") {
  // Copilot review (PR #613): LoopKey used to expose a `color()` that packed
  // `depth` and `loop_slot` into one size_t with a 12-bit slot field, guarded
  // by nothing but a comment ("loop_slot < 4096"); the numbering site
  // (peak_profile.hpp) hands out slots from an UNBOUNDED per-space counter, so
  // a large enough schedule would have aliased two distinct loops onto one
  // color and given them one batch count. The packing is gone: LoopKey is now
  // hashed and compared as the pair it is, which has no bound at all.
  using sequant::LoopKey;

  std::unordered_map<LoopKey, int> m;
  m[LoopKey{1, 0}] = 10;
  m[LoopKey{1, 1}] = 11;  // same group, sibling slot: a DISTINCT loop
  m[LoopKey{0, 1}] = 1;
  CHECK(m.size() == 3);
  CHECK(m.at(LoopKey{1, 0}) == 10);
  CHECK(m.at(LoopKey{1, 1}) == 11);
  CHECK(m.at(LoopKey{0, 1}) == 1);

  // The pair that the old 12-bit packing aliased: depth 1 / slot 4096 packed
  // to (1 << 12) | 4096 == 8192, exactly what depth 2 / slot 0 packed to.
  // They are distinct keys now.
  LoopKey const a{1, 4096}, b{2, 0};
  CHECK(a != b);
  m[a] = 40;
  m[b] = 41;
  CHECK(m.size() == 5);
  CHECK(m.at(a) == 40);
  CHECK(m.at(b) == 41);

  // Equality is on both fields; the hash agrees with it.
  CHECK(LoopKey{3, 2} == LoopKey{3, 2});
  CHECK(std::hash<LoopKey>{}(LoopKey{3, 2}) ==
        std::hash<LoopKey>{}(LoopKey{3, 2}));
}
