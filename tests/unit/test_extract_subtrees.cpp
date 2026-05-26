#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/extract_subtrees.hpp>

#include <vector>

namespace {

using Node = sequant::EvalNode<sequant::EvalExpr>;

Node bin(std::wstring_view s) {
  using namespace sequant;
  return binarize(deserialize(s, {.def_perm_symm = Symmetry::Antisymm}));
}

// Example predicate from the header docs: capture maximal internal
// subtrees whose leaves are all non-"t" tensors.
auto non_t_intermediate_pred = [](Node const& n, bool pl, bool pr) -> bool {
  if (n.leaf()) return false;
  auto child_ok = [](Node const& c, bool pc) {
    return c.leaf() ? c->is_tensor() && c->as_tensor().label() != L"t" : pc;
  };
  return child_ok(n.left(), pl) && child_ok(n.right(), pr);
};

}  // namespace

TEST_CASE("opt::extract_subtrees", "[optimize][extract_subtrees]") {
  using namespace sequant;

  SECTION("empty forest yields empty result set") {
    std::vector<Node> forest;
    auto extracted = opt::extract_subtrees(
        forest, [](Node const&, bool, bool) { return true; });
    REQUIRE(extracted.empty());
  }

  SECTION("never-true predicate leaves the forest untouched") {
    std::vector<Node> forest;
    forest.emplace_back(bin(L"g{a1,a2;i1,i2} * t{i1,i2;a1,a2}"));
    REQUIRE_FALSE(forest[0].leaf());

    auto extracted = opt::extract_subtrees(
        forest, [](Node const&, bool, bool) { return false; });
    REQUIRE(extracted.empty());
    REQUIRE_FALSE(forest[0].leaf());
  }

  SECTION("always-true predicate captures each root, collapsing forest items") {
    std::vector<Node> forest;
    forest.emplace_back(bin(L"g{a1;i1} * t{i1;a1}"));
    forest.emplace_back(bin(L"f{a1;a2} * g{a2;a1}"));
    REQUIRE_FALSE(forest[0].leaf());
    REQUIRE_FALSE(forest[1].leaf());

    auto extracted = opt::extract_subtrees(
        forest, [](Node const&, bool, bool) { return true; });
    REQUIRE(extracted.size() == 2);
    REQUIRE(forest[0].leaf());
    REQUIRE(forest[1].leaf());
  }

  SECTION("non-t pred: pure-non-t root collapses, mixed root is pruned") {
    std::vector<Node> forest;
    // mixed: contains a "t" leaf, so the top-level pred is false and
    // the lone non-t sibling is itself a leaf (also pred-false) — nothing
    // gets extracted from this item.
    forest.emplace_back(bin(L"g{a1,a2;i1,i2} * t{i1,i2;a1,a2}"));
    // pure non-t: pred bubbles true to the root → root-capture.
    forest.emplace_back(bin(L"f{a1;a2} * g{a2;a1}"));

    auto extracted = opt::extract_subtrees(forest, non_t_intermediate_pred);

    REQUIRE(extracted.size() == 1);
    REQUIRE_FALSE(forest[0].leaf());  // mixed item retains its internal shape
    REQUIRE(forest[1].leaf());        // pure-non-t root collapsed
  }

  SECTION(
      "captured subtrees under the example pred are always internal nodes") {
    std::vector<Node> forest;
    forest.emplace_back(bin(L"f{a1;a2} * g{a2;a1}"));
    auto extracted = opt::extract_subtrees(forest, non_t_intermediate_pred);
    REQUIRE(extracted.size() == 1);
    REQUIRE_FALSE(extracted.begin()->leaf());
  }

  SECTION(
      "no-recurse: nested positives inside a captured subtree stay inlined") {
    // All-non-t 3-factor product: regardless of how binarize shapes it,
    // every internal node is pred-positive (no "t" anywhere). Only the
    // root is captured; the internal grandchild is *not* a separate entry.
    std::vector<Node> forest;
    forest.emplace_back(bin(L"g{a1;a2} * g{a2;a3} * f{a3;a4}"));
    auto extracted = opt::extract_subtrees(forest, non_t_intermediate_pred);
    REQUIRE(extracted.size() == 1);
  }
}
