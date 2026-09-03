#include <catch2/catch_test_macros.hpp>

#include <SeQuant/core/eval/backends/dryrun/result.hpp>
#include <SeQuant/core/eval/cell_registry.hpp>

using sequant::eval::CellId;
using sequant::eval::CellReadResolver;
using sequant::eval::CellRegistry;
using sequant::eval::CellScope;
using sequant::eval::CellTable;
using sequant::eval::LoopKey;
using sequant::eval::ProductionKind;
using sequant::eval::Read;
using sequant::eval::TableCell;

namespace {
// cell 0: Leaf value 0 (life 2); cell 1: Build value 1 at scope [(1,0)] sliced
// pos 0 on (1,0) (life 1); cell 2: Build value 2 at the same scope reading
// cell 0 whole and cell 1 with a declared slice {0 -> (1,0)}.
CellTable make_table() {
  CellTable t;
  TableCell leaf;
  leaf.value_id = 0;
  leaf.production.kind = ProductionKind::Leaf;
  leaf.persistent = true;
  leaf.life = 2;
  t.cells.push_back(leaf);
  TableCell b1;
  b1.value_id = 1;
  b1.production.kind = ProductionKind::Build;
  b1.scope.path = {{LoopKey{1, 0}, 0}};
  b1.sliced = {{0, LoopKey{1, 0}}};
  b1.life = 1;
  t.cells.push_back(b1);
  TableCell b2;
  b2.value_id = 2;
  b2.production.kind = ProductionKind::Build;
  b2.scope.path = {{LoopKey{1, 0}, 0}};
  b2.sliced = {{0, LoopKey{1, 0}}};
  b2.life = 0;
  t.cells.push_back(b2);
  t.reads.push_back(Read{2, 0, 0, {{0, LoopKey{1, 0}}}, {}});
  t.reads.push_back(Read{2, 1, 1, {}, {}});
  return t;
}

// A regime giving space "i" extent 8 -- large enough that slicing to [2,4) is
// a genuine narrowing, checkable through lobounds_of.
sequant::eval::dryrun::SizeRegime cell_registry_test_regime() {
  sequant::eval::dryrun::SizeRegime r;
  r.space_extent = {{L"i", 8}};
  return r;
}
}  // namespace

TEST_CASE("cell registry: set, read, life and batch clearing",
          "[cell_registry]") {
  auto const t = make_table();
  CellRegistry reg(t);
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(
      cell_registry_test_regime());
  sequant::ResultPtr r = std::make_shared<sequant::eval::dryrun::ResultDryRun>(
      sequant::container::svector<sequant::Index>{sequant::Index{L"i_1"}}, cm);
  reg.set(1, r);
  CHECK(reg.peek(1) == r);
  CHECK(reg.read(1) == r);
  CHECK_THROWS(reg.read(1));  // life exhausted
  reg.set(1, r);              // a new batch's production restores life
  reg.clear_bound_to(LoopKey{1, 0});
  CHECK_FALSE(reg.peek(1));  // bound to (1,0): cleared
  reg.set(0, r);
  reg.clear_bound_to(LoopKey{1, 0});
  CHECK(reg.peek(0) == r);  // whole leaf: untouched
}

TEST_CASE("cell read resolver: declared slice against the batch context",
          "[cell_registry]") {
  auto const t = make_table();
  CellRegistry reg(t);
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(
      cell_registry_test_regime());
  sequant::container::svector<sequant::Index> idx{sequant::Index{L"i_1"}};
  sequant::ResultPtr leaf =
      std::make_shared<sequant::eval::dryrun::ResultDryRun>(idx, cm);
  sequant::ResultPtr b1 = leaf->slice_mode(0, 0, 2);  // already a batch slice
  reg.set(0, leaf);
  reg.set(1, b1);
  // node hash == value id for 0 and 1 (the two operands this fixture reads);
  // any other hash (77 below) is not a value of the table.
  CellReadResolver res(reg, [](std::size_t h) -> std::optional<std::size_t> {
    if (h == 0 || h == 1) return h;
    return std::nullopt;
  });
  res.begin_consumer(2);
  sequant::eval::BatchContext ctx;
  ctx.push_back({sequant::Index{L"i_1"},
                 sequant::eval::DagScopeLevel{1, L"i", 0, 0, 0},
                 {2, 4},
                 std::nullopt});
  auto got0 = res.fetch(0, ctx);
  REQUIRE(got0.has_value());
  // sliced pos 0 to [2,4): extent 2 with lobound 2
  CHECK(sequant::eval::dryrun::detail::lobounds_of(**got0).at(0) == 2);
  auto got1 = res.fetch(1, ctx);
  REQUIRE(got1.has_value());
  CHECK(*got1 == b1);               // whole read: same object
  CHECK_THROWS(res.fetch(1, ctx));  // no remaining Read of value 1 for cell 2
  CHECK_FALSE(res.fetch(77, ctx).has_value());  // not a value: transient
}
