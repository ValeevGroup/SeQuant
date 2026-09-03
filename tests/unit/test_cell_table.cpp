#include <catch2/catch_test_macros.hpp>

#include <SeQuant/core/eval/cell_table.hpp>

using sequant::eval::AssembleKind;
using sequant::eval::CellScope;
using sequant::eval::CellTable;
using sequant::eval::CellViolation;
using sequant::eval::LoopKey;
using sequant::eval::ProductionKind;
using sequant::eval::Read;
using sequant::eval::TableCell;
using sequant::eval::validate_cell_table;

TEST_CASE("cell table: scope prefix and equality", "[cell_table]") {
  CellScope root;
  CellScope outer;
  outer.path.push_back({LoopKey{1, 0}, 0});
  CellScope inner = outer;
  inner.path.push_back({LoopKey{2, 1}, 0});
  CellScope inner_other_pass = outer;
  inner_other_pass.path.push_back({LoopKey{2, 1}, 1});

  CHECK(root.encloses(root));
  CHECK(root.encloses(outer));
  CHECK(root.encloses(inner));
  CHECK(outer.encloses(inner));
  CHECK_FALSE(inner.encloses(outer));
  CHECK_FALSE(inner.encloses(inner_other_pass));  // same loop, other pass
  CHECK(inner == inner);
  CHECK_FALSE(inner == inner_other_pass);
}

TEST_CASE("cell table: default cell and read", "[cell_table]") {
  TableCell c;
  c.value_id = 7;
  c.production.kind = ProductionKind::Build;
  CHECK(c.sliced.empty());
  CHECK(c.scope.path.empty());
  CHECK(c.life == 0);
  CHECK_FALSE(c.persistent);
  CHECK_FALSE(c.produce_if_absent);
  Read r{/*consumer=*/0, /*operand_value_id=*/7, /*source=*/0, {}};
  CHECK(r.slice.empty());
  CellTable t;
  t.cells.push_back(c);
  t.reads.push_back(r);
  CHECK(t.cells.size() == 1);
  CHECK(t.reads.size() == 1);
}

namespace {
// value 0 = leaf at root; value 1 = built inside loop (1,0); value 2 = built
// at root, reads 1's assembled form (cell 3).
sequant::eval::CellTable make_two_level_table() {
  using namespace sequant::eval;
  CellTable t;
  TableCell leaf;
  leaf.value_id = 0;
  leaf.production.kind = ProductionKind::Leaf;
  leaf.life = 1;
  t.cells.push_back(leaf);  // cell 0
  TableCell inner;
  inner.value_id = 1;
  inner.sliced.push_back({0, LoopKey{1, 0}});
  inner.scope.path.push_back({LoopKey{1, 0}, 0});
  inner.production.kind = ProductionKind::Build;
  inner.life = 1;            // read once, by the assemble
  t.cells.push_back(inner);  // cell 1
  TableCell assembled;
  assembled.value_id = 1;
  assembled.production.kind = ProductionKind::Assemble;
  assembled.production.assemble = AssembleKind::Scatter;
  assembled.production.source = 1;
  assembled.production.scatter_map.push_back({0, LoopKey{1, 0}});
  assembled.life = 1;
  t.cells.push_back(assembled);  // cell 2
  TableCell root_consumer;
  root_consumer.value_id = 2;
  root_consumer.production.kind = ProductionKind::Build;
  root_consumer.life = 0;
  t.cells.push_back(root_consumer);  // cell 3
  t.reads.push_back(
      Read{/*consumer=*/1, /*operand_value_id=*/0, /*source=*/0, {}});
  t.reads.push_back(
      Read{/*consumer=*/3, /*operand_value_id=*/1, /*source=*/2, {}});
  return t;
}
sequant::eval::ScopeBlock empty_root() { return sequant::eval::ScopeBlock{}; }
}  // namespace

TEST_CASE("cell table validator: a consistent two-level table is valid",
          "[cell_table]") {
  auto const t = make_two_level_table();
  auto const v = validate_cell_table(t, empty_root());
  for (auto const& x : v) UNSCOPED_INFO(x.rule << ": " << x.what);
  CHECK(v.empty());
}

TEST_CASE(
    "cell table validator: a read of a non-resident cell is a "
    "visibility violation",
    "[cell_table]") {
  auto t = make_two_level_table();
  // the root consumer reads the INNER per-batch cell (1) instead of the
  // assembled form (2): not resident at root
  t.reads[1].source = 1;
  t.cells[1].life = 2;
  t.cells[2].life = 0;
  auto const v = validate_cell_table(t, empty_root());
  REQUIRE(v.size() >= 1);
  CHECK(v.front().rule == "visibility");
}

TEST_CASE("cell table validator: life must equal the read count",
          "[cell_table]") {
  auto t = make_two_level_table();
  t.cells[2].life = 5;
  auto const v = validate_cell_table(t, empty_root());
  REQUIRE(v.size() == 1);
  CHECK(v.front().rule == "life");
}

TEST_CASE(
    "cell table validator: an assemble whose source is not sliced "
    "by the assembled instance is a chain violation",
    "[cell_table]") {
  auto t = make_two_level_table();
  t.cells[1].sliced.clear();  // source no longer sliced by (1,0)
  auto const v = validate_cell_table(t, empty_root());
  REQUIRE(v.size() >= 1);
  CHECK(v.front().rule == "chain");
}

TEST_CASE(
    "cell table validator: two Build cells of one value in one scope "
    "is a uniqueness violation",
    "[cell_table]") {
  auto t = make_two_level_table();
  t.cells.push_back(t.cells[3]);  // second root Build of value 2
  auto const v = validate_cell_table(t, empty_root());
  REQUIRE(v.size() >= 1);
  CHECK(v.back().rule == "uniqueness");
}

TEST_CASE(
    "cell table validator: operands bound to different instances of "
    "one loop group is a form violation",
    "[cell_table]") {
  using namespace sequant::eval;
  CellTable t;
  TableCell a;  // value 0, sliced by (1,0)
  a.value_id = 0;
  a.sliced.push_back({0, LoopKey{1, 0}});
  a.scope.path.push_back({LoopKey{1, 0}, 0});
  a.life = 1;
  t.cells.push_back(a);  // cell 0
  TableCell b = a;  // value 1, sliced by (1,1): another instance of group 1
  b.value_id = 1;
  b.sliced[0].second = LoopKey{1, 1};
  b.scope.path[0].first = LoopKey{1, 1};
  t.cells.push_back(b);  // cell 1
  TableCell c;           // value 2 = a * b, sliced by (1,0), built inside (1,0)
  c.value_id = 2;
  c.sliced.push_back({0, LoopKey{1, 0}});
  c.scope.path.push_back({LoopKey{1, 0}, 0});
  t.cells.push_back(c);  // cell 2
  t.reads.push_back(Read{2, 0, 0, {}});
  t.reads.push_back(Read{2, 1, 1, {}});
  auto const v = validate_cell_table(t, ScopeBlock{});
  bool form = false, visibility = false;
  for (auto const& x : v) {
    if (x.rule == "form") form = true;
    if (x.rule == "visibility") visibility = true;
  }
  CHECK(form);
  CHECK(visibility);  // cell 1 lives in (1,1), not visible inside (1,0)
}
