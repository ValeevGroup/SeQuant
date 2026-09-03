#include <catch2/catch_test_macros.hpp>

#include <SeQuant/core/eval/cell_table.hpp>

using sequant::eval::CellScope;
using sequant::eval::CellTable;
using sequant::eval::LoopKey;
using sequant::eval::ProductionKind;
using sequant::eval::Read;
using sequant::eval::TableCell;

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
