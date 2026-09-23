#include <catch2/catch_test_macros.hpp>

#include <SeQuant/core/eval/cell_table.hpp>
#include <SeQuant/core/export/memory_model.hpp>
#include <SeQuant/core/export/natural_memory_model.hpp>

using namespace sequant;
using eval::CellId;
using eval::CellTable;
using eval::ProductionKind;
using eval::Read;
using eval::TableCell;

namespace {

TableCell leaf(std::size_t value_id) {
  TableCell c;
  c.value_id = value_id;
  c.production.kind = ProductionKind::Leaf;
  return c;
}

TableCell build(std::size_t value_id) {
  TableCell c;
  c.value_id = value_id;
  c.production.kind = ProductionKind::Build;
  return c;
}

TableCell assemble(std::size_t value_id, CellId source) {
  TableCell c;
  c.value_id = value_id;
  c.production.kind = ProductionKind::Assemble;
  c.production.source = source;
  return c;
}

Read read(CellId consumer, CellId source) {
  Read r;
  r.consumer = consumer;
  r.source = source;
  return r;
}

}  // namespace

TEST_CASE("natural_memory_model", "[export]") {
  SECTION("single-consumer chain: immediately adjacent -> Stack") {
    // id0, id1: leaves; id2 = Build(A) from leaves; id3: leaf; id4 =
    // Build(root) reading A and a leaf, immediately after A.
    CellTable table;
    table.cells = {leaf(0), leaf(1), build(2), leaf(3), build(4)};
    table.reads = {read(4, 2)};

    const auto natural = natural_memory_model(table);
    REQUIRE(natural.size() == table.cells.size());
    CHECK(natural.at(2) == MemoryModel::Stack);
  }

  SECTION("value shared by two consumers -> RandomAccess") {
    // id2 = Build(A); both id4 and id6 read A.
    CellTable table;
    table.cells = {leaf(0),  leaf(1), build(2), leaf(3),
                   build(4), leaf(5), build(6)};
    table.reads = {read(4, 2), read(6, 2)};

    const auto natural = natural_memory_model(table);
    CHECK(natural.at(2) == MemoryModel::RandomAccess);
  }

  SECTION(
      "single consumer, but separated by an interposed unrelated Build -> "
      "RandomAccess") {
    // id2 = Build(A); id5 = Build(Unrelated), interposed, opened while A is
    // still open and not yet closed by id5; id7 = Build(root) reads A last.
    CellTable table;
    table.cells = {leaf(0), leaf(1),  build(2), leaf(3),
                   leaf(4), build(5), leaf(6),  build(7)};
    table.reads = {read(7, 2)};

    const auto natural = natural_memory_model(table);
    CHECK(natural.at(2) == MemoryModel::RandomAccess);
  }

  SECTION(
      "Assemble consumer draining its source, adjacent -> both remain "
      "Stack") {
    // id1 = Build(A, partial form); id2 = Assemble(A, complete form,
    // source = id1), immediately closing id1; id4 = Build(root) reads the
    // assembled cell (id2) immediately after.
    CellTable table;
    table.cells = {leaf(0), build(1), assemble(1, 1), leaf(3), build(4)};
    table.reads = {read(4, 2)};

    const auto natural = natural_memory_model(table);
    CHECK(natural.at(1) == MemoryModel::Stack);
    CHECK(natural.at(2) == MemoryModel::Stack);
  }

  SECTION("leaf reads never force anything and are never classified") {
    // A consumer reading only leaves never touches the simulated stack.
    CellTable table;
    table.cells = {leaf(0), leaf(1), build(2)};
    table.reads = {};

    const auto natural = natural_memory_model(table);
    CHECK(natural.at(2) == MemoryModel::Stack);
  }
}
