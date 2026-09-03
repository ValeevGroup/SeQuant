#include <catch2/catch_test_macros.hpp>

#include <SeQuant/core/eval/cell_table.hpp>
// for ScopeBlock, which validate_cell_table names but never dereferences
#include <SeQuant/core/eval/ordered_schedule.hpp>

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

namespace {
// value 1 = built inside loop (1,0) as a partial sum over that loop (a
// synthesized per-batch cell of a value whose production fuses an operand
// contraction with the reduction, so it carries no sliced position of its
// own); value 1's Sum-assembled form at root reduces that loop away.
sequant::eval::CellTable make_sum_table() {
  using namespace sequant::eval;
  CellTable t;
  TableCell inner;  // cell 0: value 1, partial sum over loop (1,0)
  inner.value_id = 1;
  inner.scope.path.push_back({LoopKey{1, 0}, 0});
  inner.production.kind = ProductionKind::Build;
  inner.partial_over.push_back(LoopKey{1, 0});
  inner.life = 1;
  t.cells.push_back(inner);  // cell 0
  TableCell assembled;       // cell 1: Sum-assemble at root, reduces (1,0)
  assembled.value_id = 1;
  assembled.production.kind = ProductionKind::Assemble;
  assembled.production.assemble = AssembleKind::Sum;
  assembled.production.source = 0;
  assembled.life = 1;
  t.cells.push_back(assembled);  // cell 1
  TableCell root_consumer;       // cell 2: root, reads the assembled form
  root_consumer.value_id = 2;
  root_consumer.production.kind = ProductionKind::Build;
  root_consumer.life = 0;
  t.cells.push_back(root_consumer);  // cell 2
  t.reads.push_back(
      Read{/*consumer=*/2, /*operand_value_id=*/1, /*source=*/1, {}});
  return t;
}
}  // namespace

TEST_CASE(
    "cell table validator: a sum assemble requires its source's "
    "partial_over to record the closing instance, and the same value_id",
    "[cell_table]") {
  using namespace sequant::eval;
  {
    auto const t = make_sum_table();
    auto const v = validate_cell_table(t, ScopeBlock{});
    for (auto const& x : v) UNSCOPED_INFO(x.rule << ": " << x.what);
    CHECK(v.empty());
  }
  {
    // source's partial_over does not record the closing instance (1,0).
    auto t = make_sum_table();
    t.cells[0].partial_over.clear();
    auto const v = validate_cell_table(t, ScopeBlock{});
    REQUIRE(v.size() == 1);
    CHECK(v.front().rule == "chain");
  }
  {
    // source is a cell of a different value.
    auto t = make_sum_table();
    t.cells[0].value_id = 99;
    auto const v = validate_cell_table(t, ScopeBlock{});
    bool chain = false;
    for (auto const& x : v) chain = chain || x.rule == "chain";
    CHECK(chain);
  }
}

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
  std::size_t n_form = 0;
  for (auto const& x : v) {
    if (x.rule == "form") {
      form = true;
      ++n_form;
    }
    if (x.rule == "visibility") visibility = true;
  }
  CHECK(form);
  CHECK(n_form == 1);  // exactly the mismatched-instance fault, no
                       // whole-vs-bound double-count
  CHECK(visibility);   // cell 1 lives in (1,1), not visible inside (1,0)
}

TEST_CASE(
    "cell table validator: an explicitly invariant read is not a form "
    "mismatch, but the same read without the marker is",
    "[cell_table]") {
  using namespace sequant::eval;
  CellTable t;
  TableCell a;  // value 0, sliced by (1,0)
  a.value_id = 0;
  a.sliced.push_back({0, LoopKey{1, 0}});
  a.scope.path.push_back({LoopKey{1, 0}, 0});
  a.life = 1;
  t.cells.push_back(a);  // cell 0
  TableCell b;  // value 1, whole -- does not carry loop (1,0)'s mode at all
  b.value_id = 1;
  b.life = 1;
  t.cells.push_back(b);  // cell 1
  TableCell c;           // value 2 = a * b, sliced by (1,0), built inside (1,0)
  c.value_id = 2;
  c.sliced.push_back({0, LoopKey{1, 0}});
  c.scope.path.push_back({LoopKey{1, 0}, 0});
  t.cells.push_back(c);  // cell 2
  t.reads.push_back(Read{2, 0, 0, {}, {}});
  Read r_b{2, 1, 1, {}, {}};
  r_b.invariant_on.push_back(LoopKey{1, 0});
  t.reads.push_back(r_b);

  {
    // b's read is explicitly marked invariant on (1,0): not a mismatch even
    // though a is bound to (1,0).
    auto const v = validate_cell_table(t, ScopeBlock{});
    std::size_t n_form = 0;
    for (auto const& x : v)
      if (x.rule == "form") ++n_form;
    CHECK(n_form == 0);
  }
  // clearing the invariant marker restores the undecided-whole-vs-bound
  // mismatch.
  t.reads.back().invariant_on.clear();
  {
    auto const v = validate_cell_table(t, ScopeBlock{});
    std::size_t n_form = 0;
    for (auto const& x : v)
      if (x.rule == "form") ++n_form;
    CHECK(n_form == 1);
  }
  // a read that is BOTH bound to another instance of the group AND marked
  // invariant on this instance: still its own "bound to another instance"
  // violation (invariant does not suppress that), but the invariant record
  // still keeps it out of the whole-vs-bound double count.
  t.cells[1].sliced.push_back({0, LoopKey{1, 1}});  // b now bound to (1,1)
  t.reads.back().invariant_on.push_back(LoopKey{1, 0});
  {
    auto const v = validate_cell_table(t, ScopeBlock{});
    std::size_t n_form = 0;
    for (auto const& x : v)
      if (x.rule == "form") ++n_form;
    CHECK(n_form == 1);
  }
}

TEST_CASE(
    "cell table: a whole Assemble is resident at root; a plain Build inside "
    "a loop is confined to its own scope regardless of persistent",
    "[cell_table]") {
  using namespace sequant::eval;

  // valid: a whole Assemble (kind Assemble, bound to none of its own
  // enclosing loops) is resident EVERYWHERE, including root, no matter how
  // deeply nested its own scope is -- the runtime's close-store walk homes
  // a finished block output that far out.
  {
    CellTable t;
    TableCell inner;  // cell 0: value 1, built inside (1,0)/(2,0), a partial
                      // sum over the depth-2 instance (bound via
                      // partial_over)
    inner.value_id = 1;
    inner.scope.path.push_back({LoopKey{1, 0}, 0});
    inner.scope.path.push_back({LoopKey{2, 0}, 0});
    inner.partial_over.push_back(LoopKey{2, 0});
    inner.production.kind = ProductionKind::Build;
    inner.life = 1;            // read once, by the Assemble
    t.cells.push_back(inner);  // cell 0

    TableCell assembled;  // cell 1: Sum-assemble at depth 1, whole (no
                          // instances of its own bound) -- residency root
    assembled.value_id = 1;
    assembled.scope.path.push_back({LoopKey{1, 0}, 0});
    assembled.production.kind = ProductionKind::Assemble;
    assembled.production.assemble = AssembleKind::Sum;
    assembled.production.source = 0;
    assembled.life = 1;
    t.cells.push_back(assembled);  // cell 1

    TableCell root_consumer;  // cell 2: root, reads cell 1 (the Assemble)
                              // whole
    root_consumer.value_id = 2;
    root_consumer.production.kind = ProductionKind::Build;
    root_consumer.life = 0;
    t.cells.push_back(root_consumer);  // cell 2

    t.reads.push_back(
        Read{/*consumer=*/2, /*operand_value_id=*/1, /*source=*/1, {}, {}});

    auto const v = validate_cell_table(t, ScopeBlock{});
    for (auto const& x : v) UNSCOPED_INFO(x.rule << ": " << x.what);
    CHECK(v.empty());
  }

  // invalid: a PLAIN Build cell (not an Assemble) at depth 1, whole, read by
  // a root Build -- a step's value is stored in its own block's cache and
  // dies with that block, so it is NEVER visible at root, regardless of
  // persistent (checked both ways).
  for (bool persistent : {false, true}) {
    CellTable t2;
    TableCell built;  // cell 0: plain Build inside (1,0), whole
    built.value_id = 3;
    built.scope.path.push_back({LoopKey{1, 0}, 0});
    built.production.kind = ProductionKind::Build;
    built.persistent = persistent;
    built.life = 1;
    t2.cells.push_back(built);  // cell 0

    TableCell root_consumer2;  // cell 1: root, reads cell 0 whole
    root_consumer2.value_id = 4;
    root_consumer2.production.kind = ProductionKind::Build;
    root_consumer2.life = 0;
    t2.cells.push_back(root_consumer2);  // cell 1

    t2.reads.push_back(
        Read{/*consumer=*/1, /*operand_value_id=*/3, /*source=*/0, {}, {}});

    auto const v = validate_cell_table(t2, ScopeBlock{});
    std::size_t n_visibility = 0;
    for (auto const& x : v)
      if (x.rule == "visibility") ++n_visibility;
    CHECK(n_visibility == 1);
    REQUIRE(v.size() == 1);
    CHECK(v.front().rule == "visibility");
  }
}

TEST_CASE(
    "cell table validator: life weighs a read by the consumer's extra "
    "batches against the SOURCE'S RESIDENCY scope, not its raw scope",
    "[cell_table]") {
  using namespace sequant::eval;
  auto const n_batches_of = [](LoopKey const& k) -> std::size_t {
    return k.depth == 1 ? 3 : 1;
  };

  // source Build at root, consumer Build inside (1,0) reading it whole: the
  // consumer's own loop is not on the source's (root) residency path, so
  // the read is weighted by n_batches_of((1,0)) = 3.
  {
    CellTable t;
    TableCell src;  // cell 0: value 0, Build at root
    src.value_id = 0;
    src.life = 3;
    t.cells.push_back(src);  // cell 0
    TableCell consumer;      // cell 1: value 1, Build inside (1,0), reads src
                             // whole
    consumer.value_id = 1;
    consumer.scope.path.push_back({LoopKey{1, 0}, 0});
    t.cells.push_back(consumer);  // cell 1
    t.reads.push_back(Read{/*consumer=*/1,
                           /*operand_value_id=*/0,
                           /*source=*/0,
                           {},
                           {}});

    {
      auto const v = validate_cell_table(t, ScopeBlock{}, n_batches_of);
      for (auto const& x : v) UNSCOPED_INFO(x.rule << ": " << x.what);
      CHECK(v.empty());
    }
    t.cells[0].life = 1;
    {
      auto const v = validate_cell_table(t, ScopeBlock{}, n_batches_of);
      std::size_t n_life = 0;
      for (auto const& x : v)
        if (x.rule == "life") ++n_life;
      CHECK(n_life == 1);
    }
  }

  // a whole Assemble at depth 1 (residency root, per the kind-based rule)
  // read by a consumer at depth 1 inside a SIBLING nest (1,1): the
  // consumer's loop instance is not in the Assemble's RESIDENCY path (root),
  // even though the Assemble's own SCOPE is depth 1 -- multiplicity 3.
  {
    CellTable t;
    TableCell inner;  // cell 0: value 0, partial sum over depth-2 loop
                      // (2,0), built inside depth-1 loop (1,0)
    inner.value_id = 0;
    inner.scope.path.push_back({LoopKey{1, 0}, 0});
    inner.scope.path.push_back({LoopKey{2, 0}, 0});
    inner.partial_over.push_back(LoopKey{2, 0});
    inner.life = 1;
    t.cells.push_back(inner);  // cell 0

    TableCell assembled;  // cell 1: Sum-assemble at depth 1, whole
    assembled.value_id = 0;
    assembled.scope.path.push_back({LoopKey{1, 0}, 0});
    assembled.production.kind = ProductionKind::Assemble;
    assembled.production.assemble = AssembleKind::Sum;
    assembled.production.source = 0;
    assembled.life = 3;            // one read, multiplicity 3
    t.cells.push_back(assembled);  // cell 1

    TableCell consumer;  // cell 2: value 1, Build inside sibling (1,1)
    consumer.value_id = 1;
    consumer.scope.path.push_back({LoopKey{1, 1}, 0});
    t.cells.push_back(consumer);  // cell 2

    t.reads.push_back(Read{/*consumer=*/2,
                           /*operand_value_id=*/0,
                           /*source=*/1,
                           {},
                           {}});

    {
      auto const v = validate_cell_table(t, ScopeBlock{}, n_batches_of);
      for (auto const& x : v) UNSCOPED_INFO(x.rule << ": " << x.what);
      CHECK(v.empty());
    }
    t.cells[1].life = 1;
    {
      auto const v = validate_cell_table(t, ScopeBlock{}, n_batches_of);
      std::size_t n_life = 0;
      for (auto const& x : v)
        if (x.rule == "life") ++n_life;
      CHECK(n_life == 1);
    }
  }
}

TEST_CASE("cell table: empty table and out-of-range read ids", "[cell_table]") {
  using namespace sequant::eval;
  {
    CellTable t;
    auto const v = validate_cell_table(t, ScopeBlock{});
    CHECK(v.empty());
  }
  {
    CellTable t;
    TableCell leaf;
    leaf.value_id = 0;
    leaf.production.kind = ProductionKind::Leaf;
    t.cells.push_back(leaf);  // cell 0
    t.reads.push_back(Read{/*consumer=*/5,
                           /*operand_value_id=*/0,
                           /*source=*/0,
                           {}});
    auto const v = validate_cell_table(t, ScopeBlock{});
    CHECK(v.size() >= 1);
  }
}
