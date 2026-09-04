#include <catch2/catch_test_macros.hpp>

#include <SeQuant/core/eval/backends/dryrun/result.hpp>
#include <SeQuant/core/eval/cell_registry.hpp>

#include <cstdlib>
#include <string>

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

TEST_CASE(
    "cell registry: the read that spends a cell's last life hands over sole "
    "ownership",
    "[cell_registry]") {
  // A cell with no reader left this evaluation must not keep its value alive:
  // the reader that took it is its only holder (the scope cache lets go
  // through CacheManager::release_at, driven by the flag this read reports).
  auto const t = make_table();
  CellRegistry reg(t);
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(
      cell_registry_test_regime());
  sequant::ResultPtr r = std::make_shared<sequant::eval::dryrun::ResultDryRun>(
      sequant::container::svector<sequant::Index>{sequant::Index{L"i_1"}}, cm);

  reg.set(1, r);  // cell 1: Build, non-persistent, life 1
  bool exhausted = false;
  CHECK(reg.read(1, &exhausted) == r);
  CHECK(exhausted);
  CHECK_FALSE(reg.peek(1));  // the registry has let go
  reg.set(1, r);             // a later production restores value AND life
  CHECK(reg.peek(1) == r);

  reg.set(0, r);  // cell 0: Leaf, PERSISTENT, life 2
  bool persistent_exhausted = true;
  CHECK(reg.read(0, &persistent_exhausted) == r);
  CHECK_FALSE(persistent_exhausted);  // a persistent cell never exhausts
  CHECK(reg.peek(0) == r);
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
  CHECK(*got1 == b1);  // whole read: same object
  // cell 1's life was 1: that read was its last, so the caller must be told
  // to release every other reference to the buffer.
  CHECK(res.last_read_exhausted_source());
  CHECK_THROWS(res.fetch(1, ctx));  // no remaining Read of value 1 for cell 2
  CHECK_FALSE(res.fetch(77, ctx).has_value());  // not a value: transient
}

TEST_CASE(
    "cell read resolver: leaf first touch defers without consuming; a "
    "non-leaf miss throws",
    "[cell_registry]") {
  auto const t = make_table();
  CellRegistry reg(t);
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(
      cell_registry_test_regime());
  sequant::container::svector<sequant::Index> idx{sequant::Index{L"i_1"}};
  // Neither cell 0 (Leaf) nor cell 1 (Build) is set: both reads (see
  // make_table's doc comment: consumer 2 reads cell 0 sliced, cell 1 whole)
  // start with no current result.
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

  // (a) Leaf cell 0's first touch: nullopt, and the Read is NOT consumed --
  // record_leaf populates it, and the SAME Read is then served (and
  // consumed, sliced per the read) by a second fetch.
  CHECK_FALSE(res.fetch(0, ctx).has_value());
  CHECK(res.served() == 0);  // the deferring touch serves nothing
  auto leaf = std::make_shared<sequant::eval::dryrun::ResultDryRun>(idx, cm);
  res.record_leaf(0, leaf);
  // THE first touch's read (eval.hpp's leaf branch calls exactly this after
  // record_leaf): the SAME Read, now served and SLICED per its declaration --
  // finalizing the recorded leaf whole instead would drop the slice.
  auto got0 = res.fetch(0, ctx);
  REQUIRE(got0.has_value());
  CHECK(sequant::eval::dryrun::detail::lobounds_of(**got0).at(0) == 2);
  CHECK(sequant::eval::dryrun::detail::overrides_of(**got0).at(0) == 2);
  CHECK(res.served() == 1);         // consumed exactly once for the leg
  CHECK_THROWS(res.fetch(0, ctx));  // now consumed: no remaining Read

  // (b) Build cell 1 has no current result and is never recorded by the
  // caller (unlike a Leaf, nothing ever populates it out-of-band) -- this
  // is a genuine table/tree or recording gap, so fetch throws rather than
  // deferring.
  CHECK_THROWS(res.fetch(1, ctx));
}

namespace {
// A separate small fixture for cell_of: cell 0 is an UNBOUND Assemble
// (value 5) at scope [(1,0)] -- residency root (detail::residency_scope
// returns the empty scope for an Assemble bound to none of its enclosing
// loops), so it is visible from any scope. Cell 1 is a Build (value 6) at
// the same scope [(1,0)] -- a Build's residency is always its own scope, so
// it is visible only from scopes that scope encloses.
CellTable make_residency_table() {
  CellTable t;
  TableCell a;
  a.value_id = 5;
  a.production.kind = ProductionKind::Assemble;
  a.scope.path = {{LoopKey{1, 0}, 0}};
  t.cells.push_back(a);
  TableCell b;
  b.value_id = 6;
  b.production.kind = ProductionKind::Build;
  b.scope.path = {{LoopKey{1, 0}, 0}};
  t.cells.push_back(b);
  return t;
}
}  // namespace

TEST_CASE("cell registry: cell_of finds a value's form by residency",
          "[cell_registry]") {
  auto const t = make_residency_table();
  CellRegistry reg(t);
  CellScope deep;
  deep.path = {{LoopKey{1, 0}, 0}, {LoopKey{2, 0}, 0}};
  CellScope const root;  // empty path

  // The Assemble's residency is root: found from a nested query scope AND
  // from root itself.
  auto const a_deep = reg.cell_of(5, deep);
  REQUIRE(a_deep.has_value());
  CHECK(*a_deep == 0);
  auto const a_root = reg.cell_of(5, root);
  REQUIRE(a_root.has_value());
  CHECK(*a_root == 0);

  // The Build's residency is its own scope [(1,0)]: found from a scope it
  // encloses (the deeper query), not from root.
  auto const b_deep = reg.cell_of(6, deep);
  REQUIRE(b_deep.has_value());
  CHECK(*b_deep == 1);
  CHECK_FALSE(reg.cell_of(6, root).has_value());
}

TEST_CASE("table_read spends one declared life and reports exhaustion once",
          "[cell_registry]") {
  // The shared ownership helper: EVERY site that spends a table-declared
  // life routes through it -- CellReadResolver::fetch for a consumer's
  // operand reads, and the ordered executor's block-close handoffs for the
  // read an Assemble declares of its production.source (a close that takes
  // the per-batch value from a step's own result, a child's closed result or
  // a resident home still owes the table that read). The two must not drift:
  // TableRead::exhausted is set on the read that spends the LAST life and on
  // no other, and never for a persistent cell -- the caller (not table_read
  // itself, since Task 2) decides what to do with that flag.
  auto const t = make_table();
  CellRegistry reg(t);
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(
      cell_registry_test_regime());
  sequant::ResultPtr r = std::make_shared<sequant::eval::dryrun::ResultDryRun>(
      sequant::container::svector<sequant::Index>{sequant::Index{L"i_1"}}, cm);

  std::size_t released = 0;
  reg.set(1, r);  // cell 1: Build, non-persistent, life 1
  {
    auto const tr = sequant::eval::table_read(reg, 1);
    CHECK(tr.value == r);
    CHECK(tr.exhausted);
    if (tr.exhausted) ++released;
  }
  CHECK(released == 1);
  CHECK_FALSE(reg.peek(1));  // sole ownership handed to the reader

  reg.set(0, r);  // cell 0: Leaf, PERSISTENT, life 2
  {
    auto const tr1 = sequant::eval::table_read(reg, 0);
    auto const tr2 = sequant::eval::table_read(reg, 0);
    CHECK(tr1.value == r);
    CHECK(tr2.value == r);
    if (tr1.exhausted) ++released;
    if (tr2.exhausted) ++released;
  }
  CHECK(released == 1);  // a persistent cell never exhausts
  CHECK(reg.peek(0) == r);
}

TEST_CASE("cell registry owns results: bytes, fill-once, persistence",
          "[cell_registry]") {
  // Set (and later restore) SEQUANT_UT_STRICT_FILL_ONCE at the very start of
  // the case, before anything in this process can call strict_fill_once()
  // for the first time: that function caches its env lookup in a static
  // bool on first call, so a setenv anywhere later would be too late to
  // change what it returns for the rest of this process.
  char const* const prev_strict = std::getenv("SEQUANT_UT_STRICT_FILL_ONCE");
  std::string const prev_strict_val = prev_strict ? prev_strict : "";
  setenv("SEQUANT_UT_STRICT_FILL_ONCE", "1", 1);

  // A separate small table (not make_table()'s): cell 0 Leaf persistent;
  // cell 1 Build non-persistent bound to (1,0), life 1; cell 2 Build
  // persistent (whole -- unbound), life 2.
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
  b1.sliced = {{0, LoopKey{1, 0}}};  // bound_instances reads `sliced`, not
                                     // `scope.path` -- this is what makes
                                     // clear_bound_to(LoopKey{1, 0}) reach it
  b1.life = 1;
  t.cells.push_back(b1);
  TableCell b2;
  b2.value_id = 2;
  b2.production.kind = ProductionKind::Build;
  b2.persistent = true;
  b2.life = 2;
  t.cells.push_back(b2);

  sequant::eval::PersistentValueStore store;
  std::size_t bytes_seen = 0;
  sequant::eval::CellRegistryHooks hooks;
  hooks.persistent = &store;
  hooks.hash_of = [](std::size_t vid) { return 1000 + vid; };
  hooks.on_bytes_changed = [&](std::size_t b) { bytes_seen = b; };
  sequant::eval::CellRegistry reg(t, hooks);

  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(
      cell_registry_test_regime());
  sequant::ResultPtr r1 = std::make_shared<sequant::eval::dryrun::ResultDryRun>(
      sequant::container::svector<sequant::Index>{sequant::Index{L"i_1"}}, cm);
  sequant::ResultPtr r2 = std::make_shared<sequant::eval::dryrun::ResultDryRun>(
      sequant::container::svector<sequant::Index>{sequant::Index{L"i_2"}}, cm);

  reg.set(1, r1);
  CHECK(reg.live_bytes() == r1->size_in_bytes());
  CHECK(bytes_seen == reg.live_bytes());
  CHECK_THROWS(reg.set(1, r1));  // fill-once: cell 1 not read/cleared since
  reg.clear_bound_to(sequant::eval::LoopKey{1, 0});
  CHECK(reg.live_bytes() == 0);
  reg.set(1, r1);  // allowed again after the clearing boundary
  reg.set(2, r2);
  CHECK(store.holds(1002));  // persistent cell published
  CHECK(store.get(1002) == r2);
  bool ex = false;
  auto got = reg.read(1, &ex);
  CHECK(ex);
  CHECK(got == r1);
  CHECK_FALSE(reg.peek(1));
  CHECK(reg.drained(1));
  CHECK(reg.live_bytes() == r2->size_in_bytes());

  sequant::eval::CellRegistry reg2(t, hooks);
  reg2.seed_persistent();
  CHECK(reg2.peek(2) == r2);  // seeded from the store
  CHECK_FALSE(reg2.peek(1));

  if (prev_strict)
    setenv("SEQUANT_UT_STRICT_FILL_ONCE", prev_strict_val.c_str(), 1);
  else
    unsetenv("SEQUANT_UT_STRICT_FILL_ONCE");
}
