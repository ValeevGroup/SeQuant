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
// cell 0: a persistent Leaf, value 0; cell 1: Build value 1 at scope [(1,0)]
// sliced pos 0 on (1,0); cell 2: Build value 2 at the same scope reading
// cell 0 with a declared slice {0 -> (1,0)} and cell 1 whole.
CellTable make_table() {
  CellTable t;
  TableCell leaf;
  leaf.value_id = 0;
  leaf.production.kind = ProductionKind::Leaf;
  leaf.persistent = true;
  t.cells.push_back(leaf);
  TableCell b1;
  b1.value_id = 1;
  b1.production.kind = ProductionKind::Build;
  b1.scope.path = {{LoopKey{1, 0}, 0}};
  b1.sliced = {{0, LoopKey{1, 0}}};
  t.cells.push_back(b1);
  TableCell b2;
  b2.value_id = 2;
  b2.production.kind = ProductionKind::Build;
  b2.scope.path = {{LoopKey{1, 0}, 0}};
  b2.sliced = {{0, LoopKey{1, 0}}};
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

TEST_CASE("cell registry: set, read, release and batch clearing",
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
  CHECK(reg.read(1) == r);  // a read does not use the value up
  CHECK(reg.live_bytes() == r->size_in_bytes());
  reg.release(1);  // the schedule says nothing reads it again
  CHECK_FALSE(reg.peek(1));
  CHECK(reg.drained(1));
  CHECK(reg.live_bytes() == 0);
  CHECK_THROWS(reg.read(1));  // nothing to read
  reg.release(1);             // releasing an empty cell changes nothing
  reg.set(1, r);              // a new batch's production
  reg.clear_bound_to(LoopKey{1, 0});
  CHECK_FALSE(reg.peek(1));  // bound to (1,0): cleared
  reg.set(0, r);
  reg.clear_bound_to(LoopKey{1, 0});
  CHECK(reg.peek(0) == r);  // whole leaf: untouched
  reg.release(0);
  CHECK(reg.peek(0) == r);  // persistent: never released
}

TEST_CASE("last-read-hands-over-ownership", "[cell_registry]") {
  // The reader the schedule names as a cell's last takes the value: the
  // registry lets go, so the reader is its only holder.
  auto const t = make_table();
  CellRegistry reg(t);
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(
      cell_registry_test_regime());
  sequant::ResultPtr r = std::make_shared<sequant::eval::dryrun::ResultDryRun>(
      sequant::container::svector<sequant::Index>{sequant::Index{L"i_1"}}, cm);

  reg.set(1, r);  // cell 1: Build, non-persistent
  CHECK(reg.take(1) == r);
  CHECK_FALSE(reg.peek(1));  // the registry has let go
  CHECK(reg.drained(1));
  CHECK_THROWS(reg.take(1));
  reg.set(1, r);  // a later production is a fresh one
  CHECK(reg.peek(1) == r);

  reg.set(0, r);  // cell 0: Leaf, PERSISTENT
  CHECK(reg.take(0) == r);
  CHECK(reg.peek(0) == r);  // a persistent cell is never given away
  CHECK_FALSE(reg.drained(0));
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
  // The release plan's answer for this table: cell 2's read of cell 1 (the
  // second Read) is cell 1's last.
  sequant::container::vector<char> const last_reads{0, 1};
  res.set_last_reads(&last_reads);
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
  CHECK_FALSE(res.operand_drained(0));  // not the last read: still held
  auto got1 = res.fetch(1, ctx);
  REQUIRE(got1.has_value());
  CHECK(*got1 == b1);  // whole read: same object
  // That read was cell 1's last, so the reader took the value and the
  // registry no longer holds it -- operand_drained() reports the source cell
  // itself as drained.
  CHECK(res.operand_drained(1));
  CHECK_THROWS(res.fetch(1, ctx));  // no remaining Read of value 1 for cell 2
  CHECK_FALSE(res.fetch(77, ctx).has_value());  // not a value: transient
  // A hash the table has no value for at all is reported drained too: it is
  // a private transient no table cell could ever be sharing.
  CHECK(res.operand_drained(77));
}

TEST_CASE("resolver-leaf-first-touch", "[cell_registry]") {
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

TEST_CASE("cell registry owns results: bytes, fill-once, persistence",
          "[cell_registry]") {
  // Strictness is asked for HERE, on this registry's own hooks (see
  // CellRegistryHooks::strict_fill_once), rather than through the environment:
  // eval::strict_fill_once() latches its env lookup on the first call anywhere
  // in the process, so a setenv from inside a test case is a coin flip on the
  // order the cases ran in.
  //
  // A separate small table (not make_table()'s): cell 0 Leaf persistent;
  // cell 1 Build non-persistent bound to (1,0); cell 2 Build persistent
  // (whole -- unbound).
  CellTable t;
  TableCell leaf;
  leaf.value_id = 0;
  leaf.production.kind = ProductionKind::Leaf;
  leaf.persistent = true;
  t.cells.push_back(leaf);
  TableCell b1;
  b1.value_id = 1;
  b1.production.kind = ProductionKind::Build;
  b1.scope.path = {{LoopKey{1, 0}, 0}};
  b1.sliced = {{0, LoopKey{1, 0}}};  // bound_instances reads `sliced`, not
                                     // `scope.path` -- this is what makes
                                     // clear_bound_to(LoopKey{1, 0}) reach it
  t.cells.push_back(b1);
  TableCell b2;
  b2.value_id = 2;
  b2.production.kind = ProductionKind::Build;
  b2.persistent = true;
  t.cells.push_back(b2);

  sequant::eval::PersistentValueStore store;
  std::size_t bytes_seen = 0;
  sequant::eval::CellRegistryHooks hooks;
  hooks.persistent = &store;
  hooks.hash_of = [](std::size_t vid) { return 1000 + vid; };
  hooks.on_bytes_changed = [&](std::size_t b) { bytes_seen = b; };
  hooks.strict_fill_once = true;  // a fill-once violation must THROW here
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
  CHECK_THROWS(reg.set(1, r1));  // fill-once: cell 1 not taken/cleared since
  reg.clear_bound_to(sequant::eval::LoopKey{1, 0});
  CHECK(reg.live_bytes() == 0);
  reg.set(1, r1);  // allowed again after the clearing boundary
  reg.set(2, r2);
  CHECK(store.holds(1002));  // persistent cell published
  CHECK(store.get(1002) == r2);
  auto got = reg.take(1);
  CHECK(got == r1);
  CHECK_FALSE(reg.peek(1));
  CHECK(reg.drained(1));
  CHECK(reg.live_bytes() == r2->size_in_bytes());

  sequant::eval::CellRegistry reg2(t, hooks);
  reg2.seed_persistent();
  CHECK(reg2.peek(2) == r2);  // seeded from the store
  CHECK_FALSE(reg2.peek(1));
}

TEST_CASE("registry-clear-bound-to", "[cell_registry]") {
  // Three non-persistent Build cells at the SAME scope [(1,0)], differing
  // only in how they relate to that loop instance:
  //   cell 0 -- sliced by (1,0): EXPLICITLY bound (detail::bound_instances).
  //   cell 1 -- whole, not produce_if_absent: IMPLICITLY bound, because its
  //             own scope's innermost entry is (1,0) -- the executor re-runs
  //             whatever its tree position computes every batch of that loop.
  //   cell 2 -- whole and produce_if_absent: invariant to (1,0), produced on
  //             the first visit and REUSED by every later batch, so the
  //             implicit rule must NOT reach it.
  CellTable t;
  auto add = [&t](bool sliced, bool pia) {
    TableCell c;
    c.value_id = t.cells.size();
    c.production.kind = ProductionKind::Build;
    c.scope.path = {{LoopKey{1, 0}, 0}};
    if (sliced) c.sliced = {{0, LoopKey{1, 0}}};
    c.produce_if_absent = pia;
    t.cells.push_back(c);
  };
  add(/*sliced=*/true, /*pia=*/false);
  add(/*sliced=*/false, /*pia=*/false);
  add(/*sliced=*/false, /*pia=*/true);

  CellRegistry reg(t);
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(
      cell_registry_test_regime());
  sequant::ResultPtr r = std::make_shared<sequant::eval::dryrun::ResultDryRun>(
      sequant::container::svector<sequant::Index>{sequant::Index{L"i_1"}}, cm);
  reg.set(0, r);
  reg.set(1, r);
  reg.set(2, r);

  // A DIFFERENT loop instance's batch boundary touches none of them.
  reg.clear_bound_to(LoopKey{2, 0});
  CHECK(reg.peek(0) == r);
  CHECK(reg.peek(1) == r);
  CHECK(reg.peek(2) == r);

  reg.clear_bound_to(LoopKey{1, 0});
  CHECK_FALSE(reg.peek(0));  // explicitly bound: cleared
  CHECK_FALSE(reg.peek(1));  // implicitly bound (scope's innermost): cleared
  CHECK(reg.peek(2) == r);   // produce_if_absent: kept for reuse

  // A produce_if_absent cell IS cleared on an instance it is genuinely bound
  // to: give one a sliced position on an OUTER instance it is not homed in.
  CellTable t2;
  TableCell c;
  c.value_id = 0;
  c.production.kind = ProductionKind::Build;
  c.scope.path = {{LoopKey{1, 0}, 0}, {LoopKey{2, 0}, 0}};
  c.sliced = {{0, LoopKey{1, 0}}};  // bound to the OUTER loop
  c.produce_if_absent = true;       // invariant to its own innermost loop
  t2.cells.push_back(c);
  CellRegistry reg2(t2);
  reg2.set(0, r);
  reg2.clear_bound_to(LoopKey{2, 0});
  CHECK(reg2.peek(0) == r);  // its own innermost loop: kept
  reg2.clear_bound_to(LoopKey{1, 0});
  CHECK_FALSE(reg2.peek(0));  // an instance it IS bound to: cleared
}

TEST_CASE("resolver-unfetched-not-drained", "[cell_registry]") {
  // operand_drained() gates in-place accumulation. With no record of ever
  // having served this value, the safe answer is "not drained" (some cell may
  // still be holding it for another reader): out-of-place costs one
  // allocation, the wrong answer corrupts a shared buffer. A hash that is not
  // a table value at all is still reported drained -- nothing could share it.
  auto const t = make_table();
  CellRegistry reg(t);
  CellReadResolver rr(reg, [](std::size_t h) -> std::optional<std::size_t> {
    if (h == 100) return std::size_t{1};  // value 1 (a Build cell)
    return std::nullopt;                  // anything else: not a table value
  });
  CHECK_FALSE(rr.operand_drained(100));  // table value, never fetched here
  CHECK(rr.operand_drained(999));        // not a table value: private
}
