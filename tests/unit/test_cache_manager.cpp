#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <algorithm>
#include <iostream>
#include <optional>
#include <range/v3/view/zip.hpp>

#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

namespace {
using node_type = sequant::EvalNode<sequant::EvalExpr>;
using manager_type = sequant::CacheManager<node_type>;

// Helper to create distinct EvalNode keys from expressions
node_type make_node(std::wstring_view expr_str) {
  return sequant::binarize(sequant::deserialize<sequant::ResultExpr>(expr_str));
}

}  // namespace

TEST_CASE("cache_manager", "[cache_manager]") {
  using ranges::views::zip;
  using key_type = node_type;
  using count_type = size_t;
  using data_type = sequant::ResultPtr;

  auto eval_result = [](int x) {
    return sequant::eval_result<sequant::ResultScalar<int>>(x);
  };

  size_t constexpr n_decaying = 4;  // arbitrary

  // Create distinct nodes to use as cache keys
  auto const node0 = make_node(L"R{a1;i1} = f{a1;i1}");
  auto const node1 = make_node(L"R{a1;i1} = g{a1;i1}");
  auto const node2 = make_node(L"R{a1,a2;i1,i2} = g{a1,a2;i1,i2}");
  auto const node3 = make_node(L"R = A");

  auto const decaying_keys =
      std::array<key_type, n_decaying>{node0, node1, node2, node3};
  // decaying entries repeat more than once
  auto const decaying_repeats = std::array<count_type, n_decaying>{2, 2, 4, 3};
  // arbitrary vals corresponding to decaying keys
  auto const decaying_vals =
      std::array<data_type, n_decaying>{eval_result(10),  //
                                        eval_result(11),  //
                                        eval_result(20),  //
                                        eval_result(21)};

  // Build key-count pairs for construction
  using hasher_t = sequant::TreeNodeHasher<node_type>;
  using comp_t = sequant::TreeNodeEqualityComparator<node_type>;
  auto key_count_pairs =
      zip(decaying_keys, decaying_repeats) |
      ranges::to<std::unordered_map<node_type, size_t, hasher_t, comp_t>>;

  auto const man_const = manager_type(std::move(key_count_pairs));

  SECTION("Construction") {
    auto const& man = man_const;

    // verify all keys exist
    for (auto&& k : decaying_keys) REQUIRE(man.exists(k));

    // verifying the life count of decaying entries
    for (auto&& [k, c] : zip(decaying_keys, decaying_repeats))
      REQUIRE(man.life(k) == static_cast<int>(c));
  }

  SECTION("Data Access") {
    // need a non-const manager object
    auto man = man_const;
    // filling data
    for (auto&& [k, v] : zip(decaying_keys, decaying_vals)) {
      // NOTE: man.store() calls man.access() implicitly and
      // returns a ResultPtr
      // hence, a count of lifetime is lost right here
      REQUIRE(man.store(k, v));
    }

    // now accessing decaying entries' data from the cache (c - 1) times
    // where c is the corresponding entry's lifetime count
    for (auto&& [k, v, r] :
         zip(decaying_keys, decaying_vals, decaying_repeats)) {
      // r - 1: the lifetime count at this point
      for (auto i = r - 1; i > 1; --i) {
        auto entry = man.access(k);
        REQUIRE(entry);  // cannot be a nullptr
        REQUIRE(entry->get<int>() == v->get<int>());
        REQUIRE(man.life(k) == static_cast<int>(i - 1));
      }
    }

    // at this point all the decaying entries have only one lifetime left
    // accessing each decaying entry one more time should release
    // their *data* from the memory
    for (auto&& k : decaying_keys) {
      REQUIRE(man.life(k) == 1);
      REQUIRE(man.access(k));        // accessed once. non-null ptr returned
      REQUIRE_FALSE(man.access(k));  // nullptr returned
      REQUIRE(man.life(k) == 0);
    }

    // now we reset the decaying entries which restores thier lifetimes
    man.reset();
    for (auto&& [k, c] : zip(decaying_keys, decaying_repeats)) {
      REQUIRE(man.life(k) == static_cast<int>(c));
      REQUIRE_FALSE(man.access(k));  // nullptr to data returned
    }
  }

  SECTION("alive and entry_size_in_bytes") {
    auto man = man_const;

    // A key never registered in the manager.
    auto const stray = make_node(L"R = B");
    REQUIRE_FALSE(man.exists(stray));
    REQUIRE_FALSE(man.alive(stray));
    REQUIRE(man.entry_size_in_bytes(stray) == 0);

    // Registered but never stored: not alive, zero size.
    for (auto&& k : decaying_keys) {
      REQUIRE(man.exists(k));
      REQUIRE_FALSE(man.alive(k));
      REQUIRE(man.entry_size_in_bytes(k) == 0);
    }

    // After store(): alive, size matches the stored data.
    for (auto&& [k, v] : zip(decaying_keys, decaying_vals)) {
      REQUIRE(man.store(k, v));
      REQUIRE(man.alive(k));
      REQUIRE(man.entry_size_in_bytes(k) == v->size_in_bytes());
    }

    // Drain each entry's remaining life. store() consumed one access already,
    // so r - 1 accesses remain before data_p is moved out.
    for (auto&& [k, r] : zip(decaying_keys, decaying_repeats)) {
      for (auto i = r - 1; i > 0; --i) {
        REQUIRE(man.alive(k));  // still holds data before this access
        REQUIRE(man.access(k));
      }
      // Final access has drained life to 0 and moved data_p out.
      REQUIRE_FALSE(man.alive(k));
      REQUIRE(man.entry_size_in_bytes(k) == 0);
    }

    // reset() restores life counts; re-store, confirm alive, then reset()
    // and confirm entries are not-alive again.
    man.reset();
    for (auto&& [k, v] : zip(decaying_keys, decaying_vals))
      REQUIRE(man.store(k, v));
    for (auto&& k : decaying_keys) REQUIRE(man.alive(k));
    man.reset();
    for (auto&& k : decaying_keys) {
      REQUIRE_FALSE(man.alive(k));
      REQUIRE(man.entry_size_in_bytes(k) == 0);
    }
  }

  SECTION("for_each_key enumerates every registered key") {
    auto const& man = man_const;
    size_t count = 0;
    man.for_each_key([&](node_type const& k) {
      ++count;
      REQUIRE(man.exists(k));
    });
    REQUIRE(count == n_decaying);

    // empty manager: no keys, no invocations
    auto const empty = manager_type::empty();
    empty.for_each_key([](node_type const&) { FAIL("no keys expected"); });
  }

  SECTION("Hash collision safety") {
    // Two structurally different nodes
    auto const n1 = make_node(L"R{a1;i1} = f{a1;i1}");
    auto const n2 = make_node(L"R{a1,a2;i1,i2} = g{a1,a2;i1,i2}");

    // Build a CacheManager with force_hash_collisions = true
    // This forces all hashes to 0, so if the map relied only on hash equality,
    // both nodes would be treated as the same key.
    using collision_manager_type =
        sequant::CacheManager<node_type, /*force_hash_collisions=*/true>;
    using collision_hasher = sequant::TreeNodeHasher<node_type, true>;
    using collision_comp = sequant::TreeNodeEqualityComparator<node_type>;

    std::unordered_map<node_type, size_t, collision_hasher, collision_comp>
        collision_entries;
    collision_entries.emplace(n1, 2);
    collision_entries.emplace(n2, 3);

    auto cm = collision_manager_type(std::move(collision_entries));

    // Both should exist as separate entries
    REQUIRE(cm.exists(n1));
    REQUIRE(cm.exists(n2));
    REQUIRE(cm.max_life(n1) == 2);
    REQUIRE(cm.max_life(n2) == 3);

    // Store different data for each
    auto val1 = eval_result(42);
    auto val2 = eval_result(99);
    auto stored1 = cm.store(n1, val1);
    auto stored2 = cm.store(n2, val2);

    REQUIRE(stored1);
    REQUIRE(stored2);

    // Access returns the correct data for each (not mixed up)
    auto accessed1 = cm.access(n1);
    REQUIRE(accessed1);
    REQUIRE(accessed1->get<int>() == 42);

    auto accessed2 = cm.access(n2);
    REQUIRE(accessed2);
    REQUIRE(accessed2->get<int>() == 99);
  }
}

TEST_CASE("cache manager scope chain fall-through", "[cache_manager]") {
  using hasher_t = sequant::TreeNodeHasher<node_type>;
  using comp_t = sequant::TreeNodeEqualityComparator<node_type>;
  auto eval_result = [](int x) {
    return sequant::eval_result<sequant::ResultScalar<int>>(x);
  };

  auto const X = make_node(L"R{a1;i1} = f{a1;i1}");
  auto const unregistered = make_node(L"R{a1;i1} = g{a1;i1}");

  // parent registers node X with use-count 3 (NP); child is empty with parent
  // set. (Use-count 3, not 2: store() itself performs an implicit access --
  // see the NOTE in the "Data Access" section above -- consuming one life
  // before the fall-through access under test consumes a second, landing the
  // post-fall-through life count at 1 for a clean assertion.)
  std::unordered_map<node_type, size_t, hasher_t, comp_t> counts;
  counts.emplace(X, 3);
  auto parent = manager_type(std::move(counts));
  auto child = manager_type::empty();
  child.set_parent(&parent);

  (void)parent.store(X, eval_result(42));
  REQUIRE(parent.alive(X));

  // child has no local entry for X, but access() falls through to the
  // parent...
  auto found = child.access(X);
  REQUIRE(found != nullptr);
  REQUIRE(found->get<int>() == 42);
  // ...and that fall-through read decayed the parent's NP entry (3 -> 2 from
  // store()'s own implicit access, then 2 -> 1 from the fall-through read).
  REQUIRE(parent.life(X) == 1);

  // absent-everywhere returns null; a null parent (default) is a no-op.
  REQUIRE(child.access(unregistered) == nullptr);
  auto lone = manager_type::empty();
  REQUIRE(lone.access(X) == nullptr);
}

// B1 hardening target: the batched evaluator's per-iteration ("inner") cache
// is built from a per-iteration node-repeat scan that, for a hoisted
// loop-invariant node, may register the *same* canonical key locally too
// (e.g. it also happens to recur within one iteration's own terms) even
// though the value is only ever produced once, up at the ("outer") ancestor
// level, and never stored locally. access() must not let a locally-registered
// -but-never-stored entry shadow the parent: on a genuine local miss (no data
// held at this level, whether or not the key is registered here) it must
// still walk up the chain. reset() at the inner level must, symmetrically,
// touch only its own bookkeeping and never reach into the parent.
TEST_CASE("cache manager scope chain read-through survives local registration",
          "[cache_manager]") {
  using hasher_t = sequant::TreeNodeHasher<node_type>;
  using comp_t = sequant::TreeNodeEqualityComparator<node_type>;
  auto eval_result = [](int x) {
    return sequant::eval_result<sequant::ResultScalar<int>>(x);
  };

  auto const X = make_node(L"R{a1;i1} = f{a1;i1}");
  auto const Y = make_node(L"R{a1;i1} = g{a1;i1}");

  // outer registers and stores X (the hoisted, loop-invariant node).
  std::unordered_map<node_type, size_t, hasher_t, comp_t> outer_counts;
  outer_counts.emplace(X, 3);  // 3, not 2: see the NOTE in the prior test --
                               // store()'s implicit access (3 -> 2), the
                               // walk-up access under test (2 -> 1), and the
                               // post-reset ancestor-untouched check (1 -> 0)
                               // each consume one life.
  auto outer = manager_type(std::move(outer_counts));

  // inner ALSO registers X locally (e.g. its own per-iteration scan sees it
  // recur too) but never calls inner.store(X, ...) -- the value is meant to
  // come from the ancestor. inner separately registers and stores Y, a node
  // local to this scope only, to confirm reset() clears local data.
  std::unordered_map<node_type, size_t, hasher_t, comp_t> inner_counts;
  inner_counts.emplace(X, 5);
  inner_counts.emplace(Y, 2);
  auto inner = manager_type(std::move(inner_counts));
  inner.set_parent(&outer);

  (void)outer.store(X, eval_result(42));
  REQUIRE(outer.alive(X));
  REQUIRE(inner.exists(X));       // registered locally...
  REQUIRE_FALSE(inner.alive(X));  // ...but never stored locally.

  (void)inner.store(Y, eval_result(7));
  REQUIRE(inner.alive(Y));

  // Local miss on X (registered but no local data) must still walk up to the
  // ancestor's stored value, not short-circuit to null.
  auto found = inner.access(X);
  REQUIRE(found != nullptr);
  REQUIRE(found->get<int>() == 42);

  // inner.reset() clears only inner's own bookkeeping: Y's data is dropped,
  // but the ancestor's X entry (never touched by inner.reset()) survives with
  // its post-walk-up life count intact.
  inner.reset();
  REQUIRE_FALSE(inner.alive(Y));
  REQUIRE(inner.life(Y) == 2);
  REQUIRE_FALSE(inner.alive(X));  // was never alive locally to begin with

  auto still_there = outer.access(X);
  REQUIRE(still_there != nullptr);
  REQUIRE(still_there->get<int>() == 42);
}

// F2 baseline: pin the per-loop scope invariant on the already-working
// CONTRACTED path (scope_level >= 0, i.e. an actual nested-loop chain of
// three or more levels), not just the single parent/child link exercised
// above. A node stored at an inner loop's ancestor (level k) must survive
// resets of loops nested *inside* it (level > k) and be cleared only by a
// reset() at its own level.
TEST_CASE("cache manager scope chain per-loop storage survives inner resets",
          "[cache_manager]") {
  using hasher_t = sequant::TreeNodeHasher<node_type>;
  using comp_t = sequant::TreeNodeEqualityComparator<node_type>;
  auto eval_result = [](int x) {
    return sequant::eval_result<sequant::ResultScalar<int>>(x);
  };

  auto const X = make_node(L"R{a1;i1} = f{a1;i1}");
  auto const Y = make_node(L"R{a1;i1} = g{a1;i1}");

  // outer -> mid -> inner: a three-level loop-nest chain. X is registered at
  // every level (mirroring how a per-loop node-repeat scan may see the same
  // canonical key recur at each nesting level -- see the read-through test
  // above) but only ever *stored* at mid (the "level-1" scope).
  std::unordered_map<node_type, size_t, hasher_t, comp_t> outer_counts;
  outer_counts.emplace(X, 3);
  auto outer = manager_type(std::move(outer_counts));

  // 3, not 2: store()'s own implicit access (3 -> 2) and inner's fall-through
  // walk-up access under test (2 -> 1) each consume one life, landing at 1
  // for a clean post-walk-up alive() check.
  std::unordered_map<node_type, size_t, hasher_t, comp_t> mid_counts;
  mid_counts.emplace(X, 3);
  auto mid = manager_type(std::move(mid_counts));
  mid.set_parent(&outer);

  std::unordered_map<node_type, size_t, hasher_t, comp_t> inner_counts;
  inner_counts.emplace(X, 5);
  inner_counts.emplace(Y, 2);
  auto inner = manager_type(std::move(inner_counts));
  inner.set_parent(&mid);

  (void)mid.store(X, eval_result(42));
  REQUIRE(mid.alive(X));

  // inner has no local data for X (only a local registration); access()
  // must walk up through mid (a hit) rather than stop at outer or null.
  auto found = inner.access(X);
  REQUIRE(found != nullptr);
  REQUIRE(found->get<int>() == 42);
  REQUIRE(mid.life(X) == 1);

  // a distinct key Y, local to inner only, to confirm inner.reset() clears
  // its own data while leaving ancestor levels untouched.
  (void)inner.store(Y, eval_result(7));
  REQUIRE(inner.alive(Y));

  // inner.reset() (the innermost loop's per-iteration reset) must clear only
  // inner's own bookkeeping: Y is dropped, but X -- stored one level up at
  // mid -- must survive untouched.
  inner.reset();
  REQUIRE_FALSE(inner.alive(Y));
  REQUIRE(inner.life(Y) == 2);
  REQUIRE(mid.alive(X));
  REQUIRE(mid.life(X) == 1);

  // mid.reset() (its own loop-level reset) clears X, which lives at mid's
  // own level.
  mid.reset();
  REQUIRE_FALSE(mid.alive(X));
}

// access_at() surfaces the value's LIFETIME SCOPE as a hop distance: the number
// of parent (scope-chain) links crossed to reach the scope that actually holds
// the data. This is what the Enter-stage slice-on-use needs -- a value fetched
// `hops` scopes up does not have those `hops` innermost loops' slices baked in,
// so exactly those loops must be sliced on use. Chain outer -> mid -> inner,
// store X only in outer, and read it from each level.
TEST_CASE("cache manager access_at surfaces the lifetime scope hop distance",
          "[cache_manager]") {
  using hasher_t = sequant::TreeNodeHasher<node_type>;
  using comp_t = sequant::TreeNodeEqualityComparator<node_type>;
  auto eval_result = [](int x) {
    return sequant::eval_result<sequant::ResultScalar<int>>(x);
  };

  auto const X = make_node(L"R{a1;i1} = f{a1;i1}");

  // X is registered (and stored) ONLY at outer, so mid and inner miss locally
  // and fall through. A generous use-count keeps X alive across the reads
  // below (each access_at decays outer's entry once).
  std::unordered_map<node_type, size_t, hasher_t, comp_t> outer_counts;
  outer_counts.emplace(X, 10);
  auto outer = manager_type(std::move(outer_counts));
  auto mid = manager_type::empty();
  mid.set_parent(&outer);
  auto inner = manager_type::empty();
  inner.set_parent(&mid);

  (void)outer.store(X, eval_result(42));
  REQUIRE(outer.alive(X));

  // outer: local hit, zero hops.
  auto const a_outer = outer.access_at(X);
  REQUIRE(a_outer.ptr != nullptr);
  REQUIRE(a_outer.ptr->get<int>() == 42);
  REQUIRE(a_outer.hops == 0);

  // mid: one link up (mid -> outer).
  auto const a_mid = mid.access_at(X);
  REQUIRE(a_mid.ptr != nullptr);
  REQUIRE(a_mid.ptr->get<int>() == 42);
  REQUIRE(a_mid.hops == 1);

  // inner: two links up (inner -> mid -> outer).
  auto const a_inner = inner.access_at(X);
  REQUIRE(a_inner.ptr != nullptr);
  REQUIRE(a_inner.ptr->get<int>() == 42);
  REQUIRE(a_inner.hops == 2);

  // The ptr-only access() forwarder still resolves the same value (the path
  // the three non-batched callers keep using).
  REQUIRE(inner.access(X) != nullptr);
}

// access_at_hops() mirrors access_at() but fetches at an EXACT hop count
// (walk hops parent links, then ONE local entry::access()) instead of
// searching the whole chain -- what a router-directed read (Phase 2 T3) uses
// once it already knows the target scope (e.g. from a HomeTarget's
// home_depth). Same outer -> mid -> inner chain, X stored only in outer.
TEST_CASE("cache manager access_at_hops fetches at an exact hop count",
          "[cache_manager]") {
  using hasher_t = sequant::TreeNodeHasher<node_type>;
  using comp_t = sequant::TreeNodeEqualityComparator<node_type>;
  auto eval_result = [](int x) {
    return sequant::eval_result<sequant::ResultScalar<int>>(x);
  };

  auto const X = make_node(L"R{a1;i1} = f{a1;i1}");
  // Y is registered (so it exists in outer's map) but never stored -- a
  // scope that does not currently hold it.
  auto const Y = make_node(L"R{a1,a2;i1,i2} = g{a1,a2;i1,i2}");

  std::unordered_map<node_type, size_t, hasher_t, comp_t> outer_counts;
  outer_counts.emplace(X, 10);
  outer_counts.emplace(Y, 10);
  auto outer = manager_type(std::move(outer_counts));
  auto mid = manager_type::empty();
  mid.set_parent(&outer);
  auto inner = manager_type::empty();
  inner.set_parent(&mid);

  (void)outer.store(X, eval_result(42));  // store()'s internal access(): -1
  REQUIRE(outer.life(X) == 9);

  // inner -> outer is 2 hops (inner -> mid -> outer).
  auto const via_hops = inner.access_at_hops(X, 2);
  REQUIRE(via_hops != nullptr);
  CHECK(via_hops->get<int>() == 42);
  // Exactly one further decay -- the same single lifetime step access_at()
  // takes.
  CHECK(outer.life(X) == 8);

  auto const via_access_at = inner.access_at(X).ptr;
  REQUIRE(via_access_at != nullptr);
  CHECK(via_access_at->get<int>() == 42);
  CHECK(outer.life(X) == 7);

  // Both accessor paths resolve to the SAME buffer (pointer identity).
  CHECK(via_hops.get() == via_access_at.get());

  // A hop count naming a scope that does not currently hold the value (Y was
  // registered but never stored) returns null.
  CHECK(inner.access_at_hops(Y, 2) == nullptr);

  // A hop count past the chain end (walked off the root) returns null.
  CHECK(inner.access_at_hops(X, 5) == nullptr);
}

TEST_CASE("cache_manager_persistent", "[cache_manager]") {
  using hasher_t = sequant::TreeNodeHasher<node_type>;
  using comp_t = sequant::TreeNodeEqualityComparator<node_type>;
  auto eval_result = [](int x) {
    return sequant::eval_result<sequant::ResultScalar<int>>(x);
  };

  auto const np = make_node(L"R{a1;i1} = f{a1;i1}");  // non-persistent
  auto const p = make_node(L"R{a1;i1} = g{a1;i1}");   // persistent

  std::unordered_map<node_type, size_t, hasher_t, comp_t> counts;
  counts.emplace(np, 2);  // NP, drained after 2 accesses
  counts.emplace(p, 1);   // P, count irrelevant (registered once)

  comp_t eq;
  auto is_persistent = [&p, &eq](node_type const& k) { return eq(k, p); };
  auto man = manager_type(std::move(counts), is_persistent);

  // A persistent entry is never drained: arbitrarily many accesses all return
  // the stored data (unlike an NP entry, whose data is released after its
  // max_life-th access).
  man.store(p, eval_result(20));
  for (int i = 0; i < 10; ++i) {
    auto r = man.access(p);
    REQUIRE(r);
    REQUIRE(r->get<int>() == 20);
  }
  REQUIRE(man.alive(p));

  // reset() clears the non-persistent entry but keeps the persistent one, so
  // the latter's data survives across evaluations (e.g. CC iterations).
  man.store(np, eval_result(10));
  man.reset();
  REQUIRE(man.access(np) == nullptr);  // NP cleared by reset
  auto rp = man.access(p);             // P survives reset
  REQUIRE(rp);
  REQUIRE(rp->get<int>() == 20);
  REQUIRE(man.alive(p));
}

// Task 9 regression tripwire: after the combined single-DAG evaluation
// (every value built exactly once per evaluation), a NON-persistent cache
// entry re-store()'d with no intervening reset() means a duplicate producer
// survived -- a bug (see entry::store() / stored_this_eval_ in
// cache_manager.hpp).
//
// NOTE on build config: SEQUANT_ASSERT compiles to a no-op unless
// SEQUANT_ASSERT_ENABLED is #defined (not the case for this project's
// default Release config: SEQUANT_ASSERT_BEHAVIOR=IGNORE there; Debug
// defaults to ABORT, which is also not catchable). So the PRIMARY assertion
// here is the FLAG STATE itself (store() sets it, reset() clears it), which
// is meaningful in every build config. The actual throw is exercised only
// when this build happens to be configured with
// -DSEQUANT_ASSERT_BEHAVIOR=THROW (checked at runtime via
// sequant::assert_behavior(), following the pattern in test_macros.cpp).
TEST_CASE("cache_manager restore tripwire", "[cache_manager]") {
  using hasher_t = sequant::TreeNodeHasher<node_type>;
  using comp_t = sequant::TreeNodeEqualityComparator<node_type>;
  auto eval_result = [](int x) {
    return sequant::eval_result<sequant::ResultScalar<int>>(x);
  };

  auto const np = make_node(L"R{a1;i1} = f{a1;i1}");  // non-persistent
  auto const p = make_node(L"R{a1;i1} = g{a1;i1}");   // persistent

  std::unordered_map<node_type, size_t, hasher_t, comp_t> counts;
  counts.emplace(np, 3);  // enough life to survive several stores/accesses
  counts.emplace(p, 1);

  comp_t eq;
  auto is_persistent = [&p, &eq](node_type const& k) { return eq(k, p); };
  auto man = manager_type(std::move(counts), is_persistent);

  SECTION("flag state: set on store, cleared on reset (non-persistent)") {
    REQUIRE_FALSE(man.stored_this_eval(np));

    (void)man.store(np, eval_result(1));
    REQUIRE(man.stored_this_eval(np));

    // reset() clears the flag -- a subsequent store() is the legitimate
    // (non-flagged) case.
    man.reset();
    REQUIRE_FALSE(man.stored_this_eval(np));

    (void)man.store(np, eval_result(2));
    REQUIRE(man.stored_this_eval(np));
  }

  SECTION("persistent entries never set (or need) the flag") {
    REQUIRE_FALSE(man.stored_this_eval(p));

    // A persistent entry legitimately re-stores across batch replays with no
    // intervening reset() -- this must never flag, and the flag stays false
    // throughout (it is not even consulted for persistent entries).
    (void)man.store(p, eval_result(10));
    REQUIRE_FALSE(man.stored_this_eval(p));
    (void)man.store(p, eval_result(20));  // re-store, no reset(): legitimate
    REQUIRE_FALSE(man.stored_this_eval(p));
    (void)man.store(p, eval_result(30));
    REQUIRE_FALSE(man.stored_this_eval(p));
  }

  SECTION("assert-enabled + THROW: re-store without reset() throws") {
    if (sequant::assert_behavior() != sequant::AssertBehavior::Throw) {
      // Default build config (Release: IGNORE: no-op; Debug: ABORT: not
      // catchable) -- nothing to observe via REQUIRE_THROWS here. The flag-
      // state sections above already exercise the guard's logic.
      return;
    }
    (void)man.store(np,
                    eval_result(1));  // 1st store since construction/reset: OK
    REQUIRE_THROWS(man.store(np, eval_result(2)));  // 2nd, no reset(): flagged
  }
}

TEST_CASE("cache_manager_volatility_frontier", "[cache_manager]") {
  // R = f * g * t : f,g constant (NV), t the amplitude (intrinsically V).
  auto const node = make_node(L"R{a1;i1} = f{a1;a2} * g{a2;a3} * t{a3;i1}");

  auto is_volatile = [](node_type const& n) {
    return n.leaf() && n->is_tensor() && n->as_tensor().label() == L"t";
  };
  std::function<bool(node_type const&)> has_t =
      [&](node_type const& n) -> bool {
    return n.leaf() ? is_volatile(n) : (has_t(n.left()) || has_t(n.right()));
  };

  auto man = sequant::cache_manager(std::array{node}, is_volatile);

  // For every internal node, the factory's persistence must match the rule:
  // persistent  <=>  node is NV  AND  its parent is V.
  bool saw_persistent = false;
  std::function<void(node_type const&, bool)> check =
      [&](node_type const& n, bool parent_volatile) {
        if (n.leaf()) return;
        bool const node_volatile = has_t(n);
        bool const want_p = !node_volatile && parent_volatile;
        REQUIRE(man.persistent(n) == want_p);
        if (want_p) saw_persistent = true;
        check(n.left(), node_volatile);
        check(n.right(), node_volatile);
      };
  check(node, /*parent_volatile=*/false);  // root has no (volatile) consumer
  REQUIRE(saw_persistent);  // the NV product feeding the volatile root is P
}

TEST_CASE("cache_manager_footprint_gate", "[cache_manager]") {
  // R = f * g * t : the NV product (f*g) feeds the volatile root, so without a
  // footprint gate it is cached as a persistent (cross-iteration) entry.
  auto const node = make_node(L"R{a1;i1} = f{a1;a2} * g{a2;a3} * t{a3;i1}");
  auto is_volatile = [](node_type const& n) {
    return n.leaf() && n->is_tensor() && n->as_tensor().label() == L"t";
  };
  // footprint proxy: number of result indices (the NV frontier I{a1;a3} has 2).
  auto footprint_of = [](node_type const& n) -> double {
    return n.leaf() ? 0. : static_cast<double>(n->canon_indices().size());
  };

  std::function<int(node_type const&, manager_type const&)> count_persistent =
      [&](node_type const& n, manager_type const& m) -> int {
    if (n.leaf()) return 0;
    return (m.persistent(n) ? 1 : 0) + count_persistent(n.left(), m) +
           count_persistent(n.right(), m);
  };

  // no gate (max_footprint == 0): the NV/V frontier is cached and persistent.
  auto man0 = sequant::cache_manager(std::array{node}, is_volatile);
  REQUIRE(count_persistent(node, man0) >= 1);

  // a threshold above the frontier footprint leaves caching unchanged.
  auto man_hi = sequant::cache_manager(std::array{node}, is_volatile,
                                       /*min_repeats=*/2, footprint_of, 10.);
  REQUIRE(count_persistent(node, man_hi) == count_persistent(node, man0));

  // a threshold below the 2-index frontier footprint evicts it: it is not
  // cached at all (no persistent entry; absent from the cache map).
  auto man_lo = sequant::cache_manager(std::array{node}, is_volatile,
                                       /*min_repeats=*/2, footprint_of, 1.5);
  REQUIRE(count_persistent(node, man_lo) == 0);
}

TEST_CASE("cache_manager_batch_axis_veto", "[cache_manager]") {
  // R = f * g * t : the NV product (f*g) = I{a1;a3} feeds the volatile root, so
  // by default it is cached as a persistent (cross-iteration) entry. a3 is free
  // in the frontier's own result but contracted away at the root. The
  // batch-variant veto (cache_manager.hpp) refuses a node from the run-scope
  // cache iff its cross-occurrence lifetime mask (stamp_lifetime_masks) is
  // non-empty -- some External batch mode slices it in every occurrence, so its
  // value is batch-variant.
  //
  // An External entry on a mode that is FREE on the node's OWN result is folded
  // into that node's own sliced_modes by stamp_lifetime_masks (a node's own
  // External stamp on its own slot IS part of its mask; see lifetime_mask.hpp),
  // so the mask is non-empty and the veto fires: the frontier is a
  // *consistently-sliced* external node (its value differs per external block),
  // which eval.hpp places/slices at its external loop as a hoist target -- NOT
  // a whole run-scope cache entry. Only a genuinely block-agnostic gC -- whose
  // External mode is DEMOTED to empty by the cross-occurrence meet across
  // proto-incompatible occurrences -- stays all-full and cacheable at run
  // scope; a single-occurrence node cannot reproduce that demotion. No
  // annotation at all leaves the frontier all-full, hence cached and
  // persistent.
  using sequant::BatchModeType;

  auto is_volatile = [](node_type const& n) {
    return n.leaf() && n->is_tensor() && n->as_tensor().label() == L"t";
  };
  auto in = [](auto const& v, sequant::Index const& x) {
    return std::find(v.begin(), v.end(), x) != v.end();
  };
  // a3: in the frontier's result indices, absent from the root's (contracted
  // away by the time the root is formed).
  auto find_a3 = [&](node_type const& node) -> std::optional<sequant::Index> {
    auto const& root_ix = node->canon_indices();
    auto const& frontier_ix = node.left()->canon_indices();
    std::optional<sequant::Index> a3;
    for (auto const& ix : frontier_ix)
      if (!in(root_ix, ix)) a3 = ix;
    return a3;
  };

  // baseline: no node_slice_mask annotation -> all-full mask, veto inert; the
  // NV/V frontier is cached and persistent.
  {
    auto node = make_node(L"R{a1;i1} = f{a1;a2} * g{a2;a3} * t{a3;i1}");
    auto const a3 = find_a3(node);
    REQUIRE(a3);
    auto man = sequant::cache_manager(std::array{node}, is_volatile,
                                      /*min_repeats=*/2);
    REQUIRE(man.exists(node.left()));
    REQUIRE(man.persistent(node.left()));
  }

  // External + free on the frontier's own result, single occurrence -> the mode
  // survives the cross-occurrence meet and lands in the frontier's own
  // sliced_modes, so mask_all_full() is false and the veto fires: this is a
  // consistently-sliced external node (a hoist target sliced at its external
  // loop, not a whole run-scope entry), so it is VETOED from the run-scope
  // cache
  // -- not registered in the cache map, not persistent.
  {
    auto node = make_node(L"R{a1;i1} = f{a1;a2} * g{a2;a3} * t{a3;i1}");
    auto const a3 = find_a3(node);
    REQUIRE(a3);
    node.left()->set_node_slice_mask({{*a3, BatchModeType::External}});
    node.left()->set_batch_loops_opened_here({{*a3, BatchModeType::External}});
    auto man = sequant::cache_manager(std::array{node}, is_volatile,
                                      /*min_repeats=*/2);
    REQUIRE_FALSE(man.exists(node.left()));
    REQUIRE_FALSE(man.persistent(node.left()));
  }
}

TEST_CASE("cache_manager residency", "[cache_manager]") {
  using hasher_t = sequant::TreeNodeHasher<node_type>;
  using comp_t = sequant::TreeNodeEqualityComparator<node_type>;

  auto eval_result = [](int x) {
    return sequant::eval_result<sequant::ResultScalar<int>>(x);
  };

  auto const k0 = make_node(L"R{a1;i1} = f{a1;i1}");
  auto const k1 = make_node(L"R{a1;i1} = g{a1;i1}");

  // Build a manager with two keys
  std::unordered_map<node_type, size_t, hasher_t, comp_t> key_count_pairs;
  key_count_pairs.emplace(k0, 2);
  key_count_pairs.emplace(k1, 2);
  auto man = manager_type(std::move(key_count_pairs));

  // 1. Before any store(), current_residency() == 0
  REQUIRE(man.current_residency() == 0);

  // 2. After storing data into two alive keys k0, k1,
  //    current_residency() == entry_size_in_bytes(k0) + entry_size_in_bytes(k1)
  auto val0 = eval_result(42);
  auto val1 = eval_result(99);
  (void)man.store(k0, val0);
  (void)man.store(k1, val1);

  size_t expected_size =
      man.entry_size_in_bytes(k0) + man.entry_size_in_bytes(k1);
  REQUIRE(man.current_residency() == expected_size);

  // 3. current_residency() drops to 0 once both entries are drained
  //    (their life reaches 0 and data is released).
  // store() itself performs one implicit access, so one access remains.
  REQUIRE(man.access(k0));
  REQUIRE(man.access(k1));
  // Both should now be drained
  REQUIRE_FALSE(man.alive(k0));
  REQUIRE_FALSE(man.alive(k1));
  REQUIRE(man.current_residency() == 0);

  // 4. For a standalone manager (no parent),
  //    chain_residency() == current_residency()
  // First, populate the cache again
  man.reset();
  auto val0_2 = eval_result(111);
  auto val1_2 = eval_result(222);
  (void)man.store(k0, val0_2);
  (void)man.store(k1, val1_2);

  REQUIRE(man.chain_residency() == man.current_residency());

  // Drain and confirm chain_residency() also drops to 0
  REQUIRE(man.access(k0));
  REQUIRE(man.access(k1));
  REQUIRE(man.chain_residency() == 0);
}

// chain_holds(): the peak trace's de-alias check. It reports, by POINTER
// IDENTITY across the whole scope chain, whether a buffer is already held by
// some cache -- so the per-op hwmark counts that buffer once (via the cache
// residency) instead of also adding it as an operand. A DISTINCT buffer, which
// is what a sliced / permuted / phase-shifted read of a cached value produces,
// is not held and must be counted. This is the I1 fix: alive() alone is
// local-only and would miss an ancestor-resident operand read full.
TEST_CASE("cache_manager chain_holds pointer identity", "[cache_manager]") {
  using hasher_t = sequant::TreeNodeHasher<node_type>;
  using comp_t = sequant::TreeNodeEqualityComparator<node_type>;
  auto eval_result = [](int x) {
    return sequant::eval_result<sequant::ResultScalar<int>>(x);
  };

  auto const X = make_node(L"R{a1;i1} = f{a1;i1}");  // buffer lives at parent
  auto const Y = make_node(L"R{a1;i1} = g{a1;i1}");  // buffer lives at child

  std::unordered_map<node_type, size_t, hasher_t, comp_t> parent_counts;
  parent_counts.emplace(X, 2);  // store()'s implicit access -> life 1, alive
  auto parent = manager_type(std::move(parent_counts));

  std::unordered_map<node_type, size_t, hasher_t, comp_t> child_counts;
  child_counts.emplace(Y, 2);
  auto child = manager_type(std::move(child_counts));
  child.set_parent(&parent);

  auto vX = eval_result(42);     // the buffer the parent will hold
  auto vY = eval_result(7);      // the buffer the child will hold
  auto other = eval_result(99);  // a distinct buffer held by nobody
  (void)parent.store(X, vX);
  (void)child.store(Y, vY);

  // Held up the chain (parent) and locally (child): pointer identity hits, so
  // an operand aliasing either is skipped (already in chain_residency()).
  REQUIRE(child.chain_holds(vX));
  REQUIRE(child.chain_holds(vY));
  // A distinct buffer is NOT held -- the sliced/permuted-read case that must be
  // counted, not skipped. (Same integer VALUE, different buffer.)
  REQUIRE_FALSE(child.chain_holds(other));
  REQUIRE_FALSE(child.chain_holds(eval_result(42)));
  // The chain only looks UP: the parent does not see the child's entry.
  REQUIRE_FALSE(parent.chain_holds(vY));
  // Null is never held.
  REQUIRE_FALSE(child.chain_holds(nullptr));

  // Once the parent's entry is drained, its buffer is no longer held anywhere
  // on the chain.
  REQUIRE(parent.access(X) != nullptr);  // life 1 -> 0, released
  REQUIRE_FALSE(parent.alive(X));
  REQUIRE_FALSE(child.chain_holds(vX));
}
