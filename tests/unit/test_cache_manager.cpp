#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/io/shorthands.hpp>
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
  // in the frontier's result but contracted away at the root. Declaring a3 a
  // batchable axis means the runtime slices the frontier over it and the
  // optimizer prices it sliced, so the cache must NOT materialize it whole: the
  // veto drops it from the cache. i1 -- carried by the root, not the frontier
  // -- is the control: declaring it batchable must not touch the frontier.
  auto const node = make_node(L"R{a1;i1} = f{a1;a2} * g{a2;a3} * t{a3;i1}");
  auto is_volatile = [](node_type const& n) {
    return n.leaf() && n->is_tensor() && n->as_tensor().label() == L"t";
  };

  // a3: in the frontier's result indices, absent from the root's (contracted).
  // i1: in the root's result indices, absent from the frontier's.
  auto const& root_ix = node->canon_indices();
  auto const& frontier_ix = node.left()->canon_indices();  // the (f*g) product
  auto in = [](auto const& v, sequant::Index const& x) {
    return std::find(v.begin(), v.end(), x) != v.end();
  };
  std::optional<sequant::Index> a3, i1;
  for (auto const& ix : frontier_ix)
    if (!in(root_ix, ix)) a3 = ix;
  for (auto const& ix : root_ix)
    if (!in(frontier_ix, ix)) i1 = ix;
  REQUIRE(a3);
  REQUIRE(i1);

  std::function<int(node_type const&, manager_type const&)> count_persistent =
      [&](node_type const& n, manager_type const& m) -> int {
    if (n.leaf()) return 0;
    return (m.persistent(n) ? 1 : 0) + count_persistent(n.left(), m) +
           count_persistent(n.right(), m);
  };

  // baseline: no batchable axis -> the NV/V frontier is cached and persistent.
  auto man0 = sequant::cache_manager(std::array{node}, is_volatile);
  REQUIRE(count_persistent(node, man0) >= 1);

  // control: i1 is not carried by the frontier, so the veto leaves it cached.
  auto man_nomatch = sequant::cache_manager(
      std::array{node}, is_volatile, /*min_repeats=*/2,
      sequant::zero_footprint{}, /*max_footprint=*/0.,
      [&](sequant::Index const& ix) { return ix == *i1; });
  REQUIRE(count_persistent(node, man_nomatch) == count_persistent(node, man0));

  // a3 is carried free by the frontier -> the veto drops it: not cached at all.
  auto man_match = sequant::cache_manager(
      std::array{node}, is_volatile, /*min_repeats=*/2,
      sequant::zero_footprint{}, /*max_footprint=*/0.,
      [&](sequant::Index const& ix) { return ix == *a3; });
  REQUIRE(count_persistent(node, man_match) == 0);
}
