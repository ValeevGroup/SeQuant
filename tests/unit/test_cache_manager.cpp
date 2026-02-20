#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <iostream>
#include <range/v3/view.hpp>

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
