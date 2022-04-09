#include <SeQuant/domain/eval/cache_manager.hpp>
#include <iostream>
#include <range/v3/view.hpp>

#include "catch.hpp"

using data_type = int;

namespace sequant::eval {
struct TestCacheManager {};

template <>
template <>
struct CacheManager<data_type>::template access_by<TestCacheManager>{
  auto const& map(CacheManager<data_type> const& man) const {
    return man.cache_map_;
  }
};
} //

TEST_CASE("TEST_CACHE_MANAGER", "[cache_manager]") {
  using ranges::views::concat;
  using ranges::views::zip;
  using manager_type = sequant::eval::CacheManager<data_type>;
  using key_type = manager_type::key_t;
  using count_type = manager_type::count_t;
  using tester_type = manager_type::access_by<sequant::eval::TestCacheManager>;

  auto const tester = tester_type{};

  size_t constexpr n_persistent = 4;  // arbitrary
  size_t constexpr n_decaying = 4;    // arbitrary

  // arbitrary decaying keys
  auto const decaying_keys =
      std::array<key_type, n_decaying>{100, 110, 200, 210};
  // decaying entries repeat more than once
  auto const decaying_repeats = std::array<count_type, n_decaying>{2, 2, 4, 3};
  // arbitrary vals corresponding to decaying keys
  auto const decaying_vals = std::array<data_type, n_decaying>{10, 11, 20, 21};

  // arbitrary persistent vals, and keys not present in decaying_keys
  auto const persistent_keys =
      std::array<key_type, n_persistent>{111, 222, 333, 444};
  auto const persistent_vals =
      std::array<data_type, n_persistent>{11, 22, 33, 44};

  auto const man_const =
      manager_type(zip(decaying_keys, decaying_repeats), persistent_keys);

  SECTION("Construction") {
    auto const& man = man_const;
    auto const& map = tester.map(man); // access private map object

    REQUIRE(map.size() == n_persistent + n_decaying);

    // verifying the life count of decaying entries
    for (auto&& [k, c] : zip(decaying_keys, decaying_repeats))
      REQUIRE(map.find(k)->second.life_count() == c);
  }

  SECTION("Data Access") {
    // need a non-const manager object
    auto man = man_const;
    auto const& map = tester.map(man); // access private map object
    // filling data
    auto const kvs = zip(concat(decaying_keys, persistent_keys),
                         concat(decaying_vals, persistent_vals))
                     | ranges::to<sequant::container::map<key_type,data_type>>;
    for (auto&& [k, v] : kvs) {
      // NOTE: man.store() calls man.access() implicitly and
      // returns a shared_ptr to data
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
        REQUIRE(entry);          // optional<shared_ptr<..>>
        REQUIRE(entry.value());  // shared_ptr<..>
        REQUIRE(*entry.value() == v);
        auto iter = map.find(k);
        REQUIRE(iter != map.end());
        REQUIRE(iter->second.life_count() == i - 1);
      }
    }

    // at this point all the decaying entries have only one lifetime left
    // accessing each decaying entry one more time should release
    // their *data* from the memory
    for (auto&& k : decaying_keys) {
      auto iter = map.find(k);
      REQUIRE(iter->second.life_count() == 1);
      REQUIRE(man.access(k).value());  // accessed once. non-null ptr returned
      REQUIRE_FALSE(man.access(k).value());  // nullptr returned
      REQUIRE(iter->second.life_count() == 0);
    }

    // meanwhile, the persistent entries are all intact
    for (auto&& [k, v] : zip(persistent_keys, persistent_vals))
      REQUIRE(*man.access(k).value() == v);

    // now we reset the decaying entries which restores thier lifetimes
    man.reset_decaying();
    for (auto&& [k, c] : zip(decaying_keys, decaying_repeats))
      REQUIRE(map.find(k)->second.life_count() == c);

    // now we reset all entries
    man.reset_all();
    for (auto&& k : concat(decaying_keys, persistent_keys))
      REQUIRE_FALSE(man.access(k).value());  // nullptr to data returned
  }
}