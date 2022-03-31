#include <SeQuant/domain/eval/cache_manager.hpp>
#include <iostream>
#include <range/v3/view.hpp>

#include "catch.hpp"

TEST_CASE("TEST_CACHE_MANAGER", "[cache_manager]") {
  using ranges::views::concat;
  using ranges::views::zip;
  using data_type = int;
  using manager_type = sequant::eval::CacheManager<data_type>;
  using key_type = manager_type::key_t;
  using count_type = manager_type::count_t;

  auto const n_persistent = 4; // arbitrary
  auto const n_decaying = 4; // arbitrary
  // arbitrary keys and vals
  auto const decaying_keys = std::array<key_type, n_persistent>{100, 110, 200, 210};
  // decaying entries repeat more than once
  auto const decaying_repeats = std::array<count_type, n_persistent>{2, 2, 4, 3};
  auto const decaying_vals = std::array<data_type, n_persistent>{10, 11, 20, 21};

  // arbitrary vals and keys not present in decaying_keys
  auto const persistent_keys = std::array<key_type, n_decaying>{111, 222, 333, 444};
  auto const persistent_vals = std::array<data_type, n_decaying>{11, 22, 33, 44};

  auto man =
      manager_type(zip(decaying_keys, decaying_repeats), persistent_keys);

  // filling data
  for (auto&& [k,v]:
       zip(concat(decaying_keys, persistent_keys),
           concat(decaying_vals, persistent_vals))) {
    // NOTE: man.store() calls man.access() implicitly and
    // returns a shared_ptr to data
    // hence, a count of lifetime is lost right here
    REQUIRE(man.store(k, v));
  }

  auto const man_copy = man;

  SECTION("Construction") {
    REQUIRE(man.cache_map().size() == n_persistent + n_decaying);
    for (auto&& k: concat(decaying_keys, persistent_keys))
      REQUIRE(man.cache_map().find(k) != man.cache_map().end());
    for (auto&& [k, c]: zip(decaying_keys, decaying_repeats))
      // (c - 1) because a lifetime count is lost by implicit access
      // during storing
      REQUIRE(man.cache_map().find(k)->second.life_count() == c-1);
  }

  SECTION("Data Access") {
    // restore the cache manager in with full lifetimes and data
    man = man_copy;
    for (auto&& [k, v, r] :
         zip(decaying_keys, decaying_vals, decaying_repeats)) {
      // r - 1: the lifetime count at this point
      for (auto i = r - 1; i > 1; --i) {
        auto entry = man.access(k);
        REQUIRE(entry);          // optional<shared_ptr<..>>
        REQUIRE(entry.value());  // shared_ptr<..>
        REQUIRE(*entry.value() == v);
        auto iter = man.cache_map().find(k);
        REQUIRE(iter != man.cache_map().end());
        REQUIRE(iter->second.life_count() == i - 1);
      }
    }
    // at this point all the decaying entries have only one lifetime left
    // accessing each decaying entry one more time should release
    // their *data* from the memory
    for (auto&& k: decaying_keys) {
      auto iter = man.cache_map().find(k);
      REQUIRE(iter->second.life_count() == 1);
      REQUIRE(man.access(k).value()); // accessed once. non-null ptr returned
      REQUIRE_FALSE(man.access(k).value()); // nullptr returned
      REQUIRE(iter->second.life_count() == 0);
    }

    // meanwhile, the persistent entries are all intact
    for (auto&& [k,v]: zip(persistent_keys, persistent_vals))
      REQUIRE(*man.access(k).value() == v);

    // now we reset the decaying entries which restores thier lifetimes
    man.reset_decaying();
    for (auto&& [k,c]: zip(decaying_keys, decaying_repeats))
      REQUIRE(man.cache_map().find(k)->second.life_count() == c);

    // now we reset all entries
    man.reset_all();
    for (auto&& k: concat(decaying_keys, persistent_keys))
      REQUIRE_FALSE(man.access(k).value()); // nullptr to data returned
  }
}
