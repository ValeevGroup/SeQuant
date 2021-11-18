#ifndef SEQUANT_EVAL_CACHE_MANAGER_HPP
#define SEQUANT_EVAL_CACHE_MANAGER_HPP

#include <SeQuant/core/container.hpp>
#include <memory>

namespace sequant::eval {

template <typename Data>
class CacheManager {
 public:
  using key_t = size_t;
  using count_t = size_t;
  using ptr_t = std::shared_ptr<Data>;

 private:
  enum struct Lifetime { Persistent, Decaying };

  template <typename D>
  class entry {
   private:
    using ptr_t = typename CacheManager<D>::ptr_t;

    Lifetime life_t;

    count_t max_life;

    count_t life_c;

    ptr_t data_p;

   public:
    entry() noexcept
        : life_t{Lifetime::Persistent},  //
          max_life{0},                   //
          life_c{0},
          data_p{nullptr} {}

    entry(count_t count) noexcept
        : life_t{Lifetime::Decaying},  //
          max_life{count},             //
          life_c{count},               //
          data_p{nullptr} {}

    ptr_t access() noexcept {
      if (!data_p) return data_p;

      return decay() == 0 ? std::move(data_p) : data_p;
    }

    void store(D data) noexcept {
      data_p = std::make_shared<D>(std::move(data));
    }

    void reset(bool decaying_only) noexcept {
      if ((decaying_only && (life_t == Lifetime::Decaying)) || !decaying_only) {
        life_c = max_life;
        data_p = nullptr;
      }
    }

   private:
    [[nodiscard]] int decay() noexcept {
      return life_t == Lifetime::Persistent ? -1 : (life_c > 0 ? --life_c : 0);
    }

  };  // entry

  ptr_t store(entry<Data> &entry, Data data) {
    entry.store(std::move(data));
    return entry.access();
  }

  container::map<key_t, entry<Data>> cache_map;

 public:
  template <typename Iterable1 = container::map<key_t, count_t>,
            typename Iterable2 = container::svector<key_t>>
  CacheManager(Iterable1 &&decaying, Iterable2 &&persistent = {}) {
    for (auto &&[k, c] : decaying)
      cache_map.try_emplace(k, entry<Data>{static_cast<count_t>(c)});

    for (auto &&k : persistent) cache_map.try_emplace(k, entry<Data>{});
  }

  void reset_all() {
    for (auto &&[k, v] : cache_map) v.reset(false);
  }

  void reset_decaying() {
    for (auto &&[k, v] : cache_map) v.reset(true);
  }

  std::optional<ptr_t> access(key_t key) noexcept {
    if (auto &&found = cache_map.find(key); found != cache_map.end())
      return found->second.access();

    return std::nullopt;
  }

  ptr_t store(key_t key, Data data) {
    if (auto &&found = cache_map.find(key); found != cache_map.end())
      return (store(found->second, std::move(data)));
    return std::make_shared<Data>(std::move(data));
  }

};  // CacheManager

}  // namespace

#endif  // SEQUANT_EVAL_CACHE_MANAGER_HPP
