#ifndef SEQUANT_EVAL_CACHE_MANAGER_HPP
#define SEQUANT_EVAL_CACHE_MANAGER_HPP

#include <SeQuant/core/container.hpp>
#include <memory>

namespace sequant::eval {

///
/// This class implements a cache manager useful for the cases when the number
/// of times the cached objects will be accessed.
///
template <typename Data>
class CacheManager {
 public:
  using ptr_type = std::shared_ptr<Data>;
  using cached_type = Data;
  using key_type = size_t;

 private:
  template <typename D>
  class entry {
   private:
    using ptr_type = typename CacheManager<D>::ptr_type;

    size_t max_life;

    size_t life_c;

    ptr_type data_p;

   public:
    entry(size_t count) noexcept
        : max_life{count},  //
          life_c{count},    //
          data_p{nullptr} {}

    ptr_type access() noexcept {
      if (!data_p) return nullptr;

      return decay() == 0 ? std::move(data_p) : data_p;
    }

    void store(D data) noexcept {
      data_p = std::make_shared<D>(std::move(data));
    }

    void reset() noexcept {
      life_c = max_life;
      data_p = nullptr;
    }

    [[nodiscard]] size_t life_count() const noexcept { return life_c; }

    [[nodiscard]] size_t max_life_count() const noexcept { return max_life; }

#ifndef NDEBUG
    ///
    /// When data outlives life-time, it's a zombie entry: a bad state of
    /// cache_manager. Use it for debugging.
    ///
    [[nodiscard]] bool zombie() const noexcept { return life_c == 0 && data_p; }
#endif

   private:
    [[nodiscard]] int decay() noexcept { return life_c > 0 ? --life_c : 0; }

  };  // entry

  ptr_type store(entry<Data> &entry, Data data) noexcept {
    entry.store(std::move(data));
    return entry.access();
  }

  container::map<key_type, entry<Data>> cache_map_;

 public:
  template <typename Iterable1 = container::map<size_t, size_t>>
  explicit CacheManager(Iterable1 &&decaying) noexcept {
    for (auto &&[k, c] : decaying) cache_map_.try_emplace(k, entry<Data>{c});
  }

  ///
  /// Resets all cached data.
  ///
  void reset() noexcept {
    for (auto &&[k, v] : cache_map_) v.reset();
  }

  ///
  /// @brief Access cached data.
  ///
  /// @param key The that identifies the cached data.
  /// @retur shared_ptr to the cached data if it exists, nullptr otherwise.
  ptr_type access(key_type key) noexcept {
    if (auto &&found = cache_map_.find(key); found != cache_map_.end())
      return found->second.access();

    return nullptr;
  }

  ///
  /// @param key The key to identify the cached data.
  /// @param data The data to be cached.
  /// \return Pointer to the stored data. Implictly accesses the stored data,
  ///         hence, decays the lifetime if the key accesses a decaying cache
  ///         entry. Passing @c key that was not present during construction of
  ///         this CacheManager object, stores nothing, but still returns a
  ///         valid pointer to @c data.
  [[nodiscard]] ptr_type store(key_type key, Data data) noexcept {
    if (auto &&found = cache_map_.find(key); found != cache_map_.end())
      return store(found->second, std::move(data));
    return std::make_shared<Data>(std::move(data));
  }

  [[nodiscard]] bool exists(key_type key) const noexcept {
    return cache_map_.find(key) != cache_map_.end();
  }

#ifndef NDEBUG
  ///
  /// When data outlives life-time, it's a zombie entry: a bad state of
  /// cache_manager. Use it for debugging.
  ///
  [[nodiscard]] bool zombie(key_type key) const noexcept {
    auto iter = cache_map_.find(key);
    return iter != cache_map_.end() && !iter->second.zombie();
  }
#endif

  /// if the key exists in the database, return the current lifetime count of
  /// the cached data otherwise return -1
  [[nodiscard]] int life(key_type key) const noexcept {
    auto iter = cache_map_.find(key);
    auto end = cache_map_.end();
    return iter == end ? -1 : iter->second.life_count();
  }

  /// if the key exists in the database, return the maximum lifetime count of
  /// the cached data that implies the maximum number of accesses allowed for
  /// this key before the cache is released. This value was set by the c'tor.
  [[nodiscard]] int max_life(key_type key) const noexcept {
    auto iter = cache_map_.find(key);
    auto end = cache_map_.end();
    return iter == end ? -1 : iter->second.max_life_count();
  }

  ///
  /// Get an empty cache manager.
  ///
  static CacheManager<Data> empty() noexcept { return CacheManager<Data>{{}}; }

  // for unit testing
  template <typename T>
  struct access_by;
  template <typename T>
  friend struct access_by;

};  // CacheManager

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_CACHE_MANAGER_HPP
