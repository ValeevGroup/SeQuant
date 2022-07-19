#ifndef SEQUANT_EVAL_CACHE_MANAGER_HPP
#define SEQUANT_EVAL_CACHE_MANAGER_HPP

#include <SeQuant/core/container.hpp>
#include <memory>
#include <optional>

namespace sequant::eval {

///
/// This class implements a cache manager useful for the cases when the number
/// of times the cached objects will be accessed, and/or they need to be cached
/// for indefinitely many accesses.
///
/// \tparam Data Type of the data to be cached.
/// @details
///          - It is a wrapper around a map structure.
///          - The keys of the map, which are provided by the user at
///            construction, are used to store and access data.
///          - The mapped values(MVs) of the map are a special objects with
///            the following properties:
///             - MV has a store method that takes a key and the cached data.
///               Storing a data implicitly accesses it after storing. This a
///               design decision made based on the use case of this class. So,
///               decaying entries start decaying right from the storing. See
///               below for more on decaying entries.
///             - MV has an access method that returns a pointer to the cached
///               data.
///             - MV has a tag to identify itself as having either a persistent
///               or a decaying lifetime.
///             - If MV has a decaying lifetime:
///               - It takes a count at construction aka max_life.
///               - The current lifetime count is set to max_life
///                 at construction.
///               - Accessing the cached data dcreases the MV's current
///                 lifetime count.
///               - Once the current lifetime reaches zero, the pointer to the
///                 cached data becomes equivalent to a nullptr thereby freeing
///                 the cache memory.
///               - It has a reset method that restores the current lifetime to
///                 the max_life. Reseting also resets the cached data pointer.
///               - Calling store method with zero lifetime count in the MV
///                 will store nothing.
///             - If MV has a persistent lifetime:
///               - It takes no param at construction.
///               - Allows indefinitely many accesses to the cached data.
///               - Calling reset method resets the cached data pointer.
///
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

      // decay() < 0 implies persistent lifetime
      // decay() >= 0 implies decaying lifetime
      return decay() == 0 ? std::move(data_p) : data_p;
    }

    void store(D data) {
      data_p = std::make_shared<D>(std::move(data));
    }

    void reset(bool decaying_only) noexcept {
      if ((decaying_only && (life_t == Lifetime::Decaying)) || !decaying_only) {
        life_c = max_life;
        data_p = nullptr;
      }
    }

    [[nodiscard]] count_t life_count() const noexcept { return life_c; }

   private:
    ///
    /// @details life_c == 0 for objects with Lifetime::Decaying implies full
    ///          decay. They don't decay beyond zero.
    /// \return If object has persistent lifetime return -1 else
    ///         decrement lifetime, if it hasn't fully decayed and return life_c
    [[nodiscard]] int decay() noexcept {
      return life_t == Lifetime::Persistent ? -1 : (life_c > 0 ? --life_c : 0);
    }

  };  // entry

  ptr_t store(entry<Data> &entry, Data data) {
    entry.store(std::move(data));
    return entry.access();
  }

  container::map<key_t, entry<Data>> cache_map_;

 public:
  ///
  /// @brief Construct a cache manager.
  ///        CacheManger<>::key_t type keys are expected for construction.
  ///
  /// @param decaying A map-like iterable that maps the keys to the maximum
  ///                 number of times the associated data should be accessed.
  /// @param persistent An iterable of keys to the data that are to be accessed
  ///                   an indefinitely many times.
  /// @note  Repeating keys in @c decaying and @c persistent leads to an
  ///        undefined behavior.
  template <typename Iterable1 = container::map<key_t, count_t>,
            typename Iterable2 = container::svector<key_t>>
  CacheManager(Iterable1 &&decaying, Iterable2 &&persistent) {
    for (auto &&[k, c] : decaying)
      cache_map_.try_emplace(k, entry<Data>{static_cast<count_t>(c)});

    for (auto &&k : persistent) cache_map_.try_emplace(k, entry<Data>{});
  }

  ///
  /// Resets all cached data.
  ///
  void reset_all() {
    for (auto &&[k, v] : cache_map_) v.reset(false);
  }

  ///
  /// Only resets decaying cached data, which restores their lifetimes to the
  /// values they were constructed with.
  void reset_decaying() {
    for (auto &&[k, v] : cache_map_) v.reset(true);
  }

  ///
  /// @brief Access cached data.
  ///
  /// @param key The that identifies the cached data.
  /// @return Optional object to the pointer to the cached data. Only if @c key
  ///         doesn't exist in the cache database, nullopt is returned.
  ///         In other words if @c key was not passed during construction, the
  ///         return value is a std::nullopt object.
  std::optional<ptr_t> access(key_t key) noexcept {
    if (auto &&found = cache_map_.find(key); found != cache_map_.end())
      return found->second.access();

    return std::nullopt;
  }

  ///
  /// @param key The key to identify the cached data.
  /// @param data The data to be cached.
  /// \return Pointer to the stored data. Implictly accesses the stored data,
  ///         hence, decays the lifetime if the key accesses a decaying cache
  ///         entry. Passing @c key that was not present during construction of
  ///         this CacheManager object, stores nothing, but still returns a
  ///         valid pointer to @c data.
  ptr_t store(key_t key, Data data) {
    if (auto &&found = cache_map_.find(key); found != cache_map_.end())
      return store(found->second, std::move(data));
    return std::make_shared<Data>(std::move(data));
  }

  // for unit testing
  template <typename T> struct access_by;
  template <typename T> friend struct access_by;

};  // CacheManager

}  // namespace

#endif  // SEQUANT_EVAL_CACHE_MANAGER_HPP
