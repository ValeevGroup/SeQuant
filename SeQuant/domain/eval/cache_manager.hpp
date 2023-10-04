#ifndef SEQUANT_EVAL_CACHE_MANAGER_HPP
#define SEQUANT_EVAL_CACHE_MANAGER_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/domain/eval/eval_result.hpp>

#include <memory>
#include <range/v3/view.hpp>

namespace sequant {

///
/// This class implements a cache manager useful for the cases when the number
/// of times the cached objects will be accessed is known.
///
class CacheManager {
 public:
  using key_type = size_t;

 private:
  class entry {
   private:
    size_t max_life;

    size_t life_c;

    ERPtr data_p;

   public:
    explicit entry(size_t count) noexcept;

    [[nodiscard]] ERPtr access() noexcept;

    void store(ERPtr data) noexcept;

    void reset() noexcept;

    [[nodiscard]] size_t life_count() const noexcept;

    [[nodiscard]] size_t max_life_count() const noexcept;

   private:
    [[nodiscard]] int decay() noexcept;

  };  // entry

  static ERPtr store(entry &entry, ERPtr data) noexcept;

  container::map<key_type, entry> cache_map_;

 public:
  template <typename Iterable1 = container::map<size_t, size_t>>
  explicit CacheManager(Iterable1 &&decaying) noexcept {
    for (auto &&[k, c] : decaying) cache_map_.try_emplace(k, entry{c});
  }

  ///
  /// Resets all cached data.
  ///
  void reset() noexcept;

  ///
  /// @brief Access cached data.
  ///
  /// @param key The that identifies the cached data.
  /// @return ERPtr to EvalResult
  ERPtr access(key_type key) noexcept;

  ///
  /// @param key The key to identify the cached data.
  /// @param data The data to be cached.
  /// \return Pointer to the stored data. Implictly accesses the stored data,
  ///         hence, decays the lifetime if the key accesses a decaying cache
  ///         entry. Passing @c key that was not present during construction of
  ///         this CacheManager object, stores nothing, but still returns a
  ///         valid pointer to @c data.
  [[nodiscard]] ERPtr store(key_type key, ERPtr data) noexcept;

  ///
  /// \brief Check if the key exists in the database: does not check if cache
  ///        exists
  ///
  [[nodiscard]] bool exists(key_type key) const noexcept;

  /// if the key exists in the database, return the current lifetime count of
  /// the cached data otherwise return -1
  [[nodiscard]] int life(key_type key) const noexcept;

  /// if the key exists in the database, return the maximum lifetime count of
  /// the cached data that implies the maximum number of accesses allowed for
  /// this key before the cache is released. This value was set by the c'tor.
  [[nodiscard]] int max_life(key_type key) const noexcept;

  ///
  /// \return A set of all the keys present in the cache manager.
  ///
  [[nodiscard]] container::set<size_t> keys() const noexcept;

  ///
  /// Get an empty cache manager.
  ///
  static CacheManager empty() noexcept;

  // for unit testing
  template <typename T>
  struct access_by;
  template <typename T>
  friend struct access_by;

};  // CacheManager

}  // namespace sequant

#endif  // SEQUANT_EVAL_CACHE_MANAGER_HPP
