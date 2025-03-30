#ifndef SEQUANT_EVAL_CACHE_MANAGER_HPP
#define SEQUANT_EVAL_CACHE_MANAGER_HPP

#include <SeQuant/core/asy_cost.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/domain/eval/eval_fwd.hpp>

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

    ResultPtr data_p;

   public:
    explicit entry(size_t count) noexcept;

    [[nodiscard]] ResultPtr access() noexcept;

    void store(ResultPtr data) noexcept;

    void reset() noexcept;

    [[nodiscard]] size_t life_count() const noexcept;

    [[nodiscard]] size_t max_life_count() const noexcept;

    [[nodiscard]] size_t size_in_bytes() const noexcept;

   private:
    [[nodiscard]] int decay() noexcept;

  };  // entry

  static ResultPtr store(entry& entry, ResultPtr data) noexcept;

  container::map<key_type, entry> cache_map_;

 public:
  template <typename Iterable1 = container::map<size_t, size_t>>
  explicit CacheManager(Iterable1&& decaying) noexcept {
    for (auto&& [k, c] : decaying) cache_map_.try_emplace(k, entry{c});
  }

  ///
  /// Resets all cached data.
  ///
  void reset() noexcept;

  ///
  /// @brief Access cached data.
  ///
  /// @param key The that identifies the cached data.
  /// @return ResultPtr to Result
  ResultPtr access(key_type key) noexcept;

  ///
  /// @param key The key to identify the cached data.
  /// @param data The data to be cached.
  /// \return Pointer to the stored data. Implictly accesses the stored data,
  ///         hence, decays the lifetime if the key accesses a decaying cache
  ///         entry. Passing @c key that was not present during construction of
  ///         this CacheManager object, stores nothing, but still returns a
  ///         valid pointer to @c data.
  [[nodiscard]] ResultPtr store(key_type key, ResultPtr data) noexcept;

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
  /// \return The number of entries with life_count greater than zero.
  ///
  [[nodiscard]] size_t alive_count() const noexcept;

  ///
  /// \return Returns the sum of `Result::size_in_bytes` of alive entries.
  ///
  [[nodiscard]] size_t size_in_bytes() const noexcept;

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

///
/// \brief Make a cache manager from an iterable of evaluable nodes.
///
/// \param nodes An iterable of eval nodes.
///
/// \param min_repeats Minimum number of repeats for a node to be cached. By
///                    default anything repeated twice or more is cached.
///
/// \return A cache manager.
///
/// \see CacheManager
///
CacheManager cache_manager(meta::eval_node_range auto const& nodes,
                           size_t min_repeats = 2) noexcept {
  auto imed_counts = container::map<size_t, size_t>{};

  // visits a node and check if its hash value exists in imed_counts map
  // if it does increase the count and return false (to signal stop visiting
  // children nodes) otherwise returns true.
  auto imed_visitor = [&imed_counts](auto&& n) -> bool {
    auto&& end = imed_counts.end();
    auto&& h = hash::value(*n);
    if (auto&& found = imed_counts.find(h); found != end) {
      ++found->second;
      return false;
    } else
      imed_counts.emplace(h, 1);
    return true;
  };  // imed_visitor

  // visit imeds
  ranges::for_each(nodes, [&imed_visitor](auto&& tree) {
    tree.visit_internal(imed_visitor);
  });

  // remove less repeating imeds
  auto less_repeating = [min_repeats](auto&& pair) {
    return pair.second < min_repeats;
  };
  ranges::actions::remove_if(imed_counts, less_repeating);

  return CacheManager{imed_counts};
}

///
/// \brief Estimates the peak memory required to hold the intermediates that
///        repeat when a Sum is evaluated term by term.
/// \note Reordering the terms in a Sum affects the peak cache memory.
///
/// \param expr A Sum whose terms will be evaluated by reusing intermediates.
/// \return AsyCost object that represents the memory in terms of powers of
///         active occupied and active unoccupied index extents of stored
///         tensor.
///
AsyCost peak_cache(Sum const& expr);

}  // namespace sequant

#endif  // SEQUANT_EVAL_CACHE_MANAGER_HPP
