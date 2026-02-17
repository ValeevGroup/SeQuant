#ifndef SEQUANT_EVAL_CACHE_MANAGER_HPP
#define SEQUANT_EVAL_CACHE_MANAGER_HPP

#include <SeQuant/core/asy_cost.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/eval/eval_node_compare.hpp>
#include <SeQuant/core/eval/fwd.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/expr.hpp>

#include <memory>
#include <range/v3/view.hpp>
#include <unordered_map>

namespace sequant {

///
/// This class implements a cache manager useful for the cases when the number
/// of times the cached objects will be accessed is known.
///
/// \tparam TreeNode The evaluation tree node type used as the cache key.
/// \tparam force_hash_collisions If true, forces all hash values to 0 (for
///         testing collision safety).
///
template <typename TreeNode, bool force_hash_collisions>
class CacheManager {
 public:
  using key_type = TreeNode;

 private:
  using hasher_type = TreeNodeHasher<TreeNode, force_hash_collisions>;
  using comparator_type = TreeNodeEqualityComparator<TreeNode>;

  class entry {
   private:
    size_t max_life;

    size_t life_c;

    ResultPtr data_p;

   public:
    explicit entry(size_t count) noexcept
        : max_life{count}, life_c{count}, data_p{nullptr} {}

    [[nodiscard]] ResultPtr access() noexcept {
      if (!data_p) return nullptr;
      return decay() == 0 ? std::move(data_p) : data_p;
    }

    void store(ResultPtr&& data) noexcept { data_p = std::move(data); }

    void reset() noexcept {
      life_c = max_life;
      data_p = nullptr;
    }

    [[nodiscard]] size_t life_count() const noexcept { return life_c; }

    [[nodiscard]] size_t max_life_count() const noexcept { return max_life; }

    [[nodiscard]] size_t size_in_bytes() const noexcept {
      return data_p ? data_p->size_in_bytes() : 0;
    }

    [[nodiscard]] bool alive() const noexcept { return data_p ? true : false; }

   private:
    [[nodiscard]] int decay() noexcept {
      return life_c > 0 ? static_cast<int>(--life_c) : 0;
    }

  };  // entry

  static ResultPtr store(entry& ent, ResultPtr&& data) noexcept {
    ent.store(std::move(data));
    return ent.access();
  }

  std::unordered_map<TreeNode, entry, hasher_type, comparator_type> cache_map_;

 public:
  template <typename Iterable1>
  explicit CacheManager(Iterable1&& decaying) noexcept {
    for (auto&& [k, c] : decaying) cache_map_.try_emplace(k, entry{c});
  }

  ///
  /// Resets all cached data.
  ///
  void reset() noexcept {
    for (auto&& [k, v] : cache_map_) v.reset();
  }

  ///
  /// @brief Access cached data.
  ///
  /// @param key The key that identifies the cached data.
  /// @return ResultPtr to Result
  ResultPtr access(key_type const& key) noexcept {
    if (auto found = cache_map_.find(key); found != cache_map_.end())
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
  [[nodiscard]] ResultPtr store(key_type const& key, ResultPtr data) noexcept {
    if (auto found = cache_map_.find(key); found != cache_map_.end())
      return store(found->second, std::move(data));
    return data;
  }

  ///
  /// \brief Check if the key exists in the database: does not check if cache
  ///        exists
  ///
  [[nodiscard]] bool exists(key_type const& key) const noexcept {
    return cache_map_.find(key) != cache_map_.end();
  }

  /// if the key exists in the database, return the current lifetime count of
  /// the cached data otherwise return -1
  [[nodiscard]] int life(key_type const& key) const noexcept {
    auto iter = cache_map_.find(key);
    auto end = cache_map_.end();
    return iter == end ? -1 : static_cast<int>(iter->second.life_count());
  }

  /// if the key exists in the database, return the maximum lifetime count of
  /// the cached data that implies the maximum number of accesses allowed for
  /// this key before the cache is released. This value was set by the c'tor.
  [[nodiscard]] int max_life(key_type const& key) const noexcept {
    auto iter = cache_map_.find(key);
    auto end = cache_map_.end();
    return iter == end ? -1 : static_cast<int>(iter->second.max_life_count());
  }

  ///
  /// \return The number of entries with life_count greater than zero.
  ///
  [[nodiscard]] size_t alive_count() const noexcept {
    using ranges::views::filter;
    using ranges::views::transform;
    using ranges::views::values;
    return ranges::accumulate(cache_map_                            //
                                  | values                          //
                                  | filter(&entry::alive)           //
                                  | transform(&entry::life_count),  //
                              size_t{0});
  }

  ///
  /// \return Returns the sum of `Result::size_in_bytes` of alive entries.
  ///
  [[nodiscard]] size_t size_in_bytes() const noexcept {
    using ranges::views::transform;
    using ranges::views::values;
    return ranges::accumulate(
        cache_map_ | values | transform(&entry::size_in_bytes), size_t{0});
  }

  ///
  /// Get an empty cache manager.
  ///
  static CacheManager empty() noexcept {
    using map_type =
        std::unordered_map<TreeNode, size_t, hasher_type, comparator_type>;
    return CacheManager{map_type{}};
  }

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
template <bool force_hash_collisions = false>
auto cache_manager(meta::eval_node_range auto const& nodes,
                   size_t min_repeats = 2) noexcept {
  using TreeNode =
      std::ranges::range_value_t<std::remove_cvref_t<decltype(nodes)>>;
  using Hasher = TreeNodeHasher<TreeNode, force_hash_collisions>;
  using Comp = TreeNodeEqualityComparator<TreeNode>;

  // Phase 1: Scan with pointer-based map (low memory)
  std::unordered_map<const TreeNode*, size_t, Hasher, Comp> imed_counts;

  auto imed_visitor = [&imed_counts](auto&& n) -> bool {
    if (auto found = imed_counts.find(&n); found != imed_counts.end()) {
      ++found->second;
      return false;
    }
    imed_counts.emplace(&n, 1);
    return true;
  };

  ranges::for_each(nodes, [&imed_visitor](auto&& tree) {
    tree.visit_internal(imed_visitor);
  });

  // Phase 2: Copy repeated entries (node by value)
  std::unordered_map<TreeNode, size_t, Hasher, Comp> filtered;
  for (auto&& [ptr, count] : imed_counts) {
    if (count >= min_repeats) filtered.emplace(*ptr, count);
  }

  return CacheManager<TreeNode, force_hash_collisions>{std::move(filtered)};
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
