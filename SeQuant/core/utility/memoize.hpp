#ifndef SEQUANT_CORE_UTILITY_MEMOIZE_HPP
#define SEQUANT_CORE_UTILITY_MEMOIZE_HPP

#include <SeQuant/core/container.hpp>

#include <condition_variable>
#include <mutex>
#include <optional>
#include <utility>

namespace sequant::detail {

/// Thread-safe memoization helper for expensive computations.
/// Uses a placeholder pattern: inserts nullopt first to signal "computing",
/// then fills in the value and notifies waiters.
/// @tparam Key The cache key type
/// @tparam Value The cached value type
/// @tparam ComputeFn Callable returning Value
/// @param cache The cache map (must use std::optional<Value> as values)
/// @param mutex The mutex protecting the cache
/// @param cv The condition variable for waiting
/// @param key The key to look up/store
/// @param compute_fn A callable that computes the value if not cached
/// @return Reference to the cached value
template <typename Key, typename Value, typename ComputeFn>
const Value& memoize(container::map<Key, std::optional<Value>>& cache,
                     std::mutex& mutex, std::condition_variable& cv,
                     const Key& key, ComputeFn&& compute_fn) {
  {
    std::unique_lock<std::mutex> lock(mutex);
    using cache_t = container::map<Key, std::optional<Value>>;
    typename cache_t::iterator it;
    bool inserted;
    // N.B. Use std::tie instead of structured binding to allow lambda capture
    // (C++17 does not support capturing structured bindings in lambdas)
    std::tie(it, inserted) = cache.try_emplace(key, std::nullopt);
    if (!inserted) {
      cv.wait(lock, [&it] { return it->second.has_value(); });
      return it->second.value();
    }
  }

  Value result = compute_fn();

  {
    std::lock_guard<std::mutex> lock(mutex);
    cache[key] = std::move(result);
    cv.notify_all();
    return cache[key].value();
  }
}

}  // namespace sequant::detail

#endif  // SEQUANT_CORE_UTILITY_MEMOIZE_HPP
