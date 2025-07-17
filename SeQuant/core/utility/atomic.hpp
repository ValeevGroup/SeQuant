#ifndef SEQUANT_UTILITY_ATOMIC_HPP
#define SEQUANT_UTILITY_ATOMIC_HPP

#include <atomic>
#include <concepts>

namespace sequant {

/// Takes in the given atomic counter, fetches the current value and increments
/// it by the specified amount. This is done in a thread-safe manner that
/// prevents data races and hence it is guaranteed that the returned value will
/// be unique across all threads of the program.
template <std::integral T>
T fetch_and_increment(std::atomic<T> &atomic, T increment = 1) {
  T value = atomic.load();

  // This ensures that we are updating atomic with the incremented value while
  // also ensuring that no other thread has been fetching the same value as we
  // did. In other words: this ensures that we don't run into an inconsistent
  // state due to a data race.
  while (!atomic.compare_exchange_weak(value, value + increment)) {
  }

  return value;
}

}  // namespace sequant

#endif  // SEQUANT_UTILITY_ATOMIC_HPP
