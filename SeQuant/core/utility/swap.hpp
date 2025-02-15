//
// Created by Eduard Valeyev on 2/14/25.
//

#ifndef SEQUANT_CORE_UTILITY_SWAPPABLE_HPP
#define SEQUANT_CORE_UTILITY_SWAPPABLE_HPP

#include <atomic>

namespace sequant {

/// use this as a swap-countable deep wrapper for types whose swaps are not
/// counted
template <typename T>
class SwapCountable;

/// use this as a swap-countable shallow wrapper for types whose swaps are not
/// counted
template <typename T>
class SwapCountableRef;

/// counted swap for types whose default swap is not counted
template <typename T>
inline void counted_swap(T& a, T& b);

namespace detail {

/// use this for implementing custom swap that counts swaps by default
template <typename T>
inline void count_swap();

/// atomic counter, used by swap overloads for SwapCountable and
/// SwapCountableRef
template <typename T>
struct SwapCounter {
  SwapCounter() : even_num_of_swaps_(true) {}
  static SwapCounter& thread_instance() {
    static thread_local SwapCounter instance_{};
    return instance_;
  }

  bool even_num_of_swaps() const { return even_num_of_swaps_; }
  void reset() { even_num_of_swaps_ = true; }

 private:
  std::atomic<bool> even_num_of_swaps_;
  void toggle() { even_num_of_swaps_ = !even_num_of_swaps_; }

  friend class SwapCountable<T>;
  friend class SwapCountableRef<T>;
  friend inline void counted_swap<T>(T& a, T& b);
  friend inline void count_swap<T>();
};

template <typename T>
inline void count_swap() {
  detail::SwapCounter<T>::thread_instance().toggle();
}

}  // namespace detail

/// Wraps `T` to make its swap countable
template <typename T>
class SwapCountable {
 public:
  template <typename U>
  explicit SwapCountable(U&& v) : value_(std::forward<U>(v)) {}

 private:
  T value_;

  friend inline void swap(SwapCountable& a, SwapCountable& b) {
    using std::swap;
    swap(a.value_, b.value_);
    detail::SwapCounter<T>::thread_instance().toggle();
  }

  friend inline bool operator<(const SwapCountable& a, const SwapCountable& b) {
    return a.value_ < b.value_;
  }
};

/// Wraps `T&` to make its swap countable
template <typename T>
class SwapCountableRef {
 public:
  explicit SwapCountableRef(T& ref) : ref_(ref) {}

 private:
  T& ref_;

  // NB swapping const wrappers swaps the payload
  friend inline void swap(const SwapCountableRef& a,
                          const SwapCountableRef& b) {
    using std::swap;
    swap(a.ref_, b.ref_);
    detail::SwapCounter<T>::thread_instance().toggle();
  }

  friend inline bool operator<(const SwapCountableRef& a,
                               const SwapCountableRef& b) {
    return a.ref_ < b.ref_;
  }
};

template <typename T>
void reset_ts_swap_counter() {
  detail::SwapCounter<T>::thread_instance().reset();
}

template <typename T>
bool ts_swap_counter_is_even() {
  return detail::SwapCounter<T>::thread_instance().even_num_of_swaps();
}

template <typename T>
inline void counted_swap(T& a, T& b) {
  using std::swap;
  swap(a, b);
  detail::SwapCounter<T>::thread_instance().toggle();
}

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_SWAPPABLE_HPP
