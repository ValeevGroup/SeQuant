//
// Created by Eduard Valeyev on 5/17/23.
//

#ifndef SEQUANT_CORE_MATH_HPP
#define SEQUANT_CORE_MATH_HPP

#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <atomic>
#include <cstddef>

namespace sequant {

constexpr std::size_t pow2(std::size_t n) {
  SEQUANT_ASSERT(n <= 63);
  return 1ul << n;
}

/// Returns the factorial of `n`

/// @param n the argument; there is no practical limit on the supported
/// magnitude of `n`
/// @return `n!`
/// @internal uses pretabulated factorials for up to `n=20`, and memoized values
/// for `n>20`.
/// @note this function is reentrant
inline sequant::intmax_t factorial(std::size_t n) {
  using result_t = sequant::intmax_t;
  constexpr std::size_t n_max_precomputed =
      20;  // 20! is the largest factorial that fits into uint64_t
  const result_t values[n_max_precomputed + 1] = {1ull,
                                                  1ull,
                                                  2ull,
                                                  6ull,
                                                  24ull,
                                                  120ull,
                                                  720ull,
                                                  5040ull,
                                                  40320ull,
                                                  362880ull,
                                                  3628800ull,
                                                  39916800ull,
                                                  479001600ull,
                                                  6227020800ull,
                                                  87178291200ull,
                                                  1307674368000ull,
                                                  20922789888000ull,
                                                  355687428096000ull,
                                                  6402373705728000ull,
                                                  121645100408832000ull,
                                                  2432902008176640000ull};
  if (n <= n_max_precomputed)
    return values[n];
  else {
    // this memoizes values for n>n_max_precomputed
#if defined(__cpp_lib_atomic_shared_ptr) && \
    __cpp_lib_atomic_shared_ptr >= 201711L
    using ValHolder =
        std::atomic<std::shared_ptr<std::vector<sequant::intmax_t>>>;
    auto load = [](auto &&holder) { return holder.load(); };
    auto store = [](auto &&holder, ValHolder::value_type arg) {
      holder.store(std::move(arg));
    };
#else
    using ValHolder = std::shared_ptr<std::vector<sequant::intmax_t>>;
    auto load = [](auto &&holder) { return std::atomic_load(&holder); };
    auto store = [](auto &&holder, ValHolder arg) {
      std::atomic_store(&holder, std::move(arg));
    };
#endif
    static ValHolder memvals;
    // used to serialize access to memvals
    static std::mutex memvals_mutex;

    const std::size_t n_in_memvals = n - n_max_precomputed - 1;
    const auto memvals_handle = load(memvals);
    if (memvals_handle && memvals_handle->size() > n_in_memvals) {
      return (*memvals_handle)[n_in_memvals];
    } else {
      std::lock_guard<std::mutex> lock(memvals_mutex);
      // memvals might have been updated while we were waiting for the lock
      // hence reload the state
      const auto memvals_handle = load(memvals);
      if (memvals_handle && memvals_handle->size() > n_in_memvals) {
        return (*memvals_handle)[n_in_memvals];
      } else {
        auto new_memvals = std::make_shared<std::vector<sequant::intmax_t>>(
            memvals_handle
                ? *memvals_handle
                : std::vector<sequant::intmax_t>{});  // copy old memvals
        new_memvals->reserve(n_in_memvals + 1);
        for (auto i = new_memvals->size(); i <= n_in_memvals; ++i)
          new_memvals->emplace_back(
              (i == 0 ? values[n_max_precomputed] : new_memvals->back()) *
              (i + n_max_precomputed + 1));
        store(memvals, new_memvals);
        return (*new_memvals)[n_in_memvals];
      }
    }
  }
}

}  // namespace sequant

#endif  // SEQUANT_CORE_MATH_HPP
