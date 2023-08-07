//
// Created by Eduard Valeyev on 8/1/23.
//

#ifndef SEQUANT_CORE_UTILITY_NODISCARD_HPP
#define SEQUANT_CORE_UTILITY_NODISCARD_HPP

#include <utility>

namespace sequant::detail {

/// can be used to wrap a callable and make its call operator [[nodiscard]]
/// @note this will not be needed in C++23
/// @internal see https://stackoverflow.com/a/43687812
template <typename F>
struct NoDiscard {
  F f;
  NoDiscard(const F& ff) : f(ff) {}
  NoDiscard(F&& ff) : f(std::move(ff)) {}

  template <typename... T>
  [[nodiscard]] constexpr auto operator()(T&&... t) const
      noexcept(noexcept(f(std::forward<T>(t)...))) {
    return f(std::forward<T>(t)...);
  }
};

}  // namespace sequant::detail

#endif  // SEQUANT_CORE_UTILITY_NODISCARD_HPP
