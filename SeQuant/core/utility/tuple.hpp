//
// Created by Eduard Valeyev on 1/29/25.
//

#ifndef SEQUANT_CORE_UTILITY_TUPLE_HPP
#define SEQUANT_CORE_UTILITY_TUPLE_HPP

#include <functional>
#include <tuple>

namespace sequant::detail {

/// compares Kth elements of two tuples
template <std::size_t K>
struct tuple_less {
  template <typename... Ts>
  std::enable_if_t<std::less{}(K, sizeof...(Ts)), bool> operator()(
      const std::tuple<Ts...> &a, const std::tuple<Ts...> &b) const {
    return std::get<K>(a) < std::get<K>(b);
  }
};

}  // namespace sequant::detail

#endif  // SEQUANT_CORE_UTILITY_STRONG_HPP
