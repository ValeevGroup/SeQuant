//
// Created by Eduard Valeyev on 7/18/23.
//

#ifndef SEQUANT_CORE_BITSET_HPP
#define SEQUANT_CORE_BITSET_HPP

#include <SeQuant/core/fwd.hpp>

#include <concepts>

namespace sequant {
template <typename T>
concept convertible_to_bitset =
    requires { static_cast<bitset_t>(std::declval<T>()); };

template <convertible_to_bitset T>
bitset_t to_bitset(T e) {
  return static_cast<bitset_t>(e);
}
}  // namespace sequant

#endif  // SEQUANT_CORE_BITSET_HPP
