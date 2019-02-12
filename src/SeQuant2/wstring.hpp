//
// Created by Eduard Valeyev on 2019-02-11.
//

#ifndef SEQUANT2_WSTRING_HPP
#define SEQUANT2_WSTRING_HPP

#include <type_traits>

namespace sequant2 {

/// Converts integral type to its std::wstring representation
template <typename T>
std::enable_if_t<std::is_integral_v<std::decay_t<T>>, std::wstring> to_wstring(
    T&& t) {
  return std::to_wstring(t);
}

/// Converts real type to its std::wstring representation, converting to integer
/// if possible
template <typename T>
std::enable_if_t<std::is_floating_point_v<std::decay_t<T>>, std::wstring>
to_wstring(T&& t) {
  if (std::floor(t) == t)
    return std::to_wstring(static_cast<long long>(t));
  else
    return std::to_wstring(t);
}

}  // namespace sequant2

#endif  // SEQUANT2_WSTRING_HPP
