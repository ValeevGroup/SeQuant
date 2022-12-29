//
// Created by Eduard Valeyev on 2019-02-11.
//

#ifndef SEQUANT_WSTRING_HPP
#define SEQUANT_WSTRING_HPP

#include <cmath>
#include <string>
#include <type_traits>

namespace sequant {

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

/// @brief (potentially) narrowing character converter.
///
/// Converts a UTF-8 encoded std::wstring to a UTF-8 encoded std::string
std::string to_string(const std::wstring& wstr_utf8);

}  // namespace sequant

#endif  // SEQUANT_WSTRING_HPP
