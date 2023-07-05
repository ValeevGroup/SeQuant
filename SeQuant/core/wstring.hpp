//
// Created by Eduard Valeyev on 2019-02-11.
//

#ifndef SEQUANT_WSTRING_HPP
#define SEQUANT_WSTRING_HPP

#include <cmath>
#include <string>
#include <string_view>
#include <type_traits>

#include <boost/locale/encoding_utf.hpp>

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
/// Converts a UTF-8 encoded std::basic_string_view<Char> to a UTF-8 encoded
/// std::string
template <typename Char, typename Traits>
std::string to_string(
    const std::basic_string_view<Char, Traits>& str_utf8_view) {
  using boost::locale::conv::utf_to_utf;
  return utf_to_utf<char>(str_utf8_view.data(),
                          str_utf8_view.data() + str_utf8_view.size());
}

/// @brief (potentially) narrowing character converter.
///
/// Converts a UTF-8 encoded std::basic_string_view<Char> to a UTF-8 encoded
/// std::string
template <typename Char, typename Traits, typename Allocator>
std::string to_string(
    const std::basic_string<Char, Traits, Allocator>& str_utf8) {
  using boost::locale::conv::utf_to_utf;
  return utf_to_utf<char>(str_utf8.data(), str_utf8.data() + str_utf8.size());
}

/// @brief (potentially) narrowing character converter.
///
/// Converts a UTF-8 encoded std::basic_string_view<Char> to a UTF-8 encoded
/// std::wstring
template <typename SourceChar, typename SourceTraits>
std::wstring to_wstring(
    const std::basic_string_view<SourceChar, SourceTraits>& str_utf8_view) {
  using boost::locale::conv::utf_to_utf;
  return utf_to_utf<wchar_t>(str_utf8_view.data(),
                             str_utf8_view.data() + str_utf8_view.size());
}

/// @brief (potentially) narrowing character converter.
///
/// Converts a UTF-8 encoded std::basic_string_view<Char> to a UTF-8 encoded
/// std::wstring
template <typename SourceChar, typename SourceTraits, typename SourceAllocator>
std::wstring to_wstring(const std::basic_string<SourceChar, SourceTraits,
                                                SourceAllocator>& str_utf8) {
  using boost::locale::conv::utf_to_utf;
  return utf_to_utf<wchar_t>(str_utf8.data(),
                             str_utf8.data() + str_utf8.size());
}

}  // namespace sequant

#endif  // SEQUANT_WSTRING_HPP
