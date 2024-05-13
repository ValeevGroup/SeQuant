//
// Created by Eduard Valeyev on 2019-02-11.
//

#ifndef SEQUANT_WSTRING_HPP
#define SEQUANT_WSTRING_HPP

#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

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
/// Converts a UTF-8 encoded string (C or C++) to a UTF-8 encoded
/// std::string
template <typename S, typename = meta::EnableIfAnyStringConvertible<S>>
std::string to_string(S&& str_utf8) {
  auto str_utf8_view = to_basic_string_view(std::forward<S>(str_utf8));
  using boost::locale::conv::utf_to_utf;
  return utf_to_utf<char>(str_utf8_view.data(),
                          str_utf8_view.data() + str_utf8_view.size());
}

/// Optimized to_string for std::string
inline std::string to_string(std::string&& str_utf8) {
  return std::move(str_utf8);
}

/// @brief (potentially) narrowing character converter.
///
/// Converts a UTF-8 encoded std::basic_string_view<Char> to a UTF-8 encoded
/// std::wstring
template <typename S, typename = meta::EnableIfAnyStringConvertible<S>>
std::wstring to_wstring(S&& str_utf8) {
  auto str_utf8_view = to_basic_string_view(std::forward<S>(str_utf8));
  using boost::locale::conv::utf_to_utf;
  return utf_to_utf<wchar_t>(str_utf8_view.data(),
                             str_utf8_view.data() + str_utf8_view.size());
}

/// Optimized to_wstring for std::wstring
inline std::wstring to_wstring(std::wstring&& str_utf8) {
  return std::move(str_utf8);
}

#if __cplusplus >= 202002L

namespace detail {

/// selects the string literal matching the code unit given by @p Char
template <typename Char, std::size_t NChar, std::size_t NWChar,
          std::size_t NChar8, std::size_t NChar16, std::size_t NChar32>
constexpr decltype(auto) select_string_literal(
    const char (&str)[NChar], const wchar_t (&wstr)[NWChar],
    const char8_t (&u8str)[NChar8], const char16_t (&u16str)[NChar16],
    const char32_t (&u32str)[NChar32]) {
  if constexpr (std::is_same_v<Char, char>)
    return str;
  else if constexpr (std::is_same_v<Char, wchar_t>)
    return wstr;
  else if constexpr (std::is_same_v<Char, char8_t>)
    return u8str;
  else if constexpr (std::is_same_v<Char, char16_t>)
    return u16str;
  else if constexpr (std::is_same_v<Char, char32_t>)
    return u32str;
  else
    abort();
}
}  // namespace detail

/// @brief Converts a base string literal to the desired code unit type
/// @param code_unit The desired code unit type; either `char` or `wchar_t`
/// @param char_string_literal `char`-based string literal
#define SQ_STRLIT(code_unit, char_string_literal)                  \
  detail::select_string_literal<code_unit>(                        \
      char_string_literal, SEQUANT_CONCAT(L, char_string_literal), \
      SEQUANT_CONCAT(u8, char_string_literal),                     \
      SEQUANT_CONCAT(u, char_string_literal),                      \
      SEQUANT_CONCAT(U, char_string_literal))

#else  // __cplusplus < 202002L

namespace detail {

/// selects the string literal matching the code unit given by @p Char
template <typename Char, std::size_t NChar, std::size_t NWChar>
constexpr decltype(auto) select_string_literal(const char (&str)[NChar],
                                               const wchar_t (&wstr)[NWChar]) {
  if constexpr (std::is_same_v<Char, char>)
    return str;
  else if constexpr (std::is_same_v<Char, wchar_t>)
    return wstr;
  else
    abort();
}

}  // namespace detail

/// @brief Converts a base string literal to the desired code unit type
/// @param code_unit The desired code unit type; either `char` or `wchar_t`
/// @param char_string_literal `char`-based string literal
#define SQ_STRLIT(code_unit, char_string_literal) \
  detail::select_string_literal<code_unit>(       \
      char_string_literal, SEQUANT_CONCAT(L, char_string_literal))

#endif  // __cplusplus < 202002L

}  // namespace sequant

#endif  // SEQUANT_WSTRING_HPP
