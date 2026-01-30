//
// Created by Robert Adam on 2023-10-10
//

#ifndef SEQUANT_CORE_UTILITY_STRING_HPP
#define SEQUANT_CORE_UTILITY_STRING_HPP

#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <string>
#include <string_view>

namespace sequant {

namespace meta {

template <typename T>
constexpr inline bool is_string_or_view_v =
    std::is_same_v<std::remove_cvref_t<T>, std::string> ||
    std::is_same_v<std::remove_cvref_t<T>, std::string_view>;

template <typename T>
constexpr inline bool is_wstring_or_view_v =
    std::is_same_v<std::remove_cvref_t<T>, std::wstring> ||
    std::is_same_v<std::remove_cvref_t<T>, std::wstring_view>;

template <typename T>
constexpr inline bool is_string_convertible_v =
    is_string_or_view_v<T> || std::is_same_v<std::remove_cvref_t<T>, char[]> ||
    std::is_same_v<std::remove_all_extents_t<std::remove_cvref_t<T>>, char> ||
    std::is_same_v<std::remove_cvref_t<T>, char *> ||
    std::is_same_v<std::remove_cvref_t<T>, const char *> ||
    std::is_same_v<std::remove_cvref_t<T>, char>;

template <typename T>
constexpr inline bool is_wstring_convertible_v =
    is_wstring_or_view_v<T> ||
    std::is_same_v<std::remove_cvref_t<T>, wchar_t[]> ||
    std::is_same_v<std::remove_all_extents_t<std::remove_cvref_t<T>>,
                   wchar_t> ||
    std::is_same_v<std::remove_cvref_t<T>, wchar_t *> ||
    std::is_same_v<std::remove_cvref_t<T>, const wchar_t *> ||
    std::is_same_v<std::remove_cvref_t<T>, wchar_t>;

template <typename T>
constexpr inline bool is_basic_string_convertible_v =
    is_string_convertible_v<T> || is_wstring_convertible_v<T>;

template <typename T, typename Enabler = void>
struct char_type;

template <typename T>
struct char_type<T, std::enable_if_t<is_basic_string_convertible_v<T>>> {
  using type = std::conditional_t<is_string_convertible_v<T>, char, wchar_t>;
};

template <typename T>
using char_t = typename char_type<T>::type;

}  // namespace meta

template <typename T>
concept basic_string_convertible =
    meta::is_basic_string_convertible_v<std::remove_cvref_t<T>>;

/// Converts the given wide-string to a UTF-8 encoded narrow string
std::string toUtf8(std::wstring_view str);

/// Pass-through (only converts to std::string)
std::string toUtf8(std::string_view str);

inline std::string toUtf8(char c) { return toUtf8(std::string_view(&c, 1)); }

inline std::string toUtf8(wchar_t c) {
  return toUtf8(std::wstring_view(&c, 1));
}

/// Converts the given UTF-8 encoded narrow-string to a UTF-16 encoded
/// wide-string
std::wstring toUtf16(std::string_view str);

/// Pass-through (only converts to std::wstring)
std::wstring toUtf16(std::wstring_view str);

inline std::wstring toUtf16(char c) { return toUtf16(std::string_view(&c, 1)); }

inline std::wstring toUtf16(wchar_t c) {
  return toUtf16(std::wstring_view(&c, 1));
}

template <basic_string_convertible S>
std::basic_string_view<meta::char_t<S>> to_basic_string_view(S &&str) {
  if constexpr (meta::is_char_v<std::remove_reference_t<S>>)
    return {&str, 1};
  else
    return str;
}

/// Converts integral type to its std::wstring representation
std::wstring to_wstring(meta::integral auto t) { return std::to_wstring(t); }

/// Converts real type to its std::wstring representation, converting to integer
/// if possible
std::wstring to_wstring(std::floating_point auto t) {
  if (std::floor(t) == t)
    return std::to_wstring(static_cast<long long>(t));
  else
    return std::to_wstring(t);
}

#if 0
/// @brief (potentially) narrowing character converter.
///
/// Converts a UTF-8 encoded string (C or C++) to a UTF-8 encoded
/// std::string
template <basic_string_convertible S>
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
template <basic_string_convertible S>
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
#endif

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

  SEQUANT_ABORT("Unhandled character type");
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

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_STRING_HPP
