//
// Created by Robert Adam on 2023-10-10
//

#ifndef SEQUANT_CORE_UTILITY_STRING_HPP
#define SEQUANT_CORE_UTILITY_STRING_HPP

#include <SeQuant/core/meta.hpp>
#include <string>
#include <string_view>

namespace sequant {

namespace meta {

template <typename T>
constexpr inline bool is_string_or_view_v =
    std::is_same_v<remove_cvref_t<T>, std::string> ||
    std::is_same_v<remove_cvref_t<T>, std::string_view>;

template <typename T>
constexpr inline bool is_wstring_or_view_v =
    std::is_same_v<remove_cvref_t<T>, std::wstring> ||
    std::is_same_v<remove_cvref_t<T>, std::wstring_view>;

template <typename T>
constexpr inline bool is_string_convertible_v =
    is_string_or_view_v<T> || std::is_same_v<remove_cvref_t<T>, char[]> ||
    std::is_same_v<remove_cvref_t<T>, const char[]> ||
    std::is_same_v<remove_cvref_t<T>, char *> ||
    std::is_same_v<remove_cvref_t<T>, const char *> ||
    std::is_same_v<remove_cvref_t<T>, char>;

template <typename T>
constexpr inline bool is_wstring_convertible_v =
    is_wstring_or_view_v<T> || std::is_same_v<remove_cvref_t<T>, wchar_t[]> ||
    std::is_same_v<remove_cvref_t<T>, const wchar_t[]> ||
    std::is_same_v<remove_cvref_t<T>, wchar_t *> ||
    std::is_same_v<remove_cvref_t<T>, const wchar_t *> ||
    std::is_same_v<remove_cvref_t<T>, wchar_t>;

}  // namespace meta

/// Converts the given wide-string to a UTF-8 encoded narrow string
std::string toUtf8(std::wstring_view str);

/// Converts the given UTF-8 encoded narrow-string to a UTF-16 encoded
/// wide-string
std::wstring toUtf16(std::string_view str);

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_STRING_HPP
