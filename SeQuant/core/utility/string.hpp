//
// Created by Robert Adam on 2023-10-10
//

#ifndef SEQUANT_CORE_UTILITY_STRING_HPP
#define SEQUANT_CORE_UTILITY_STRING_HPP

#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <cmath>
#include <functional>
#include <string>
#include <string_view>

// All strings in SeQuant are UTF-8 encoded narrow strings, which requires
// narrow string literals to be UTF-8 encoded as well
static_assert(sizeof("α") == 3,
              "SeQuant requires the execution character set to be UTF-8");

namespace sequant {

namespace meta {

template <typename T>
constexpr inline bool is_string_or_view_v =
    std::is_same_v<std::remove_cvref_t<T>, std::string> ||
    std::is_same_v<std::remove_cvref_t<T>, std::string_view>;

template <typename T>
constexpr inline bool is_string_convertible_v =
    is_string_or_view_v<T> || std::is_same_v<std::remove_cvref_t<T>, char[]> ||
    std::is_same_v<std::remove_all_extents_t<std::remove_cvref_t<T>>, char> ||
    std::is_same_v<std::remove_cvref_t<T>, char *> ||
    std::is_same_v<std::remove_cvref_t<T>, const char *> ||
    std::is_same_v<std::remove_cvref_t<T>, char>;

}  // namespace meta

template <typename T>
concept string_convertible =
    meta::is_string_convertible_v<std::remove_cvref_t<T>>;

/// @return true if @p str is a valid UTF-8 encoded string
bool is_valid_utf8(std::string_view str);

/// Decodes the given UTF-8 encoded string into its code points
/// @pre @p str is valid UTF-8
std::u32string to_utf32(std::string_view str);

/// Encodes the given code points as UTF-8
std::string to_utf8(std::u32string_view str);

/// Encodes the given code point as UTF-8
std::string to_utf8(char32_t code_point);

/// @return the UTF-8 encoded first code point of @p str (empty if @p str is
/// empty)
/// @pre @p str is valid UTF-8
std::string_view first_code_point(std::string_view str);

/// @return the UTF-8 encoded last code point of @p str (empty if @p str is
/// empty)
/// @pre @p str is valid UTF-8
std::string_view last_code_point(std::string_view str);

/// Translates a regular expression (ECMAScript syntax) written in terms of
/// (Unicode) characters into an equivalent one that can be used with
/// `std::regex` (or `boost::regex`) to match UTF-8 encoded strings. In the
/// returned pattern, `.`, negated bracket expressions and `\D`, `\W` and `\S`
/// match entire code points (rather than individual bytes), quantifiers
/// apply to entire (multi-byte) characters and bracket expressions may contain
/// non-ASCII characters.
/// @param pattern UTF-8 encoded regular expression
/// @return the translated regular expression
/// @throw Exception if @p pattern contains ranges of non-ASCII
/// characters or negated character class escapes inside bracket expressions
/// (neither is supported)
std::string utf8_regex_pattern(std::string_view pattern);

/// @return true if @p c is an ASCII digit; unlike std::isdigit, this is safe to
/// use on the individual code units of a UTF-8 encoded string
constexpr bool is_ascii_digit(char c) { return c >= '0' && c <= '9'; }

/// @return true if @p c is an ASCII letter; unlike std::isalpha, this is
/// locale-independent and safe to use on the individual code units of a UTF-8
/// encoded string
constexpr bool is_ascii_letter(char c) {
  return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}

template <string_convertible S>
std::string_view to_string_view(S &&str) {
  if constexpr (meta::is_char_v<std::remove_reference_t<S>>)
    return {&str, 1};
  else
    return str;
}

/// Joins the elements of a range into a string
/// @param rng the range to join
/// @param separator inserted between consecutive elements
/// @param proj maps each element to something appendable to a `std::string`
/// @return the projected elements of @p rng, separated by @p separator
template <typename Range, typename Proj = std::identity>
std::string join_strings(Range &&rng, std::string_view separator,
                         Proj proj = {}) {
  std::string result;
  bool first = true;
  for (auto &&elem : rng) {
    if (!first) result += separator;
    first = false;
    result += std::invoke(proj, elem);
  }
  return result;
}

/// Converts integral type to its std::string representation
std::string to_string(meta::integral auto t) { return std::to_string(t); }

/// Converts real type to its std::string representation, converting to integer
/// if possible
std::string to_string(std::floating_point auto t) {
  if (std::floor(t) == t)
    return std::to_string(static_cast<long long>(t));
  else
    return std::to_string(t);
}

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_STRING_HPP
