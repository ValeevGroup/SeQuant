#ifndef SEQUANT_UTILITY_CONVERSION_HPP
#define SEQUANT_UTILITY_CONVERSION_HPP

#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <algorithm>
#include <charconv>
#include <concepts>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>
#include <typeinfo>

namespace sequant {

class ConversionException : public Exception {
 public:
  using Exception::Exception;
};

namespace detail {
template <typename T, typename Arg>
  requires(std::integral<T> || std::floating_point<T>)
T string_to_impl(std::string_view str, Arg &&arg) {
  T val = 0;

  std::from_chars_result result = std::from_chars(
      str.data(), str.data() + str.size(), val, std::forward<Arg>(arg));

  // Note: Beware that typeid(T).name() yields something implementation defined,
  // which may or may not be useful for a human
  if (result.ec != std::errc{}) {
    switch (result.ec) {
      case std::errc::invalid_argument:
        throw ConversionException("'" + std::string(str) +
                                  "' is not a valid '" + typeid(T).name() +
                                  "'");
      case std::errc::result_out_of_range:
        throw ConversionException("'" + std::string(str) +
                                  "' is out of range for type '" +
                                  typeid(T).name() + "'");
      default:
        throw ConversionException("Unexpected conversion error");
    }
  }

  if (result.ptr != str.data() + str.size()) {
    // Incomplete parse
    throw ConversionException("'" + std::string(str) +
                              "' could not be fully parsed as a '" +
                              typeid(T).name() + "'");
  }

  return val;
}

}  // namespace detail

template <typename T>
concept from_chars_supports =
    requires(const char *c, T &v) { std::from_chars(c, c + 1, v); };

/// Converts the provided string to the desired integral type.
/// Contrary to standard functions, this function does perform explicit
/// error checking and will throw a ConversionException if the provided
/// string can't be parsed in its entirety (no partial conversions, no
/// whitespace skipping!). The function will also throw if the parsed
/// value can't be represented within the range of T.
///
/// @tparam T The integral type one wishes to obtain
/// @param str The string to convert
/// @param base The base used to represent the value in its string format
/// @returns The parsed value
///
/// @throws ConversionException in case the conversion is unsuccessful or
/// the parsed value can't be represented as a T.
template <std::integral T>
T string_to(std::string_view str, int base = 10) {
  static_assert(from_chars_supports<T>,
                "Your C++ standard library is missing a std::from_chars "
                "implementation for this integral type");
  return detail::string_to_impl<T>(str, base);
}

template <std::integral T>
T string_to(std::wstring_view str, int base = 10) {
  if (std::ranges::any_of(str, [](wchar_t c) {
        return c > std::numeric_limits<char>::max();
      })) {
    throw ConversionException(
        "'" + toUtf8(str) +
        "' contains characters that can't possibly represent a numeric value");
  }

  std::string buff;
  buff.resize(str.size());
  std::ranges::transform(str, buff.begin(),
                         [](wchar_t c) { return static_cast<char>(c); });

  return string_to<T>(buff, base);
}

/// Converts the provided string to the desired floating-point type.
/// Contrary to standard functions, this function does perform explicit
/// error checking and will throw a ConversionException if the provided
/// string can't be parsed in its entirety (no partial conversions, no
/// whitespace skipping!). The function will also throw if the parsed
/// value can't be represented within the range of T.
///
/// @tparam T The float-point type one wishes to obtain
/// @param str The string to convert
/// @param base The base used to represent the value in its string format
/// @returns The parsed value
///
/// @throws ConversionException in case the conversion is unsuccessful or
/// the parsed value can't be represented as a T.
template <std::floating_point T>
T string_to(std::string_view str,
            std::chars_format fmt = std::chars_format::general) {
#ifndef __APPLE__
  static_assert(from_chars_supports<T>,
                "Your C++ standard library is missing a std::from_chars "
                "implementation for this floating point type");

  return detail::string_to_impl<T>(str, fmt);
#else
  // For some reason it seems that Apple is unable to supply an implementation
  // of std::from_chars so we need to work around its (potential) absence
  if constexpr (from_chars_supports<T>) {
    // In case they update their standard lib…
    return detail::string_to_impl<T>(str, fmt);
  }

  // Workaround implementation that doesn't do any error checking - not great
  // but better than not being able to use this function at all
  std::stringstream stream{std::string(str)};
  T val = 0;
  stream >> val;

  return val;
#endif
}

template <std::floating_point T>
T string_to(std::wstring_view str,
            std::chars_format fmt = std::chars_format::general) {
  if (std::ranges::any_of(str, [](wchar_t c) {
        return c > std::numeric_limits<char>::max();
      })) {
    throw ConversionException(
        "'" + toUtf8(str) +
        "' contains characters that can't possibly represent a numeric value");
  }

  std::string buff;
  buff.resize(str.size());
  std::ranges::transform(str, buff.begin(),
                         [](wchar_t c) { return static_cast<char>(c); });

  return string_to<T>(buff, fmt);
}

/// Converts the provided string to an exact rational number. Accepted are an
/// optionally signed integer ("-12"), decimal ("3.25", ".5") or fraction
/// ("-1/3") of arbitrary size, which is converted without rounding (e.g.
/// "0.1" yields 1/10). As for the other overloads, the string must be parsable
/// in its entirety (no whitespace skipping).
///
/// @param str The string to convert
/// @returns The parsed value
///
/// @throws ConversionException in case the conversion is unsuccessful
template <std::same_as<rational> T>
T string_to(std::string_view str) {
  auto fail = [&str](const std::string &reason) -> ConversionException {
    return ConversionException("'" + std::string(str) +
                               "' is not a valid rational number: " + reason);
  };
  auto is_digits = [](std::string_view digits) {
    return std::ranges::all_of(digits,
                               [](char c) { return c >= '0' && c <= '9'; });
  };

  std::string_view rest = str;
  const bool negative = !rest.empty() && rest.front() == '-';
  if (!rest.empty() && (rest.front() == '-' || rest.front() == '+'))
    rest.remove_prefix(1);

  rational value;
  if (const auto slash = rest.find('/'); slash != std::string_view::npos) {
    const std::string_view num = rest.substr(0, slash);
    const std::string_view den = rest.substr(slash + 1);
    if (num.empty() || den.empty() || !is_digits(num) || !is_digits(den))
      throw fail("expected digits on both sides of '/'");
    const intmax_t denominator{std::string(den)};
    if (denominator == 0) throw fail("zero denominator");
    value = rational(intmax_t{std::string(num)}, denominator);
  } else {
    const auto dot = rest.find('.');
    const std::string_view whole = rest.substr(0, dot);
    const std::string_view fraction = dot == std::string_view::npos
                                          ? std::string_view{}
                                          : rest.substr(dot + 1);
    if (whole.empty() && fraction.empty()) throw fail("expected digits");
    if (!is_digits(whole) || !is_digits(fraction) ||
        (dot != std::string_view::npos && fraction.empty()))
      throw fail("expected an integer or decimal number");
    const std::string digits = std::string(whole) + std::string(fraction);
    intmax_t scale = 1;
    for (std::size_t i = 0; i < fraction.size(); ++i) scale *= 10;
    value = rational(intmax_t{digits}, scale);
  }

  return negative ? rational(-value) : value;
}

}  // namespace sequant

#endif
