#ifndef SEQUANT_UTILITY_CONVERSION_HPP
#define SEQUANT_UTILITY_CONVERSION_HPP

#include <SeQuant/core/utility/exception.hpp>

#include <charconv>
#include <concepts>
#include <string_view>
#include <system_error>
#include <typeinfo>

namespace sequant {

class ConversionException : public Exception {
 public:
  using Exception::Exception;
};

namespace {
template <typename T, typename Arg>
  requires(std::integral<T> || std::floating_point<T>)
T string_to_impl(std::string_view str, Arg &&arg) {
  T val = 0;

  std::from_chars_result result =
      std::from_chars(str.begin(), str.end(), val, std::forward<Arg>(arg));

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

  if (result.ptr != str.end()) {
    // Incomplete parse
    throw ConversionException("'" + std::string(str) +
                              "' could not be fully parsed as a '" +
                              typeid(T).name() + "'");
  }

  return val;
}

}  // namespace

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
  return string_to_impl<T>(str, base);
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
  return string_to_impl<T>(str, fmt);
}

}  // namespace sequant

#endif
