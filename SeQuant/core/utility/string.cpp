#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <iterator>
#include <stdexcept>
#include <string>

#include <utf8/checked.h>
#include <utf8/unchecked.h>

namespace sequant {

bool is_valid_utf8(std::string_view str) { return utf8::is_valid(str); }

void validate_utf8(std::string_view str) {
  if (!is_valid_utf8(str)) {
    throw std::invalid_argument(
        "String is not valid UTF-8 (invalid byte at "
        "offset " +
        std::to_string(utf8::find_invalid(str)) + ")");
  }
}

std::u32string to_utf32(std::string_view str) {
  SEQUANT_ASSERT(is_valid_utf8(str));
  std::u32string result;
  utf8::unchecked::utf8to32(str.begin(), str.end(), std::back_inserter(result));
  return result;
}

std::string to_utf8(std::u32string_view str) {
  std::string result;
  utf8::utf32to8(str.begin(), str.end(), std::back_inserter(result));
  return result;
}

std::string to_utf8(char32_t code_point) {
  return to_utf8(std::u32string_view(&code_point, 1));
}

std::string_view first_code_point(std::string_view str) {
  SEQUANT_ASSERT(is_valid_utf8(str));
  if (str.empty()) return str;
  auto it = str.begin();
  utf8::unchecked::next(it);
  return str.substr(0, std::distance(str.begin(), it));
}

std::string_view last_code_point(std::string_view str) {
  SEQUANT_ASSERT(is_valid_utf8(str));
  if (str.empty()) return str;
  auto it = str.end();
  utf8::unchecked::prior(it);
  return str.substr(std::distance(str.begin(), it));
}

}  // namespace sequant
