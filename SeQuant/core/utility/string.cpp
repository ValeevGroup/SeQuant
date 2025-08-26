#include <SeQuant/core/utility/string.hpp>

#include <algorithm>
#include <cassert>
#include <limits>
#include <string>

#include <utfcpp-src/source/utf8.h>

namespace sequant {

std::string toUtf8(const std::wstring_view str) {
  assert(std::none_of(str.begin(), str.end(), [](wchar_t c) {
    return c > std::numeric_limits<char16_t>::max();
  }));

  std::u16string u16str(str.begin(), str.end());

  auto stream = utf8::utf16to8(u16str);

  return std::string(stream.begin(), stream.end());
}

std::wstring toUtf16(const std::string_view str) {
  auto stream = utf8::utf8to16(str);

  static_assert(sizeof(wchar_t) >= sizeof(char16_t));

  return std::wstring(stream.begin(), stream.end());
}

}  // namespace sequant
