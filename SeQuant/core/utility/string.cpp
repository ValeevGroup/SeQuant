#include <SeQuant/core/utility/string.hpp>

#include <iterator>
#include <string>

#include <utf8/checked.h>

namespace sequant {

// a std::wstring holds UTF-16 where wchar_t is 16 bits wide (Windows) and
// UTF-32 where it is 32 bits wide; only the two differ outside the BMP
static_assert(sizeof(wchar_t) == sizeof(char16_t) ||
              sizeof(wchar_t) == sizeof(char32_t));

std::string toUtf8(const std::wstring_view str) {
  std::string result;
  if constexpr (sizeof(wchar_t) == sizeof(char16_t)) {
    utf8::utf16to8(str.begin(), str.end(), std::back_inserter(result));
  } else {
    utf8::utf32to8(str.begin(), str.end(), std::back_inserter(result));
  }
  return result;
}

std::string toUtf8(std::string_view str) { return std::string(str); }

std::wstring toUtf16(const std::string_view str) {
  std::wstring result;
  if constexpr (sizeof(wchar_t) == sizeof(char16_t)) {
    utf8::utf8to16(str.begin(), str.end(), std::back_inserter(result));
  } else {
    utf8::utf8to32(str.begin(), str.end(), std::back_inserter(result));
  }
  return result;
}

std::wstring toUtf16(std::wstring_view str) { return std::wstring(str); }

}  // namespace sequant
