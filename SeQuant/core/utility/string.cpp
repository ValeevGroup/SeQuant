#include "string.hpp"

#include <codecvt>
#include <locale>

namespace sequant {

std::string toUtf8(const std::wstring_view str) {
  std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;

  return converter.to_bytes(str.begin(), str.end());
}

std::wstring toUtf16(const std::string_view str) {
  std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;

  return converter.from_bytes(str.begin(), str.end());
}

}  // namespace sequant
