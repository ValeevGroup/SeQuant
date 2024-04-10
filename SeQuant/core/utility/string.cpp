#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <codecvt>
#include <locale>

namespace sequant {

SEQUANT_PRAGMA_CLANG(diagnostic push)
SEQUANT_PRAGMA_CLANG(diagnostic ignored "-Wdeprecated-declarations")
SEQUANT_PRAGMA_GCC(diagnostic push)
SEQUANT_PRAGMA_GCC(diagnostic ignored "-Wdeprecated-declarations")

std::string toUtf8(const std::wstring_view str) {
  std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;

  return converter.to_bytes(str.begin(), str.end());
}

std::wstring toUtf16(const std::string_view str) {
  std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;

  return converter.from_bytes(str.begin(), str.end());
}

SEQUANT_PRAGMA_CLANG(diagnostic pop)
SEQUANT_PRAGMA_GCC(diagnostic pop)

}  // namespace sequant
